#!/usr/bin/env python3
"""
wbeT Mutation Tracker: Detects Ogawa to Inaba serotype transitions.
Inactivation of wbeT (often via S158P, frameshifts, or truncations) 
converts the Ogawa serotype to Inaba.
"""
import json
import os
import subprocess
import gzip
from pathlib import Path
from Bio import SeqIO


def check_vcf_heterozygosity(vcf_path):
    """
    Scans VCF for heterozygous calls in the wbeT gene (CP003069.1:2687324-2688226).
    Returns {"status": "MIXED/PURE", "details": ...}
    """
    if not vcf_path or not os.path.exists(vcf_path):
        return None

    print(f"  🔍 Scanning VCF for Heterozygosity (Hikojima Check): {vcf_path}")
    
    # wbeT Coordinates (2010EL-1786)
    CHROM = "CP003069.1"
    START = 2687324
    END = 2688226
    
    mixed_sites = []
    
    try:
        # Use zgrep to quickly filter for the region to avoid full parse
        # Note: This assumes VCF is sorted and indexed, or at least small enough.
        # Since we use parsing, let's just interpret lines.
        
        # Open GZIP or Plain
        if vcf_path.endswith(".gz"):
            opener = gzip.open(vcf_path, 'rt')
        else:
            opener = open(vcf_path, 'r')
            
        with opener as f:
            for line in f:
                if line.startswith("#"): continue
                
                parts = line.split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                
                if chrom != CHROM: continue
                if pos < START: continue
                if pos > END: break # Sorted VCF assumption
                
                # We are in wbeT gene
                info_field = parts[7]
                format_field = parts[8]
                sample_field = parts[9]
                
                # Look for Allele Frequency (AF) in INFO or AD in FORMAT
                # GATK/LoFreq usually provides AF or DP4
                
                af = None
                
                # Strategy 1: AF in INFO (LoFreq style)
                if "AF=" in info_field:
                    for tag in info_field.split(';'):
                        if tag.startswith("AF="):
                            af = float(tag.split('=')[1])
                            break
                            
                # Strategy 2: AD in FORMAT (GATK style)
                # FORMAT: GT:PL:AD
                elif "AD" in format_field:
                    fmt_parts = format_field.split(':')
                    try:
                        ad_idx = fmt_parts.index("AD")
                        sample_parts = sample_field.split(':')
                        ad_str = sample_parts[ad_idx] # e.g. "10,5" (Ref, Alt)
                        counts = [int(x) for x in ad_str.split(',')]
                        total_depth = sum(counts)
                        if total_depth > 0:
                            af = counts[1] / total_depth
                    except:
                        pass
                
                # Check for Mixed Status (20% - 80%)
                if af is not None:
                    if 0.20 <= af <= 0.80:
                        mixed_sites.append(f"Pos {pos}: AF={af:.2f}")
                        
    except Exception as e:
        print(f"  Warning: VCF parsing failed: {e}")
        return None

    if mixed_sites:
        print(f"  🚨 HETEROGENEITY DETECTED in wbeT: {len(mixed_sites)} sites mixed.")
        print(f"     Example: {mixed_sites[0]}")
        return {
            "status": "MIXED",
            "sites": mixed_sites,
            "count": len(mixed_sites)
        }
    else:
        print("  ✓ VCF confirms clonal purity (No heterozygous wbeT sites).")
        return {"status": "PURE"}

def track_wbeT(consensus_fasta, reference_fasta, vcf_path=None):
    print("🧬 wbeT MUTATION TRACKER (Ogawa-Inaba Surveillance)")
    print("=" * 70)
    
    # Check if wbeT reference exists
    wbeT_ref = "data/references/wbeT_ref.fasta"
    if not os.path.exists(wbeT_ref):
        print(f"  Warning: wbeT reference not found at {wbeT_ref}")
        return {"status": "UNKNOWN", "error": "Reference missing"}

    print(f"Analyzing sample: {consensus_fasta}")
    
    # Use blastn to find wbeT in the sample and get the subject sequence
    # outfmt 15 is JSON, but outfmt 6 with 'sseq' is easy to parse
    cmd = [
        "blastn",
        "-query", wbeT_ref,
        "-subject", consensus_fasta,
        "-outfmt", "6 qseqid sseqid pident length qlen slen qstart qend sstart send sseq",
        "-perc_identity", "80"
    ]
    
    serotype_result = {
        "serotype_status": "Unknown", 
        "wbeT_status": "UNKNOWN",
        "details": {}
    }
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        hits = result.stdout.strip().split("\n")
        
        if not hits or not hits[0]:
            print("  ⚠️  wbeT gene NOT FOUND in sample. Possible major deletion or non-O1 lineage.")
            return {"serotype_status": "Non-O1/Deleted", "wbeT_status": "ABSENT"}

        # Analyze the top hit
        parts = hits[0].split("\t")
        pident = float(parts[2])
        length = int(parts[3])
        qlen = int(parts[4])
        slen = int(parts[5])
        qstart = int(parts[6])
        sseq = parts[10].replace("-", "") # Remove gaps for translation

        coverage = length / qlen
        print(f"  wbeT detected: {pident}% identity, {coverage:.1%} coverage (qstart: {qstart})")

        # 🚀 ADVANCED DETECTION: Translation check
        # wbeT in reformulated reference starts at index 0 (1st base)
        CDS_START_IN_QUERY = 1
        
        # Calculate frame offset to stay in sync with CDS
        if qstart <= CDS_START_IN_QUERY:
            frame_offset = CDS_START_IN_QUERY - qstart
            # sseq_to_translate starts at the ATG
            sseq_to_translate = sseq[frame_offset:]
            aa_start_offset = 0
        else:
            # Hit starts after the ATG
            # Find the offset to the next codon start
            frame_offset = (3 - (qstart - CDS_START_IN_QUERY) % 3) % 3
            sseq_to_translate = sseq[frame_offset:]
            # We are starting 'N' amino acids into the protein
            aa_start_offset = (qstart - CDS_START_IN_QUERY + frame_offset) // 3

        from Bio.Seq import Seq
        seq_obj = Seq(sseq_to_translate)
        protein = seq_obj.translate(to_stop=False)
        
        # Check for premature stops (truncation)
        # Note: A real O1 wbeT has a stop at the end (840 bp -> 280 aa + stop)
        stop_pos = protein.find("*")
        is_truncated = False
        # If stop is found before the expected end (approx 280 aa)
        if stop_pos != -1 and (aa_start_offset + stop_pos) < 275:
            is_truncated = True
            print(f"  🚨 DISRUPTION DETECTED: Premature stop codon at amino acid {aa_start_offset + stop_pos + 1}!")
            
        # Check for length variation (Indels)
        # HC795 has a frameshift Indel that causes Inaba serotype
        # CRITICAL REFINEMENT: Only trigger truncation-Inaba if the subject is a full genome/assembly
        # (slen > 100kb). Subsampled reads (slen small) will naturally have low coverage.
        if length < (qlen * 0.95) and slen > 100000:
             print(f"  ⚠️  Significant sequence truncation/deletion detected in assembly/consensus.")
             is_truncated = True 
        elif length < (qlen * 0.95):
             print(f"  ℹ️  Low coverage in subsampled reads. Truncation check inconclusive.")

        # Logic for Inaba
        # 1. Truncation (anywhere significant)
        # 2. S158P mutation (position 158)
        # 3. Large deletion
        
        inaba_detected = is_truncated
        
        # S158P check (if sequence covers it)
        # Pos 158 (1-indexed) corresponds to aa_index 157
        target_aa_index = 157 - aa_start_offset
        if 0 <= target_aa_index < len(protein):
            actual_aa = protein[target_aa_index]
            if actual_aa == "P":
                print("  🚨 DISRUPTION DETECTED: S158P (Ser158Pro) mutation found!")
                inaba_detected = True
            elif actual_aa == "S":
                print("  ✓ Functional Serine found at position 158 (Ogawa-like)")

        if inaba_detected:
            print("  🎯 SEROTYPE PREDICTION: O1 Inaba (wbeT Inactive)")
            serotype_result = {
                "serotype_status": "O1 Inaba",
                "wbeT_status": "INACTIVE/DISRUPTED",
                "identity": pident,
                "coverage": coverage,
                "protein_length": len(protein) + aa_start_offset,
                "premature_stop": is_truncated,
                "s158p": inaba_detected and not is_truncated, # Simplified
                "alert_triggered": True
            }
        else:
            print("  ✓ wbeT appears functional. Serotype: O1 Ogawa")
            serotype_result = {
                "serotype_status": "O1 Ogawa",
                "wbeT_status": "FUNCTIONAL",
                "identity": pident,
                "coverage": coverage,
                "alert_triggered": False
            }

        # 🧬 NEW: VCF Cross-Check for Hikojima (Mixed)
        if vcf_path:
             vcf_check = check_vcf_heterozygosity(vcf_path)
             if vcf_check and vcf_check["status"] == "MIXED":
                  print(f"  🚨 OVERRIDE: Hikojima / Mixed Serotype Detected via VCF!")
                  serotype_result["serotype_status"] += " / Hikojima (Mixed)"
                  serotype_result["mixed_evidence"] = vcf_check["sites"]
                  serotype_result["alert_triggered"] = True

        return serotype_result

    except Exception as e:
        print(f"  Error tracking wbeT: {e}")
        return {"status": "ERROR", "error": str(e)}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("consensus")
    parser.add_argument("--ref", default="data/references/2010EL-1786.fasta")
    parser.add_argument("--vcf", help="Optional VCF for heterozygosity check")
    args = parser.parse_args()
    
    result = track_wbeT(args.consensus, args.ref, args.vcf)
    print("\nJSON Summary:")
    print(json.dumps(result, indent=2))
