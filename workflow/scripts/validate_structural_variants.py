#!/usr/bin/env python3
import argparse
import json
import subprocess
import os
import tempfile
from Bio import SeqIO
import pysam

def discover_region_via_blast(query_name, target_fasta):
    """
    Search for the locus using BLAST to find coordinates in the current reference.
    """
    ref_loci_path = "data/references/reference_loci.fasta"
    if not os.path.exists(ref_loci_path):
        return None

    # 1. Extract the query sequence from reference_loci.fasta
    query_seq = None
    for record in SeqIO.parse(ref_loci_path, "fasta"):
        if record.id == query_name:
            query_seq = str(record.seq)
            break
    
    if not query_seq:
        return None

    # 2. BLAST against target
    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w") as tmp_query:
        tmp_query.write(f">{query_name}\n{query_seq}\n")
        tmp_query.flush()

        cmd = [
            "blastn",
            "-query", tmp_query.name,
            "-subject", str(target_fasta),
            "-outfmt", "6 sseqid sstart send pident length",
            "-perc_identity", "70",
            "-word_size", "7",
            "-max_target_seqs", "1"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = result.stdout.strip().split("\n")
            if not hits or not hits[0]:
                return None
            
            # 3. Extract the hit coordinates
            parts = hits[0].split("\t")
            sseqid, sstart, send = parts[0], int(parts[1]), int(parts[2])
            
            # Ensure sstart < send
            start, end = min(sstart, send), max(sstart, send)
            return f"{sseqid}:{start}-{end}"
        except Exception:
            return None

def get_region_coverage(bam_path, region):
    """Calculate average coverage for a region."""
    if not region:
        return 0.0
    cmd = ["samtools", "depth", "-a", "-r", region, bam_path]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        depths = [int(line.split('\t')[2]) for line in result.stdout.splitlines() if len(line.split('\t')) >= 3]
        return sum(depths) / len(depths) if depths else 0.0
    except subprocess.CalledProcessError:
        return 0.0

def check_wbeT_snp_manual(bam_path, chrom, pos):
    """Check for C->T at specific position using mpileup."""
    try:
        # pos is 1-based coordinate in the discovered chrom
        cmd = f"samtools mpileup -r {chrom}:{pos}-{pos} {bam_path} 2>/dev/null"
        output = subprocess.check_output(cmd, shell=True, text=True).strip()
        if not output:
            return None
        
        # Format: CHR POS REF DEPTH BASES QUAL
        parts = output.split('\t')
        if len(parts) >= 5:
            bases = parts[4].upper()
            depth = int(parts[3])
            # Count T's (Inaba) vs ,/. or C (Ogawa)
            count_t = bases.count('T')
            if depth > 0 and (count_t / depth) > 0.5:
                return f"Inaba (T at {pos}, quality {parts[5]})"
        return None
    except Exception as e:
        print(f"Error in manual SNP check: {e}")
        return None

def call_wbeT_variants(bam_path, reference, region):
    if not region:
        return []
    try:
        cmd = (
            f"bcftools mpileup -Ou -f {reference} -r {region} {bam_path} 2>/dev/null | "
            f"bcftools call -mv -Ou 2>/dev/null | "
            "bcftools query -f '%POS:%REF>%ALT\n' 2>/dev/null"
        )
        output = subprocess.check_output(cmd, shell=True, text=True).strip()
        return output.splitlines() if output else []
    except Exception:
        return []

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-bam", required=True, help="High-precision alignment (Illumina) for SNPs")
    parser.add_argument("--structural-bam", required=True, help="High-sensitivity alignment (Minimap2) for SVs")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--collection-year", type=int, default=2021)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Regions (Reference CP003069.1)
    DEFAULT_REGIONS = {
        "CTX_Phage": "CP003069.1:1041000-1047000", # Toxin
        "SXT_Element": "CP003069.1:98000-170000",   # MDR ICE
        "wbeT_Gene":  "CP003069.1:2678186-2678980"  # Critical for Inaba
    }
    
    # Gene names for BLAST discovery
    BLAST_QUERIES = {
        "CTX_Phage": "ctxB", # use ctxB as proxy for CTX phage
        "SXT_Element": "sxtMO", # proxy for SXT
        "wbeT_Gene": "wbeT"
    }

    report = {
        "sample": os.path.basename(args.primary_bam).replace("_aligned.sorted.bam", ""),
        "structural_variants": {},
        "inaba_status": "Unknown",
        "alerts": []
    }

    # Verify if default regions exist in BAM
    samfile = pysam.AlignmentFile(args.primary_bam, "rb")
    bam_refs = samfile.references
    samfile.close()
    
    active_regions = {}
    for key, region in DEFAULT_REGIONS.items():
        chrom = region.split(":")[0]
        if chrom in bam_refs:
            active_regions[key] = region
        else:
            print(f"   ⚠️  Default region {key} ({chrom}) not in BAM. Attempting BLAST discovery...")
            discovered = discover_region_via_blast(BLAST_QUERIES[key], args.reference)
            if discovered:
                print(f"   🎯 Discovered {key} at {discovered}")
                active_regions[key] = discovered
            else:
                print(f"   ❌ Could not locate {key} in reference.")
                active_regions[key] = None

    # 1. Check wbeT (Inaba Switch) - Primary BAM (BWA)
    wbeT_region = active_regions["wbeT_Gene"]
    wbeT_cov = get_region_coverage(args.primary_bam, wbeT_region) if wbeT_region else 0.0
    
    # RESCUE LOGIC TRIGGER
    used_rescue = False
    active_bam = args.primary_bam
    
    if wbeT_region and wbeT_cov < 5.0:
        # Check Secondary BAM (Minimap2)
        wbeT_cov_secondary = get_region_coverage(args.structural_bam, wbeT_region)
        if wbeT_cov_secondary > 5.0:
            print(f"Rescue Triggered for wbeT: Low primary coverage ({wbeT_cov}x), using structural alignment ({wbeT_cov_secondary}x)")
            active_bam = args.structural_bam
            wbeT_cov = wbeT_cov_secondary
            used_rescue = True
            report["alerts"].append("Targeted Rescue Alignment (Minimap2) used for wbeT analysis.")
    
    report["structural_variants"]["wbeT_coverage"] = round(wbeT_cov, 2)
    report["structural_variants"]["wbeT_rescue_used"] = used_rescue

    if wbeT_region and wbeT_cov > 5:
        # 1a. Direct SNP Check (Robust)
        # Offset calculation if it was CP003069.1:2678186-2678980 and SNP at 2678724
        # Offset is 2678724 - 2678186 = 538
        chrom_new, coords_new = wbeT_region.split(":")
        start_new = int(coords_new.split("-")[0])
        snp_pos = start_new + 538
        
        inaba_snp = check_wbeT_snp_manual(active_bam, chrom_new, snp_pos)
        
        # 1b. General Variant Calling
        variants = call_wbeT_variants(active_bam, args.reference, wbeT_region)
        
        if inaba_snp:
            report["inaba_status"] = "Inaba (Confirmed SNP)"
            report["structural_variants"]["inaba_snp_detail"] = inaba_snp
            report["alerts"].append(f"Inaba Switch Variant Detected (Q121*): {inaba_snp}")
        elif variants:
            report["inaba_status"] = "Inaba (Candidate Variants)"
            report["structural_variants"]["wbeT_variants"] = variants
            report["alerts"].append(f"wbeT Variants Found: {variants}")
        else:
            report["inaba_status"] = "Ogawa (Functional wbeT / Reversion Candidate)"
            report["alerts"].append("wbeT Gene Intact: Confirmed Ogawa Serotype.")
            report["structural_variants"]["ogawa_reversion"] = True
    else:
        report["inaba_status"] = "Inaba (Likely Deletion/Divergence)"
        if not wbeT_region:
             report["alerts"].append("wbeT Gene Not Found in genome.")
        else:
             report["alerts"].append("wbeT Gene Low Coverage even after rescue.")

    # 2. Check CTX (Toxin) - Use STRUCTURAL BAM (Sensitivity)
    ctx_region = active_regions["CTX_Phage"]
    ctx_cov = get_region_coverage(args.structural_bam, ctx_region) if ctx_region else 0.0
    report["structural_variants"]["CTX_coverage"] = round(ctx_cov, 2)
    
    if ctx_cov > 5:
        report["toxigenic"] = True
    else:
        report["toxigenic"] = False
        if not ctx_region:
             report["alerts"].append("CTX Phage Locus Not Found.")
        else:
             report["alerts"].append("CTX Phage Missing (Non-Toxigenic?)")

    # 2b. Check SXT (MDR Element) - Use STRUCTURAL BAM
    sxt_region = active_regions["SXT_Element"]
    sxt_cov = get_region_coverage(args.structural_bam, sxt_region) if sxt_region else 0.0
    report["structural_variants"]["SXT_coverage"] = round(sxt_cov, 2)

    # Write Output
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)

    print(f"Validation Report Saved: {args.output}")

if __name__ == "__main__":
    main()