#!/usr/bin/env python3
"""
Verify Coverage Integrity & Library Strategy Gatekeeper
-----------------------------------------------------
1. Checks for "True Deletion" vs "Missing Data" using global/target depth ratio.
2. Detects LIBRARY STRATEGY (WGS vs RNA-Seq) using Coefficient of Variation (CV).

Forensic Logic:
- WGS: Flat coverage (Low CV < 1.0)
- RNA-Seq: Spiky coverage (High CV > 2.0)
- Mixed: High CV + High Base Coverage
"""

import sys
import argparse
import pysam
import json
import numpy as np
import os
import tempfile
import subprocess
from Bio import SeqIO

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

def calculate_coverage_metrics(bam_path, target_region):
    """
    Calculate coverage stats for Target Region and Whole Genome.
    Returns: Global Depth, Target Depth, CV (Coefficient of Variation)
    """
    samfile = pysam.AlignmentFile(bam_path, "rb")
    
    # 1. Target Coverage (e.g., wbeT)
    # Parse region string "Chrom:Start-End"
    try:
        chrom, coords = target_region.split(":")
        start, end = map(int, coords.split("-"))
        
        target_depths = []
        # Check if contig exists in BAM
        if chrom in samfile.references:
            for pileup in samfile.pileup(chrom, start, end, truncate=True):
                target_depths.append(pileup.nsegments)
        else:
            print(f"   ⚠️  Target chrom {chrom} not found in BAM. Reference mismatch?")
            return 0, 0, 0, 0

        # Handle zero coverage
        if not target_depths:
            avg_target_depth = 0
            min_target_depth = 0
        else:
            avg_target_depth = np.mean(target_depths)
            min_target_depth = np.min(target_depths)

        # 2. Global Coverage & Uniformity (CV)
        # Use the first contig for global estimate
        ref_name = samfile.references[0]
        chrom_len = samfile.get_reference_length(ref_name)
        sample_depths = []
        
        for pileup in samfile.pileup(ref_name, 0, chrom_len, step=100):
            sample_depths.append(pileup.nsegments)
            
        if not sample_depths:
            global_depth = 0
            cv = 0
        else:
            global_depth = np.mean(sample_depths)
            std_dev = np.std(sample_depths)
            cv = std_dev / global_depth if global_depth > 0 else 0

        return avg_target_depth, min_target_depth, global_depth, cv
    finally:
        samfile.close()

def main():
    parser = argparse.ArgumentParser(description="Verify Coverage Integrity & Library Strategy")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--target-region", required=True, help="Target region (e.g., CP003069.1:2678187-2678980)")
    parser.add_argument("--target-name", required=True, help="Name of target gene (e.g., wbeT)")
    parser.add_argument("--reference", help="Reference FASTA for BLAST discovery fallback")
    
    args = parser.parse_args()
    
    print(f"🔍 Analyzing Coverage Integrity: {args.bam}")
    
    target_region = args.target_region
    
    # Check if target chrom exists in BAM
    samfile = pysam.AlignmentFile(args.bam, "rb")
    chrom = target_region.split(":")[0]
    if chrom not in samfile.references and args.reference:
        print(f"   ⚠️  Default target {chrom} not in BAM. Attempting BLAST discovery...")
        new_region = discover_region_via_blast(args.target_name, args.reference)
        if new_region:
            print(f"   🎯 Discovered {args.target_name} at {new_region}")
            target_region = new_region
        else:
            print(f"   ❌ BLAST discovery failed for {args.target_name}")
    samfile.close()
            
    print(f"   Target: {args.target_name} ({target_region})")
    
    try:
        avg_target, min_target, global_depth, cv = calculate_coverage_metrics(args.bam, target_region)
        
        # --- LOGIC GATES ---
        
        # 1. Library Strategy Gatekeeper (CV Check)
        if cv < 1.2:
            lib_strategy = "WGS (Genomic DNA)"
            lib_status = "PASS"
            lib_warning = None
        elif cv < 2.5:
            lib_strategy = "AMBIGUOUS / MIXED"
            lib_status = "WARNING"
            lib_warning = "High coverage variance detected. Possible contaminated WGS or Mixed DNA/RNA."
        else:
            lib_strategy = "RNA-Seq (Transcriptomic)"
            lib_status = "FAIL"
            lib_warning = "Extreme coverage spikes detected. Incompatible with Variant Calling."

        # 2. Deletion Integrity (Ratio Check)
        if global_depth < 5:
            status = "LOW_COVERAGE (Global)"
            confidence = "LOW"
            verdict = "INCONCLUSIVE"
        else:
            ratio = avg_target / global_depth
            
            if ratio < 0.1:
                status = "DELETED (True Deletion)"
                confidence = "HIGH"
                verdict = "CONFIRMED_ABSENT"
            elif ratio < 0.5:
                 status = "PARTIAL / MOSAIC"
                 confidence = "MEDIUM"
                 verdict = "AMBIGUOUS"
            else:
                status = "PRESENT"
                confidence = "HIGH"
                verdict = "CONFIRMED_PRESENT"
        
        # Check valid coverage at target
        target_completeness = bool(min_target > 0)
        
        results = {
            "target": args.target_name,
            "region": target_region,
            "metrics": {
                "global_depth": round(global_depth, 2),
                "target_depth": round(avg_target, 2),
                "coverage_ratio": round(avg_target / global_depth if global_depth > 0 else 0, 2),
                "coefficient_of_variation": round(cv, 2)
            },
            "integrity": {
                "status": status,
                "confidence": confidence,
                "verdict": verdict,
                "target_completeness": target_completeness
            },
            "library_strategy": {
                "detected_type": lib_strategy,
                "status": lib_status,
                "warning": lib_warning
            }
        }
        
        print(f"   📊 Global Depth: {global_depth:.2f}x")
        print(f"   🎯 Target Depth: {avg_target:.2f}x")
        print(f"   📉 Coeff. Variation: {cv:.2f} ({lib_strategy})")
        print(f"   ✅ Verdict: {status} ({confidence})")

        with open(args.output, 'w') as f:
            json.dump(results, f, indent=2)
            
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"🚨 Error calculating coverage: {e}")
        # Fail gracefully JSON
        err = {"error": str(e), "status": "ERROR"}
        with open(args.output, 'w') as f:
            json.dump(err, f)
        sys.exit(1)

if __name__ == "__main__":
    main()