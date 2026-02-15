#!/usr/bin/env python3
"""
Resolve haplotypes for heterogeneous surveillance loci.
Extracts reads for mixed genes and attempts to reconstruct Major/Minor alleles.
"""
import pysam
import json
import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter

def get_consensus_from_reads(reads, ref_start, ref_end):
    """Simple consensus from a list of AlignedSegment objects"""
    if not reads: return "N" * (ref_end - ref_start)
    
    # Initialize counts
    counts = {pos: Counter() for pos in range(ref_start, ref_end)}
    
    for read in reads:
        pairs = read.get_aligned_pairs(matches_only=True)
        read_seq = read.query_sequence
        for q_pos, r_pos in pairs:
            if ref_start <= r_pos < ref_end:
                base = read_seq[q_pos]
                counts[r_pos][base] += 1
                
    consensus = []
    for pos in range(ref_start, ref_end):
        if not counts[pos]:
            consensus.append("N")
        else:
            # Pick most common base
            consensus.append(counts[pos].most_common(1)[0][0])
            
    return "".join(consensus)

def resolve_haplotypes(bam_file, snp_report_file, output_file):
    try:
        with open(snp_report_file) as f:
            report = json.load(f)
    except Exception:
        # If report doesn't exist or is empty
        with open(output_file, 'w') as f: json.dump({}, f)
        return

    alerts = report.get('heterogeneity_alerts', [])
    if not alerts:
        with open(output_file, 'w') as f: json.dump({"status": "NO_ALERTS"}, f)
        return

    bam = pysam.AlignmentFile(bam_file, "rb")
    results = {}

    print(f"Resolving haplotypes for {len(alerts)} heterogeneous loci...")

    for alert in alerts:
        gene = alert.get('gene')
        chrom = alert.get('chrom')
        snp_pos = alert.get('pos')
        
        if not chrom or not gene or snp_pos is None: continue
        
        # Define window: +/- 50bp around SNP for context
        start = max(0, snp_pos - 50)
        end = snp_pos + 51
        
        # Group reads by base at SNP position
        groups = defaultdict(list)
        
        try:
            for read in bam.fetch(chrom, snp_pos, snp_pos+1):
                 pairs = read.get_aligned_pairs(matches_only=True)
                 for q_pos, r_pos in pairs:
                     if r_pos == snp_pos:
                         base = read.query_sequence[q_pos]
                         groups[base].append(read)
                         break
        except Exception as e:
            print(f"Error fetching reads for {gene}: {e}")
            continue
        
        # Sort groups by size
        sorted_groups = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)
        
        if len(sorted_groups) < 2:
            results[gene] = {"status": "SINGLE_ALLELE_FOUND_IN_WINDOW"}
            continue
            
        major_base, major_reads = sorted_groups[0]
        minor_base, minor_reads = sorted_groups[1]
        
        # Generate consensus for the window for each group
        seq_major = get_consensus_from_reads(major_reads, start, end)
        seq_minor = get_consensus_from_reads(minor_reads, start, end)
        
        results[gene] = {
            "chrom": chrom,
            "snp_pos": snp_pos,
            "window": f"{start}-{end}",
            "haplotypes": {
                "major": {
                    "base": major_base,
                    "count": len(major_reads),
                    "sequence": seq_major
                },
                "minor": {
                    "base": minor_base,
                    "count": len(minor_reads),
                    "sequence": seq_minor
                }
            }
        }
        
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
        
    print(f"✅ Haplotype resolution complete. Saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--snp-report", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    
    resolve_haplotypes(args.bam, args.snp_report, args.output)
