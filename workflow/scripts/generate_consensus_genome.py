#!/usr/bin/env python3
"""
Generate consensus genome from aligned BAM with quality-aware base calling.
Supports multi-contig references (essential for fragmented regional archetypes).
"""

import pysam
import numpy as np
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

def generate_consensus_genome(bam_file, reference_fasta, output_dir, min_depth=1, min_quality=20, sample_id="sample"):
    """
    Generate consensus genome with quality awareness across all reference contigs.
    """
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"📖 Generating consensus genome from: {bam_file}")
    print(f"📋 Reference: {reference_fasta}")
    print(f"🎯 Parameters: min_depth={min_depth}X, min_quality=Q{min_quality}")
    print()
    
    # Load reference genome
    print("Loading reference genome...")
    ref_records = list(SeqIO.parse(reference_fasta, "fasta"))
    ref_seqs = {record.id: str(record.seq).upper() for record in ref_records}
    
    total_ref_length = sum(len(s) for s in ref_seqs.values())
    print(f"✓ Reference loaded: {total_ref_length:,} bp across {len(ref_records)} contigs")
    
    # Initialize consensus, depth and stats per contig
    consensus_dict = {}
    depth_dict = {}
    heterogeneous_positions = {}
    coverage_regions = []
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Process each contig
    total_processed_bases = 0
    
    for record in ref_records:
        contig_id = record.id
        contig_len = len(record.seq)
        
        # Initialize for this contig
        contig_consensus = list(str(record.seq).upper())
        contig_depth = np.zeros(contig_len, dtype=int)
        
        print(f"  Processing contig {contig_id} ({contig_len:,} bp)...")
        
        # Use pileup for this specific contig
        # Check if contig is in BAM
        if contig_id not in bam.references:
            print(f"    ⚠️  Contig {contig_id} not found in BAM. Marking as N.")
            contig_consensus = ['N'] * contig_len
            consensus_dict[contig_id] = contig_consensus
            depth_dict[contig_id] = contig_depth
            continue

        for pileup_column in bam.pileup(contig_id, ignore_index=False, ignore_orphans=False, min_base_quality=0):
            pos = pileup_column.pos
            if pos >= contig_len:
                continue
            
            base_tally = defaultdict(int)
            valid_bases = 0
            
            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip:
                    continue
                
                # Check for query_position (None if read has deletion at this pos)
                if pileup_read.query_position is None:
                    continue
                    
                base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                qual = pileup_read.alignment.query_qualities[pileup_read.query_position]
                
                if qual >= min_quality:
                    base_tally[base] += 1
                    valid_bases += 1
            
            contig_depth[pos] = valid_bases
            
            if valid_bases > 0:
                sorted_tally = sorted(base_tally.items(), key=lambda x: x[1], reverse=True)
                most_common_base = sorted_tally[0][0]
                
                # Heterogeneity tracking
                minor_alleles = [
                    {"base": b, "freq": round(count / valid_bases, 3)} 
                    for b, count in sorted_tally[1:] 
                    if (count / valid_bases) >= 0.02
                ]
                
                if minor_alleles:
                    heterogeneous_positions[f"{contig_id}:{pos}"] = {
                        "contig": contig_id,
                        "pos": int(pos),
                        "dominant": most_common_base,
                        "minor": minor_alleles,
                        "depth": int(valid_bases)
                    }
                
                if valid_bases >= min_depth:
                    contig_consensus[pos] = most_common_base
                else:
                    contig_consensus[pos] = 'N'
            else:
                contig_consensus[pos] = 'N'
                
            total_processed_bases += 1
            
        # Analyze coverage regions for this contig
        current_region_start = None
        for i in range(contig_len):
            if contig_depth[i] >= min_depth:
                if current_region_start is None:
                    current_region_start = i
            else:
                if current_region_start is not None:
                    coverage_regions.append({
                        'contig': contig_id,
                        'start': current_region_start,
                        'end': i - 1,
                        'length': i - current_region_start,
                        'avg_depth': np.mean(contig_depth[current_region_start:i])
                    })
                    current_region_start = None
        if current_region_start is not None:
            coverage_regions.append({
                'contig': contig_id,
                'start': current_region_start,
                'end': contig_len - 1,
                'length': contig_len - current_region_start,
                'avg_depth': np.mean(contig_depth[current_region_start:])
            })
            
        consensus_dict[contig_id] = contig_consensus
        depth_dict[contig_id] = contig_depth

    bam.close()
    
    # Aggregate Stats
    all_depths = np.concatenate([d for d in depth_dict.values()])
    covered_positions = np.sum(all_depths >= 1)
    confident_positions = np.sum(all_depths >= min_depth)
    
    stats = {
        'reference_length': total_ref_length,
        'total_aligned_positions': int(covered_positions),
        'confident_positions': int(confident_positions),
        'coverage_percentage': float(covered_positions / total_ref_length * 100) if total_ref_length > 0 else 0,
        'mean_depth': float(np.mean(all_depths)) if len(all_depths) > 0 else 0,
        'median_depth': float(np.median(all_depths)) if len(all_depths) > 0 else 0,
        'max_depth': int(np.max(all_depths)) if len(all_depths) > 0 else 0,
        'heterogeneous_sites_count': len(heterogeneous_positions)
    }
    
    # Write consensus FASTA
    consensus_file = output_dir / f"{sample_id}_consensus.fasta"
    with open(consensus_file, 'w') as f:
        for contig_id, seq_list in consensus_dict.items():
            f.write(f">{sample_id}_{contig_id}_consensus\n")
            consensus_seq = ''.join(seq_list)
            for i in range(0, len(consensus_seq), 80):
                f.write(consensus_seq[i:i+80] + '\n')
    
    # Write depth profile (condensed)
    summary_file = output_dir / "consensus_generation_summary.json"
    summary = {
        'sample': sample_id,
        'statistics': stats,
        'coverage_regions': coverage_regions[:100], # Cap for JSON size
        'output_files': {
            'consensus_fasta': str(consensus_file)
        }
    }
    
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
        
    print(f"\n✅ Consensus generation complete! {total_ref_length:,} bp processed.")
    return summary

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--min-depth", type=int, default=1)
    parser.add_argument("--min-qual", type=int, default=20)
    
    args = parser.parse_args()
    generate_consensus_genome(args.bam, args.ref, args.outdir, args.min_depth, args.min_qual, args.sample)