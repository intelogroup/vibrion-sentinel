#!/usr/bin/env python3
"""
Generate consensus genome from aligned BAM with quality-aware base calling.
Produces Strict (Majority) and IUPAC (Ambiguous) consensus FASTA files.
"""

import pysam
import numpy as np
import json
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO

IUPAC_CODES = {
    frozenset(['A']): 'A', frozenset(['C']): 'C', frozenset(['G']): 'G', frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R', frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S', frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K', frozenset(['A', 'C']): 'M',
    frozenset(['C', 'G', 'T']): 'B', frozenset(['A', 'G', 'T']): 'D',
    frozenset(['A', 'C', 'T']): 'H', frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

def get_iupac(base_counts, depth, threshold=0.2):
    """Return IUPAC code for bases above threshold frequency."""
    alleles = []
    for base, count in base_counts.items():
        if count / depth >= threshold:
            alleles.append(base)
    if not alleles: return 'N'
    return IUPAC_CODES.get(frozenset(alleles), 'N')

def generate_consensus_genome(bam_file, reference_fasta, output_dir, min_depth=1, min_quality=20, sample_id="sample"):
    """
    Generate consensus genome with quality awareness across all reference contigs.
    """
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"📖 Generating consensus genome from: {bam_file}")
    print(f"📋 Reference: {reference_fasta}")
    print(f"🎯 Parameters: min_depth={min_depth}X, min_quality=Q{min_quality}")
    
    # Load reference genome
    ref_records = list(SeqIO.parse(reference_fasta, "fasta"))
    ref_seqs = {record.id: str(record.seq).upper() for record in ref_records}
    total_ref_length = sum(len(s) for s in ref_seqs.values())
    
    consensus_dict = {}
    consensus_iupac_dict = {}
    depth_dict = {}
    heterogeneous_positions = {}
    coverage_regions = []
    
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    for record in ref_records:
        contig_id = record.id
        contig_len = len(record.seq)
        
        # Initialize
        contig_consensus = list(str(record.seq).upper())
        contig_consensus_iupac = list(str(record.seq).upper())
        contig_depth = np.zeros(contig_len, dtype=int)
        
        print(f"  Processing contig {contig_id} ({contig_len:,} bp)...")
        
        if contig_id not in bam.references:
            print(f"    ⚠️  Contig {contig_id} not found in BAM. Marking as N.")
            contig_consensus = ['N'] * contig_len
            contig_consensus_iupac = ['N'] * contig_len
            consensus_dict[contig_id] = contig_consensus
            consensus_iupac_dict[contig_id] = contig_consensus_iupac
            depth_dict[contig_id] = contig_depth
            continue

        for pileup_column in bam.pileup(contig_id, ignore_index=False, ignore_orphans=False, min_base_quality=0):
            pos = pileup_column.pos
            if pos >= contig_len: continue
            
            base_tally = defaultdict(int)
            valid_bases = 0
            
            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del or pileup_read.is_refskip: continue
                if pileup_read.query_position is None: continue
                    
                base = pileup_read.alignment.query_sequence[pileup_read.query_position].upper()
                qual = pileup_read.alignment.query_qualities[pileup_read.query_position]
                
                if qual >= min_quality:
                    base_tally[base] += 1
                    valid_bases += 1
            
            contig_depth[pos] = valid_bases
            
            if valid_bases > 0:
                # Majority Consensus
                sorted_tally = sorted(base_tally.items(), key=lambda x: x[1], reverse=True)
                most_common_base = sorted_tally[0][0]
                
                if valid_bases >= min_depth:
                    contig_consensus[pos] = most_common_base
                    # IUPAC Consensus (Ambiguous if mixed > 20%)
                    contig_consensus_iupac[pos] = get_iupac(base_tally, valid_bases, threshold=0.2)
                else:
                    contig_consensus[pos] = 'N'
                    contig_consensus_iupac[pos] = 'N'
                
                # Heterogeneity stats
                minor_alleles = [
                    {"base": b, "freq": round(count / valid_bases, 3)} 
                    for b, count in sorted_tally[1:] 
                    if (count / valid_bases) >= 0.05
                ]
                if minor_alleles:
                    heterogeneous_positions[f"{contig_id}:{pos}"] = {
                        "contig": contig_id, "pos": int(pos), "depth": int(valid_bases),
                        "dominant": most_common_base, "minor": minor_alleles
                    }
            else:
                contig_consensus[pos] = 'N'
                contig_consensus_iupac[pos] = 'N'
            
        consensus_dict[contig_id] = contig_consensus
        consensus_iupac_dict[contig_id] = contig_consensus_iupac
        depth_dict[contig_id] = contig_depth

        # Coverage regions logic (omitted for brevity, handled by depth stats)

    bam.close()
    
    # Stats aggregation
    all_depths = np.concatenate([d for d in depth_dict.values()]) if depth_dict else []
    stats = {
        'reference_length': total_ref_length,
        'coverage_percentage': float(np.sum(all_depths >= 1) / total_ref_length * 100) if total_ref_length > 0 else 0,
        'mean_depth': float(np.mean(all_depths)) if len(all_depths) > 0 else 0,
        'heterogeneous_sites_count': len(heterogeneous_positions),
        'confident_positions': int(np.sum(all_depths >= min_depth))
    }
    
    # Write Strict Consensus
    consensus_file = output_dir / f"{sample_id}_consensus.fasta"
    with open(consensus_file, 'w') as f:
        for contig_id, seq_list in consensus_dict.items():
            f.write(f">{sample_id}_{contig_id}_consensus\n")
            s = ''.join(seq_list)
            for i in range(0, len(s), 80): f.write(s[i:i+80] + '\n')

    # Write IUPAC Consensus
    iupac_file = output_dir / f"{sample_id}_consensus_iupac.fasta"
    with open(iupac_file, 'w') as f:
        for contig_id, seq_list in consensus_iupac_dict.items():
            f.write(f">{sample_id}_{contig_id}_consensus_iupac\n")
            s = ''.join(seq_list)
            for i in range(0, len(s), 80): f.write(s[i:i+80] + '\n')
            
    print(f"✅ Generated Strict Consensus: {consensus_file}")
    print(f"✅ Generated IUPAC Consensus: {iupac_file}")

    # Summary
    summary_file = output_dir / "consensus_generation_summary.json"
    summary = {
        'sample': sample_id,
        'statistics': stats,
        'output_files': {
            'consensus_fasta': str(consensus_file),
            'iupac_fasta': str(iupac_file)
        }
    }
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2)
        
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
