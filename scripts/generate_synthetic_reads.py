#!/usr/bin/env python3
"""
Vibrion Sentinel: Synthetic Read Generator
Generates perfect coverage reads from reference genome for pipeline validation.
"""

import random
import gzip
from Bio import SeqIO
from pathlib import Path

def generate_synthetic_reads(ref_path, output_path, num_reads=50000, read_length=150):
    """
    Generate synthetic reads from reference genome.
    
    Args:
        ref_path: Path to reference FASTA
        output_path: Output FASTQ path
        num_reads: Number of reads to generate
        read_length: Length of each read (bp)
    """
    print(f"🧬 Generating {num_reads} synthetic reads from {ref_path}")
    
    # Load reference
    record = next(SeqIO.parse(ref_path, "fasta"))
    seq = str(record.seq).upper()
    genome_length = len(seq)
    
    print(f"   Reference: {record.id} ({genome_length:,} bp)")
    
    # Calculate expected coverage
    total_bases = num_reads * read_length
    coverage = total_bases / genome_length
    print(f"   Target coverage: {coverage:.1f}x")
    
    # Generate reads
    reads = []
    for i in range(num_reads):
        # Random start position
        start = random.randint(0, genome_length - read_length)
        
        # Extract read sequence
        read_seq = seq[start:start + read_length]
        
        # Skip if contains N's (gaps in reference)
        if 'N' in read_seq:
            continue
        
        # Perfect quality scores (Phred 40)
        qual = "I" * len(read_seq)
        
        # FASTQ format
        reads.append(f"@SYNTHETIC_READ_{i} start={start}\n{read_seq}\n+\n{qual}")
    
    # Write output (gzipped if requested)
    output_path = Path(output_path)
    
    if output_path.suffix == '.gz':
        with gzip.open(output_path, 'wt') as f:
            f.write("\n".join(reads) + "\n")
    else:
        with open(output_path, 'w') as f:
            f.write("\n".join(reads) + "\n")
    
    print(f"✅ Created {len(reads)} reads")
    print(f"   Output: {output_path}")
    print(f"   Expected coverage: {len(reads) * read_length / genome_length:.1f}x")
    return len(reads)

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Generate synthetic reads for pipeline validation")
    parser.add_argument("--reference", required=True, help="Reference genome FASTA")
    parser.add_argument("--output", required=True, help="Output FASTQ file (.fastq.gz)")
    parser.add_argument("--reads", type=int, default=50000, help="Number of reads (default: 50000)")
    parser.add_argument("--length", type=int, default=150, help="Read length (default: 150bp)")
    
    args = parser.parse_args()
    
    generate_synthetic_reads(
        ref_path=args.reference,
        output_path=args.output,
        num_reads=args.reads,
        read_length=args.length
    )
