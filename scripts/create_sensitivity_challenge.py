#!/usr/bin/env python3
"""
Create Sensitivity Challenge Sample (Needle in a Haystack)
Mixes simulated reads from a target genome into a background FASTQ file.
"""

import argparse
import gzip
import random
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

def simulate_reads(fasta_path, num_reads, read_len=150, error_rate=0.001):
    """Simulate single-end Illumina reads from a FASTA reference."""
    print(f"   🧬 Simulating {num_reads} reads from {fasta_path}...")
    
    # Load genome
    genome_seqs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome_seqs.append(str(record.seq))
    
    if not genome_seqs:
        raise ValueError("No sequences found in reference FASTA")
    
    genome = "".join(genome_seqs)
    genome_len = len(genome)
    
    reads = []
    qual_char = chr(33 + 30) # High quality score (approx Q30)
    
    for i in range(num_reads):
        # Pick random start position
        start = random.randint(0, genome_len - read_len - 1)
        seq = list(genome[start : start + read_len].upper())
        
        # Introduce errors
        for j in range(len(seq)):
            if random.random() < error_rate:
                seq[j] = random.choice(['A', 'C', 'G', 'T'])
        
        # Reverse complement? (50% chance)
        seq_str = "".join(seq)
        if random.random() < 0.5:
            seq_str = str(Seq(seq_str).reverse_complement())
            
        header = f"@SIMULATED_READ_{i} length={read_len}"
        reads.append(f"{header}\n{seq_str}\n+\n{qual_char * read_len}\n")
        
    return reads

def load_fastq_reads(fastq_path, max_reads=None):
    """Load reads from a gzipped FASTQ file."""
    print(f"   📂 Loading background reads from {fastq_path}...")
    reads = []
    with gzip.open(fastq_path, "rt") as f:
        while True:
            header = f.readline()
            if not header: break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            reads.append(f"{header}{seq}{plus}{qual}")
            if max_reads and len(reads) >= max_reads:
                break
    return reads

def main():
    parser = argparse.ArgumentParser(description="Create Sensitivity Challenge Sample")
    parser.add_argument("--target", required=True, help="Target reference FASTA (Vibrio)")
    parser.add_argument("--background", required=True, help="Background FASTQ (E. coli)")
    parser.add_argument("--purity", type=float, default=0.01, help="Target purity (0.0 - 1.0)")
    parser.add_argument("--total-reads", type=int, default=10000, help="Total output reads")
    parser.add_argument("--output", required=True, help="Output FASTQ.gz path")
    args = parser.parse_args()

    print(f"🧪 creating Sensitivity Challenge Sample")
    print(f"   Target Purity: {args.purity:.1%}")
    print(f"   Total Reads: {args.total_reads}")

    # Calculate read counts
    n_target = int(args.total_reads * args.purity)
    n_background = args.total_reads - n_target
    
    print(f"   → Target Reads: {n_target}")
    print(f"   → Background Reads: {n_background}")

    # 1. Simulate Target Reads
    target_reads = simulate_reads(args.target, n_target)

    # 2. Load Background Reads
    bg_reads = load_fastq_reads(args.background, max_reads=n_background * 2) # Load extra to sample from
    if len(bg_reads) < n_background:
        print(f"   ⚠️  Warning: Not enough background reads! Repeating to fill.")
        while len(bg_reads) < n_background:
            bg_reads += bg_reads
    
    sampled_bg_reads = random.sample(bg_reads, n_background)

    # 3. Mix and Shuffle
    mixed_reads = target_reads + sampled_bg_reads
    random.shuffle(mixed_reads)

    # 4. Write Output
    print(f"   💾 Writing to {args.output}...")
    with gzip.open(args.output, "wt") as f:
        for read in mixed_reads:
            f.write(read)

    print("✅ Done.")

if __name__ == "__main__":
    main()
