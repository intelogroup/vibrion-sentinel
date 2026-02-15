#!/usr/bin/env python3
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import gzip
import os

def generate_reads(fasta_path, coverage, read_len=150, error_rate=0.01):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    total_len = sum(len(r.seq) for r in records)
    num_reads = int((total_len * coverage) / read_len)
    
    reads = []
    bases = ['A', 'C', 'G', 'T']
    
    print(f"Generating {num_reads} reads from {fasta_path}...")
    
    for _ in range(num_reads):
        # Pick random record
        rec = random.choice(records)
        if len(rec.seq) <= read_len:
            start = 0
            seq = str(rec.seq)
        else:
            start = random.randint(0, len(rec.seq) - read_len)
            seq = str(rec.seq[start:start+read_len])
        
        # Add random errors
        seq_list = list(seq)
        for i in range(len(seq_list)):
            if random.random() < error_rate:
                seq_list[i] = random.choice(bases)
        
        read_seq = "".join(seq_list)
        reads.append(read_seq)
        
    return reads

def main():
    o1_ref = "data/references/2010EL-1786.fasta"
    o139_cluster = "data/serogroup_reference/clusters/LC594838.1_O139.fasta"
    output_fastq = "data/raw_reads/Haiti_O1_O139_Mixed.fastq.gz"
    
    # 1. Generate O1 reads (50x)
    o1_reads = generate_reads(o1_ref, 25, read_len=150)
    
    # 2. Generate O139 cluster reads (high coverage to ensure detection)
    # Since it's just the cluster, we'll generate enough to match the O1 signal
    o139_reads = generate_reads(o139_cluster, 500, read_len=150)
    
    # 3. Mix and write to GZipped FASTQ
    all_reads = o1_reads + o139_reads
    random.shuffle(all_reads)
    
    print(f"Writing {len(all_reads)} reads to {output_fastq}...")
    with gzip.open(output_fastq, "wt") as f:
        for i, seq in enumerate(all_reads):
            f.write(f"@READ_{i}\n{seq}\n+\n{'I'*len(seq)}\n")
            
    print("✅ Co-infection dataset created.")

if __name__ == "__main__":
    main()
