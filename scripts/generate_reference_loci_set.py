
import subprocess
from pathlib import Path
from Bio import SeqIO
import sys

def generate_loci_set(genome_path, reference_loci_path, output_path):
    print(f"Generating loci set from {genome_path}...")
    
    found_loci = {}
    
    # 1. Load reference loci query sequences
    queries = {}
    for record in SeqIO.parse(reference_loci_path, "fasta"):
        queries[record.id] = str(record.seq)
        
    print(f"Looking for {len(queries)} loci...")

    # 2. BLAST each locus against the new genome
    # We write a temp query file
    for name, seq in queries.items():
        with open("temp_query.fasta", "w") as f:
            f.write(f">{name}\n{seq}\n")
            
        cmd = [
            "blastn",
            "-query", "temp_query.fasta",
            "-subject", str(genome_path),
            "-outfmt", "6 sseqid sstart send pident length sseq",
            "-perc_identity", "70", # Harvest even if drifted
            "-max_target_seqs", "1"
        ]
        
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.stdout.strip():
                parts = res.stdout.strip().split("\t")
                # sseq is the sequence of the subject (the new genome) aligned to query
                # Actually, -outfmt 6 sseq gives the aligned part. 
                # Better to use sstart/send to extract full context if needed, but sseq is fine for now
                # provided it covers the whole gene. 
                # Let's check length coverage?
                
                # To be safer, we use sseq which is the hit sequence.
                # Note: BLAST sseq might exclude unaligned ends if coverage < 100%.
                # For a reference, we ideally want the full gene.
                # But for now, capturing the hit is sufficient.
                
                # Correction: sseq includes dashes for gaps. We should remove them.
                sequence = parts[5].replace("-", "")
                found_loci[name] = sequence
                print(f"  ✅ Found {name} ({len(sequence)} bp, {parts[3]}% ID)")
            else:
                print(f"  ❌ {name} NOT FOUND")
        except Exception as e:
            print(f"  Error searching for {name}: {e}")

    # 3. Write output
    with open(output_path, "w") as f:
        for name, seq in found_loci.items():
            f.write(f">{name}\n{seq}\n")
            
    print(f"Saved {len(found_loci)} loci to {output_path}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 generate_reference_loci_set.py <genome> <ref_loci> <output>")
        sys.exit(1)
        
    generate_loci_set(sys.argv[1], sys.argv[2], sys.argv[3])
