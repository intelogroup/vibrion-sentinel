
import random
import gzip
from pathlib import Path

def read_fasta(fasta_path):
    seqs = {}
    name = None
    seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name:
            seqs[name] = "".join(seq)
    return seqs

def reverse_complement(seq):
    compl = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(compl.get(b, 'N') for b in reversed(seq))

def simulate_reads(genome_path, output_path, num_reads=250000, read_len=150, error_rate=0.01):
    genome = read_fasta(genome_path)
    contigs = list(genome.keys())
    
    print(f"Simulating {num_reads} reads from {len(contigs)} contigs...")
    
    with gzip.open(output_path, 'wt') as f_out:
        for i in range(num_reads):
            # Pick random contig weighted by length? No, just random.
            # Ideally weighted, but simple random is fine for coverage spread
            contig = random.choice(contigs)
            seq = genome[contig]
            
            if len(seq) < read_len + 50:
                continue
                
            start = random.randint(0, len(seq) - read_len)
            fragment = seq[start:start+read_len].upper()
            
            # Add errors
            read_seq = []
            qual = []
            for base in fragment:
                if random.random() < error_rate:
                    read_seq.append(random.choice(['A', 'C', 'G', 'T'])) 
                    qual.append('(') # Low quality
                else:
                    read_seq.append(base)
                    qual.append('I') # High quality (Phred 40)
            
            final_seq = "".join(read_seq)
            final_qual = "".join(qual)
            
            # Write R1 (Pairing is overkill for this logic check, SE is easier but pipeline handles PE)
            # Pipeline is configured for PE usually? "reads: ..._1.fastq.gz".
            # The pipeline detects PE if _1 and _2 exist. Or SE if just _1 is detected?
            # Let's generate SE file named _1.fastq.gz.
            
            f_out.write(f"@SIM_{i} 1/1\n{final_seq}\n+\n{final_qual}\n")

if __name__ == "__main__":
    ref_path = Path("data/references/Haiti_2022_Resurgence.fasta")
    out_path = Path("data/raw_reads/Haiti_Resurgence_2022.fastq.gz")
    
    if not ref_path.exists():
        print(f"Error: {ref_path} not found.")
        exit(1)
        
    simulate_reads(ref_path, out_path)
