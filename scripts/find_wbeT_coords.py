from Bio import SeqIO
from Bio.Seq import Seq

ref_path = "data/references/2010EL-1786.fasta"
# wbeT sequence snippet (N16961 VC0241)
target_seq = "ATGAGTTTATTTATTGC"

def find_coords():
    for record in SeqIO.parse(ref_path, "fasta"):
        idx = record.seq.find(target_seq)
        if idx != -1:
            print(f"FOUND in {record.id} at index {idx} (1-based: {idx+1})")
            return
        
        # Try reverse complement
        idx_rc = record.seq.find(Seq(target_seq).reverse_complement())
        if idx_rc != -1:
            print(f"FOUND (RC) in {record.id} at index {idx_rc}")
            return

if __name__ == "__main__":
    find_coords()
