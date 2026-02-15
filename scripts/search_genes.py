import sys
import re

def search_gene(fasta_path, start_seq):
    with open(fasta_path, 'r') as f:
        content = f.read()
    
    # Remove headers and newlines
    chroms = re.split(r'>', content)[1:]
    for chrom_data in chroms:
        lines = chrom_data.split('\n')
        header = lines[0]
        seq = ''.join(lines[1:])
        
        matches = [m.start() for m in re.finditer(start_seq, seq)]
        if matches:
            print(f"Header: {header}")
            for m in matches:
                print(f"  Match at {m+1}")

if __name__ == "__main__":
    fasta = "data/references/2010EL-1786.fasta"
    gene_starts = {
        "ctxB": "ATGATTAAATTAAA",
        "tcpA": "ATGTTCAAAATTCAAA",
        "hapR": "ATGAAAAAATTAGGAGTAAC",
        "gyrA": "ATGAGCGTGGTTGTAA",
        "parC": "ATGAGCGAAATTATTT",
        "wbeT": "ATGAAAATATTGTTTCCAG",
        "katB": "ATGTCAAACCAAAACGAT",
        "ahpC": "ATGTCATTAATTAAACC",
        "scrA": "ATGTCAAACTTTATTAGCATT",
        "lip": "ATGAAGAAATTAATCATTTTG",
        "rbmA": "ATGAAAAAACTATTACTTGCCG"
    }
    
    for name, start in gene_starts.items():
        print(f"Searching for {name} ({start})...")
        search_gene(fasta, start)
