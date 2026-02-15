from Bio import SeqIO
import sys

def compare_fasta(file1, file2):
    try:
        dict1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta"))
        dict2 = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))
    except Exception as e:
        print(f"Error reading FASTA files: {e}")
        return

    common_loci = set(dict1.keys()) & set(dict2.keys())
    
    print(f"{'Locus':<15} | {'Length (Ref1)':<15} | {'Length (Ref2)':<15} | {'Status'}")
    print("-" * 60)
    
    identical_count = 0
    different_count = 0
    
    for locus in sorted(common_loci):
        seq1 = str(dict1[locus].seq).upper()
        seq2 = str(dict2[locus].seq).upper()
        
        status = "IDENTICAL" if seq1 == seq2 else "**DIFFERENT**"
        if seq1 == seq2:
            identical_count += 1
        else:
            different_count += 1
            
        print(f"{locus:<15} | {len(seq1):<15} | {len(seq2):<15} | {status}")

    print("-" * 60)
    print(f"Summary: {identical_count} Identical, {different_count} Different")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 compare_loci_sequences.py <fasta1> <fasta2>")
        sys.exit(1)
    compare_fasta(sys.argv[1], sys.argv[2])
