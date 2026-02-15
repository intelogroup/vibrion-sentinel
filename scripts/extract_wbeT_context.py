
from Bio import SeqIO

# Load wbeT
record = None
with open("data/references/reference_loci.fasta") as f:
    for r in SeqIO.parse(f, "fasta"):
        if r.id == "wbeT":
            record = r
            break

if record:
    seq = str(record.seq)
    # 1-based pos 793 means 0-based index 792
    # Let's extract 792 +/- 15bp
    snp_pos = 792
    context = seq[snp_pos-15 : snp_pos+16]
    print(f"Full Length: {len(seq)}")
    print(f"Context (Ogawa G at {snp_pos+1}): {context}")
    print(f"Base at 793: {seq[snp_pos]}")
else:
    print("wbeT not found")
