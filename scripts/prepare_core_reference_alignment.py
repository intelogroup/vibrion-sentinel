import os
import sys
from pathlib import Path
from Bio import SeqIO

# Add workflow scripts to path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))
from utils_bio import extract_loci_from_bed, run_mafft_alignment

def main():
    bed_path = "data/references/surveillance_loci.bed"
    ref_dir = Path("data/production_references")
    out_dir = Path("data/core_alignment")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Define Archetypes
    archetypes = [
        ("Haiti_2010", Path("data/references/2010EL-1786.fasta")),
        ("Bengal_1993", ref_dir / "Bengal_1993_Combined.fasta"),
        ("Classical_569B", ref_dir / "Classical_569B_Combined.fasta"),
        ("Environmental_NOVC", ref_dir / "Environmental_NOVC_Contig1.fasta")
    ]
    
    records = []
    for name, path in archetypes:
        rec = extract_loci_from_bed(str(path), bed_path, name)
        if rec:
            records.append(rec)
            
    if not records:
        print("Error: No records extracted.")
        return
        
    # 2. Save Unaligned
    unaligned_path = out_dir / "unaligned_core_references.fasta"
    SeqIO.write(records, unaligned_path, "fasta")
    print(f"Saved unaligned core sequences to {unaligned_path}")
    
    # 3. Align (MAFFT)
    aligned_path = out_dir / "reference_core_alignment.fasta"
    
    print(f"Aligning with MAFFT...")
    if run_mafft_alignment(str(unaligned_path), str(aligned_path), threads=4):
        print(f"✅ Pre-computed core alignment saved to {aligned_path}")
    else:
        print("Failed to align")

if __name__ == "__main__":
    main()
