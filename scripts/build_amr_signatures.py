#!/usr/bin/env python3
"""
Build individual Sourmash signatures for each AMR gene in CARD database.
This enables containment-based detection of specific genes.
"""

import subprocess
from pathlib import Path
from Bio import SeqIO

def main():
    # Paths
    card_fasta = Path("data/amr_signatures/nucleotide_fasta_protein_homolog_model.fasta")
    output_dir = Path("data/amr_signatures/individual_genes")
    output_dir.mkdir(exist_ok=True)
    
    print(f"📦 Building individual AMR gene signatures from {card_fasta}")
    print(f"📁 Output directory: {output_dir}")
    
    # Parse FASTA and create individual signatures
    count = 0
    for record in SeqIO.parse(card_fasta, "fasta"):
        # Extract gene name from header
        # Example: >gb|AY234334.1|+|0-846|ARO:3000600|Erm(34) [Alkalihalobacillus clausii]
        header = record.description
        
        # Try to extract gene name
        gene_name = None
        if "|" in header:
            parts = header.split("|")
            if len(parts) >= 5:
                # Get the part after ARO number
                gene_part = parts[5] if len(parts) > 5 else parts[4]
                # Clean up
                gene_name = gene_part.split("[")[0].strip()
        
        if not gene_name:
            gene_name = record.id
        
        # Sanitize filename
        safe_name = gene_name.replace("/", "_").replace("(", "").replace(")", "").replace(" ", "_")
        
        # Write individual FASTA
        temp_fasta = output_dir / f"{safe_name}.fasta"
        with open(temp_fasta, "w") as f:
            f.write(f">{gene_name}\n{record.seq}\n")
        
        # Build signature
        sig_file = output_dir / f"{safe_name}.sig"
        cmd = [
            "sourmash", "sketch", "dna",
            "-p", "k=31,scaled=1000",
            str(temp_fasta),
            "-o", str(sig_file),
            "--name", gene_name
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            count += 1
            if count % 500 == 0:
                print(f"  ✓ Processed {count} genes...")
        
        # Clean up temp FASTA
        temp_fasta.unlink()
    
    print(f"\n✅ Built {count} individual AMR gene signatures")
    print(f"📊 Combining into searchable database...")
    
    # Combine all signatures into a single .zip for faster searching
    all_sigs = list(output_dir.glob("*.sig"))
    if all_sigs:
        cmd = [
            "sourmash", "sig", "cat",
            *[str(s) for s in all_sigs],
            "-o", "data/amr_signatures/card_amr_genes.zip"
        ]
        subprocess.run(cmd, check=True)
        print(f"✅ Created searchable database: card_amr_genes.zip")

if __name__ == "__main__":
    main()
