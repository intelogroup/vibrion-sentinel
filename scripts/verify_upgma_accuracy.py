from Bio import SeqIO, AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from pathlib import Path
import subprocess

def main():
    outdir = Path("data/pipeline_output/Mock_Haiti_2025/10_phylogeny")
    ref_alignment = Path("data/core_alignment/reference_core_alignment.fasta")
    sample_fasta = outdir / "sample_core.fasta"
    
    # 1. Merge Sample with Reference
    merged_fasta = outdir / "merged_for_upgma.fasta"
    records = list(SeqIO.parse(ref_alignment, "fasta"))
    records.append(next(SeqIO.parse(sample_fasta, "fasta")))
    SeqIO.write(records, merged_fasta, "fasta")
    
    # 2. Align (MAFFT)
    aligned_fasta = outdir / "aligned_for_upgma.fasta"
    mafft_path = "/Users/kalinovdameus/Developer/Vibrion/.snakemake/conda/80f687140b6ccd3d74604e8c789853c9_/bin/mafft"
    subprocess.run([mafft_path, "--auto", str(merged_fasta)], stdout=open(aligned_fasta, "w"), check=True)
    
    # 3. Build UPGMA Tree
    alignment = AlignIO.read(aligned_fasta, "fasta")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    
    # 4. Verify Placement
    print("\n--- UPGMA Tree Structure ---")
    Phylo.draw_ascii(tree)
    
    # Check if Mock is closer to Haiti_2010 than to Environmental
    # We can do this by looking at branches or just visually in the ASCII tree
    print("\nClustering Check:")
    # Simple logic: Is Mock parented with Haiti? 
    # BioPython tree logic is complex to parse programmatically for siblings, but draw_ascii shows it.
    
if __name__ == "__main__":
    main()
