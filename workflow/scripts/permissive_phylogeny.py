#!/usr/bin/env python3
"""
Vibrion Sentinel: Permissive Phylogeny (MIT/BSD Stack)
Legal Shield compliant tree building.
Stack: Biopython (BSD) + MAFFT (BSD).
Replaces: Augur (AGPL) + IQ-TREE (GPL).
"""

import argparse
import subprocess
import sys
from pathlib import Path
from Bio import SeqIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def run_mafft(input_fasta, output_fasta, threads=2):
    """
    Run MAFFT alignment via subprocess (Legal Gap).
    """
    cmd = ["mafft", "--retree", "2", "--maxiterate", "0", "--thread", str(threads), input_fasta]
    print(f"Running MAFFT: {' '.join(cmd)}")
    
    with open(output_fasta, "w") as out:
        try:
            subprocess.run(cmd, stdout=out, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"MAFFT Error: {e.stderr.decode()}")
            sys.exit(1)

def build_tree(aligned_fasta, output_newick, output_xml):
    """
    Build Neighbor-Joining tree using Biopython (Pure Python).
    """
    print("Reading alignment...")
    aln = SeqIO.parse(aligned_fasta, "fasta")
    # Convert iterator to MultipleSeqAlignment object? 
    # Bio.Phylo.TreeConstruction expects a MultipleSeqAlignment
    from Bio.Align import MultipleSeqAlignment
    msa = MultipleSeqAlignment(list(aln))
    
    print(f"Building Distance Matrix ({len(msa)} sequences)...")
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)
    
    print("Constructing NJ Tree...")
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # Save Outputs
    print(f"Saving Newick: {output_newick}")
    Phylo.write(tree, output_newick, "newick")
    
    print(f"Saving PhyloXML: {output_xml}")
    Phylo.write(tree, output_xml, "phyloxml")
    
    return tree

def draw_tree_png(tree, output_png):
    """
    Draw tree to PNG using Matplotlib (Agg backend).
    """
    try:
        import matplotlib
        matplotlib.use('Agg') # Headless mode
        import matplotlib.pyplot as plt
        
        print(f"Drawing tree to PNG: {output_png}")
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(1, 1, 1)
        
        # Draw tree
        Phylo.draw(tree, do_show=False, axes=ax)
        
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("✅ PNG saved.")
    except ImportError:
        print("⚠️ Matplotlib not found. Skipping PNG generation.")
    except Exception as e:
        print(f"❌ Failed to draw PNG: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", required=True, help="Sample consensus FASTA")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--outdir", required=True, help="Output directory for all phylo artifacts")
    parser.add_argument("--threads", default=4, help="Threads for MAFFT")
    args = parser.parse_args()
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    combined_fasta = outdir / "combined_genomes.fasta"
    aligned_fasta = outdir / "msa.fasta"
    tree_newick = outdir / "tree.nwk"
    tree_xml = outdir / "tree.xml"
    tree_png = outdir / "tree.png"
    
    # 1. Combine Sequences
    print("Step 1: Combining Genomes...")
    # Load Sample
    sample_seq = str(next(SeqIO.parse(args.consensus, "fasta")).seq)
    
    # Load References (Mocking the list for simplicity, in real implementation we scan a dir)
    # References should be at 'data/production_references' or 'data/references'
    # We will grab known references.
    
    ref_paths = [
        ("Haiti_2010", Path("data/references/2010EL-1786.fasta")),
        # If Bengal exists
        ("Bengal_1993", Path("data/production_references/Bengal_1993_Combined.fasta")),
         ("Environmental", Path("data/production_references/Environmental_NOVC_Contig1.fasta"))
    ]
    
    with open(combined_fasta, "w") as f:
        # Write Sample
        f.write(f">{args.sample_id}\n{sample_seq}\n")
        
        # Write Refs
        for name, path in ref_paths:
            if path.exists():
                # Take first sequence
                ref_rec = next(SeqIO.parse(path, "fasta"))
                f.write(f">{name}\n{ref_rec.seq}\n")
            else:
                print(f"Warning: Ref {name} missing at {path}")
                
    # 2. Align (MAFFT)
    print("Step 2: Aligning...")
    run_mafft(str(combined_fasta), str(aligned_fasta), args.threads)
    
    # 3. Build Tree (Bio.Phylo)
    print("Step 3: Building Tree...")
    tree = build_tree(str(aligned_fasta), str(tree_newick), str(tree_xml))
    
    # 4. Draw PNG
    draw_tree_png(tree, str(tree_png))
    
    print("✅ Permissive Phylogeny Complete.")

if __name__ == "__main__":
    main()
