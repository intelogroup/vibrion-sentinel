#!/usr/bin/env python3
"""
Vibrion Sentinel: CholeraSeq Modular Tree Builder
Professional-grade phylogenetic inference using pre-computed core alignments.
Implementation: IQ-TREE2 (+ASC) + TreeTime (Temporal Scaling).
"""

import argparse
import subprocess
import json
import sys
import shutil
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import shared utilities
sys.path.insert(0, str(Path(__file__).parent))
from utils_bio import (
    extract_loci_from_bed,
    check_sequence_quality,
    get_mafft_path,
    filter_empty_sequences,
    build_upgma_tree,
    render_tree_png
)

def run_cmd(cmd, check=True):
    """Run shell command with logging."""
    print(f"Executing: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, check=check, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        if check: sys.exit(1)
        return e

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--loci", help="Path to pre-extracted loci FASTA (optional)")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--threads", default=4, type=int)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    ref_alignment = Path("data/core_alignment/reference_core_alignment.fasta")
    bed_path = Path("data/references/surveillance_loci.bed")
    
    # Get MAFFT path from utility
    try:
        mafft_path = get_mafft_path()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    # Prefer iqtree2 if available, otherwise fall back to iqtree
    iqtree_path = shutil.which("iqtree2") or shutil.which("iqtree")
    if not iqtree_path:
        print("Error: IQ-TREE not found in PATH. Install iqtree/iqtree2 in the analysis environment.")
        sys.exit(1)

    def create_stub_outputs(message):
        """Generates placeholder files to prevent pipeline crash."""
        print(f"⚠️  {message}")
        print("Generating stub tree files to allow pipeline continuation...")
        
        final_tree_path = outdir / "tree.nwk"
        with open(final_tree_path, 'w') as f:
            f.write("();\n")
        
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            
            fig = plt.figure(figsize=(8, 6))
            plt.text(0.5, 0.5, f"Phylogeny Unavailable\n{message}", 
                    ha='center', va='center', fontsize=12)
            plt.axis('off')
            plt.savefig(outdir / "tree.png", dpi=150, bbox_inches='tight')
            plt.close(fig)
        except Exception as e:
            print(f"Warning: Could not create stub PNG: {e}")
            # Ensure file exists even if plotting fails
            (outdir / "tree.png").touch()
            
        print("✅ CholeraSeq Tree Module Complete (Stub Mode).")
        sys.exit(0)

    # 1. Extract Sample Loci
    print("Step 1: Extracting Sample Core Loci...")
    sample_fasta = outdir / "sample_core.fasta"
    
    if args.loci and Path(args.loci).exists():
        print(f"   Using pre-extracted loci from {args.loci}")
        # Concatenate multi-record FASTA into single core sequence
        loci_records = list(SeqIO.parse(args.loci, "fasta"))
        if not loci_records:
            create_stub_outputs("No sequences found in loci file (Low Coverage)")
            
        combined_seq = "".join(str(rec.seq) for rec in loci_records)
        sample_rec = SeqRecord(Seq(combined_seq), id=args.sample_id, description=f"Core loci from {args.loci}")
        SeqIO.write(sample_rec, sample_fasta, "fasta")
    else:
        sample_rec = extract_loci_from_bed(args.consensus, str(bed_path), args.sample_id)
        if not sample_rec:
            create_stub_outputs("Failed to extract sample loci (Low Coverage)")
        SeqIO.write(sample_rec, sample_fasta, "fasta")
    
    # 1b. Check sequence quality
    is_valid, qc_message = check_sequence_quality(sample_rec, ref_alignment)
    if not is_valid:
        create_stub_outputs(f"Sample sequence quality insufficient: {qc_message}")

    # 2. Filter empty sequences from reference alignment (IQ-TREE fix)
    print("Step 2a: Filtering empty sequences from reference alignment...")
    filtered_ref = outdir / "reference_filtered.fasta"
    kept = filter_empty_sequences(str(ref_alignment), str(filtered_ref), min_bases=1)
    print(f"   Kept {kept} sequences with real bases")
    
    # 2b. Append to Reference Alignment (MAFFT --add)
    print("Step 2b: Appending Sample to Reference Matrix...")
    merged_alignment = outdir / "merged_core_alignment.fasta"
    with open(merged_alignment, "w") as out_f:
        run_mafft = [mafft_path, "--add", str(sample_fasta), "--reorder", "--thread", str(args.threads), str(filtered_ref)]
        subprocess.run(run_mafft, stdout=out_f, check=True)

    # 3. Infer ML Tree (IQ-TREE2)
    print("Step 3: Inferring ML Tree (IQ-TREE2)...")
    # Using GTR+G (ASC is only for SNP-only alignments, which this is not)
    iqtree_cmd = [
        iqtree_path,
        "-s", str(merged_alignment),
        "-m", "GTR+G",
        "-nt", str(args.threads),
        "-pre", str(outdir / "iqtree_output"),
        "-redo"
    ]
    run_cmd(iqtree_cmd)
    
    # 4. Temporal Scaling (TreeTime)
    print("Step 4: Temporal Scaling (TreeTime)...")
    # Prepare dates file
    dates_path = outdir / "dates.csv"
    archetype_dates = json.load(open("data/metadata/archetype_dates.json"))
    with open(dates_path, "w") as f:
        f.write("name,date\n")
        for name, date in archetype_dates.items():
            f.write(f"{name},{date}\n")
        f.write(f"{args.sample_id},{datetime.now().year + datetime.now().month/12.0:.2f}\n")
    
    treetime_out = outdir / "treetime_output"
    treetime_cmd = [
        "treetime",
        "--tree", str(outdir / "iqtree_output.treefile"),
        "--aln", str(merged_alignment),
        "--dates", str(dates_path),
        "--outdir", str(treetime_out)
    ]
    # Note: TreeTime might fail if divergence is too low/insufficient signal, 
    # or if the executable is missing.
    try:
        run_cmd(treetime_cmd, check=False)
    except FileNotFoundError:
        print("⚠️  Warning: treetime executable not found. Skipping temporal scaling.")
    except Exception as e:
        print(f"⚠️  Warning: TreeTime failed: {e}")

    # 5. Final Outputs & Visualization
    print("Step 5: Finalizing Artifacts...")
    # Newick
    final_tree_path = outdir / "tree.nwk"
    iqtree_tree = outdir / "iqtree_output.treefile"
    if iqtree_tree.exists():
        shutil.copy(iqtree_tree, final_tree_path)
    
    # PNG using shared utility
    render_tree_png(
        str(final_tree_path),
        str(outdir / "tree.png"),
        title=f"CholeraSeq Modular Phylogeny: {args.sample_id}"
    )

    print("✅ CholeraSeq Tree Module Complete.")

if __name__ == "__main__":
    main()
