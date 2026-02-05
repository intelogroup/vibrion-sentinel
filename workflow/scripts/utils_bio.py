#!/usr/bin/env python3
"""
Shared bioinformatics utilities for Vibrion Sentinel.
Consolidates common operations: loci extraction, MAFFT alignment, tree building.
"""

import shutil
import subprocess
from pathlib import Path
from typing import Optional
from Bio import SeqIO, Phylo, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def get_mafft_path() -> str:
    """Get MAFFT executable path (conda env or system)."""
    system_mafft = shutil.which("mafft")
    if system_mafft:
        return system_mafft
    
    raise FileNotFoundError("MAFFT not found in conda env or system PATH")


def extract_loci_from_bed(fasta_path: str, bed_path: str, seq_id: str) -> Optional[SeqRecord]:
    """
    Extract and concatenate genomic loci defined in a BED file.
    
    Args:
        fasta_path: Path to input FASTA file
        bed_path: Path to BED file with coordinates
        seq_id: ID for the output sequence record
    
    Returns:
        SeqRecord with concatenated loci, or None if extraction fails
    """
    fasta_path = Path(fasta_path)
    bed_path = Path(bed_path)
    
    if not fasta_path.exists():
        print(f"Warning: {fasta_path} not found.")
        return None
    
    if not bed_path.exists():
        print(f"Warning: {bed_path} not found.")
        return None
    
    try:
        # Read BED file
        with open(bed_path) as f:
            bed_lines = [line for line in f if not line.startswith("#") and line.strip()]
        
        # Parse reference sequences into a dictionary
        seq_dict = {record.id: record.seq for record in SeqIO.parse(fasta_path, "fasta")}
        
        # Extract and concatenate loci
        combined_seq = ""
        for line in bed_lines:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            
            if chrom in seq_dict:
                combined_seq += str(seq_dict[chrom][start:end])
            else:
                # Fallback: if chrom not found by ID, try index if ID looks like 1 or 2
                print(f"Warning: Chromosome {chrom} not found in {fasta_path}. Skipping locus.")
        
        if not combined_seq:
            print(f"Error: No loci could be extracted from {fasta_path}")
            return None

        return SeqRecord(Seq(combined_seq), id=seq_id, description=f"Core loci from {fasta_path.name}")
    
    except Exception as e:
        print(f"Error extracting loci from {fasta_path}: {e}")
        return None


def run_mafft_alignment(input_fasta: str, output_fasta: str, threads: int = 4, 
                        add_to_alignment: Optional[str] = None) -> bool:
    """
    Run MAFFT alignment.
    
    Args:
        input_fasta: Input sequences to align
        output_fasta: Output aligned sequences
        threads: Number of threads
        add_to_alignment: If provided, add input_fasta to this existing alignment
    
    Returns:
        True if successful, False otherwise
    """
    mafft_path = get_mafft_path()
    
    try:
        if add_to_alignment:
            # Add mode: append new sequence to existing alignment
            cmd = [mafft_path, "--add", str(input_fasta), "--reorder", 
                   "--thread", str(threads), str(add_to_alignment)]
        else:
            # De novo alignment
            cmd = [mafft_path, "--auto", "--thread", str(threads), str(input_fasta)]
        
        with open(output_fasta, "w") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, 
                                  check=True, text=True)
        return True
    
    except subprocess.CalledProcessError as e:
        print(f"MAFFT failed: {e.stderr}")
        return False
    except Exception as e:
        print(f"MAFFT error: {e}")
        return False


def build_upgma_tree(alignment_path: str, output_tree: str, method: str = "upgma") -> bool:
    """
    Build phylogenetic tree using distance-based method (UPGMA or NJ).
    
    Args:
        alignment_path: Input aligned FASTA
        output_tree: Output newick tree file
        method: "upgma" or "nj" (neighbor-joining)
    
    Returns:
        True if successful, False otherwise
    """
    try:
        alignment = AlignIO.read(alignment_path, "fasta")
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        
        if method == "upgma":
            tree = constructor.upgma(dm)
        else:
            tree = constructor.nj(dm)
        
        Phylo.write(tree, output_tree, "newick")
        return True
    
    except Exception as e:
        print(f"Tree building failed: {e}")
        return False


def render_tree_png(tree_path: str, output_png: str, title: str = "Phylogenetic Tree",
                   figsize: tuple = (12, 10), dpi: int = 300) -> bool:
    """
    Render phylogenetic tree as PNG image.
    
    Args:
        tree_path: Input newick tree file
        output_png: Output PNG file
        title: Plot title
        figsize: Figure size (width, height)
        dpi: Image resolution
    
    Returns:
        True if successful, False otherwise
    """
    try:
        tree = Phylo.read(tree_path, "newick")
        tree.ladderize()
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, do_show=False, axes=ax)
        plt.title(title)
        plt.savefig(output_png, dpi=dpi, bbox_inches='tight')
        plt.close(fig)
        return True
    
    except Exception as e:
        print(f"Tree rendering failed: {e}")
        return False


def filter_empty_sequences(input_fasta: str, output_fasta: str, min_bases: int = 1) -> int:
    """
    Filter out sequences with too few real bases (excluding N and gaps).
    
    Args:
        input_fasta: Input FASTA file
        output_fasta: Output filtered FASTA
        min_bases: Minimum number of real bases required
    
    Returns:
        Number of sequences kept
    """
    kept = 0
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            seq_str = str(record.seq).replace("-", "").replace("N", "").replace("n", "")
            if len(seq_str) >= min_bases:
                SeqIO.write(record, out_f, "fasta")
                kept += 1
            else:
                print(f"   Filtered out empty sequence: {record.id}")
    return kept


def check_sequence_quality(seq_record: SeqRecord, ref_alignment_path: Optional[str] = None,
                          max_ambiguous_pct: float = 50.0) -> tuple[bool, str]:
    """
    Check if sequence is suitable for phylogenetic analysis.
    
    Args:
        seq_record: Sequence to check
        ref_alignment_path: Optional reference alignment to check for duplicates
        max_ambiguous_pct: Maximum percentage of N/gaps allowed
    
    Returns:
        (is_valid, message) tuple
    """
    seq_str = str(seq_record.seq).upper()
    total_len = len(seq_str)
    
    if total_len == 0:
        return False, "Empty sequence"
    
    # Count ambiguous bases
    n_count = seq_str.count('N')
    gap_count = seq_str.count('-')
    ambiguous_pct = (n_count + gap_count) / total_len * 100
    
    if ambiguous_pct > max_ambiguous_pct:
        return False, f"Too many ambiguous bases ({ambiguous_pct:.1f}%)"
    
    # Check for identity with reference (zero coverage indicator)
    if ref_alignment_path and Path(ref_alignment_path).exists():
        try:
            ref_records = list(SeqIO.parse(ref_alignment_path, "fasta"))
            sample_nogaps = seq_str.replace('-', '')
            for ref_rec in ref_records:
                ref_nogaps = str(ref_rec.seq).upper().replace('-', '')
                if sample_nogaps == ref_nogaps:
                    return False, f"Identical to reference {ref_rec.id} (zero coverage)"
        except Exception:
            pass  # Continue if reference check fails
    
    return True, "OK"
