#!/usr/bin/env python3
"""
Medaka Polishing Wrapper (Nanopore)
Performs 2-round neural network polishing for homopolymer correction
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_medaka_round(draft_fasta: Path, reads_fastq: Path, output_dir: Path, 
                     model: str, round_num: int) -> Path:
    """
    Run a single round of Medaka polishing.
    
    Args:
        draft_fasta: Input draft consensus
        reads_fastq: Original reads for polishing
        output_dir: Output directory
        model: Medaka model (e.g., r1041_e82_400bps_hac_v4.2.0)
        round_num: Round number (1 or 2)
    
    Returns:
        Path to polished consensus
    """
    round_dir = output_dir / f"round_{round_num}"
    round_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"\n🔬 Medaka Round {round_num}")
    print(f"   Model: {model}")
    print(f"   Input: {draft_fasta.name}")
    
    # Medaka consensus command
    cmd = [
        "medaka_consensus",
        "-i", str(reads_fastq),
        "-d", str(draft_fasta),
        "-o", str(round_dir),
        "-m", model,
        "-t", "4"  # threads
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        polished = round_dir / "consensus.fasta"
        
        if not polished.exists():
            raise FileNotFoundError(f"Medaka output not found: {polished}")
        
        print(f"   ✅ Complete: {polished}")
        return polished
    
    except subprocess.CalledProcessError as e:
        print(f"   ❌ Medaka failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Medaka polishing (2 rounds)")
    parser.add_argument("--draft", required=True, help="Draft consensus FASTA")
    parser.add_argument("--reads", required=True, help="Original reads FASTQ")
    parser.add_argument("--output", required=True, help="Output polished FASTA")
    parser.add_argument("--outdir", required=True, help="Working directory")
    parser.add_argument("--model", default="r1041_e82_400bps_hac_v4.2.0",
                       help="Medaka model")
    parser.add_argument("--rounds", type=int, default=2, help="Polish rounds")
    
    args = parser.parse_args()
    
    draft_fasta = Path(args.draft)
    reads_fastq = Path(args.reads)
    output_fasta = Path(args.output)
    work_dir = Path(args.outdir)
    
    print("=" * 60)
    print("🧬 Medaka Neural Polishing (Nanopore)")
    print("=" * 60)
    
    current_draft = draft_fasta
    
    # Run polishing rounds
    for round_num in range(1, args.rounds + 1):
        polished = run_medaka_round(
            current_draft, reads_fastq, work_dir, 
            args.model, round_num
        )
        current_draft = polished
    
    # Copy final polished consensus to output
    import shutil
    shutil.copy(current_draft, output_fasta)
    
    print("\n✅ Polishing complete!")
    print(f"   Final: {output_fasta}")


if __name__ == "__main__":
    main()
