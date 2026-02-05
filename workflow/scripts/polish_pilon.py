#!/usr/bin/env python3
"""
Pilon Polishing Wrapper (Illumina)
Performs high-precision base correction using k-mer pileups
"""

import argparse
import subprocess
import sys
from pathlib import Path


def run_pilon(draft_fasta: Path, bam_file: Path, output_fasta: Path, 
              work_dir: Path) -> None:
    """
    Run Pilon polishing.
    
    Args:
        draft_fasta: Input draft consensus
        bam_file: Aligned reads BAM
        output_fasta: Output polished consensus
        work_dir: Working directory
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n🔬 Pilon Polishing (Illumina)")
    print(f"   Draft: {draft_fasta.name}")
    print(f"   BAM: {bam_file.name}")
    
    # Pilon requires indexed BAM
    print("   Indexing BAM...")
    subprocess.run(["samtools", "index", str(bam_file)], check=True)
    
    # Run Pilon
    pilon_output_prefix = work_dir / "pilon"
    
    cmd = [
        "pilon",
        "--genome", str(draft_fasta),
        "--frags", str(bam_file),
        "--output", "pilon",
        "--outdir", str(work_dir),
        "--fix", "snps,indels",  # Skip gaps/local (requires paired reads)
        "--changes",     # Output change log
        "--vcf"          # Output VCF of changes
    ]
    
    import os
    env = os.environ.copy()
    env["_JAVA_OPTIONS"] = "-Xmx4G"
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, env=env)
        
        # Pilon outputs as <prefix>.fasta
        pilon_fasta = work_dir / "pilon.fasta"
        
        if not pilon_fasta.exists():
            print(f"   ⚠️  Pilon stdout: {result.stdout}", file=sys.stderr)
            print(f"   ⚠️  Pilon stderr: {result.stderr}", file=sys.stderr)
            raise FileNotFoundError(f"Pilon output not found: {pilon_fasta}")
        
        # Strip '_pilon' suffix from headers for downstream compatibility
        print("   Cleaning fasta headers...")
        with open(pilon_fasta, 'r') as infile:
            with open(output_fasta, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        outfile.write(line.replace('_pilon', ''))
                    else:
                        outfile.write(line)
        
        print(f"   ✅ Complete: {output_fasta}")
        
        # Report changes
        changes_file = work_dir / "pilon.changes"
        if changes_file.exists():
            with open(changes_file) as f:
                changes = f.readlines()
            print(f"   📊 Changes made: {len(changes)}")
    
    
    except subprocess.CalledProcessError as e:
        print("   ❌ Pilon failed!", file=sys.stderr)
        print(f"   STDOUT: {e.stdout}", file=sys.stderr)
        print(f"   STDERR: {e.stderr}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"   ❌ Error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Pilon polishing (Illumina)")
    parser.add_argument("--draft", required=True, help="Draft consensus FASTA")
    parser.add_argument("--bam", required=True, help="Aligned reads BAM")
    parser.add_argument("--output", required=True, help="Output polished FASTA")
    parser.add_argument("--outdir", required=True, help="Working directory")
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("🧬 Pilon Polishing (Illumina)")
    print("=" * 60)
    
    run_pilon(
        Path(args.draft),
        Path(args.bam),
        Path(args.output),
        Path(args.outdir)
    )
    
    print("\n✅ Polishing complete!")


if __name__ == "__main__":
    main()
