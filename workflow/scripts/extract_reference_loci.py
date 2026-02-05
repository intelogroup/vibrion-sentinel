#!/usr/bin/env python3
"""
Extract reference baseline loci from 2010EL-1786.fasta (Haiti 2010 strain).
These serve as the baseline for delta-anomaly calculations in Evo2 analysis.
"""

import sys
import subprocess
from pathlib import Path


def extract_locus_from_reference(reference_path, chrom, start, end):
    """Extract sequence from reference FASTA using samtools faidx."""
    try:
        cmd = ['samtools', 'faidx', str(reference_path), f"{chrom}:{start+1}-{end}"]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        lines = result.stdout.strip().split('\n')
        sequence = ''.join(lines[1:])  # Skip header line
        
        return sequence.upper()
    except subprocess.CalledProcessError as e:
        print(f"Error extracting {chrom}:{start}-{end}: {e}", file=sys.stderr)
        return None


def main():
    # Parse Snakemake inputs
    reference_path = Path(snakemake.input.reference) # noqa: F821
    bed_path = Path(snakemake.input.bed) # noqa: F821
    output_fasta = Path(snakemake.output.reference_loci) # noqa: F821
    sample_id = snakemake.params.sample # noqa: F821
    
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    
    # Index reference if needed
    if not Path(str(reference_path) + ".fai").exists():
        subprocess.run(['samtools', 'faidx', str(reference_path)], check=True)
    
    print(f"📍 Extracting reference baseline loci from {reference_path.name}")
    
    loci_extracted = 0
    
    with open(output_fasta, 'w') as out:
        with open(bed_path) as bed:
            for line in bed:
                if line.startswith('#'):
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                
                sequence = extract_locus_from_reference(reference_path, chrom, start, end)
                
                if sequence:
                    out.write(f">{name}\n")
                    for i in range(0, len(sequence), 80):
                        out.write(sequence[i:i+80] + '\n')
                    
                    loci_extracted += 1
    
    print(f"   ✅ Extracted {loci_extracted} reference loci for baseline comparison")


if __name__ == "__main__":
    main()
