#!/usr/bin/env python3
"""
Map surveillance loci to the provided reference genome using BLASTN.
Generates a JSON map of coordinates for downstream tools (call_variants, etc).
"""
import subprocess
import argparse
import json
import sys
from pathlib import Path

def map_loci(reference_fasta, loci_fasta, output_json):
    """
    Map surveillance loci to the provided reference genome using BLASTN.
    """
    # Verify inputs
    if not Path(reference_fasta).exists():
        print(f"Reference not found: {reference_fasta}", file=sys.stderr)
        sys.exit(1)
    if not Path(loci_fasta).exists():
        print(f"Loci FASTA not found: {loci_fasta}", file=sys.stderr)
        sys.exit(1)

    cmd = [
        "blastn",
        "-subject", str(reference_fasta),
        "-query", str(loci_fasta),
        "-outfmt", "6 qseqid sseqid sstart send pident length qlen slen",
        "-perc_identity", "90", # Stricter identity for surveillance (prevents pseudogene hits)
        "-evalue", "1e-10", # Filter weak hits
        "-max_target_seqs", "1" 
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLASTn failed: {e.stderr}", file=sys.stderr)
        sys.exit(1)
        
    mapping = {}
    
    for line in result.stdout.strip().split('\n'):
        if not line: continue
        parts = line.split('\t')
        gene = parts[0]
        chrom = parts[1]
        sstart = int(parts[2])
        send = int(parts[3])
        pident = float(parts[4])
        length = int(parts[5])
        qlen = int(parts[6])
        
        # Handle coordinates (1-based to 0-based) and strand
        if sstart < send:
            start = sstart - 1
            end = send
            strand = "+"
        else:
            start = send - 1
            end = sstart
            strand = "-"
            
        # Filter: Must cover at least 90% of the query gene
        # Ensures we are tracking the FULL functional gene
        if length / qlen < 0.9:
            print(f"Skipping {gene}: Low coverage ({length/qlen:.2%})")
            continue
            
        mapping[gene] = {
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "identity": pident,
            "description": f"Mapped via BLAST ({pident:.1f}%)"
        }
        
    with open(output_json, 'w') as f:
        json.dump(mapping, f, indent=2)
        
    print(f"Mapped {len(mapping)} loci to {reference_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True, help="Target reference FASTA")
    parser.add_argument("--loci", required=True, help="Query loci FASTA")
    parser.add_argument("--output", required=True, help="Output JSON mapping")
    args = parser.parse_args()
    
    map_loci(args.reference, args.loci, args.output)
