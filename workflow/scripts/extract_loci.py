#!/usr/bin/env python3
"""
Extract surveillance loci sequences from reference genome.
Simple extraction of reference sequences at BED coordinates.
"""

import os
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO

def find_locus_via_blast(query_name, target_fasta):
    """
    Fallback: If faidx fails (e.g. ID mismatch), search for the locus using BLAST.
    Uses data/references/reference_loci.fasta as the sequence source.
    """
    ref_loci_path = "data/references/reference_loci.fasta"
    if not os.path.exists(ref_loci_path):
        return None

    # 1. Extract the query sequence from reference_loci.fasta
    query_seq = None
    for record in SeqIO.parse(ref_loci_path, "fasta"):
        if record.id == query_name:
            query_seq = str(record.seq)
            break
    
    if not query_seq:
        return None

    # 2. BLAST against target
    with tempfile.NamedTemporaryFile(suffix=".fasta", mode="w") as tmp_query:
        tmp_query.write(f">{query_name}\n{query_seq}\n")
        tmp_query.flush()

        cmd = [
            "blastn",
            "-query", tmp_query.name,
            "-subject", str(target_fasta),
            "-outfmt", "6 sseqid sstart send pident length",
            "-perc_identity", "70", # Relaxed from 80 for "Deep Drift"
            "-word_size", "7",      # More sensitive seeding
            "-max_target_seqs", "1"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = result.stdout.strip().split("\n")
            if not hits or not hits[0]:
                return None
            
            # 3. Extract the hit sequence
            parts = hits[0].split("\t")
            sseqid, sstart, send = parts[0], int(parts[1]), int(parts[2])
            
            # Ensure sstart < send for faidx
            start, end = min(sstart, send), max(sstart, send)
            print(f"  🔍 Discovered {query_name} in consensus via BLAST at {sseqid}:{start}-{end}")
            
            return extract_locus_from_reference(target_fasta, sseqid, start-1, end)
        except Exception as e:
            print(f"  Warning: BLAST discovery failed for {query_name}: {e}")
            return None

def extract_locus_from_reference(reference_path, chrom, start, end):
    """
    Extract sequence from reference FASTA using samtools faidx.
    """
    # Fix for coordinate systems (1-based vs 0-based)
    coord_str = f"{chrom}:{start+1}-{end}"
    cmd = ['samtools', 'faidx', str(reference_path), coord_str]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split('\n')
        if len(lines) < 2:
            return None
        sequence = ''.join(lines[1:])
        return sequence.upper()
    except subprocess.CalledProcessError:
        return None


def main():
    # Parse Snakemake inputs
    reference_path = Path(snakemake.input.reference) # noqa: F821
    bed_path = Path(snakemake.input.bed) # noqa: F821
    output_fasta = Path(snakemake.output.loci_fasta) # noqa: F821
    log_file = Path(snakemake.log[0]) # noqa: F821
    
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Index reference if missing or stale
    fai_path = Path(str(reference_path) + ".fai")
    if not fai_path.exists() or fai_path.stat().st_mtime < reference_path.stat().st_mtime:
        print(f"   Indexing reference: {reference_path}")
        subprocess.run(['samtools', 'faidx', str(reference_path)], check=True)
    
    with open(log_file, 'w') as log:
        log.write(f"Extracting loci from reference: {reference_path}\n")
        log.write(f"BED file: {bed_path}\n\n")
        
        loci_extracted = 0
        
        with open(output_fasta, 'w') as out:
            with open(bed_path) as bed:
                for line in bed:
                    # Skip comments
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue
                    
                    chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
                    
                    log.write(f"Extracting {name} ({chrom}:{start}-{end})\n")
                    
                    sequence = extract_locus_from_reference(
                        reference_path, chrom, start, end
                    )
                    
                    # DISCOVERY FALLBACK: If coordinate extraction failed, try BLAST search
                    if not sequence:
                        log.write(f"  ⚠️  Coordinates {chrom}:{start}-{end} failed. Attempting BLAST discovery...\n")
                        sequence = find_locus_via_blast(name, reference_path)
                    
                    if sequence:
                        # Write FASTA entry
                        out.write(f">{name}\n")
                        # Wrap at 80 characters
                        for i in range(0, len(sequence), 80):
                            out.write(sequence[i:i+80] + '\n')
                        
                        loci_extracted += 1
                        log.write(f"  Extracted {len(sequence)} bp\n")
                    else:
                        log.write("  ERROR: Failed to extract sequence\n")
        
        log.write(f"\nTotal loci extracted: {loci_extracted}\n")
        print(f"Extracted {loci_extracted} loci to {output_fasta}")


if __name__ == "__main__":
    main()
