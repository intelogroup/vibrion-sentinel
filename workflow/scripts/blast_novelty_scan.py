#!/usr/bin/env python3
"""
BLAST Novelty Scan (Tier 3 Pre-Filter: Stranger Detection)

Scans consensus genome against reference pangenome to find:
- Novel sequences not present in reference (potential mobile elements, recombination)
- Low-identity regions (possible ancestral imports)

License: Uses NCBI BLAST+ (Public Domain)
"""

import subprocess
import os
import json
from Bio import SeqIO

def check_blast_available():
    """Check if BLAST+ is installed."""
    try:
        result = subprocess.run(["blastn", "-version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def create_blast_db(reference_fasta, db_prefix):
    """Create a BLAST database from reference genome."""
    cmd = [
        "makeblastdb",
        "-in", reference_fasta,
        "-dbtype", "nucl",
        "-out", db_prefix
    ]
    subprocess.run(cmd, check=True)

def run_blast(query_fasta, db_prefix, output_file, threads=4):
    """Run BLAST alignment of consensus against reference."""
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_prefix,
        "-out", output_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen",
        "-num_threads", str(threads),
        "-evalue", "1e-10"
    ]
    subprocess.run(cmd, check=True)

def find_novel_regions(blast_output, query_fasta, min_identity=80, min_gap_length=500):
    """
    Identify regions in query that have low/no coverage in reference.
    
    Returns:
        List of (contig, start, end, type) tuples for 'stranger' regions
    """
    # Parse BLAST hits
    hits = []
    if os.path.exists(blast_output) and os.path.getsize(blast_output) > 0:
        with open(blast_output, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 13:
                    qseqid = parts[0]
                    pident = float(parts[2])
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    qlen = int(parts[12])
                    
                    if pident >= min_identity:
                        hits.append({
                            'contig': qseqid,
                            'start': min(qstart, qend),
                            'end': max(qstart, qend),
                            'pident': pident,
                            'qlen': qlen
                        })
    
    # Get query contig lengths
    contig_lengths = {}
    for record in SeqIO.parse(query_fasta, "fasta"):
        contig_lengths[record.id] = len(record.seq)
    
    # Find gaps in coverage (novel regions)
    strangers = []
    
    for contig, length in contig_lengths.items():
        # Get all hits for this contig
        contig_hits = [h for h in hits if h['contig'] == contig]
        
        # Sort by start position
        contig_hits.sort(key=lambda x: x['start'])
        
        # Find gaps
        covered = set()
        for hit in contig_hits:
            covered.update(range(hit['start'], hit['end'] + 1))
        
        # Find uncovered stretches
        uncovered = []
        start = None
        for pos in range(1, length + 1):
            if pos not in covered:
                if start is None:
                    start = pos
            else:
                if start is not None:
                    if pos - start >= min_gap_length:
                        uncovered.append((start, pos - 1))
                    start = None
        # Handle trailing gap
        if start is not None and length - start >= min_gap_length:
            uncovered.append((start, length))
        
        for gap_start, gap_end in uncovered:
            strangers.append({
                'contig': contig,
                'start': gap_start,
                'end': gap_end,
                'length': gap_end - gap_start + 1,
                'type': 'novel_region'
            })
    
    # Also check for low-identity hits (potential recombination)
    for hit in hits:
        if hit['pident'] < 95 and (hit['end'] - hit['start']) > min_gap_length:
            strangers.append({
                'contig': hit['contig'],
                'start': hit['start'],
                'end': hit['end'],
                'length': hit['end'] - hit['start'],
                'type': 'low_identity',
                'identity': hit['pident']
            })
    
    return strangers

def extract_stranger_sequences(query_fasta, strangers, output_fasta):
    """Extract stranger sequences to FASTA for Evo2 analysis."""
    sequences = {record.id: record for record in SeqIO.parse(query_fasta, "fasta")}
    
    with open(output_fasta, 'w') as f:
        for i, stranger in enumerate(strangers):
            contig = stranger['contig']
            if contig in sequences:
                seq = sequences[contig].seq[stranger['start']-1:stranger['end']]
                f.write(f">{contig}_stranger_{i+1}_{stranger['start']}-{stranger['end']}_{stranger['type']}\n")
                f.write(str(seq) + "\n")

def main():
    # Snakemake integration
    consensus = snakemake.input.consensus # noqa: F821
    reference = snakemake.input.reference # noqa: F821
    
    output_bed = snakemake.output.stranger_bed # noqa: F821
    output_fasta = snakemake.output.stranger_fasta # noqa: F821
    stats_file = snakemake.output.stats # noqa: F821
    log_file = snakemake.log[0] # noqa: F821
    
    threads = snakemake.threads # noqa: F821
    min_gap = snakemake.params.get("min_gap_length", 500) # noqa: F821
    min_identity = snakemake.params.get("min_identity", 80) # noqa: F821
    
    output_dir = os.path.dirname(output_bed)
    os.makedirs(output_dir, exist_ok=True)
    
    # Paths
    db_prefix = os.path.join(output_dir, "ref_db")
    blast_output = os.path.join(output_dir, "blast_results.txt")
    
    if not check_blast_available():
        raise RuntimeError("BLAST+ not found. Install with: conda install -c bioconda blast")
    
    try:
        # 1. Create BLAST database from reference
        print(f"Creating BLAST database from {reference}...")
        create_blast_db(reference, db_prefix)
        
        # 2. Run BLAST
        print("Running BLAST alignment...")
        run_blast(consensus, db_prefix, blast_output, threads)
        
        # 3. Find novel regions
        print("Analyzing for novel/stranger regions...")
        strangers = find_novel_regions(blast_output, consensus, min_identity, min_gap)
        
        # 4. Write BED file
        with open(output_bed, 'w') as f:
            for s in strangers:
                f.write(f"{s['contig']}\t{s['start']}\t{s['end']}\t{s['type']}\n")
        
        # 5. Extract sequences
        extract_stranger_sequences(consensus, strangers, output_fasta)
        
        # 6. Stats
        stats = {
            "method": "blast_novelty_scan",
            "novel_regions_found": len([s for s in strangers if s['type'] == 'novel_region']),
            "low_identity_regions": len([s for s in strangers if s['type'] == 'low_identity']),
            "total_stranger_bp": sum(s['length'] for s in strangers),
            "strangers": strangers[:20]  # Top 20 for report
        }
        
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        print(f"Found {len(strangers)} stranger regions.")
        
    except Exception as e:
        print(f"Error in BLAST novelty scan: {e}")
        raise

if __name__ == "__main__":
    main()
