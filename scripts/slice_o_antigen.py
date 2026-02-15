#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def run_blast(query, subject):
    """Run blastn of query against subject and return top hit coords."""
    output_fmt = "6"
    cmd = [
        "blastn",
        "-query", query,
        "-subject", subject,
        "-outfmt", output_fmt
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = result.stdout.strip().split("\n")
        if not lines or not lines[0]:
            return None
        
        # Take the top hit
        parts = lines[0].split("\t")
        # Field indices for -outfmt 6:
        # 0:qseqid, 1:sseqid, 2:pident, 3:length, 4:mismatch, 5:gaps, 6:qstart, 7:qstop, 8:sstart, 9:sstop, 10:evalue, 11:bitscore
        return {
            "pident": float(parts[2]),
            "length": int(parts[3]),
            "sstart": int(parts[8]),
            "sstop": int(parts[9]),
            "contig": parts[1]
        }
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e.stderr}")
        return None

def main():
    parser = argparse.ArgumentParser(description="Slice O-antigen region (gmhD to rjg)")
    parser.add_argument("input", help="Target genome FASTA")
    parser.add_argument("--gmhD", default="data/serogroup_reference/seeds/gmhD_seed.fasta", help="gmhD seed")
    parser.add_argument("--rjg", default="data/serogroup_reference/seeds/rjg_seed.fasta", help="rjg seed")
    parser.add_argument("--output", help="Output sliced FASTA")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found")
        sys.exit(1)
        
    print(f"Slicing O-antigen from {args.input}...")
    
    gmhD_hit = run_blast(args.gmhD, args.input)
    rjg_hit = run_blast(args.rjg, args.input)
    
    if not gmhD_hit or not rjg_hit:
        print("Error: Could not find one or both anchor genes (gmhD, rjg).")
        if not gmhD_hit: print("  gmhD: NOT FOUND")
        if not rjg_hit: print("  rjg: NOT FOUND")
        sys.exit(1)
        
    if gmhD_hit["contig"] != rjg_hit["contig"]:
        print(f"Error: Anchor genes found on different contigs: {gmhD_hit['contig']} and {rjg_hit['contig']}")
        sys.exit(1)
        
    contig_id = gmhD_hit["contig"]
    start = min(gmhD_hit["sstart"], gmhD_hit["sstop"], rjg_hit["sstart"], rjg_hit["sstop"])
    end = max(gmhD_hit["sstart"], gmhD_hit["sstop"], rjg_hit["sstart"], rjg_hit["sstop"])
    
    # Pad some length (optional, but paper says gmhD to rjg is the region)
    # We'll take exactly that range
    
    print(f"O-AGC found on {contig_id} at {start}-{end} (length: {end-start} bp)")
    
    # Load genome and extract
    genome = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    if contig_id not in genome:
        print(f"Error: Contig {contig_id} not found in genome dict")
        sys.exit(1)
        
    sub_seq = genome[contig_id].seq[start-1:end]
    
    sliced_record = SeqRecord(
        sub_seq,
        id=f"{contig_id}_OAGC",
        description=f"Sliced O-antigen biosynthesis gene cluster ({start}-{end})"
    )
    
    out_path = args.output if args.output else f"{os.path.splitext(args.input)[0]}_oagc.fasta"
    SeqIO.write(sliced_record, out_path, "fasta")
    print(f"Saved sliced region to {out_path}")

if __name__ == "__main__":
    main()
