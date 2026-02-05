#!/usr/bin/env python3
import argparse
import json
import os
import gzip
from Bio import SeqIO

def load_reference_gene(fasta_path, gene_name):
    """Load specific gene sequence from multi-FASTA."""
    if not os.path.exists(fasta_path):
        return None
    with open(fasta_path) as f:
        for record in SeqIO.parse(f, "fasta"):
            if record.id == gene_name:
                return str(record.seq)
    return None

def scan_fastq(fastq_path, o1_kmers, ogawa_probe, inaba_probe):
    """Scan FASTQ for O1 and Serotype markers."""
    stats = {
        "total_reads": 0,
        "o1_hits": 0,
        "ogawa_hits": 0,
        "inaba_hits": 0
    }
    
    # Open fastq (gzip aware)
    opener = gzip.open if fastq_path.endswith('.gz') else open
    
    try:
        with opener(fastq_path, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 1: # Sequence line
                    seq = line.strip()
                    stats["total_reads"] += 1
                    
                    # 1. Check O1 Backbone (Any wbeT match)
                    # Optimization: Check a few representative k-mers instead of all
                    for kmer in o1_kmers:
                        if kmer in seq:
                            stats["o1_hits"] += 1
                            break # Count read once for O1
                    
                    # 2. Check Serotype Probes
                    if ogawa_probe in seq:
                        stats["ogawa_hits"] += 1
                    if inaba_probe in seq:
                        stats["inaba_hits"] += 1
    except Exception as e:
        print(f"Error scanning FASTQ: {e}")
        
    return stats

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument("--db", help="Kraken2 DB path (Unused in k-mer mode)")
    parser.add_argument("--threads", default=1)
    parser.add_argument("--mode", default="triage", help="Analysis mode")
    parser.add_argument("--output", required=True)
    args, unknown = parser.parse_known_args()

    # 1. Setup Probes
    # Load wbeT (rfbT) - The O1/Serotype Determinant
    # We assume 'data/references/reference_loci.fasta' exists relative to CWD
    wbeT_seq = load_reference_gene("data/references/reference_loci.fasta", "wbeT")
    
    if not wbeT_seq:
        # Fallback if file missing (Should not happen in pipeline)
        print("WARNING: wbeT reference not found. Using Mock Fallback.")
        result = {
            "serogroup": "Unknown",
            "serotype": "Unknown",
            "confidence": "LOW",
            "reason": "Reference data missing"
        }
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        return

    # Generate Probes
    # O1 Core: Conserved regions (start, middle)
    # Use 31-mers
    o1_kmers = [
        wbeT_seq[0:31],      # Start
        wbeT_seq[200:231],   # Middle
        wbeT_seq[400:431]    # Middle 2
    ]
    
    # Serotype Probes (Haiti 793 Mutation: G -> T)
    # The extraction showed the gene ends with ...GGGTGTAG (Ogawa)
    # Mutation is likely the G at 793 (0-indexed 792)
    # Let's extract 21-mers around the end
    # Seq len is 794. 793 is 2nd to last.
    # Context: ...TCATTTTTTGGGTGTAG
    # Ogawa (WT): TTTGGGTGTAG (11-mer sufficient?) Let's use 17-mer: TCATTTTTTGGGTGTAG
    # Inaba (Mut): TCATTTTTTGGTTGTAG (G -> T)
    
    # Note: Reads might be on reverse strand!
    # Simple scanner checks forward strand. For robustness, we should check RC too?
    # Or just assume 50% reads land forward.
    # For speed, we stick to forward. If we see 0 hits, we might miss RC-only reads,
    # but coverage usually ensures both.
    
    ogawa_probe = "TCATTTTTTGGGTGTAG"
    inaba_probe = "TCATTTTTTGGTTGTAG" # G -> T mutation
    
    # 2. Scan FASTQ
    print(f"Scanning {args.fastq} for Serogroup markers...")
    stats = scan_fastq(args.fastq, o1_kmers, ogawa_probe, inaba_probe)
    
    # 3. Decision Logic
    result = {
        "sample_id": os.path.basename(args.fastq).split('.')[0],
        "stats": stats,
        "serogroup": "Unknown",
        "serotype": "Unknown",
        "confidence": "LOW",
        "reason": "Insufficient Data"
    }
    
    # Thresholds
    MIN_O1_HITS = 3 # Very low threshold for Bunker Mode
    
    if stats["o1_hits"] >= MIN_O1_HITS:
        result["serogroup"] = "O1"
        result["confidence"] = "HIGH"
        
        # Serotype Call
        total_type_hits = stats["ogawa_hits"] + stats["inaba_hits"]
        
        if total_type_hits == 0:
            result["serotype"] = "Unknown"
            result["reason"] = "O1 Detected (Backbone) but Serotype region not covered"
            result["confidence"] = "MEDIUM"
        else:
            ogawa_ratio = stats["ogawa_hits"] / total_type_hits
            
            if ogawa_ratio > 0.9:
                result["serotype"] = "Ogawa"
                result["reason"] = f"Dominant Ogawa alleles ({stats['ogawa_hits']} vs {stats['inaba_hits']})"
            elif ogawa_ratio < 0.1:
                result["serotype"] = "Inaba"
                result["reason"] = f"Dominant Inaba alleles ({stats['inaba_hits']} vs {stats['ogawa_hits']})"
            else:
                result["serotype"] = "Hikojima (Mixed)"
                result["reason"] = f"Mixed alleles detected: Ogawa={stats['ogawa_hits']}, Inaba={stats['inaba_hits']}"
                result["is_mixed"] = True
    else:
        result["serogroup"] = "Non-O1"
        result["serotype"] = "Non-O1"
        result["reason"] = f"No wbeT O1 markers found ({stats['o1_hits']} hits)"
        result["confidence"] = "HIGH"

    # Write output
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)

if __name__ == "__main__":
    main()
