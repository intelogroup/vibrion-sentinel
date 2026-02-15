import os
import sys
from pathlib import Path
from Bio import SeqIO

def calculate_n50(lengths):
    """Calculate N50 and L50 of a list of contig lengths."""
    lengths.sort(reverse=True)
    total_len = sum(lengths)
    running_sum = 0
    l50 = 0
    n50 = 0
    for i, length in enumerate(lengths):
        running_sum += length
        if running_sum >= total_len / 2:
            n50 = length
            l50 = i + 1
            break
    return n50, l50

def check_fasta(fasta_path):
    """Check assembly quality metrics (N50, L50, etc.)."""
    try:
        # Load and filter out tiny contigs (<200bp) which are usually artifacts
        all_records = list(SeqIO.parse(fasta_path, "fasta"))
        records = [r for r in all_records if len(r.seq) >= 200]
        
        if not records:
            return False, "No contigs found after 200bp filter"
        
        lengths = [len(r.seq) for r in records]
        total_len = sum(lengths)
        n50, l50 = calculate_n50(lengths)
        num_contigs = len(records)
        rejected = len(all_records) - len(records)
        
        if total_len < 3000000: # Vibrio cholerae is ~4MB
            return False, f"Genome size too small ({total_len} bp)"
            
        continuity_status = "STABLE" if n50 > 50000 else "FRAGMENTED"
        msg = (f"Size: {total_len/1e6:.2f}MB | N50: {n50:,} | L50: {l50} | "
               f"Contigs: {num_contigs} (Filtered {rejected}) | Status: {continuity_status}")
        
        # Industry standard: <500 contigs is great, <1000 is acceptable for surveillance
        is_valid = total_len > 3000000 and num_contigs < 1000
        return is_valid, msg
    except Exception as e:
        return False, str(e)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 scripts/check_references.py <ref_dir>")
        sys.exit(1)
        
    ref_dir = Path(sys.argv[1])
    if not ref_dir.exists():
        print(f"Directory {ref_dir} does not exist.")
        sys.exit(1)
        
    fastas = list(ref_dir.glob("*.fasta"))
    print(f"Checking {len(fastas)} references in {ref_dir}...")
    
    results = []
    for fasta in fastas:
        valid, msg = check_fasta(fasta)
        status = "✅" if valid else "❌"
        print(f"{status} {fasta.name}: {msg}")
        results.append(valid)
        
    if all(results) and len(results) > 0:
        print("\nAll references look good!")
        sys.exit(0)
    else:
        print("\nSome references failed validation.")
        sys.exit(1)

if __name__ == "__main__":
    main()
