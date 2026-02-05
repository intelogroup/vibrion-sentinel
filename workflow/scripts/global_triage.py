import os
import sys
import json
import argparse
import tempfile
from pathlib import Path
import subprocess
import csv

def run_sourmash_compare(query_sig, ref_sigs, output_json):
    """Compare query signature against a list of reference signatures."""
    results = []
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for ref_id, ref_sig in ref_sigs.items():
            search_csv = os.path.join(tmpdir, f"search_{ref_id}.csv")
            
            # Use 'sourmash search' for similarity and containment
            # -o for CSV output in newer sourmash
            cmd = [
                "sourmash", "search", query_sig, ref_sig, 
                "--threshold", "0.001", 
                "-o", search_csv
            ]
            try:
                subprocess.run(cmd, capture_output=True, check=True)
                
                if os.path.exists(search_csv):
                    with open(search_csv) as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            # similarity,name,filename,md5,query_filename,query_name,query_md5,containment,max_containment
                            results.append({
                                "lineage": ref_id,
                                "similarity": float(row.get("similarity", 0.0)),
                                "containment": float(row.get("containment", 0.0))
                            })
                            # We only need the best match per reference file for this loop
                            break
            except subprocess.CalledProcessError as e:
                print(f"Warning: sourmash search failed for {ref_id}: {e.stderr.decode()}", file=sys.stderr)
    
    # Sort by similarity
    results = sorted(results, key=lambda x: x['similarity'], reverse=True)
    
    with open(output_json, 'w') as f:
        json.dump({
            "global_matches": results, 
            "status": "SUCCESS",
            "count": len(results)
        }, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description="Vibrion Global Triage")
    parser.add_argument("--sample-fasta", required=True)
    parser.add_argument("--ref-dir", required=True)
    parser.add_argument("--lineages", nargs="+", required=True)
    parser.add_argument("--ksize", type=int, default=31)
    parser.add_argument("--scaled", type=int, default=1000)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ref_dir = Path(args.ref_dir)
    
    with tempfile.NamedTemporaryFile(suffix=".sig", delete=False) as tmp_query:
        query_sig = tmp_query.name
    
    try:
        # Sketch query
        subprocess.run([
            "sourmash", "sketch", "dna", args.sample_fasta, 
            "-o", query_sig, "--param-string", f"k={args.ksize},scaled={args.scaled}"
        ], check=True, capture_output=True)
        
        ref_sigs = {}
        for lineage in args.lineages:
            # Look for .sig first, then .fasta
            ref_sig = ref_dir / f"{lineage}.sig"
            ref_fasta = ref_dir / f"{lineage}.fasta"
            
            if ref_sig.exists():
                ref_sigs[lineage] = str(ref_sig)
            elif ref_fasta.exists():
                # Sketch ref if fasta exists but sig doesn't
                try:
                    subprocess.run([
                        "sourmash", "sketch", "dna", str(ref_fasta), 
                        "-o", str(ref_sig), "--param-string", f"k={args.ksize},scaled={args.scaled}"
                    ], check=True, capture_output=True)
                    ref_sigs[lineage] = str(ref_sig)
                except subprocess.CalledProcessError as e:
                    print(f"Warning: Failed to sketch {lineage}: {e.stderr.decode()}", file=sys.stderr)
                    continue
                
        run_sourmash_compare(query_sig, ref_sigs, args.output)
    finally:
        if os.path.exists(query_sig):
            os.remove(query_sig)

if __name__ == "__main__":
    main()