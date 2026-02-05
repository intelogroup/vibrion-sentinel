#!/usr/bin/env python3
import json
import subprocess
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-sig", required=True)
    parser.add_argument("--ref-dir", default="data/references")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # List of reference signatures to compare against
    ref_sigs = [
        "N16961.fasta.sig",
        "2010EL-1786.fasta.sig",
        "Haiti_2022_Resurgence.fasta.sig",
        "Yemen_2017.fasta.sig",
        "Malawi_2023.fasta.sig",
        "India_Wave3.fasta.sig"
    ]

    valid_sigs = []
    for sig in ref_sigs:
        path = os.path.join(args.ref_dir, sig)
        if os.path.exists(path):
            valid_sigs.append(path)

    if not valid_sigs:
        # Fallback to standard Haiti reference
        with open(args.output, "w") as f:
            json.dump({"best_reference": "data/references/2010EL-1786.fasta", "distances": {}}, f)
        return

    # Run sourmash compare
    results = {}
    for sig_path in valid_sigs:
        ref_name = os.path.basename(sig_path).replace(".fasta.sig", "")
        try:
            cmd = ["sourmash", "compare", args.sample_sig, sig_path, "--containment"]
            # sourmash compare outputs a matrix, we just want the similarity
            # Actually, sourmash comparison is easier with 'gather' for proportions, 
            # but for distance between two signatures, we use 'compare' or custom python.
            
            # Simple python API usage would be better but let's stick to CLI for robustness
            res = subprocess.run(cmd, capture_output=True, text=True)
            # Matrix output is usually:
            # 0-sample.sig [1.0 0.9]
            # 1-ref.sig    [0.9 1.0]
            # We parse the first line
            for line in res.stdout.split("\n"):
                if "0-" in line:
                    parts = line.strip().split("[")[1].split("]")[0].split()
                    if len(parts) >= 2:
                        similarity = float(parts[1])
                        results[ref_name] = similarity
        except Exception as e:
            print(f"Error comparing against {ref_name}: {e}")

    # Determine best reference
    best_ref = "2010EL-1786" # Default
    if results:
        best_ref = max(results, key=results.get)
        # TIE BREAKER: If differences are < 0.001, prioritize the sample's likely region if found in name
        max_sim = results[best_ref]
        sample_name = os.path.basename(args.sample_sig).lower()
        
        for name, sim in results.items():
            if abs(sim - max_sim) < 0.005: # Narrow margin
                if name.lower() in sample_name:
                    best_ref = name
                    break

    # Construct the path
    ref_path = f"data/references/{best_ref}.fasta"
    
    output_data = {
        "best_reference": ref_path,
        "best_match": best_ref,
        "similarity": results.get(best_ref, 0),
        "all_similarities": results
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=2)

if __name__ == "__main__":
    main()
