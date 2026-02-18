#!/usr/bin/env python3
import json
import subprocess
import argparse
import os
import csv
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-sig", required=True)
    parser.add_argument("--ref-dir", default="data/references")
    parser.add_argument("--output", required=True)
    parser.add_argument("--serogroup", help="Path to serogroup report for tie-breaking")
    parser.add_argument("--min-evidence", type=int, default=500, help="Minimum hashes for confidence")
    args = parser.parse_args()

    # List of reference signatures to compare against
    ref_sigs = [
        "N16961.fasta.sig",
        "2010EL-1786.fasta.sig",
        "Haiti_2022_Resurgence.fasta.sig",
        "Yemen_2017.fasta.sig",
        "Malawi_2023.fasta.sig",
        "India_Wave3.fasta.sig",
        "O139_MO10.fasta.sig",
        "O139_SG24.fasta.sig",
        "Inaba_A487.fasta.sig"
    ]

    valid_sigs = []
    for sig in ref_sigs:
        path = os.path.join(args.ref_dir, sig)
        if os.path.exists(path):
            valid_sigs.append(path)

    if not valid_sigs:
        # Fallback to standard Haiti reference
        with open(args.output, "w") as f:
            json.dump({"best_reference": "data/references/2010EL-1786.fasta", "best_match": "2010EL-1786", "similarity": 0, "distances": {}}, f)
        return

    # Run sourmash search for better precision and robustness
    results = {}
    csv_out = args.output + ".search.csv"
    try:
        # Use --containment to handle low-coverage sample reads in full reference
        cmd = ["sourmash", "search", args.sample_sig, *valid_sigs, "--containment", "--output", csv_out, "--threshold", "0.0"]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if os.path.exists(csv_out):
            with open(csv_out, "r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    # 'similarity' in search CSV is containment
                    sim = float(row["similarity"])
                    name = os.path.basename(row["filename"]).replace(".fasta.sig", "")
                    results[name] = sim
            os.remove(csv_out)
    except Exception as e:
        print(f"Error during sourmash search: {e}", file=sys.stderr)

    # If search fails, try to populate results with 0 for all valid sigs
    if not results:
        for sig_path in valid_sigs:
            name = os.path.basename(sig_path).replace(".fasta.sig", "")
            results[name] = 0.0

    # Load serogroup info for tie-breaking
    serogroup = "Unknown"
    if args.serogroup and os.path.exists(args.serogroup):
        try:
            with open(args.serogroup, "r") as f:
                sg_data = json.load(f)
                serogroup = sg_data.get("serogroup", "Unknown")
        except:
            pass

    # Determine best reference
    best_ref = "2010EL-1786" # Default
    if results:
        max_sim = max(results.values())
        # Use a wider margin for low-coverage samples (1%)
        margin = 0.01 if max_sim < 0.2 else 0.005
        top_candidates = [name for name, sim in results.items() if abs(sim - max_sim) < margin]
        
        # Tie breaking logic
        if len(top_candidates) > 1:
            # 1. If Non-O1, prefer O139 explicitly
            if "Non-O1" in serogroup or "O139" in serogroup:
                o139_candidates = [c for c in top_candidates if "O139" in c]
                if o139_candidates:
                    best_ref = o139_candidates[0]
                else:
                    best_ref = top_candidates[0]
            # 2. If O1 or Unknown, prefer most recent outbreak matching the name if possible
            else:
                sample_name = os.path.basename(args.sample_sig).lower()
                named_matches = [c for c in top_candidates if c.lower() in sample_name]
                if named_matches:
                    best_ref = named_matches[0]
                else:
                    outbreak_candidates = [c for c in top_candidates if "Resurgence" in c or "2022" in c]
                    if outbreak_candidates:
                        best_ref = outbreak_candidates[0]
                    else:
                        best_ref = top_candidates[0]
        else:
            best_ref = max(results, key=results.get)

    # Construct the path
    # Special case: Inaba_A487 has no assembly, map to N16961 instead
    if best_ref == "Inaba_A487":
        ref_path = "data/references/N16961.fasta"
    else:
        ref_path = f"data/references/{best_ref}.fasta"
    
    output_data = {
        "best_reference": ref_path,
        "best_match": best_ref,
        "similarity": results.get(best_ref, 0),
        "serogroup_context": serogroup,
        "all_similarities": results
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=2)

if __name__ == "__main__":
    main()
