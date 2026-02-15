#!/usr/bin/env python3
import argparse
import subprocess
import os
import sys
import json

def main():
    parser = argparse.ArgumentParser(description="Universal V. cholerae Serogroup Detector")
    parser.add_argument("input", help="Path to input genome FASTA (assembly)")
    parser.add_argument("--db", default="data/kraken2_serogroup", help="Path to Kraken2 serogroup DB")
    parser.add_argument("--output-dir", default="data/pipeline_output", help="Output directory")
    parser.add_argument("--keep-slice", action="store_true", help="Keep the sliced O-AGC FASTA")
    
    args = parser.parse_args()
    
    # Ensure relative paths are handled if needed, but we use absolute paths where possible
    script_dir = os.path.dirname(os.path.abspath(__file__))
    workspace_dir = os.path.dirname(script_dir)
    
    db_path = os.path.join(workspace_dir, args.db)
    slicer_script = os.path.join(script_dir, "slice_o_antigen.py")
    detector_script = os.path.join(script_dir, "detect_serogroup.py")
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    base_name = os.path.basename(args.input).split(".")[0]
    slice_path = os.path.join(args.output_dir, f"{base_name}_o_antigen.fasta")
    
    print(f"--- Step 1: Slicing O-antigen region ---")
    slice_cmd = [
        "python3", slicer_script,
        args.input,
        "--output", slice_path
    ]
    
    try:
        subprocess.run(slice_cmd, check=True)
    except subprocess.CalledProcessError:
        print("Slicing failed. Attempting detection on full genome (less specific).")
        slice_path = args.input
        
    print(f"\n--- Step 2: Running Serogroup Identification ---")
    detect_cmd = [
        "python3", detector_script,
        slice_path,
        "--db", db_path
    ]
    
    try:
        subprocess.run(detect_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Detection failed: {e}")
        sys.exit(1)
        
    if not args.keep_slice and slice_path != args.input:
        if os.path.exists(slice_path):
            os.remove(slice_path)

if __name__ == "__main__":
    main()
