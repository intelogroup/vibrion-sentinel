#!/usr/bin/env python3
import os
import sys
import subprocess
import argparse

# Default Validation Candidates
DEFAULT_SAMPLES = {
    "SRR23509888": "Inaba_Switch_2016",
    "SRR23509871": "Haiti_Endemic_2022",
    "SRR22265446": "Stranger_NonO1_Env"
}

OUTPUT_DIR = "data/raw_reads"

def download_and_subsample(srr, name, max_reads=600000):
    print(f"[{name}] Processing {srr}...")
    
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    target = os.path.join(OUTPUT_DIR, f"{srr}.fastq.gz")
    
    if os.path.exists(target):
        print(f"  -> File already exists: {target}")
        return

    # Check if sra-tools is available
    try:
        subprocess.run(["fastq-dump", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("  [ERROR] fastq-dump not found. Please install sra-tools.")
        return

    cmd = [
        "fastq-dump",
        "--split-3", 
        "-X", str(max_reads),
        "--outdir", OUTPUT_DIR,
        "--gzip",
        srr
    ]
    
    print(f"  -> Executing: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
        
        # fastq-dump produces _1.fastq.gz, _2.fastq.gz for paired end
        r1 = os.path.join(OUTPUT_DIR, f"{srr}_1.fastq.gz")
        r2 = os.path.join(OUTPUT_DIR, f"{srr}_2.fastq.gz")
        single = os.path.join(OUTPUT_DIR, f"{srr}.fastq.gz")
        
        if os.path.exists(r1):
            os.rename(r1, target)
            if os.path.exists(r2):
                os.remove(r2) # Keep it simple with R1 for fast validation
            print(f"  -> Success: {target}")
        elif os.path.exists(single):
            print(f"  -> Success: {target}")
        else:
            # Maybe it produced .fastq.gz directly
            pass
            
    except subprocess.CalledProcessError as e:
        print(f"  -> FAILED to download {srr}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Prepare validation data for Vibrion Sentinel")
    parser.add_argument("--accession", help="Specific SRA accession to download")
    parser.add_argument("--name", help="Name for the specific accession")
    parser.add_argument("--subsample", type=int, default=1000000, help="Number of reads to subsample")
    parser.add_argument("--all", action="store_true", help="Download all default validation samples")
    
    args = parser.parse_args()
    
    if args.all:
        for srr, name in DEFAULT_SAMPLES.items():
            download_and_subsample(srr, name, args.subsample)
    elif args.accession:
        name = args.name if args.name else args.accession
        download_and_subsample(args.accession, name, args.subsample)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
