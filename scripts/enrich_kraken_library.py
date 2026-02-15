
import os
import subprocess
from pathlib import Path

VIBRION_ROOT = "/Users/kalinovdameus/Developer/Vibrion"
LIBRARY_PATH = f"{VIBRION_ROOT}/data/kraken2_library"
TAXONOMY_PATH = f"{VIBRION_ROOT}/data/kraken2_haiti_custom/taxonomy"

def run_cmd(cmd):
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return result.stdout.strip()

def main():
    # 1. Decompress existing library genomes
    print("--- Decompressing library genomes ---")
    for root, dirs, files in os.walk(LIBRARY_PATH):
        for file in files:
            if file.endswith(".gz"):
                run_cmd(f"gunzip -f {os.path.join(root, file)}")

    # 2. Fix taxonomy
    print("\n--- Fixing taxonomy ---")
    os.makedirs(TAXONOMY_PATH, exist_ok=True)
    standard_db_path = f"{VIBRION_ROOT}/data/kraken2_standard_8gb"
    for f in ["names.dmp", "nodes.dmp"]:
        src = f"{standard_db_path}/{f}"
        if os.path.exists(src):
            run_cmd(f"cp {src} {TAXONOMY_PATH}/{f}")
            print(f"Copied {f} to {TAXONOMY_PATH}")
        else:
            print(f"Warning: {src} not found")

    # 3. Add Nepal 2010 if we have it in reference_genomes
    nepal_src = f"{VIBRION_ROOT}/data/reference_genomes/Nepal_2010_genomic.fna"
    nepal_dest = f"{LIBRARY_PATH}/haiti_2010_lineage/Nepal_2010.fasta"
    if os.path.exists(nepal_src):
        run_cmd(f"cp {nepal_src} {nepal_dest}")
        print("Added Nepal 2010 to library")
    elif os.path.exists(nepal_src + ".gz"):
        run_cmd(f"gunzip -c {nepal_src}.gz > {nepal_dest}")
        print("Added Nepal 2010 (from gz) to library")

    # 4. List current library
    print("\n--- Current Library Contents ---")
    run_cmd(f"ls -R {LIBRARY_PATH}")

if __name__ == "__main__":
    main()
