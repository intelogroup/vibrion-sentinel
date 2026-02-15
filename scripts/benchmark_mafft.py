import time
import subprocess
from pathlib import Path

def benchmark_mafft(cmd, label):
    print(f"--- Benchmarking: {label} ---")
    start = time.time()
    try:
        subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception as e:
        print(f"Error: {e}")
        return None
    end = time.time()
    elapsed = end - start
    print(f"Elapsed Time: {elapsed:.2f} seconds\n")
    return elapsed

def main():
    mafft_path = "/Users/kalinovdameus/Developer/Vibrion/.snakemake/conda/80f687140b6ccd3d74604e8c789853c9_/bin/mafft"
    sample_loci = "data/pipeline_output/Bangladesh_O139_Challenge/10_phylogeny/sample_core.fasta"
    ref_loci = "data/core_alignment/reference_core_alignment.fasta"
    
    full_sample = "data/pipeline_output/Bangladesh_O139_Challenge/09_consensus/Bangladesh_O139_Challenge_consensus.fasta"
    full_ref = "data/references/2010EL-1786.fasta"

    # 1. Core-Alignment Mode (Fast)
    core_cmd = f"{mafft_path} --add {sample_loci} --reorder --thread 4 {ref_loci} > /dev/null"
    t1 = benchmark_mafft(core_cmd, "Core-Alignment Mode (Loci only)")

    # 2. Reference-Only Mode (Full Genome - Slow)
    # Merge first as mafft expects a single file
    merged_full = "data/core_alignment/merged_full_benchmark.fasta"
    with open(merged_full, "w") as out:
        for p in [full_sample, full_ref]:
            out.write(open(p).read())
    
    slow_cmd = f"{mafft_path} --auto --thread 4 {merged_full} > /dev/null"
    t2 = benchmark_mafft(slow_cmd, "Reference-Only Mode (Full Genome)")

    if t2 and t1:
        delta = t2 - t1
        print(f"RESULT: Time Delta = {delta:.2f}s ({delta/60:.2f} mins)")
        if t2 > 1800: # 30 mins
             print("CRITICAL: Reference-Only Mode exceeds 30 minute threshold!")
        else:
             print("INFO: Reference-Only Mode is within safety limits for 4Mb genomes.")

if __name__ == "__main__":
    main()
