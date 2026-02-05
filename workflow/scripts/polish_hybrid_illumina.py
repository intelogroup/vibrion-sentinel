#!/usr/bin/env python3
"""
Hybrid Polishing Wrapper (Illumina) - v2.1
Integrates Polypolish (for repeats) + Pilon (for variants)
Author: Vibrion Sentinel Agent
"""

import argparse
import subprocess
import sys
import shutil
from pathlib import Path
import os

def run_command(cmd, desc, env=None):
    print(f"   ⏳ {desc}...")
    try:
        subprocess.run(cmd, check=True, shell=False, env=env, text=True, capture_output=True)
        print(f"   ✅ {desc} Complete")
    except subprocess.CalledProcessError as e:
        print(f"   ❌ {desc} Failed!", file=sys.stderr)
        print(f"   CMD: {' '.join(cmd)}", file=sys.stderr)
        print(f"   STDOUT: {e.stdout}", file=sys.stderr)
        print(f"   STDERR: {e.stderr}", file=sys.stderr)
        raise e

def run_hybrid_polishing(draft_fasta: Path, reads_r1: Path, output_fasta: Path, work_dir: Path, threads: int = 4):
    work_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n🧬 Hybrid Polishing Pipeline v2.1 (Polypolish -> Pilon)")
    print(f"   Draft: {draft_fasta.name}")
    print(f"   R1: {reads_r1.name}")
    
    # Check for Paired Read
    # Safely find R2 in the same directory as R1
    r1_dir = reads_r1.parent
    r1_name = reads_r1.name
    
    # Strategy 1: Replace _1 with _2 in filename only
    r2_name = None
    if "_1_clean" in r1_name:
        # Handle the case where sample name is SRR..._1 and output is SRR..._1_clean
        # and pair is named SRR..._1_2_clean (as per current Snakefile logic)
        r2_name = r1_name.replace("_1_clean", "_1_2_clean")
    elif "_1" in r1_name:
        # rfind/rsplit to only replace the last occurrence in the filename
        parts = r1_name.rsplit("_1", 1)
        if len(parts) == 2:
            r2_name = "_2".join(parts)
    elif "_R1" in r1_name:
        parts = r1_name.rsplit("_R1", 1)
        if len(parts) == 2:
            r2_name = "_R2".join(parts)
            
    if r2_name:
        reads_r2 = r1_dir / r2_name
    else:
        reads_r2 = reads_r1 # Will fail the exists check
    
    is_paired = reads_r2.exists() and reads_r2 != reads_r1
    
    if is_paired:
        print(f"   R2: {reads_r2.name} (Paired-End Found ✅)")
    else:
        # Strategy 2: Look for ANY file ending in _2_clean or _R2_clean in the same dir
        candidates = list(r1_dir.glob("*_2_clean.fastq.gz")) + list(r1_dir.glob("*_R2_clean.fastq.gz"))
        if candidates and candidates[0] != reads_r1:
            reads_r2 = candidates[0]
            is_paired = True
            print(f"   R2: {reads_r2.name} (Paired-End Found via Glob ✅)")

    if not is_paired:
        print(f"   R2: Not Found (searched for {r2_name if r2_name else 'N/A'} in {r1_dir}) ⚠️")

        print("   ⚠️ Polypolish SKIPPED (Requires Paired-End)")
        
    # Variables definition
    current_reference = draft_fasta

    
    # ---------------------------------------------------------
    # ROUND 1: Polypolish (Repeat Fixing) - ONLY IF PAIRED
    # ---------------------------------------------------------
    if is_paired:
        poly_work = work_dir / "polypolish"
        poly_work.mkdir(exist_ok=True)
        
        # 1. Index Draft
        # bwa index might create files next to draft, need to ensure writable or symlink
        # Safest to copy draft to workdir to avoid permission issues/clutter
        local_draft = poly_work / "draft_input.fasta"
        shutil.copy(current_reference, local_draft)
        current_reference = local_draft
        
        run_command(["bwa", "index", str(current_reference)], "Indexing Draft for Polypolish")
        
        # 2. Align R1 and R2 separately
        sam1 = poly_work / "alignments_1.sam"
        sam2 = poly_work / "alignments_2.sam"
        
        print("   ⏳ Aligning R1 for Polypolish...")
        with open(sam1, "w") as f:
            subprocess.run(["bwa", "mem", "-t", str(threads), "-a", str(current_reference), str(reads_r1)], check=True, stdout=f)
        
        print("   ⏳ Aligning R2 for Polypolish...")
        with open(sam2, "w") as f:
            subprocess.run(["bwa", "mem", "-t", str(threads), "-a", str(current_reference), str(reads_r2)], check=True, stdout=f)
            
        # 3. Filter alignments
        filtered1 = poly_work / "filtered_1.sam"
        filtered2 = poly_work / "filtered_2.sam"
        run_command([
            "polypolish", "filter", 
            "--in1", str(sam1), "--in2", str(sam2), 
            "--out1", str(filtered1), "--out2", str(filtered2)
        ], "Filtering alignments")
        
        # 4. Run Polypolish Core
        polypolished_fasta = poly_work / "polypolished.fasta"
        print("   ⏳ Running Polypolish core...")
        with open(polypolished_fasta, "w") as f:
            subprocess.run(
                ["polypolish", "polish", str(current_reference), str(filtered1), str(filtered2)],
                check=True, stdout=f
            )
            
        print("   ✅ Polypolish Complete")
        
        # Verify output exists and is not empty
        if polypolished_fasta.exists() and polypolished_fasta.stat().st_size > 0:
            current_reference = polypolished_fasta
            print("   ✅ Switching Reference -> Polypolished Fasta")
        else:
            print("   ⚠️ Polypolish failed to generate output, falling back to draft.")

    # ---------------------------------------------------------
    # ROUND 2: Pilon (Variant/Indel Fixing)
    # ---------------------------------------------------------
    pilon_work = work_dir / "pilon"
    pilon_work.mkdir(exist_ok=True)
    
    # 1. Index Current Reference (Draft or Polypolished)
    # Again, copy to pilon dir to be safe
    final_input_ref = pilon_work / "pilon_input_ref.fasta"
    shutil.copy(current_reference, final_input_ref)
    current_reference = final_input_ref
    
    run_command(["bwa", "index", str(current_reference)], "Indexing Reference for Pilon")
    
    # 2. Re-Align Reads (Standard Mapping)
    sorted_bam = pilon_work / "aligned.sorted.bam"
    
    # Stream bwa -> samtools sort -> bam to save space/time
    print("   ⏳ Aligning for Pilon (Standard)...")
    bwa_cmd = ["bwa", "mem", "-t", str(threads), str(current_reference), str(reads_r1)]
    if is_paired:
        bwa_cmd.append(str(reads_r2))
        
    sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), "-"]
    
    p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.run(sort_cmd, stdin=p1.stdout, check=True)
    p1.wait()
    
    # 3. Index BAM
    run_command(["samtools", "index", str(sorted_bam)], "Indexing BAM")
    
    # 4. Run Pilon
    print("   ⏳ Running Pilon...")
    java_opts = os.environ.copy()
    java_opts["_JAVA_OPTIONS"] = "-Xmx8G" # Give it memory
    
    pilon_out_prefix = "pilon"
    
    pilon_cmd = [
        "pilon",
        "--genome", str(current_reference),
        "--frags", str(sorted_bam),
        "--output", pilon_out_prefix,
        "--outdir", str(pilon_work),
        "--fix", "snps,indels",
        "--changes",
        "--vcf",
        "--mindepth", "5"
    ]
    
    # Pilon often writes log to stderr
    run_command(pilon_cmd, "Run Pilon Core", env=java_opts)
    
    pilon_fasta = pilon_work / f"{pilon_out_prefix}.fasta"
    
    if not pilon_fasta.exists():
        raise FileNotFoundError(f"Pilon output missing: {pilon_fasta}")
        
    # 5. Clean Headers
    print("   🧹 Cleaning Headers...")
    with open(pilon_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                outfile.write(line.replace('_pilon', ''))
            else:
                outfile.write(line)
                
    print(f"   🎉 Hybrid Polishing Completed -> {output_fasta}")


def main():
    parser = argparse.ArgumentParser(description="Hybrid Polishing (Polypolish+Pilon)")
    parser.add_argument("--draft", required=True, help="Draft Fasta")
    parser.add_argument("--reads", required=True, help="R1 Reads")
    parser.add_argument("--output", required=True, help="Output Fasta")
    parser.add_argument("--outdir", required=True, help="Work Dir")
    parser.add_argument("--threads", type=int, default=4)
    
    args = parser.parse_args()
    
    run_hybrid_polishing(Path(args.draft), Path(args.reads), Path(args.output), Path(args.outdir), args.threads)

if __name__ == "__main__":
    main()
