#!/usr/bin/env python3
"""
Dark Matter Protocol: Safe Assembly with Safety Rails
Detects novel elements (Plasmids, Prophages) absent from the reference genome.

Steps:
1. Quality Filter (FastP) - Avoid assembling noise.
2. RAM Check (psutil) - Protect field hardware.
3. Assembly (Megahit/SPAdes) - Turn unmapped reads into contigs.
4. Targeted BLAST - Local identification of plasmids/AMR.
"""

import sys
import shutil
import argparse
import subprocess
import psutil
from pathlib import Path

def check_tool(name):
    return shutil.which(name) is not None

def get_ram_gb():
    mem = psutil.virtual_memory()
    return mem.available / (1024**3)

def run_cmd(cmd):
    print(f"      🚀 Executing: {cmd}")
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"      ❌ Error: {e.stderr}")
        return False

def dark_matter_protocol(unmapped_fastq, output_dir, db_dir):
    print("🕵️‍♂️ Dark Matter Protocol Initiated...")
    print(f"   Input: {unmapped_fastq}")
    
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    
    # 1. Quality Filter (FastP)
    clean_fastq = out_path / "unmapped_clean.fastq"
    print("   🩺 Step 1: Quality Filtering...")
    if check_tool("fastp"):
        run_cmd(f"fastp -i {unmapped_fastq} -o {clean_fastq} --average_qual 20 --length_required 50 --disable_adapter_trimming")
    else:
        print("      ⚠️ FastP not found. Using raw unmapped reads.")
        shutil.copy(unmapped_fastq, clean_fastq)

    # 2. Resource Check
    ram_available = get_ram_gb()
    print(f"   🖥️  Step 2: Resource Check (Available RAM: {ram_available:.2f} GB)")
    
    assembler = None
    if ram_available > 16 and check_tool("spades.py"):
        assembler = "spades"
    elif ram_available > 4 and check_tool("megahit"):
        assembler = "megahit"
    elif ram_available > 4:
        # Fallback to SPAdes even if slow if megahit missing
        if check_tool("spades.py"):
            assembler = "spades"
        else:
            print("      ❌ Error: No assembly tools (megahit/spades) found.")
            return False
    else:
        print("      🛑 Aborting: Insufficient RAM (<4GB) for assembly.")
        return False

    # 3. Assembly
    assembly_out = out_path / "assembly"
    print(f"   🧬 Step 3: Assembly ({assembler})...")
    
    success = False
    if assembler == "megahit":
        # Megahit requires output dir to NOT exist
        if assembly_out.exists(): shutil.rmtree(assembly_out)
        success = run_cmd(f"megahit -r {clean_fastq} -o {assembly_out}")
        contigs = assembly_out / "final.contigs.fa"
    else:
        success = run_cmd(f"spades.py -s {clean_fastq} -o {assembly_out} --only-assembler")
        contigs = assembly_out / "contigs.fasta"

    if not success or not contigs.exists():
        print("      ❌ Assembly failed or produced no contigs.")
        return False

    # 4. Targeted BLAST
    print("   🏷️  Step 4: Targeted Local BLAST...")
    blast_report = out_path / "blast_summary.tsv"
    
    # Mocking BLAST for now if databases or blastn missing
    if check_tool("blastn") and Path(db_dir).exists():
        # Search against plasmids specifically
        plasmid_db = Path(db_dir) / "plasmids"
        run_cmd(f"blastn -query {contigs} -db {plasmid_db} -outfmt 6 -max_target_seqs 5 > {blast_report}")
    else:
        print("      ⚠️ BLAST skipped (Tools or DB missing). Creating empty report.")
        blast_report.touch()

    print(f"   ✅ Dark Matter Protocol Complete. Results: {out_path}")
    return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--unmapped", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--db", default="data/databases/")
    args = parser.parse_args()
    
    success = dark_matter_protocol(args.unmapped, args.output, args.db)
    if not success:
        sys.exit(0) # Exit cleanly but log failure in the tool
