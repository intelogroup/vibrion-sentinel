#!/usr/bin/env python3
"""
Create a Mixed Ogawa/Inaba Sample for Testing
This script merges reads from known Ogawa and Inaba samples to test mixed serotype detection
"""
import subprocess
import sys
from pathlib import Path

def create_mixed_sample():
    """
    Mix Haiti 2022 Ogawa (SRR22265437) with Inaba (SRR8364255) to simulate co-infection
    """
    print("🧬 Creating Mixed Ogawa/Inaba Test Sample")
    print("=" * 70)
    
    # Known samples
    ogawa_r1 = Path("data/raw_reads/SRR22265437_1.fastq.gz")
    ogawa_r2 = Path("data/raw_reads/SRR22265437_2.fastq.gz")
    inaba_r1 = Path("data/raw_reads/PRJNA510624/SRR8364255_1.fastq.gz")
    inaba_r2 = Path("data/raw_reads/PRJNA510624/SRR8364255_2.fastq.gz")
    
    # Output files
    output_r1 = Path("data/raw_reads/mixed_ogawa_inaba_1.fastq.gz")
    output_r2 = Path("data/raw_reads/mixed_ogawa_inaba_2.fastq.gz")
    
    # Check input files exist
    for f in [ogawa_r1, ogawa_r2, inaba_r1, inaba_r2]:
        if not f.exists():
            print(f"❌ ERROR: Input file not found: {f}")
            sys.exit(1)
    
    print(f"✓ Ogawa sample (Haiti 2022): SRR22265437")
    print(f"✓ Inaba sample (Inaba founder): SRR8364255")
    print()
    
    # Strategy: Take first 200k lines (50k reads) from each to simulate equal co-infection
    # 4 lines per read in FASTQ format
    
    print("📊 Step 1: Subsampling Ogawa reads (50k reads = 200k lines)...")
    subprocess.run(
        f"gunzip -c {ogawa_r1} | head -n 200000 | gzip > /tmp/ogawa_sub_1.fastq.gz",
        shell=True, check=True
    )
    subprocess.run(
        f"gunzip -c {ogawa_r2} | head -n 200000 | gzip > /tmp/ogawa_sub_2.fastq.gz",
        shell=True, check=True
    )
    
    print("📊 Step 2: Subsampling Inaba reads (50k reads = 200k lines)...")
    subprocess.run(
        f"gunzip -c {inaba_r1} | head -n 200000 | gzip > /tmp/inaba_sub_1.fastq.gz",
        shell=True, check=True
    )
    subprocess.run(
        f"gunzip -c {inaba_r2} | head -n 200000 | gzip > /tmp/inaba_sub_2.fastq.gz",
        shell=True, check=True
    )
    
    print("🔀 Step 3: Merging reads...")
    subprocess.run(
        f"cat /tmp/ogawa_sub_1.fastq.gz /tmp/inaba_sub_1.fastq.gz > {output_r1}",
        shell=True, check=True
    )
    subprocess.run(
        f"cat /tmp/ogawa_sub_2.fastq.gz /tmp/inaba_sub_2.fastq.gz > {output_r2}",
        shell=True, check=True
    )
    
    print("🧹 Step 4: Cleanup...")
    subprocess.run("rm /tmp/*_sub_*.fastq.gz", shell=True)
    
    print()
    print("✅ SUCCESS: Mixed Ogawa/Inaba sample created!")
    print(f"   R1: {output_r1}")
    print(f"   R2: {output_r2}")
    print()
    print("📋 Next: Run pipeline with this sample to test mixed serotype detection")
    print("   snakemake --configfile workflow/mixed_test_config.yaml --cores 4")

if __name__ == "__main__":
    create_mixed_sample()
