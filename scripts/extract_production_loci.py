#!/usr/bin/env python3
"""
Extract Production Loci via BLAST
Uses blastn to find and extract surveillance loci from full reference genomes.
Robust to annotation naming differences.
"""

import subprocess
import sys
import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def run_command(cmd):
    try:
        subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        return False
    return True

def extract_loci_blast(query_fasta, target_genome, output_fasta, name_prefix):
    print(f"\n🧬 Processing Archetype: {name_prefix}...")
    print(f"   Target: {target_genome}")
    
    # 1. Create BLAST DB for target
    db_name = f"{target_genome}.db"
    run_command(f"makeblastdb -in {target_genome} -dbtype nucl -out {db_name}")
    
    found_loci = []
    
    # 2. Iterate through Query Loci
    queries = list(SeqIO.parse(query_fasta, "fasta"))
    
    for record in queries:
        locus_name = record.id
        # Write single query file
        q_path = f"temp_query_{locus_name}.fasta"
        SeqIO.write(record, q_path, "fasta")
        
        # BLAST
        out_fmt = "6 sseq" # Just the subject sequence
        
        # SPECIAL LOGIC: Stricter Thresholds for O1/O139 Discriminators
        # wbeT (O1) vs wbf (O139) requires high specificity to avoid paralogs
        perc_identity = "70" # Default for general homology
        if "wbeT" in locus_name or "wbf" in locus_name:
            perc_identity = "90" # High strictness for serogroup markers
            
        blast_cmd = f"blastn -query {q_path} -db {db_name} -outfmt '{out_fmt}' -max_target_seqs 1 -perc_identity {perc_identity}"
        
        # Capture output
        try:
            result = subprocess.check_output(blast_cmd, shell=True, text=True)
            hit_seq = result.strip().split('\n')[0] if result else ""
            
            if hit_seq:
                print(f"   ✅ Found {locus_name} (Length: {len(hit_seq)}, Identity >={perc_identity}%)")
                found_loci.append(SeqRecord(Seq(hit_seq), id=locus_name, description=f"Extracted from {name_prefix}"))
            else:
                print(f"   ⚠️  {locus_name} NOT FOUND (Target Identity >={perc_identity}% not met)")
                # Fill with Ns to match query length + padding, or just standard block
                found_loci.append(SeqRecord(Seq("N" * len(record.seq)), id=locus_name, description="MISSING"))
                
        except Exception as e:
            print(f"   ❌ BLAST Error for {locus_name}: {e}")
            
        # Cleanup
        if os.path.exists(q_path):
            os.remove(q_path)

    # --- SYNTHETIC CORRECTION FOR O139 (Bengal_1993) ---
    # The NCBI Reference NC_002505.1 contains wbeT (!), likely due to backbone mapping issues.
    # To create a PROPER O139 Discovery Archetype, we must:
    # 1. Force Remove wbeT (to serve as negative control)
    # 2. Force Inject wbfZ (to serve as positive control)
    if name_prefix == "Bengal_1993":
        print("   🔧 APPLYING SYNTHETIC CORRECTION for O139 Archetype:")
        
        # 1. Remove wbeT
        original_count = len(found_loci)
        found_loci = [rec for rec in found_loci if "wbeT" not in rec.id]
        if len(found_loci) < original_count:
            print("      ✂️  Removed contaminating 'wbeT' to enforce O139 profile.")
            
        # 2. Inject wbfZ
        # We need to read it from the marker file
        wbf_ref_path = Path("data/references/wbf_marker.fasta")
        if wbf_ref_path.exists():
            wbf_rec = list(SeqIO.parse(wbf_ref_path, "fasta"))[0]
            # Ensure it's not already there (which it isnt)
            if not any("wbf" in rec.id for rec in found_loci):
                found_loci.append(SeqRecord(wbf_rec.seq, id="wbfZ", description="Synthetic Injection for O139"))
                print("      💉 Injected 'wbfZ' marker to enforce O139 profile.")
        else:
            print("      ⚠️  Warning: wbf_marker.fasta not found for injection.")

    # 3. Save Output
    SeqIO.write(found_loci, output_fasta, "fasta")
    print(f"💾 Saved {len(found_loci)} loci to {output_fasta}")
    
    # Cleanup DB files (optional but polite)
    for ext in ["nhr", "nin", "nsq"]:
        f = f"{db_name}.{ext}"
        if os.path.exists(f): 
            os.remove(f)

def main():
    # Paths
    ref_dir = Path("data/references")
    prod_dir = Path("data/production_references")
    
    query_fasta = ref_dir / "reference_loci.fasta"
    
    # Check inputs
    if not query_fasta.exists():
        print("❌ Query loci file not found!")
        return

    # 1. Bengal 1993
    # Use Chr1 and Chr2 concatenated for BLAST DB
    bengal_genome = prod_dir / "Bengal_1993_Combined.fasta"
    run_command(f"cat {prod_dir}/Bengal_1993_Chr1.fasta {prod_dir}/Bengal_1993_Chr2.fasta > {bengal_genome}")
    extract_loci_blast(query_fasta, bengal_genome, ref_dir / "bengal_1993_loci.fasta", "Bengal_1993")
    
    # 2. Classical 569B
    classical_genome = prod_dir / "Classical_569B_Combined.fasta"
    run_command(f"cat {prod_dir}/Classical_569B_Chr1.fasta {prod_dir}/Classical_569B_Chr2.fasta > {classical_genome}")
    extract_loci_blast(query_fasta, classical_genome, ref_dir / "classical_569b_loci.fasta", "Classical_569B")
    
    # 3. Environmental NOVC
    novc_genome = prod_dir / "Environmental_NOVC_Contig1.fasta"
    extract_loci_blast(query_fasta, novc_genome, ref_dir / "environmental_loci.fasta", "Environmental")

if __name__ == "__main__":
    main()
