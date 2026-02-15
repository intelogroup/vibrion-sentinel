#!/usr/bin/env python3
import os
import sys
import subprocess
import json
from Bio import SeqIO
from Bio.Seq import Seq

# Config
REF_FASTA = "data/references/2010EL-1786.fasta"
SAMPLES = ["SRR23509888", "SRR23509871", "SRR22265446"]
SAMTOOLS_BIN = "/Users/kalinovdameus/Developer/Vibrion/.snakemake/conda/ff8213fe4b44fd04c84ecf895c15e02a_/bin/samtools"

# O191 / V. mimicus and NOVC Virulence markers (from public DBs or target seqs)
# Cholix toxin (chxA): DQ840121.1
# T3SS-1 (vopF): CP000626.1 region
VIRulence_targets = {
    "chxA": "ATGAAAAAATATTTTATTTTTGCA", # Cholix toxin start
    "vopF": "ATGTCTAATATTAATTCTTTT",   # T3SS vopF start
    "vopM": "ATGAGTAATATTAATTCTTTT"    # T3SS vopM start
}

def get_base_consensus(bam, chrom, pos):
    """Get consensus base at position from BAM."""
    cmd = f"{SAMTOOLS_BIN} mpileup -r {chrom}:{pos}-{pos} {bam} 2>/dev/null"
    try:
        out = subprocess.check_output(cmd, shell=True, text=True).strip()
        if not out: return "N"
        p = out.split('\t')
        if len(p) < 5: return "N"
        ref_base = p[2].upper()
        b = p[4].upper()
        matches = b.count('.') + b.count(',')
        counts = {'A': b.count('A'), 'C': b.count('C'), 'G': b.count('G'), 'T': b.count('T'), ref_base: matches}
        return max(counts, key=counts.get) if sum(counts.values()) > 0 else "N"
    except: return "N"

def run_blast_search(fastq, query_name, query_seq):
    """Search for a virulence gene in FASTQ using blastn."""
    q_file = f"temp_{query_name}.fasta"
    with open(q_file, "w") as f:
        f.write(f">{query_name}\n{query_seq}\n")
    
    # We blast against the consensus genome if possible, or just the reads
    # For speed, we'll check if a consensus exists
    base = fastq.split('_vibrio')[0]
    consensus = base.replace("03_vibrio", "07_consensus") + "_consensus.fasta"
    
    if os.path.exists(consensus):
        db = consensus
    else:
        # Fallback to a small subset of reads
        return "Consensus missing, cannot BLAST"

    cmd = ["blastn", "-query", q_file, "-subject", db, "-outfmt", "6 pident length qlen"]
    try:
        out = subprocess.check_output(cmd, text=True).strip()
        os.remove(q_file)
        if not out: return "ABSENT"
        return f"PRESENT ({out.splitlines()[0]})"
    except:
        if os.path.exists(q_file): os.remove(q_file)
        return "ERROR"

def main():
    results = {}
    print("Loading Reference...")
    ref = {r.id: r for r in SeqIO.parse(REF_FASTA, "fasta")}
    
    for sample in SAMPLES:
        print(f"Processing {sample}...")
        results[sample] = {}
        bam = f"data/pipeline_output/validation_run/{sample}/04_alignment/{sample}_aligned.sorted.bam"
        fastq = f"data/pipeline_output/validation_run/{sample}/03_vibrio/{sample}_vibrio_complete.fastq.gz"
        
        # 1. wbeT Inaba check (codon 180)
        # Codon start: 2678186 + (180-1)*3 = 2678186 + 537 = 2678723.
        # Genomic bases: 2678723, 2678724, 2678725
        # Index: 2678722, 2678723, 2678724
        codon_ref = str(ref["CP003069.1"].seq[2678722:2678725])
        codon_sample = "".join([get_base_consensus(bam, "CP003069.1", p) for p in [2678723, 2678724, 2678725]])
        results[sample]["wbeT_codon_180"] = {
            "ref": codon_ref,
            "sample": codon_sample,
            "aa": f"{str(Seq(codon_ref).translate())} -> {str(Seq(codon_sample).translate()) if 'N' not in codon_sample else 'N'}"
        }

        # 2. gyrA/parC
        # gyrA 83: 807201 + (83-1)*3 = 807201 + 246 = 807447.
        # Index: 807446, 807447, 807448
        # parC 85: On - strand. Starts at 2075748. 2075748 - (85-1)*3 = 2075748 - 252 = 2075496.
        # Genomic bases for mpileup: 2075496, 2075495, 2075494. (RC them)
        gyrA_ref = str(ref["CP003069.1"].seq[807446:807449])
        gyrA_sample = "".join([get_base_consensus(bam, "CP003069.1", p) for p in [807447, 807448, 807449]])
        
        parC_ref_g = str(ref["CP003069.1"].seq[2075493:2075496])
        parC_ref = str(Seq(parC_ref_g).reverse_complement())
        parC_sample_g = "".join([get_base_consensus(bam, "CP003069.1", p) for p in [2075494, 2075495, 2075496]])
        parC_sample = str(Seq(parC_sample_g).reverse_complement())
        
        results[sample]["resistance"] = {
            "gyrA_83_aa": f"{str(Seq(gyrA_ref).translate())} -> {str(Seq(gyrA_sample).translate()) if 'N' not in gyrA_sample else 'N'}",
            "parC_85_aa": f"{str(Seq(parC_ref).translate())} -> {str(Seq(parC_sample).translate()) if 'N' not in parC_sample else 'N'}"
        }

        # 3. NOVC Virulence (Stranger only)
        if sample == "SRR22265446":
            for gene, seq in VIRulence_targets.items():
                 results[sample][f"virulence_{gene}"] = run_blast_search(fastq, gene, seq)

        # 4. hapR integrity
        # hapR start codon: 1086835. Genomic 1086835-1086837.
        hapR_ref = str(ref["CP003069.1"].seq[1086834:1086837])
        hapR_sample = "".join([get_base_consensus(bam, "CP003069.1", p) for p in [1086835, 1086836, 1086837]])
        results[sample]["hapR_start"] =f"{hapR_ref} -> {hapR_sample}"

    print(json.dumps(results, indent=2))

if __name__ == "__main__":
    main()
