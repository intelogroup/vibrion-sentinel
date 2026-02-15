#!/usr/bin/env python3
"""
Fetch Real Archetype Loci from NCBI
Production script to replace simulated reference data with real genomic sequences.
Uses Entrez API to fetch specific loci coordinates for:
1. Bengal 1993 (O139) - Strain MO10
2. Classical 569B - Strain 569B
3. Envrionmental - Strain AM-19226
"""

from Bio import Entrez, SeqIO
from pathlib import Path
import os
import sys
import time

# Set your email for Entrez
Entrez.email = "sentinel_pipeline@vibrion.org"

# Target Archetypes (Genbank Accessions)
ARCHETYPES = {
    "Bengal_1993": {
        "accession": "NC_002505.1", # Chr1 (MO10)
        "accession_chr2": "NC_002506.1", # Chr2
        "desc": "Vibrio cholerae O139 strain MO10"
    },
    "Classical_569B": {
        "accession": "NC_012668.1", # Chr1
        "accession_chr2": "NC_012667.1", # Chr2
        "desc": "Vibrio cholerae O1 biovar Classical strain 569B"
    },
    "Environmental_NOVC": {
        "accession": "NZ_AAWD01000001.1", # Contig from AM-19226 WGS
        "desc": "Vibrio cholerae strain AM-19226 (Non-O1/Non-O139)"
    }
}

# Loci Map (Approximation of Haiti coordinates mapped to these genomes)
# In a real production system, we would BLAST the Haiti locus to find the exact coordinates.
# For this script, we will fetch the FULL GENOME and then Extract matching genes by Name/BLAST 
# to ensure we get the ortholog, not just the coordinate.
# BUT, downloading full genomes takes time. 
# Optimized Strategy: Use efetch with seq_start/stop? No, we don't know coords.
# Strategy: Download Feature Table, find gene names, get coords, then fetch.

def fetch_locus_sequence(accession, gene_name, flank=200):
    """
    Search for a gene in a genbank record and extract its sequence + flanking regions.
    """
    try:
        print(f"  🔍 Searching for gene '{gene_name}' in {accession}...")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        
        found_seq = None
        for feature in record.features:
            if feature.type == "gene" or feature.type == "CDS":
                q_gene = feature.qualifiers.get("gene", [""])[0]
                q_locus = feature.qualifiers.get("locus_tag", [""])[0]
                
                if gene_name.lower() in q_gene.lower() or gene_name.lower() in q_locus.lower():
                    # biomol extraction
                    start = max(0, int(feature.location.start) - flank)
                    end = min(len(record.seq), int(feature.location.end) + flank)
                    found_seq = record.seq[start:end]
                    print(f"    ✅ Found {gene_name} ({len(found_seq)}bp)")
                    return str(found_seq)
        
        print(f"    ⚠️  Gene '{gene_name}' not found in {accession}")
        return None
        
    except Exception as e:
        print(f"    ❌ Error fetching {gene_name}: {e}")
        return None

def fetch_archetype_dataset(name, info, locus_list):
    """Generate the multi-fasta for a specific archetype"""
    print(f"\n🧬 Building Archetype: {name} ({info['desc']})")
    
    loci_seqs = {}
    
    # Locus-to-Gene Mapping (Robust Aliases)
    # Maps internal locus names to list of potential NCBI gene/CDS matches
    gene_map = {
        "wbeT": ["wbeT", "rfbT", "manC", "VC0240"], # O139 antigen markers
        "ctxB": ["ctxB", "VC1456", "VC_RS07050", "cholera enterotoxin subunit B"],
        "tcpA": ["tcpA", "VC0828", "VC_RS04060", "toxin-coregulated pilus major subunit"],
        "gyrA": ["gyrA", "VC1258", "VC_RS06085", "DNA gyrase subunit A"],
        "parE": ["parE", "VC2434", "VC_RS11680", "DNA topoisomerase IV subunit B"],
        "rstR": ["rstR", "VC1455", "VC_RS07045"],
        "rtxA": ["rtxA", "VC1451", "VC_RS07025", "repeats-in-toxin"]
    }
    
    for locus in locus_list:
        targets = gene_map.get(locus, [locus])
        found = False
        
        for gene_target in targets:
            print(f"    Trying alias: {gene_target}")
            seq = fetch_locus_sequence(info["accession"], gene_target)
            
            # Try Chr2 if not found
            if not seq and "accession_chr2" in info:
                 seq = fetch_locus_sequence(info["accession_chr2"], gene_target)
            
            if seq:
                loci_seqs[locus] = seq
                found = True
                break
             
        if not found:
            # Fallback: Create placeholder N-string to avoid crashing 40B
            print(f"    ⚠️  Creating N-gap for missing locus {locus} (all aliases failed)")
            loci_seqs[locus] = "N" * 500
            
        time.sleep(0.5) # Courtesy limit
        
    return loci_seqs

def save_fasta(loci, path):
    with open(path, "w") as f:
        for name, seq in loci.items():
            f.write(f">{name}\n{seq}\n")
    print(f"💾 Saved to {path}")

def main():
    root = Path("data/references")
    # Identify loci to fetch from the existing simulated files or hardcoded list
    target_loci = ["wbeT", "ctxB", "tcpA", "gyrA", "parE", "rstR", "rtxA"]
    
    # 1. Bengal 1993
    bengal = fetch_archetype_dataset("Bengal_1993", ARCHETYPES["Bengal_1993"], target_loci)
    save_fasta(bengal, root / "bengal_1993_loci.fasta")
    
    # 2. Classical 569B
    classical = fetch_archetype_dataset("Classical_569B", ARCHETYPES["Classical_569B"], target_loci)
    save_fasta(classical, root / "classical_569b_loci.fasta")
    
    # 3. Environmental (AM-19226)
    novc = fetch_archetype_dataset("Environmental_NOVC", ARCHETYPES["Environmental_NOVC"], target_loci)
    save_fasta(novc, root / "environmental_loci.fasta")

if __name__ == "__main__":
    main()
