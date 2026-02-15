#!/usr/bin/env python3
"""
Simulate future genomic drift (2027) and novel toxin HGT.
Tests if the Triage system can detect 'unseen' threats.
"""

import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def simulate_drift(input_fasta, drift_rate=0.001):
    """Introduce random SNPs to simulate 1-2 years of future evolution."""
    record = list(SeqIO.parse(input_fasta, "fasta"))[0]
    seq_list = list(str(record.seq))
    bases = ['A', 'C', 'G', 'T']
    
    mutations = 0
    for i in range(len(seq_list)):
        if random.random() < drift_rate:
            original = seq_list[i]
            possible = [b for b in bases if b != original]
            seq_list[i] = random.choice(possible)
            mutations += 1
            
    print(f"🧬 Introduced {mutations} random SNPs (Drift Rate: {drift_rate})")
    record.seq = Seq("".join(seq_list))
    return record

def inject_synthetic_toxin(record, locus_name="ctxB_2027_omega"):
    """Replace a core toxin with a highly mutated/novel variant."""
    # Find ctxB coordinates (approximate for 2010EL ref)
    # CP003069.1:1041239-1041612
    start, end = 1041239, 1041612
    
    # Generate random sequence of same length
    novel_toxin = "".join(random.choice(['A', 'C', 'G', 'T']) for _ in range(end - start))
    
    seq_list = list(str(record.seq))
    seq_list[start:end] = list(novel_toxin)
    
    print(f"⚠️  Injected NOVEL Synthetic Toxin '{locus_name}' at {start}-{end}")
    record.seq = Seq("".join(seq_list))
    return record

if __name__ == "__main__":
    input_ref = "data/references/2010EL-1786.fasta"
    output_dir = Path("data/raw_reads/simulated_2027")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Create Baseline Drift (2027 Routine)
    drifted_record = simulate_drift(input_ref, drift_rate=0.0005)
    SeqIO.write(drifted_record, output_dir / "Haiti_2027_Routine.fasta", "fasta")
    
    # 2. Create Threat Scenario (2027 Super-Strain)
    threat_record = inject_synthetic_toxin(drifted_record)
    SeqIO.write(threat_record, output_dir / "Haiti_2027_SuperStrain.fasta", "fasta")
    
    print(f"\n✅ Simulation complete. Files ready in {output_dir}")
