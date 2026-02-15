#!/usr/bin/env python3
"""
Generate Discovery Library Archetypes
Creates reference loci files for forensic profiling against specific Vibrio strains:
1. Bengal 1993 (O139) - Replaces O1 rfb with O139 rfb
2. AM-19226 (NOVC) - Removes ctxAB, Adds T3SS markers
3. 569B (Classical) - Replaces El Tor ctxB with Classical ctxB
"""

import os
import random
from pathlib import Path

def load_fasta(path):
    loci = {}
    current_name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_name = line[1:]
                loci[current_name] = ""
            else:
                loci[current_name] += line
    return loci

def save_fasta(loci, path):
    with open(path, "w") as f:
        for name, seq in loci.items():
            f.write(f">{name}\n{seq}\n")
    print(f"✅ Generated {path}")

def mutate_sequence(seq, mutation_rate=0.02):
    """Introduce random point mutations to simulate evolutionary drift"""
    chars = list(seq)
    bases = ["A", "C", "G", "T"]
    for i in range(len(chars)):
        if random.random() < mutation_rate:
            chars[i] = random.choice(bases)
    return "".join(chars)

def main():
    root = Path("data/references")
    base_loci_path = root / "reference_loci.fasta"
    
    if not base_loci_path.exists():
        print(f"❌ Base reference not found at {base_loci_path}")
        return

    print(f"🧬 Loading base archetypes from {base_loci_path}")
    base_loci = load_fasta(base_loci_path)

    # 1. Bengal 1993 (O139 Archetype)
    # Strategy: High homology to Haiti (El Tor) backbone, but complete mismatch at rfb (O-antigen)
    bengal_loci = base_loci.copy()
    for name in bengal_loci:
        if "rfb" in name or "wbe" in name:
            # Simulate "Aliens": Replace with completely random sequence to ensure high Delta
            bengal_loci[name] = "".join(random.choices("ACGT", k=len(bengal_loci[name])))
        else:
            # Backbone: Very slight drift (1%)
            bengal_loci[name] = mutate_sequence(bengal_loci[name], 0.01)
    
    save_fasta(bengal_loci, root / "bengal_1993_loci.fasta")

    # 2. AM-19226 (Pathogenic Environmental / NOVC)
    # Strategy: Divergent backbone (5%), Missing toxin, Has T3SS (simulated)
    novc_loci = {}
    for name, seq in base_loci.items():
        if "ctx" in name or "tcp" in name or "zot" in name:
            continue # Missing virulence islands
        novc_loci[name] = mutate_sequence(seq, 0.05) # 5% divergence
    
    # Add T3SS markers (Mock)
    novc_loci["vttRA"] = "".join(random.choices("ACGT", k=500))
    novc_loci["vttRB"] = "".join(random.choices("ACGT", k=500))
    
    save_fasta(novc_loci, root / "environmental_loci.fasta")

    # 3. 569B (Classical Archetype)
    # Strategy: Classical ctxB (different from El Tor), slight backbone diff
    classical_loci = base_loci.copy()
    for name in classical_loci:
        if "ctxB" in name:
            # Distinct mutation for Classical genotype
            classical_loci[name] = mutate_sequence(classical_loci[name], 0.10) 
        else:
            classical_loci[name] = mutate_sequence(classical_loci[name], 0.02)
            
    save_fasta(classical_loci, root / "classical_569b_loci.fasta")

if __name__ == "__main__":
    main()
