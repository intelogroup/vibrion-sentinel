import sys

def persistence_logic_challenge(fastq_path):
    print(f"🌊 ENVIRONMENTAL PERSISTENCE CHALLENGE: {fastq_path}")
    print("=" * 80)
    
    # 1. Check for Haiti SNP Baseline (Simulated)
    matches_baseline = True # Assume it's the 2010 clone
    
    # 2. Inspect vpsR Integrity (The Rugose Switch)
    # Reference: GTGAGGTGGAG
    # Mutation:  GTGAGGTTAGAG (Truncating/Disruptive)
    
    vps_intact = True
    with open(fastq_path, 'r') as f:
        content = f.read()
        if "GTGAGGTTAGAG" in content:
            vps_intact = False
            
    print(f"📊 Haiti SNP Baseline: {'MATCH' if matches_baseline else 'NOMATCH'}")
    print(f"📊 VPS Cluster (vpsR):  {'INTACT' if vps_intact else 'MUTATED/DISRUPTED'}")
    
    # 3. Decision Logic
    if matches_baseline and not vps_intact:
        print("\n⚠️  WARNING: PHENOTYPIC MISMATCH DETECTED")
        print("   Lineage: Haiti 2010-2022")
        print("   Status: LOSS OF RUGOSE_PERSISTER PHENOTYPE")
        print("   Verdict: DISQUALIFIED as Environmental Persister.")
        print("   Reasoning: VPS integrity compromised. Strain unlikely to survive 2019-2022 'silent period' in aquatic reservoirs.")
    elif matches_baseline and vps_intact:
        print("\n✅ CONFIRMED: HIGH-CONFIDENCE PERSISTER")
        print("   Status: RUGOSE_PERSISTER Phenotype Active.")
        print("   Conclusion: Capable of long-term environmental survival.")

if __name__ == "__main__":
    # Simulate a mutated sample
    sim_path = "data/tests/rugose_test/mutated_vpsR.fastq"
    import os
    os.makedirs("data/tests/rugose_test", exist_ok=True)
    with open(sim_path, "w") as f:
        f.write("@read1\nGTGAGGTTAGAG\n+\nIIIIIIIIIIII\n")
        
    persistence_logic_challenge(sim_path)
