import json
import sys
from pathlib import Path

def load_sim(path):
    with open(path) as f:
        data = json.load(f)["loci_analysis"]
    return {locus: data[locus]["cosine_similarity"] for locus in data}

def main():
    sim_2010_path = "results/hyena_test/haiti_2026_sim.json" # 2026 vs 2010
    sim_2022_path = "results/hyena_test/haiti_2026_vs_2022_sim.json" # 2026 vs 2022
    
    try:
        sim_2010 = load_sim(sim_2010_path)
        sim_2022 = load_sim(sim_2022_path)
    except FileNotFoundError as e:
        print(f"Error loading files: {e}")
        return

    print(f"{'Locus':<15} | {'vs 2010 (Ref)':<15} | {'vs 2022 (Resurg)':<15} | {'Closer To?'}")
    print("-" * 65)
    
    # Track "winning" lineage
    wins_2010 = 0
    wins_2022 = 0
    
    for locus in sim_2010:
        if locus not in sim_2022: continue
        
        s10 = sim_2010[locus]
        s22 = sim_2022[locus]
        
        diff = s22 - s10
        winner = "2022" if diff > 0.0001 else ("2010" if diff < -0.0001 else "Equal")
        
        if winner == "2022": wins_2022 += 1
        if winner == "2010": wins_2010 += 1
        
        # Highlight significant shifts
        hl = ""
        if abs(diff) > 0.01: hl = "**"
        
        print(f"{locus:<15} | {s10:<15.4f} | {s22:<15.4f} | {winner} {hl}")

    print("-" * 65)
    print(f"Summary: 2026 sample is structurally closer to:")
    print(f"  2010 Sentinel: {wins_2010} loci")
    print(f"  2022 Resurgence: {wins_2022} loci")

if __name__ == "__main__":
    main()
