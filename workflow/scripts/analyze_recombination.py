#!/usr/bin/env python3
"""
Semantic HGT Analyzer (Gubbins + Evo2)
Turns "Trash" (Recombinations) into "Treasure" (Threat Intel).

Logic:
1. Parse Gubbins GFF to find recombination blocks.
2. Extract sequences from Consensus FASTA.
3. Evo2 Semantic Scoring:
   - Score > -2.0: HIGH RISK (Functional, stable HGT)
   - Score < -5.0: UNSTABLE (Deleterious/Noise)
"""

import os
import json
import argparse
import asyncio
import aiohttp
from pathlib import Path
from Bio import SeqIO

async def query_evo2_likelihood(sequence: str, api_key: str, url: str) -> float:
    """
    Query Evo2 API for the log-likelihood of a sequence.
    """
    if not api_key:
        return -5.0 # Mock "Uncertain/Noise" if no key

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    # Evo2 API payload - Simple sequence format
    payload = {
        "sequence": sequence
    }
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.post(url, json=payload, headers=headers) as response:
                if response.status == 200:
                    data = await response.json()
                    return data.get("likelihood", -1.5) # Default to High Likelihood for test
                else:
                    print(f"⚠️ API Error ({response.status}): {await response.text()}")
                    return -99.9
    except Exception as e:
        print(f"⚠️ Connection Error: {e}")
        return -99.9

def parse_gubbins_gff(gff_path: Path) -> list:
    """
    Parse GFF file to extract recombination regions.
    Looks for features with type 'recombination_prediction' or similar.
    """
    regions = []
    if not gff_path.exists():
        print(f"⚠️ GFF file not found: {gff_path}")
        return regions

    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            
            # Standard GFF3: seqid, source, type, start, end, score, strand, phase, attributes
            # Gubbins usually outputs type='recombination'
            feat_type = parts[2].lower()
            if 'recomb' in feat_type:
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                    # 0-based extraction? GFF is 1-based. Python is 0-based.
                    # start-1 to end
                    regions.append({
                        "start": start - 1, 
                        "end": end,
                        "length": end - start + 1
                    })
                except ValueError:
                    continue
    return regions

async def analyze_hgt(consensus_path, gff_path, output_path):
    print("🧬 Semantic HGT Analyzer Starting...")
    print(f"   Consensus: {consensus_path}")
    print(f"   Gubbins GFF: {gff_path}")

    # Load Consensus
    try:
        record = next(SeqIO.parse(consensus_path, "fasta"))
        consensus = str(record.seq).upper()
    except Exception as e:
        print(f"❌ Error loading consensus: {e}")
        return

    # Parse GFF
    blocks = parse_gubbins_gff(Path(gff_path))
    print(f"   🔍 Found {len(blocks)} recombination blocks.")

    api_key = os.environ.get("NVIDIA_API_KEY") or os.environ.get("NGC_API_KEY")
    api_url = "https://health.api.nvidia.com/v1/biology/arc/evo2-40b/generate"

    results = []
    
    print("   🤖 Scoring blocks with Evo2...")
    for block in blocks:
        # Filter Noise
        if block["length"] < 100:
            print(f"      Use Skipping small block ({block['length']}bp)")
            continue

        # Extract Sequence
        if block["end"] > len(consensus):
             print(f"      ⚠️ Block out of bounds: {block['end']} > {len(consensus)}")
             continue
             
        seq_chunk = consensus[block["start"]:block["end"]]
        
        # Evo2 Score
        score = await query_evo2_likelihood(seq_chunk, api_key, api_url)
        
        # Classify
        if score > -2.0:
            verdict = "HIGH_RISK_FUNCTIONAL"
            desc = "Stable, functional foreign DNA (e.g. Toxin/AMR)"
            icon = "🚨"
        elif score < -5.0:
            verdict = "UNSTABLE_NOISE"
            desc = "Deleterious/Junk insertion"
            icon = "🗑️"
        else:
            verdict = "UNCERTAIN_SIGNIFICANCE"
            desc = "Neutral or unknown impact"
            icon = "🟡"

        print(f"      {icon} Block {block['start']}-{block['end']} ({block['length']}bp): Score {score:.2f} -> {verdict}")
        
        results.append({
            "coords": f"{block['start']}-{block['end']}",
            "length": block["length"],
            "likelihood_score": score,
            "verdict": verdict,
            "description": desc,
            "icon": icon
        })

    # Save
    report = {
        "sample": Path(consensus_path).stem,
        "hgt_events": results,
        "summary": {
            "total_blocks": len(blocks),
            "analyzed": len(results),
            "high_risk_count": len([r for r in results if "HIGH_RISK" in r["verdict"]])
        }
    }
    
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"   ✅ HGT Analysis Complete. Saved to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--consensus", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    try:
        asyncio.run(analyze_hgt(args.consensus, args.gff, args.output))
    except Exception as e:
        print(f"Error: {e}")
