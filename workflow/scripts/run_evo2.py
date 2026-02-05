#!/usr/bin/env python3
import json
import asyncio
import aiohttp
import random
import sys
import argparse
from typing import Dict
from Bio import SeqIO
import numpy as np

# Increase recursion depth for deep nesting
sys.setrecursionlimit(2000)

def log(msg: str):
    """Log to stderr so it shows up in Snakemake logs."""
    sys.stderr.write(f"{msg}\n")
    sys.stderr.flush()

async def get_evo2_metrics(
    session: aiohttp.ClientSession, 
    api_url: str, 
    api_key: str, 
    sequence: str
) -> Dict:
    """
    Call the NVIDIA BioNeMo Evo2 API with the verified 'arc/evo2-40b/generate' contract.
    """
    if not api_key:
        return {"success": False, "error": "No API key provided"}

    # Clean sequence: Upper case and ACGT only
    clean_seq = "".join(c for c in sequence.upper() if c in "ACGT")
    
    if len(clean_seq) < 10:
        return {"success": False, "error": "Sequence too short (likely unmapped/Ns)"}
    
    # NVIDIA NIM arc/evo2-40b/generate contract
    analysis_seq = clean_seq[:1000] 
    payload = {
        "sequence": analysis_seq,
        "num_tokens": 1,
        "temperature": 0.1,
        "top_k": 1,
        "enable_sampled_probs": True
    }

    headers = {
        "Authorization": f"Bearer {api_key}",
        "Content-Type": "application/json",
        "Accept": "application/json"
    }

    max_retries = 5
    for attempt in range(max_retries):
        try:
            if attempt > 0:
                await asyncio.sleep(5 * (2 ** (attempt - 1)) + random.uniform(0.1, 1.0))

            async with session.post(api_url, headers=headers, json=payload, timeout=120) as response:
                if response.status == 200:
                    result = await response.json()
                    probs = result.get("sampled_probs", [])
                    if probs:
                        confidence = sum(probs) / len(probs)
                        return {"success": True, "confidence": confidence}
                    return {"success": True, "confidence": 0.5, "note": "No probs returned"}
                
                elif response.status == 429:
                    log(f"      ⚠️  Rate Limited (429). Attempt {attempt+1}/{max_retries}")
                    continue
                else:
                    error_val = await response.text()
                    return {"success": False, "error": f"HTTP {response.status}: {error_val[:50]}"}
        except Exception as e:
            if attempt == max_retries - 1:
                return {"success": False, "error": str(e)}
    return {"success": False, "error": "Max retries exceeded"}

async def main():
    parser = argparse.ArgumentParser(description="Evo2 Locus Likelihood Analysis (Temporal Surveillance)")
    parser.add_argument("--sample-fasta", required=True, help="Path to sample loci FASTA")
    parser.add_argument("--sentinel-fasta", required=True, help="Path to 2010 sentinel loci FASTA")
    parser.add_argument("--outbreak-fasta", required=True, help="Path to 2022 outbreak loci FASTA")
    parser.add_argument("--output", required=True, help="Path to output JSON")
    parser.add_argument("--api-url", required=True, help="Evo2 API Endpoint")
    parser.add_argument("--api-key", required=True, help="Evo2 API Key")
    args = parser.parse_args()

    sample_loci = list(SeqIO.parse(args.sample_fasta, "fasta"))
    sentinel_loci = {s.id: str(s.seq) for s in SeqIO.parse(args.sentinel_fasta, "fasta")}
    outbreak_loci = {s.id: str(s.seq) for s in SeqIO.parse(args.outbreak_fasta, "fasta")}
    
    loci_results = []
    sem = asyncio.Semaphore(1) # Strict serial processing for trial keys

    async with aiohttp.ClientSession() as session:
        for sample in sample_loci:
            locus_id = sample.id
            log(f"🔍 Analyzing Locus: {locus_id}")
            
            # 1. Get Sample Confidence
            async with sem:
                sample_res = await get_evo2_metrics(session, args.api_url, args.api_key, str(sample.seq))
            
            # 2. Get Sentinel (2010) Baseline
            sentinel_seq = sentinel_loci.get(locus_id)
            sentinel_res = {"success": False}
            if sentinel_seq:
                async with sem:
                    sentinel_res = await get_evo2_metrics(session, args.api_url, args.api_key, sentinel_seq)

            # 3. Get Outbreak (2022) Baseline
            outbreak_seq = outbreak_loci.get(locus_id)
            outbreak_res = {"success": False}
            if outbreak_seq:
                async with sem:
                    outbreak_res = await get_evo2_metrics(session, args.api_url, args.api_key, outbreak_seq)

            # --- Calculate Dual Deltas ---
            sample_conf = sample_res.get("confidence", 0.0) if sample_res["success"] else 0.0
            sentinel_conf = sentinel_res.get("confidence", 0.0) if sentinel_res["success"] else 0.0
            outbreak_conf = outbreak_res.get("confidence", 0.0) if outbreak_res["success"] else 0.0
            
            sentinel_delta = abs(sample_conf - sentinel_conf)
            outbreak_delta = abs(sample_conf - outbreak_conf)

            # Threat Detection: Any shift > 0.20 from Sentinel OR > 0.10 from Outbreak
            threat_level = "low"
            if sentinel_delta > 0.20 or outbreak_delta > 0.10:
                threat_level = "high"
            elif sentinel_delta > 0.10 or outbreak_delta > 0.05:
                threat_level = "medium"

            locus_result = {
                "locus": locus_id,
                "api_status": "success" if sample_res["success"] else "error",
                "sample_confidence": round(sample_conf, 5),
                "sentinel_confidence": round(sentinel_conf, 5),
                "outbreak_confidence": round(outbreak_conf, 5),
                "delta_anomaly": round(sentinel_delta, 5), # Primary for report sorting
                "sentinel_delta": round(sentinel_delta, 5),
                "outbreak_delta": round(outbreak_delta, 5),
                "threat_level": threat_level,
                "interpretation": "Normal" if sentinel_delta < 0.05 else "Drifted",
                "reference_archetype": "Haiti_2010_Ancestor" # Default for sorting
            }
            loci_results.append(locus_result)
            log(f"   ✅ Result: {threat_level.upper()} (SentDelta: {sentinel_delta:.2f}, OutDelta: {outbreak_delta:.2f})")

    # Calculate Aggregate Profiles for Report Generator
    sentinel_deltas = [r["sentinel_delta"] for r in loci_results if r["api_status"] == "success"]
    outbreak_deltas = [r["outbreak_delta"] for r in loci_results if r["api_status"] == "success"]
    
    avg_sentinel = np.mean(sentinel_deltas) if sentinel_deltas else 0
    avg_outbreak = np.mean(outbreak_deltas) if outbreak_deltas else 0
    
    best_match = "Haiti_2010_Ancestor" if avg_sentinel <= avg_outbreak else "Haiti_2022_Resurgence"
    
    output_data = {
        "loci_analysis": loci_results,
        "best_archetype_match": best_match,
        "archetype_profiles": {
            "Haiti_2010_Ancestor": {
                "average_delta_anomaly": round(avg_sentinel, 5),
                "threat_level": "LOW" if avg_sentinel < 0.1 else "MEDIUM"
            },
            "Haiti_2022_Resurgence": {
                "average_delta_anomaly": round(avg_outbreak, 5),
                "threat_level": "LOW" if avg_outbreak < 0.1 else "MEDIUM"
            }
        },
        "summary": {
            "total_loci": len(loci_results),
            "success_rate": len(sentinel_deltas) / len(loci_results) if loci_results else 0,
            "threat_level": "LOW" if min(avg_sentinel, avg_outbreak) < 0.1 else "MEDIUM"
        }
    }

    with open(args.output, "w") as f:
        json.dump(output_data, f, indent=2)
    log(f"🎉 Analysis saved to {args.output}")

if __name__ == "__main__":
    asyncio.run(main())