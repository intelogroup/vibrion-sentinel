#!/usr/bin/env python3
"""
Tier 0 Triage: sourmash K-mer Identity Check

This script uses sourmash to calculate Jaccard similarity between a sample
and known reference strains (2010 sentinel, 2022 outbreak). If the sample
is >99.9% identical to a known strain, it's marked as "ROUTINE" and skips
all AI analysis.

This is the "thermometer check" before ordering the X-ray.
"""

import argparse
import json
import sys
from pathlib import Path

try:
    import sourmash
except ImportError:
    print("ERROR: sourmash not installed. Run: conda install -c bioconda sourmash", file=sys.stderr)
    sys.exit(1)

try:
    import screed
except ImportError:
    print("ERROR: screed not installed. Run: conda install -c bioconda screed", file=sys.stderr)
    sys.exit(1)


def sketch_fasta(fasta_path: Path, ksize: int, scaled: int) -> sourmash.MinHash:
    """Create a MinHash sketch from a FASTA file."""
    mh = sourmash.MinHash(n=0, ksize=ksize, scaled=scaled)
    
    for record in screed.open(str(fasta_path)):
        mh.add_sequence(record.sequence, force=True)
    
    return mh


def calculate_jaccard(sample_mh: sourmash.MinHash, ref_mh: sourmash.MinHash) -> float:
    """Calculate Jaccard similarity between two MinHash sketches."""
    return sample_mh.similarity(ref_mh)


def main():
    parser = argparse.ArgumentParser(description="Tier 0: sourmash K-mer Triage (Cascading Mode)")
    parser.add_argument("--sample-fasta", required=True, help="Path to sample consensus FASTA")
    parser.add_argument("--sentinel-ref", required=True, help="Path to 2010 sentinel reference")
    parser.add_argument("--outbreak-ref", required=True, help="Path to 2022 outbreak reference")
    parser.add_argument("--ksize", type=int, default=31, help="K-mer size")
    parser.add_argument("--scaled", type=int, default=1000, help="Scaled factor")
    parser.add_argument("--threshold", type=float, default=0.999, help="Jaccard threshold for 'routine'")
    parser.add_argument("--min-divergent-kmers", type=int, default=10, help="Min divergent k-mers to trigger Tier 1")
    parser.add_argument("--output", required=True, help="Output JSON path")
    
    args = parser.parse_args()
    
    print("🧬 Tier 0: sourmash K-mer Identity Check (Cascading Mode)")
    print(f"   Sample: {args.sample_fasta}")
    print(f"   K-mer size: {args.ksize}, Scaled: {args.scaled}")
    
    # Create sketches
    print("\n📊 Creating MinHash sketches...")
    sample_mh = sketch_fasta(Path(args.sample_fasta), args.ksize, args.scaled)
    sentinel_mh = sketch_fasta(Path(args.sentinel_ref), args.ksize, args.scaled)
    outbreak_mh = sketch_fasta(Path(args.outbreak_ref), args.ksize, args.scaled)
    
    # Calculate similarities
    sentinel_jaccard = calculate_jaccard(sample_mh, sentinel_mh)
    outbreak_jaccard = calculate_jaccard(sample_mh, outbreak_mh)
    
    print("\n🔍 Jaccard Similarities:")
    print(f"   vs 2010 Sentinel: {sentinel_jaccard:.4f}")
    print(f"   vs 2022 Outbreak: {outbreak_jaccard:.4f}")
    
    # Identify divergent k-mers (k-mers in sample but not in reference)
    closest_ref_mh = sentinel_mh if sentinel_jaccard > outbreak_jaccard else outbreak_mh
    closest_ref = "2010_Sentinel" if sentinel_jaccard > outbreak_jaccard else "2022_Outbreak"
    
    # Get k-mers unique to sample
    sample_hashes = set(sample_mh.hashes)
    ref_hashes = set(closest_ref_mh.hashes)
    divergent_hashes = sample_hashes - ref_hashes
    
    num_divergent = len(divergent_hashes)
    
    print(f"\n🔬 Divergent K-mers: {num_divergent}")
    
    # Determine if routine
    max_similarity = max(sentinel_jaccard, outbreak_jaccard)
    is_routine = max_similarity >= args.threshold and num_divergent < args.min_divergent_kmers
    
    # For cascading: output divergent k-mer positions
    # Note: sourmash MinHash doesn't store positions directly, but we can flag regions
    # In a real implementation, you'd use a different k-mer counter (e.g., KMC, Jellyfish)
    # For now, we'll create synthetic "regions of interest" based on divergence
    
    divergent_regions = []
    if not is_routine and num_divergent > 0:
        # Simulate regions (in production, use actual k-mer positions)
        # For now, flag the whole genome for HyenaDNA to scan
        divergent_regions.append({
            "start": 0,
            "end": 4000000,  # Approximate V. cholerae genome size
            "num_divergent_kmers": num_divergent
        })
    
    result = {
        "tier": 0,
        "method": "sourmash",
        "ksize": args.ksize,
        "scaled": args.scaled,
        "sentinel_jaccard": round(sentinel_jaccard, 4),
        "outbreak_jaccard": round(outbreak_jaccard, 4),
        "max_similarity": round(max_similarity, 4),
        "threshold": args.threshold,
        "num_divergent_kmers": num_divergent,
        "min_divergent_kmers": args.min_divergent_kmers,
        "is_routine": is_routine,
        "closest_reference": closest_ref,
        "divergent_regions": divergent_regions,
        "cascade_to_tier1": not is_routine,
        "decision": "ROUTINE - Skip AI" if is_routine else f"ANOMALY ({num_divergent} divergent k-mers) - Cascade to Tier 1"
    }
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    
    if is_routine:
        print(f"\n✅ ROUTINE: {max_similarity:.2%} match to {closest_ref}")
        print("   → Skipping AI analysis (cost savings!)")
    else:
        print(f"\n⚠️  ANOMALY DETECTED: {num_divergent} divergent k-mers")
        print("   → Cascading to Tier 1 (HyenaDNA will focus on divergent regions)")
    
    print(f"\n📄 Results written to: {args.output}")



if __name__ == "__main__":
    main()
