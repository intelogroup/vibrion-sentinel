#!/usr/bin/env python3
"""
Platform Detection for Consensus Polishing
Auto-detects Nanopore vs Illumina vs Hybrid sequencing data
"""

import gzip
from pathlib import Path
from typing import Tuple, Dict
import json


def calculate_read_stats(fastq_path: Path, sample_size: int = 10000) -> Dict:
    """
    Calculate read length statistics from FASTQ file.
    Samples first N reads for speed.
    """
    lengths = []
    quality_scores = []
    
    opener = gzip.open if str(fastq_path).endswith('.gz') else open
    
    with opener(fastq_path, 'rt') as f:
        line_count = 0
        for line in f:
            line_count += 1
            
            # FASTQ format: @header, seq, +, qual (4 lines per read)
            if line_count % 4 == 2:  # Sequence line
                lengths.append(len(line.strip()))
            elif line_count % 4 == 0:  # Quality line
                # Average quality score (Phred)
                qual_line = line.strip()
                avg_qual = sum(ord(c) - 33 for c in qual_line) / len(qual_line)
                quality_scores.append(avg_qual)
            
            # Stop after sampling
            if len(lengths) >= sample_size:
                break
    
    if not lengths:
        raise ValueError(f"No reads found in {fastq_path}")
    
    avg_length = sum(lengths) / len(lengths)
    avg_quality = sum(quality_scores) / len(quality_scores)
    max_length = max(lengths)
    
    return {
        "avg_length": avg_length,
        "avg_quality": avg_quality,
        "max_length": max_length,
        "reads_sampled": len(lengths)
    }


def detect_platform(fastq_path: Path) -> Tuple[str, Dict]:
    """
    Detect sequencing platform from read characteristics.
    
    Returns:
        platform: "NANOPORE", "ILLUMINA", or "HYBRID"
        stats: Read statistics used for detection
    """
    stats = calculate_read_stats(fastq_path)
    
    avg_length = stats["avg_length"]
    avg_quality = stats["avg_quality"]
    
    # Decision logic
    if avg_length > 500:
        # Long reads = Nanopore or PacBio
        if avg_quality < 15:
            platform = "NANOPORE"  # Lower quality, characteristic of ONT
        else:
            platform = "PACBIO"    # Higher quality HiFi reads
    elif avg_length < 350:
        platform = "ILLUMINA"      # Short, high-quality reads
    else:
        platform = "HYBRID"        # Ambiguous or mixed
    
    return platform, stats


def select_polisher(platform: str, nanopore_model: str = None) -> Dict:
    """
    Select appropriate polishing tool based on platform.
    
    Args:
        platform: Detected platform
        nanopore_model: Specific Medaka model (if known)
    
    Returns:
        Polishing configuration
    """
    if platform == "NANOPORE":
        # Default to R10.4.1 HAC model (most common in 2024+)
        model = nanopore_model or "r1041_e82_400bps_hac_v4.2.0"
        return {
            "tool": "medaka",
            "rounds": 2,
            "model": model,
            "conda_env": "polishing_nanopore"
        }
    
    elif platform == "ILLUMINA":
        return {
            "tool": "pilon",
            "rounds": 1,
            "model": None,
            "conda_env": "polishing_illumina"
        }
    
    elif platform == "PACBIO":
        return {
            "tool": "arrow",  # PacBio's native polisher
            "rounds": 1,
            "model": None,
            "conda_env": "polishing_pacbio"
        }
    
    else:  # HYBRID
        return {
            "tool": "polypolish",  # Handles mixed data
            "rounds": 1,
            "model": None,
            "conda_env": "polishing_hybrid"
        }


def main():
    """CLI entry point"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Detect sequencing platform")
    parser.add_argument("--fastq", required=True, help="Input FASTQ file")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--nanopore-model", help="Override Medaka model")
    
    args = parser.parse_args()
    
    fastq_path = Path(args.fastq)
    
    print(f"🔍 Detecting platform from {fastq_path.name}...")
    
    platform, stats = detect_platform(fastq_path)
    polisher_config = select_polisher(platform, args.nanopore_model)
    
    result = {
        "platform": platform,
        "read_stats": stats,
        "polisher_config": polisher_config,
        "recommendation": f"Use {polisher_config['tool']} with {polisher_config['rounds']} rounds"
    }
    
    # Save report
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)
    
    print(f"   ✅ Platform: {platform}")
    print(f"   📊 Avg read length: {stats['avg_length']:.0f} bp")
    print(f"   🔧 Polisher: {polisher_config['tool']}")
    print(f"   📁 Report saved: {args.output}")


if __name__ == "__main__":
    # Check if we are running inside Snakemake
    if 'snakemake' in globals():  # noqa: F821
        # Running as a Snakemake script
        fastq_path = Path(snakemake.input.fastq)  # noqa: F821
        output_report = Path(snakemake.output.report)  # noqa: F821
        nanopore_model = getattr(snakemake.params, 'nanopore_model', None)  # noqa: F821
        
        print(f"🔍 Detecting platform from {fastq_path.name}...")
        
        platform, stats = detect_platform(fastq_path)
        polisher_config = select_polisher(platform, nanopore_model)
        
        result = {
            "platform": platform,
            "read_stats": stats,
            "polisher_config": polisher_config,
            "recommendation": f"Use {polisher_config['tool']} with {polisher_config['rounds']} rounds"
        }
        
        # Save report
        with open(output_report, 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"   ✅ Platform: {platform}")
        print(f"   📊 Avg read length: {stats['avg_length']:.0f} bp")
        print(f"   🔧 Polisher: {polisher_config['tool']}")
    else:
        # Running as a standalone CLI script
        main()
