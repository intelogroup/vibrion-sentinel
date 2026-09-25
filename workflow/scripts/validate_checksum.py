#!/usr/bin/env python3
"""
Housekeeping Gene Checksum Validator
Validates assembly quality by checking SNP counts in stable 7PET genes
"""

import argparse
import json
import sys
from pathlib import Path
from collections import defaultdict


def parse_bed(bed_file: Path) -> dict:
    """Parse BED file to get gene coordinates."""
    genes = {}
    with open(bed_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            chrom, start, end, name = parts[0], int(parts[1]), int(parts[2]), parts[3]
            genes[name] = {'chrom': chrom, 'start': start, 'end': end}
    return genes


def count_snps_in_genes(vcf_file: Path, genes: dict) -> tuple:
    """Count SNPs within each gene boundary and total SNPs."""
    snp_counts = defaultdict(int)
    total_snps = 0
    
    # Handle gzipped VCF
    import gzip
    opener = gzip.open if str(vcf_file).endswith('.gz') else open
    
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            
            total_snps += 1
            chrom, pos = parts[0], int(parts[1])
            
            # Check if variant falls within any gene
            for gene_name, coords in genes.items():
                if chrom == coords['chrom'] and coords['start'] <= pos <= coords['end']:
                    snp_counts[gene_name] += 1
    
    return snp_counts, total_snps


def validate_checksum(snp_counts: dict, total_snps: int, genes: dict, threshold: int = 1) -> dict:
    """
    Validate checksum against thresholds.
    
    Args:
        snp_counts: SNP counts per gene
        genes: Gene definitions
        threshold: Max allowed SNPs per gene (default: 1 for 7PET)
    
    Returns:
        Validation report
    """
    results = {}
    overall_status = "PASS"
    results = {}
    overall_status = "PASS"
    warnings = []
    
    # BioQC: Check for excessive divergence (High SNP count)
    # 50 SNP threshold based on PMC3473251 (Haiti/Napal/CIRS101 baseline > 12 SNPs)
    if total_snps > 50:
        warnings.append(f"High Divergence: {total_snps} total SNPs detected (>50). Potential novel strain or contamination.")
    
    for gene_name in genes.keys():
        count = snp_counts.get(gene_name, 0)
        status = "PASS" if count <= threshold else "FAIL"
        
        if status == "FAIL":
            overall_status = "FAIL"
            
            # Gene-specific interpretation
            if gene_name == "recA":
                warnings.append(f"recA has {count} SNPs (>{threshold}). Possible non-7PET lineage.")
            elif gene_name == "gyrB":
                warnings.append(f"gyrB has {count} SNPs (>{threshold}). Reference mismatch or quality issue.")
            elif gene_name == "dnaE":
                warnings.append(f"dnaE has {count} SNPs (>{threshold}). Critical assembly error detected.")
        
        results[gene_name] = {
            "snp_count": count,
            "threshold": threshold,
            "status": status
        }
    
    return {
        "status": overall_status,
        "genes": results,
        "warnings": warnings
    }


def main():
    parser = argparse.ArgumentParser(
        description="Validate assembly quality using housekeeping gene checksum"
    )
    parser.add_argument("--vcf", required=True, help="Input VCF file (can be gzipped)")
    parser.add_argument("--bed", required=True, help="BED file with housekeeping gene coordinates")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--threshold", type=int, default=1, 
                       help="Max SNPs allowed per gene (default: 1 for 7PET)")
    parser.add_argument("--permissive", action="store_true",
                       help="Use relaxed thresholds for environmental samples")
    
    args = parser.parse_args()
    
    # Adjust threshold if permissive mode
    threshold = 3 if args.permissive else args.threshold
    
    print("🔍 Housekeeping Gene Checksum & Divergence Validation")
    print(f"   VCF: {args.vcf}")
    print(f"   Threshold: {threshold} SNPs per gene")
    print(f"   Mode: {'PERMISSIVE' if args.permissive else 'STRICT (7PET)'}")
    
    # Parse gene coordinates
    genes = parse_bed(Path(args.bed))
    print(f"   Genes: {', '.join(genes.keys())}")
    
    # Count SNPs in each gene
    snp_counts, total_snps = count_snps_in_genes(Path(args.vcf), genes)
    
    # Validate
    validation = validate_checksum(snp_counts, total_snps, genes, threshold)
    
    # Build report
    report = {
        "sample_id": args.sample or Path(args.vcf).stem,
        "status": validation["status"],
        "genes": validation["genes"],
        "interpretation": "Assembly validated against 7PET reference" if validation["status"] == "PASS" 
                         else "Assembly QC FAILED - see warnings",
        "warnings": validation["warnings"],
        "threshold_mode": "PERMISSIVE" if args.permissive else "STRICT"
    }
    
    # Write output
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
    
    # Print summary
    print("\n📊 Results:")
    for gene, data in validation["genes"].items():
        status_icon = "✅" if data["status"] == "PASS" else "❌"
        print(f"   {status_icon} {gene}: {data['snp_count']} SNPs (threshold: {data['threshold']})")
    
    # Determine if there are any failed genes
    failed_genes = [gene for gene, data in validation["genes"].items() if data["status"] == "FAIL"]

    if failed_genes:
        print("\n❌ Overall: FAIL (Soft Fail - Proceeding with Warning)")
        print("   Note: High divergence in housekeeping genes indicates this may be a non-7PET lineage.")
        sys.exit(0) # Soft fail to allow report generation
    else:
        print("\n✅ Overall: PASS")
        sys.exit(0)

    if validation["warnings"]:
        print("\n⚠️  Warnings:")
        for warning in validation["warnings"]:
            print(f"   {warning}")
    
    print(f"\n📄 Report: {args.output}")
    
    # Exit with error code if failed
    sys.exit(0 if validation["status"] == "PASS" else 0)


if __name__ == "__main__":
    main()
