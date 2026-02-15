#!/usr/bin/env python3
import argparse
import json
import os
import sys
from datetime import datetime

# 7PET Lineage Constants
MAX_7PET_SNPS = 37000  # Threshold for non-7PET lineages (e.g. non-O1/O139)
NORMAL_DRIFT_PER_YEAR = 4.4  # Upper bound of expected hqSNP accumulation

def calculate_velocity(snp_count, reference_year, sample_year):
    years = sample_year - reference_year
    if years <= 0:
        return 0.0
    return snp_count / years

def clonal_filter(vcf_path, reference_year=2010):
    print(f"🔬 VIBRION CLONAL FILTER & VELOCITY TRACKER")
    print("=" * 70)
    
    if not os.path.exists(vcf_path):
        print(f"Error: VCF not found at {vcf_path}")
        return None

    # Count hqSNPs
    snp_count = 0
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.split("\t")
            ref, alt = parts[3], parts[4]
            if len(ref) == 1 and len(alt) == 1:
                snp_count += 1
                
    # Sample year (try to extract from filename or use current)
    # Assume 2022 if Haiti-2022 mentioned
    sample_year = 2022 
    velocity = calculate_velocity(snp_count, reference_year, sample_year)
    
    alarm_37k = snp_count > MAX_7PET_SNPS
    velocity_anomaly = velocity > (NORMAL_DRIFT_PER_YEAR * 5) # 5x acceleration
    
    print(f"  hqSNP Count: {snp_count:,}")
    print(f"  Drift Velocity: {velocity:.2f} SNPs/year (Ref: {reference_year}, Sample: {sample_year})")
    
    if alarm_37k:
        print(f"  🚨 ALARM: SNP count ({snp_count:,}) exceeds 37k threshold!")
        print(f"     Interpretation: Sample is likely a non-7PET lineage (potential decoy).")
    
    if velocity_anomaly:
        print(f"  🚨 ALARM: Drift velocity ({velocity:.2f}) is {velocity/NORMAL_DRIFT_PER_YEAR:.1f}x normal!")
        print(f"     Interpretation: Rapid evolution or recombination event detected.")

    return {
        "snp_count": snp_count,
        "reference_year": reference_year,
        "sample_year": sample_year,
        "drift_velocity": velocity,
        "alarm_37k_triggered": alarm_37k,
        "velocity_anomaly": velocity_anomaly,
        "lineage_consistency": "LOW" if alarm_37k or velocity_anomaly else "HIGH"
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf", help="Input VCF file")
    parser.add_argument("--ref_year", type=int, default=2010)
    args = parser.parse_args()
    
    result = clonal_filter(args.vcf, args.ref_year)
    if result:
        output_path = args.vcf.replace(".vcf", "_clonal_stats.json")
        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\nClonal stats saved to {output_path}")
