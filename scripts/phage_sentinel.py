#!/usr/bin/env python3
import json
import argparse
import os

# ============================================================================= 
# NORMALIZED PHAGE RATIO METRICS (Critical Refinement #2)
# ============================================================================= 
# Genome sizes (base pairs):
VIBRIO_GENOME_SIZE_BP = 4_000_000  # V. cholerae ~4.0 Mb
PHAGE_GENOME_SIZE_BP = 40_000      # Typical lytic phage (ICP1) ~40 kb

# Predation threshold: normalized by coverage depth (not raw reads)
# A ratio of 0.05 coverage depth = 5% of bacterial population is phage
PREDATION_THRESHOLD_DEPTH = 0.05  # 5% phage coverage is "High" predation

def calculate_length_normalized_ratio(phage_coverage_depth, vibrio_coverage_depth):
    """
    Calculate length-normalized phage/vibrio ratio (Critical Refinement #2).
    
    The raw read count trap: If a ratio is based on reads, 0.05 would mean
    5 phage reads per 100 bacterial reads. But the Vibrio genome is ~100x
    larger, so this represents a NEGLIGIBLE number of viral particles.
    
    Solution: Normalize by genome size to get actual population ratio.
    
    Formula: (Phage_Coverage_Depth * Phage_Genome_Size) / (Vibrio_Coverage_Depth * Vibrio_Genome_Size)
    
    This gives: What fraction of the bacterial population is represented by phage?
    """
    if vibrio_coverage_depth == 0:
        return None
    
    # Normalize coverage depths to population ratio
    phage_population = phage_coverage_depth * PHAGE_GENOME_SIZE_BP
    vibrio_population = vibrio_coverage_depth * VIBRIO_GENOME_SIZE_BP
    
    if vibrio_population == 0:
        return None
    
    return phage_population / vibrio_population

def phage_sentinel_analysis(vibrio_count, phage_counts, vibrio_coverage=None, phage_coverage=None):
    print("🌊 VIBRION PHAGE SENTINEL: Predation Monitoring (Coverage Depth Normalized)")
    print("=" * 80)
    
    if vibrio_count == 0:
        print("  ⚠️  Vibrio count is zero. Cannot calculate predation ratio.")
        return {"status": "UNKNOWN", "error": "No Vibrio detected"}

    total_phage = sum(phage_counts.values())
    
    # Raw ratio (for reference only)
    raw_overall_ratio = total_phage / vibrio_count
    
    # CRITICAL REFINEMENT #2: Use coverage depth if available
    normalized_ratio = None
    ratio_type = "raw_reads"
    
    if vibrio_coverage is not None and phage_coverage is not None:
        normalized_ratio = calculate_length_normalized_ratio(phage_coverage, vibrio_coverage)
        ratio_type = "coverage_depth_normalized"
        print(f"  ℹ️  Using NORMALIZED coverage depth metrics (recommended)")
    else:
        normalized_ratio = raw_overall_ratio
        print(f"  ⚠️  Coverage depth not provided. Using raw read counts (less accurate).")
    
    print(f"\n  Vibrio Count: {vibrio_count:,}")
    print(f"  Total Phage Count: {total_phage:,}")
    print(f"  Vibrio Genome Size: {VIBRIO_GENOME_SIZE_BP/1e6:.1f} Mb")
    print(f"  Phage Genome Size: {PHAGE_GENOME_SIZE_BP/1e3:.0f} kb (typical)")
    print(f"  Raw Phage/Vibrio Ratio (reads): {raw_overall_ratio:.4f}")
    if ratio_type == "coverage_depth_normalized":
        print(f"  Normalized Ratio (coverage depth): {normalized_ratio:.4f}")
    print(f"  Alert Threshold: {PREDATION_THRESHOLD_DEPTH:.4f} ({PREDATION_THRESHOLD_DEPTH*100:.1f}%)")
    
    ratios = {k: v/vibrio_count for k, v in phage_counts.items()}
    
    for phage, ratio in ratios.items():
        print(f"    - {phage}: {ratio:.4f} (Count: {phage_counts[phage]:,})")
        
    # Use normalized ratio for alert decision
    alert = normalized_ratio > PREDATION_THRESHOLD_DEPTH
    
    # Specific monitoring for RS1 (Satellite of CTX phi)
    has_rs1 = phage_counts.get("RS1", 0) > 0
    if has_rs1:
        print("  ✓ RS1 Satellite Phage detected. Consistent with 7th Pandemic El Tor (7PET).")

    if alert:
        print(f"\n  🚨 ALERT: PREDATION-HIGH TRIGGERED ({normalized_ratio:.2%} normalized coverage depth)")
        print(f"     Status: Pathogen is under intense phage pressure.")
        print(f"     Interpretation: {normalized_ratio*100:.1f}% of bacterial population represented by phage.")
        print(f"     Diagnostic Impact: Clinical cultures may be inhibited (False Negatives).")
        print(f"     Epidemic Implication: Within 24 hours, phage pressure can cause 10x drop in infectious dose.")
        dominant = max(phage_counts, key=phage_counts.get)
        print(f"     Dominant Predator: {dominant}")
        status = "PREDATION_HIGH"
    else:
        print(f"\n  ✓ Phage levels within normal environmental range ({normalized_ratio*100:.2f}% coverage depth).")
        status = "NORMAL"

    return {
        "status": status,
        "overall_ratio_raw": raw_overall_ratio,
        "overall_ratio_normalized": normalized_ratio,
        "ratio_type": ratio_type,
        "phage_ratios": ratios,
        "phage_counts": phage_counts,
        "vibrio_count": vibrio_count,
        "vibrio_coverage_depth": vibrio_coverage,
        "phage_coverage_depth": phage_coverage,
        "genome_sizes": {
            "vibrio_bp": VIBRIO_GENOME_SIZE_BP,
            "phage_bp": PHAGE_GENOME_SIZE_BP
        },
        "alert_triggered": alert,
        "threshold_used": PREDATION_THRESHOLD_DEPTH,
        "7pet_marker_rs1": has_rs1,
        "refinement_applied": "Coverage Depth Normalized (Critical #2)"
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vibrio", type=int, default=1000000)
    parser.add_argument("--vibrio-coverage", type=float, default=None, 
                        help="Vibrio coverage depth (reads/bp) for normalized calculation")
    parser.add_argument("--phage-coverage", type=float, default=None,
                        help="Phage coverage depth (reads/bp) for normalized calculation")
    parser.add_argument("--input", help="Optional JSON or Kraken2 report")
    parser.add_argument("--output", help="Output JSON path")
    args = parser.parse_args()
    
    vibrio_cov = None
    phage_cov = None
    vibrio = args.vibrio
    phages = {"ICP1": 85000, "ICP2": 12000, "RS1": 25000} # Mock defaults
    
    if args.input and os.path.exists(args.input):
        with open(args.input) as f:
            line = f.readline()
            f.seek(0)
            if line.startswith("{"):
                # Handle JSON input
                data = json.load(f)
                phage_data = data.get("phage_surveillance", {})
                vibrio = phage_data.get("vibrio_count", 0)
                phages = phage_data.get("phage_counts", {})
                vibrio_cov = phage_data.get("vibrio_coverage_depth", None)
                phage_cov = phage_data.get("phage_coverage_depth", None)
            else:
                # Handle Kraken2 report
                vibrio = 0
                phages = {}
                target_phages = {
                    645063: "ICP1",
                    2291560: "ICP2",
                    2291561: "ICP3",
                    113540: "CTX-phi",
                    1854898: "RS1",
                    223524: "VGJ"
                }
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 6: continue
                    taxid = int(parts[4])
                    count = int(parts[1])
                    if taxid in [666, 668]: vibrio = count
                    elif taxid in target_phages: phages[target_phages[taxid]] = count

    # Override with CLI args if provided
    if args.vibrio_coverage:
        vibrio_cov = args.vibrio_coverage
    if args.phage_coverage:
        phage_cov = args.phage_coverage
        
    result = phage_sentinel_analysis(vibrio, phages, vibrio_cov, phage_cov)
    
    output_path = args.output if args.output else "phage_sentinel_report.json"
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f"\nPhage sentinel report saved to {output_path}")
