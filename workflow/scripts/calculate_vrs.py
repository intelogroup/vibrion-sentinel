"""
Calculate VRS (Vibrio Resurgence Score)
Integrates Evo2 Sentinel analysis with decontamination metrics
for CDC-compliant public health risk assessment
"""

import json
from pathlib import Path
from datetime import datetime

# Snakemake inputs/outputs
evo2_result_file = Path(snakemake.input.evo2_result) # noqa: F821
vibrio_stats_file = Path(snakemake.input.vibrio_stats) # noqa: F821
hostile_stats_file = Path(snakemake.input.hostile_stats) # noqa: F821
output_file = Path(snakemake.output.vrs) # noqa: F821
sample_id = snakemake.params.sample_id # noqa: F821

output_file.parent.mkdir(parents=True, exist_ok=True)

print(f"📊 VRS Calculation: {sample_id}")
print("   Integrating Evo2 Sentinel + abundance + decontamination metrics")

# Load all analysis components
with open(evo2_result_file) as f:
    evo2_data = json.load(f)

with open(vibrio_stats_file) as f:
    vibrio_stats = json.load(f)

with open(hostile_stats_file) as f:
    hostile_stats = json.load(f)

# Extract key metrics
sentinel_score = evo2_data.get("sentinel_score", 5)  # 0-10 scale
alert_level = evo2_data.get("alert_level", "advisory")
classification = evo2_data.get("classification", "Unknown")
trajectory = evo2_data.get("evolutionary_trajectory", "STABLE_ENDEMIC")
vibrio_reads = vibrio_stats.get("vibrio_reads", 0)
total_reads = vibrio_stats.get("total_reads", 1)
vibrio_percentage = (vibrio_reads / total_reads * 100) if total_reads > 0 else 0

print(f"   Sentinel Score: {sentinel_score}/10")
print(f"   Alert Level: {alert_level}")
print(f"   Vibrio Abundance: {vibrio_percentage:.1f}%")

# =============================================================================
# VRS CALCULATION ALGORITHM (0-100 scale)
# =============================================================================
# Based on three components:
# 1. Genomic threat level (Evo2 Sentinel score) - 60% weight
# 2. Vibrio abundance (environmental load) - 30% weight
# 3. Alert level multiplier - amplifies high-threat scenarios

# Component 1: Genomic Threat (0-60 points)
# Sentinel score 0-10 maps to 0-60 points
genomic_threat_score = (sentinel_score / 10) * 60

# Component 2: Abundance (0-30 points)
# Higher Vibrio percentage = higher environmental risk
# Cap at 100% for percentage, scale to 30 points
abundance_score = min(vibrio_percentage, 100) * 0.3

# Component 3: Alert Level Multiplier
alert_multipliers = {
    "critical": 1.5,      # 50% increase for critical alerts
    "automated": 1.3,     # 30% increase for automated alerts
    "advisory": 1.1,      # 10% increase for advisory
    "log_only": 1.0       # No increase for log-only
}
multiplier = alert_multipliers.get(alert_level, 1.0)

# Baseline VRS before multiplier
baseline_vrs = genomic_threat_score + abundance_score

# Final VRS with alert multiplier
vrs_raw = baseline_vrs * multiplier

# Cap at 100
vrs = min(100, round(vrs_raw, 1))

# Risk categorization thresholds
def categorize_risk(vrs_value: float) -> dict:
    """Categorize VRS into actionable risk levels"""
    if vrs_value >= 80:
        return {
            "category": "CRITICAL",
            "color": "red",
            "action": "Immediate public health intervention required",
            "response_time": "< 24 hours"
        }
    elif vrs_value >= 60:
        return {
            "category": "HIGH",
            "color": "orange",
            "action": "Enhanced surveillance and contact tracing",
            "response_time": "< 72 hours"
        }
    elif vrs_value >= 40:
        return {
            "category": "MODERATE",
            "color": "yellow",
            "action": "Routine monitoring with elevated vigilance",
            "response_time": "< 1 week"
        }
    else:
        return {
            "category": "LOW",
            "color": "green",
            "action": "Standard surveillance protocol",
            "response_time": "Routine"
        }

risk_info = categorize_risk(vrs)

print("   VRS Calculation:")
print(f"      Genomic Threat: {genomic_threat_score:.1f}/60")
print(f"      Abundance: {abundance_score:.1f}/30")
print(f"      Alert Multiplier: {multiplier}x")
print(f"      Final VRS: {vrs}/100 ({risk_info['category']})")

# Compile comprehensive result for MongoDB and API
result = {
    "sample_id": sample_id,
    "timestamp": datetime.utcnow().isoformat(),
    
    # VRS Score and Risk Assessment
    "vrs": vrs,
    "risk_category": risk_info["category"],
    "risk_color": risk_info["color"],
    "recommended_action": risk_info["action"],
    "response_time": risk_info["response_time"],
    
    # VRS Component Breakdown
    "vrs_components": {
        "genomic_threat_score": round(genomic_threat_score, 1),
        "abundance_score": round(abundance_score, 1),
        "baseline_vrs": round(baseline_vrs, 1),
        "alert_multiplier": multiplier,
        "final_vrs": vrs
    },
    
    # Evo2 Sentinel Analysis Summary
    "evo2_sentinel": {
        "score": sentinel_score,
        "alert_level": alert_level,
        "classification": classification,
        "evolutionary_trajectory": trajectory,
        "confidence": evo2_data.get("confidence", "Medium"),
        "dread_signals": evo2_data.get("dread_signals", {}),
        "forward_trajectory": evo2_data.get("forward_trajectory", {})
    },
    
    # Decontamination and Abundance Metrics
    "decontamination": {
        "total_reads": total_reads,
        "vibrio_reads": vibrio_reads,
        "vibrio_percentage": round(vibrio_percentage, 2),
        "human_reads_removed": hostile_stats.get("reads_removed", 0),
        "reads_after_cleaning": hostile_stats.get("reads_out", total_reads)
    },
    
    # Pipeline Status
    "pipeline_status": "complete",
    "pipeline_version": "snakemake.v1.0",
    "analysis_complete": True
}

print("   ✅ VRS Calculation Complete")
print(f"   🎯 Final Score: {vrs}/100 - {risk_info['category']} Risk")
print(f"   📋 Recommended Action: {risk_info['action']}")

# Write final result
with open(output_file, 'w') as f:
    json.dump(result, f, indent=2)

print(f"   📁 Result saved to {output_file.name}")
