#!/usr/bin/env python3
"""
Verification script for Metabolic Grounding (scrABC/lip cluster analysis).
Tests the system's ability to exclude V. mimicus decoys.
"""
import json
import os
import subprocess
from pathlib import Path

def verify_mimicus_exclusion():
    print("🧪 VERIFYING METABOLIC GROUNDING (V. mimicus Exclusion)")
    print("=" * 70)
    
    decoy_input = "data/decoy_validation/V_MIMICUS_DECOY_evo2_input.json"
    output_dir = "data/decoy_validation/V_MIMICUS_DECOY_results"
    
    # Ensure decoy input exists with the metabolic markers flagged
    # In a real run, the functional_annotation.py would do this
    if not os.path.exists(decoy_input):
        print(f"Creating mock decoy input: {decoy_input}")
        mock_input = {
            "metadata": {
                "sample_id": "V_MIMICUS_DECOY",
                "organism": "Vibrio species (V. mimicus candidate)",
                "serotype": "Unknown",
                "location": "Synthetic Test",
                "collection_date": "2026-01-22",
                "reference_strain": "2010EL-1786 (VC O1)",
                "coverage_depth": 50.0,
                "coverage_percentage": 95.0,
                "reads_analyzed": 100000
            },
            "genomic_data": {
                "total_variants": 35000,
                "snps": 32000,
                "high_quality_snps": 1200,
                "indels": 3000,
                "consensus_length": 4000000
            },
            "surveillance_context": {
                "known_resistance_mutations": [],
                "known_virulence_mutations": [],
                "surveillance_loci_variants": [
                    {"gene": "metabolic_scrABC", "gene_type": "metabolic", "change": "DIVERGENT"},
                    {"gene": "metabolic_lip", "gene_type": "metabolic", "change": "DIVERGENT"},
                    # Simulate 50 metabolic variants to trigger the >10 threshold
                    *[{"gene": "metabolic_locus", "gene_type": "metabolic"} for _ in range(50)]
                ]
            }
        }
        os.makedirs(os.path.dirname(decoy_input), exist_ok=True)
        with open(decoy_input, 'w') as f:
            json.dump(mock_input, f, indent=2)

    # Run the inference engine (which now has the exclusion logic)
    print("Running EVO2 Inference Engine...")
    cmd = [
        "python3", "run_evo2_inference.py",
        "--input", decoy_input,
        "--outdir", output_dir
    ]
    subprocess.run(cmd, capture_output=True, text=True)
    
    # Check the result
    result_file = Path(output_dir) / "evo2_threat_assessment.json"
    if result_file.exists():
        with open(result_file, 'r') as f:
            report = json.load(f)
            
        print("\n📈 ANALYSIS RESULT:")
        print(f"  Threat Level: {report['threat_assessment']['threat_level']}")
        print(f"  Threat Category: {report['threat_assessment']['threat_category']}")
        print(f"  Delta Score: {report['delta_anomaly_analysis']['delta_score']:.4f}")
        
        metabolic = report['delta_anomaly_analysis']['components']['metabolic_divergence']
        print(f"  Metabolic Component: {metabolic['interpretation']}")
        print(f"  Exclusion Triggered: {metabolic.get('exclusion_triggered', False)}")
        
        if report['threat_assessment']['threat_level'] == "EXCLUDE_DECOY":
            print("\n✅ SUCCESS: V. mimicus was successfully identified and excluded via metabolic grounding.")
        else:
            print("\n❌ FAILURE: Exclusion logic did not trigger as expected.")
    else:
        print(f"❌ Error: {result_file} not generated.")

if __name__ == "__main__":
    verify_mimicus_exclusion()
