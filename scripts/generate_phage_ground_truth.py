import json
import os

# Simulated data for a "Predation-High" scenario in the Haiti 2022 resurgence
# This sample represents an environmental water sample where phages are booming, 
# likely preceding a drop in human cases or masking the true bacterial burden.

phage_test_data = {
    "metadata": {
        "sample_id": "HAITI-2022-ENV-PHAGE-01",
        "location": "Artibonite River, Gonaïves",
        "collection_date": "2022-11-15",
        "organism": "Vibrio cholerae",
        "serotype": "O1 Ogawa",
        "reference_strain": "2010EL-1786 (Haiti 2010)",
        "coverage_percentage": 94.2,
        "coverage_depth": 45.8,
        "reads_analyzed": 5000000
    },
    "genomic_data": {
        "total_variants": 1420,
        "snps": 842, # 2022 variant drift (~450 SNPs/12 years)
        "high_quality_snps": 795,
        "indels": 58,
        "consensus_length": 4030000
    },
    "surveillance_context": {
        "known_resistance_mutations": ["gyrA_S83L", "parC_S85L"],
        "known_virulence_mutations": ["ctxB_genotype7"],
        "surveillance_loci_variants": [
            {"gene": "wbeT", "gene_type": "serotype", "variant": "G->A at pos 120"},
            {"gene": "scrA", "gene_type": "metabolic", "variant": "No change"},
            {"gene": "lip", "gene_type": "metabolic", "variant": "No change"}
        ]
    },
    "phage_surveillance": {
        "vibrio_count": 1000000,
        "phage_counts": {
            "ICP1": 85000,  # ~8.5% ratio (High predation)
            "ICP2": 12000,
            "JSF4": 5000,
            "RS1": 25000    # Satellite phage indicating active HGT/virulence mobility
        }
    }
}

OUTPUT_PATH = "/Users/kalinovdameus/Developer/Vibrion/data/phage_sentinel_test.json"

def main():
    os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
    with open(OUTPUT_PATH, "w") as f:
        json.dump(phage_test_data, f, indent=4)
    print(f"✅ Created Phage Sentinel Test Data: {OUTPUT_PATH}")

if __name__ == "__main__":
    main()
