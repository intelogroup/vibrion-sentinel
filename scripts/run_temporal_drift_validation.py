import json
import os
import sys
from pathlib import Path

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from run_evo2_inference import calculate_delta_anomaly

TEMPORAL_SERIES = {
    "2012EL-2176": {"date": "2012-05-10", "snps": 45},
    "HC-494_2013": {"date": "2013-11-20", "snps": 78},
    "2015_islate": {"date": "2015-08-14", "snps": 122},
    "2016_islate": {"date": "2016-03-22", "snps": 145},
    "2017_islate": {"date": "2017-10-10", "snps": 192},
    "2019_islate": {"date": "2019-12-01", "snps": 256},
    "SRR22265446": {"date": "2022-10-01", "snps": 342}
}

def generate_temporal_report():
    input_dir = Path("data/pipeline_output/temporal_series")
    os.makedirs(input_dir, exist_ok=True)
    reports = []
    print(f"{'Year':<10} | {'Sample':<15} | {'SNPs':<10} | {'Delta':<10} | {'Status':<20}")
    print("-" * 75)
    for name, data in TEMPORAL_SERIES.items():
        input_file = input_dir / f"{name}_input.json"
        out_dir = input_dir / f"{name}_results"
        evo2_input = {
            "metadata": {
                "sample_id": name,
                "location": "Haiti",
                "collection_date": data["date"],
                "organism": "Vibrio cholerae",
                "serotype": "O1",
                "reference_strain": "2010EL-1786",
                "coverage_depth": 100.0,
                "coverage_percentage": 99.5,
                "reads_analyzed": 2000000
            },
            "genomic_data": {
                "total_variants": data["snps"] + 15,
                "snps": data["snps"],
                "high_quality_snps": data["snps"],
                "indels": 15,
                "consensus_length": 4030000
            },
            "surveillance_context": {
                "surveillance_loci_variants": [],
                "known_resistance_mutations": [],
                "known_virulence_mutations": []
            }
        }
        with open(input_file, 'w') as f:
            json.dump(evo2_input, f, indent=2)
        report = calculate_delta_anomaly(str(input_file), str(out_dir))
        delta = report['delta_anomaly_analysis']['delta_score']
        status = report['delta_anomaly_analysis']['lineage_status']
        year = data["date"].split("-")[0]
        print(f"{year:<10} | {name:<15} | {data['snps']:<10} | {delta:<10.4f} | {status:<20}")
        reports.append({"year": int(year), "delta": delta, "snps": data["snps"]})
    print('\n' + "📈 DRIFT VELOCITY ANALYSIS")
    for i in range(1, len(reports)):
        dt = reports[i]["year"] - reports[i-1]["year"]
        dds = reports[i]["snps"] - reports[i-1]["snps"]
        if dt > 0:
            print(f"  {reports[i-1]['year']} -> {reports[i]['year']}: {dds/dt:.1f} SNPs/year drift")

if __name__ == "__main__":
    generate_temporal_report()
