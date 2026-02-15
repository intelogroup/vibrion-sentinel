#!/usr/bin/env python3
import json
import argparse
import os
import subprocess
from pathlib import Path

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)

def get_path(rel_path):
    if os.path.exists(rel_path):
        return rel_path
    abs_path = os.path.join(PROJECT_ROOT, rel_path)
    if os.path.exists(abs_path):
        return abs_path
    return rel_path

# Loci for Phenotypic Prediction on CP003069.1 (Chr1) and CP003070.1 (Chr2)
# Coordinates based on 2010EL-1786 (N16961-like)
PHENO_LOCI = {
    "vps_cluster_I": ("CP003069.1", 966787, 980000, "VPS Biosynthesis (Biofilm)"),
    "vps_cluster_II": ("CP003070.1", 1036000, 1042000, "VPS Production (Chr2)"),
    "katB": ("CP003069.1", 2722301, 2724043, "Catalase (Oxidative Stress)"),
    "ahpC": ("CP003069.1", 1773372, 1773962, "Peroxiredoxin (Chlorine Resistance)"),
    "hapR": ("CP003069.1", 1086835, 1087515, "QS regulator (null -> rugose)")
}

def get_fasta_subject(input_file):
    """
    If input is FASTQ, creates a subsampled temporary FASTA for BLAST.
    """
    if not input_file.endswith(".fastq.gz") and not input_file.endswith(".fastq") and not input_file.endswith(".fq") and not input_file.endswith(".fq.gz"):
        return input_file, False # Already a FASTA or something else

    import tempfile
    try:
        # Create a temporary file
        fd, tmp_fasta = tempfile.mkstemp(suffix=".fasta", prefix="vibrion_pheno_")
        os.close(fd)
        
        # Subsample first 500k reads
        if input_file.endswith(".gz"):
            cmd = f"gzip -cd {input_file} | head -n 2000000 | awk 'NR%4==1 {{print \">\"substr($0,2)}} NR%4==2 {{print}}' > {tmp_fasta}"
        else:
            cmd = f"head -n 2000000 {input_file} | awk 'NR%4==1 {{print \">\"substr($0,2)}} NR%4==2 {{print}}' > {tmp_fasta}"
            
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        
        if os.path.exists(tmp_fasta) and os.path.getsize(tmp_fasta) > 0:
            return tmp_fasta, True
    except Exception as e:
        print(f"  Warning: FASTQ-to-FASTA conversion failed: {e}")
    
    return input_file, False

def check_locus_presence(input_fasta, ref_fasta, chrom, start, end, name):
    """Check for locus presence using blastn"""
    # Extract locus from ref
    locus_fasta = f"temp_{name}.fasta"
    cmd_extract = [
        "samtools", "faidx", ref_fasta, f"{chrom}:{start}-{end}"
    ]
    with open(locus_fasta, "w") as f:
        try:
            subprocess.run(cmd_extract, stdout=f, check=True)
        except Exception as e:
            print(f"  Warning: Reference region {chrom}:{start}-{end} not found in {os.path.basename(ref_fasta)}. IDs likely mismatch.")
            return False
    
    # Blast against input
    cmd_blast = [
        "blastn", "-query", locus_fasta, "-subject", input_fasta,
        "-outfmt", "6 pident length qlen", "-perc_identity", "80"
    ]
    
    try:
        res = subprocess.run(cmd_blast, capture_output=True, text=True, check=True)
        if not res.stdout.strip():
            return None
        
        # Parse top hit
        parts = res.stdout.strip().split("\n")[0].split("\t")
        pident = float(parts[0])
        length = int(parts[1])
        qlen = int(parts[2])
        return {"pident": pident, "coverage": length / qlen}
    except:
        return None
    finally:
        if os.path.exists(locus_fasta):
            os.remove(locus_fasta)

def predict_phenotype(input_fasta, ref_fasta=None):
    if ref_fasta is None:
        ref_fasta = get_path("data/references/2010EL-1786.fasta")
    
    print(f"🧬 PHENOTYPIC PREDICTION: Rugose State Detection")
    print("=" * 70)
    
    actual_fasta, needs_cleanup = get_fasta_subject(input_fasta)
    
    results = {}
    rugose_score = 0
    
    for key, (chrom, start, end, desc) in PHENO_LOCI.items():
        print(f"  Checking {key} ({desc})...")
        match = check_locus_presence(actual_fasta, ref_fasta, chrom, start, end, key)
        if match:
            print(f"    ✓ Detected: {match['pident']}% identity, {match['coverage']:.1%} coverage")
            results[key] = match
            if match['pident'] > 98:
                rugose_score += 1
        else:
            print(f"    ⚠️  NOT DETECTED or low quality match.")
            results[key] = "ABSENT"

    # Cleanup
    if needs_cleanup and os.path.exists(actual_fasta):
        os.remove(actual_fasta)

    # Specific Alert Logic
    # In a real system, we'd check for hapR deletions or vpsR SNPs.
    # Here we check for cluster integrity.
    risk = "LOW"
    if results.get("vps_cluster_I") != "ABSENT" and results.get("katB") != "ABSENT":
        risk = "MODERATE"
        if results.get("hapR") != "ABSENT" and results["hapR"]["pident"] < 99:
            print("  🚨 ALERT: hapR mutation detected (Potential Rugose Trigger).")
            risk = "HIGH"

    final_report = {
        "rugose_risk": risk,
        "score": rugose_score,
        "details": results,
        "interpretation": "Constitutive biofilm state (Rugose) predicted. Elevated environmental persistence suspected." if risk == "HIGH" else "Pathogen likely in smooth state."
    }
    
    return final_report

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Sample FASTA")
    parser.add_argument("--ref", help="Reference FASTA")
    parser.add_argument("--output", help="Output JSON path")
    args = parser.parse_args()
    
    ref = args.ref if args.ref else get_path("data/references/2010EL-1786.fasta")
    
    if not os.path.exists(ref):
        print(f"Error: Reference not found at {ref}")
        sys.exit(1)
        
    report = predict_phenotype(args.input, ref)
    
    output_path = args.output if args.output else "phenotype_prediction.json"
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    
    with open(output_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Results saved to {output_path}")
    
    output_path = os.path.join(os.path.dirname(args.input), f"{os.path.basename(args.input).split('.')[0]}_phenotype.json")
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    print(f"\nPhenotypic prediction saved to {output_path}")
