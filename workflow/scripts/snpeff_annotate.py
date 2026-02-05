#!/usr/bin/env python3
"""
SnpEff Variant Annotation Wrapper (with BCFtools fallback)

Primary: SnpEff (LGPLv3) - Best bacterial annotation
Fallback: BCFtools csq (MIT) - Commercial-safe alternative

This script annotates VCF files with predicted variant effects.
"""

import subprocess
import os
import json
import gzip

# SnpEff configuration - MATCHES our reference genome (2010EL-1786)
SNPEFF_DATABASE = "Vibrio_cholerae_o1_str_2010el_1786"

def check_snpeff_available():
    """Check if SnpEff is installed and database is available."""
    try:
        result = subprocess.run(["snpEff", "-version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def check_bcftools_available():
    """Check if BCFtools is available for fallback."""
    try:
        result = subprocess.run(["bcftools", "--version"], capture_output=True, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def run_snpeff(input_vcf, output_vcf, stats_file, log_file):
    """Run SnpEff variant annotation."""
    cmd = [
        "snpEff", "ann",
        "-v",
        "-stats", stats_file.replace(".json", ".html"),
        "-csvStats", stats_file.replace(".json", ".csv"),
        SNPEFF_DATABASE,
        input_vcf
    ]
    
    print(f"Running SnpEff: {' '.join(cmd)}")
    
    try:
        with open(log_file, 'w') as log:
            # Run SnpEff and pipe to gzip in python
            with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=log) as proc:
                with gzip.open(output_vcf, 'wb') as out:
                    import shutil
                    if proc.stdout:
                        shutil.copyfileobj(proc.stdout, out)
                
                proc.wait()
                return proc.returncode == 0
    except Exception as e:
        print(f"Error running SnpEff: {e}")
        return False

def run_bcftools_fallback(input_vcf, output_vcf, gff_file, log_file):
    """Run BCFtools csq as fallback (MIT license)."""
    # BCFtools csq requires a GFF file for annotation
    if not os.path.exists(gff_file):
        print(f"Warning: GFF file not found at {gff_file}. Cannot run BCFtools csq.")
        return False
    
    cmd = [
        "bcftools", "csq",
        "-f", "data/references/2010EL-1786.fasta",
        "-g", gff_file,
        "-o", output_vcf,
        "-Oz",
        input_vcf
    ]
    
    print(f"Running BCFtools csq (fallback): {' '.join(cmd)}")
    
    with open(log_file, 'w') as log:
        proc = subprocess.run(cmd, stderr=log)
    
    return proc.returncode == 0

def parse_snpeff_vcf(annotated_vcf):
    """Parse SnpEff annotated VCF and extract impact summary."""
    impacts = {"HIGH": [], "MODERATE": [], "LOW": [], "MODIFIER": []}
    
    open_func = gzip.open if annotated_vcf.endswith('.gz') else open
    mode = 'rt' if annotated_vcf.endswith('.gz') else 'r'
    
    with open_func(annotated_vcf, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
            
            info = fields[7]
            # Parse ANN field from SnpEff
            if 'ANN=' in info:
                ann_part = [x for x in info.split(';') if x.startswith('ANN=')]
                if ann_part:
                    ann = ann_part[0].replace('ANN=', '')
                    annotations = ann.split(',')
                    for annot in annotations:
                        parts = annot.split('|')
                        if len(parts) >= 3:
                            impact = parts[2]  # Impact is the 3rd field
                            if impact in impacts:
                                impacts[impact].append({
                                    "position": f"{fields[0]}:{fields[1]}",
                                    "ref": fields[3],
                                    "alt": fields[4],
                                    "effect": parts[1] if len(parts) > 1 else "unknown",
                                    "gene": parts[3] if len(parts) > 3 else "unknown"
                                })
    
    return impacts

def main():
    # Snakemake integration
    input_vcf = snakemake.input.vcf # noqa: F821
    output_vcf = snakemake.output.annotated_vcf # noqa: F821
    stats_file = snakemake.output.stats # noqa: F821
    log_file = snakemake.log[0] # noqa: F821
    
    gff_file = snakemake.params.get("gff", "data/references/2010EL-1786.gff3") # noqa: F821
    
    success = False
    method = None
    
    # Try SnpEff first (best for bacteria)
    if check_snpeff_available():
        print("Using SnpEff (LGPLv3) for variant annotation...")
        success = run_snpeff(input_vcf, output_vcf, stats_file, log_file)
        method = "snpeff"
    
    # Fallback to BCFtools csq (MIT license)
    if not success and check_bcftools_available():
        print("SnpEff failed. Falling back to BCFtools csq (MIT)...")
        success = run_bcftools_fallback(input_vcf, output_vcf, gff_file, log_file)
        method = "bcftools_csq"
    
    if not success:
        raise RuntimeError("Both SnpEff and BCFtools csq failed. Cannot annotate variants.")
    
    # Parse and summarize impacts
    impacts = parse_snpeff_vcf(output_vcf)
    
    stats = {
        "method": method,
        "database": SNPEFF_DATABASE if method == "snpeff" else "GFF annotation",
        "high_impact_count": len(impacts["HIGH"]),
        "moderate_impact_count": len(impacts["MODERATE"]),
        "low_impact_count": len(impacts["LOW"]),
        "modifier_count": len(impacts["MODIFIER"]),
        "high_impact_variants": impacts["HIGH"][:10],  # Top 10 for report
        "moderate_impact_variants": impacts["MODERATE"][:10]
    }
    
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"Annotation complete. {stats['high_impact_count']} HIGH, {stats['moderate_impact_count']} MODERATE impact variants found.")

if __name__ == "__main__":
    main()
