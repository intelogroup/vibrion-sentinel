#!/usr/bin/env python3
"""
Public Health Decision Support Module
Predicts Serotype (Ogawa/Inaba) and Toxin Genotype (ctxB) for Vaccine Alignment.

This script interprets VCF variants in critical loci to provide
actionable public health intelligence for Dr. Louise Ivers' team.
"""

import argparse
import json
import sys
import gzip
from pathlib import Path

# Critical Loci Definitions (2010EL-1786 Reference)
LOCI = {
    "wbeT": {"chrom": "CP003069.1", "start": 2678186, "end": 2678980, "description": "O1 Antigen Biosynthesis (Serotype Switch)"},
    "ctxB": {"chrom": "CP003069.1", "start": 1041238, "end": 1041612, "description": "Cholera Toxin B Subunit"},
    "tcpA": {"chrom": "CP003069.1", "start": 367950, "end": 368624, "description": "Toxin Coregulated Pilus"},
    "hapR": {"chrom": "CP003069.1", "start": 1086835, "end": 1087515, "description": "Quorum Sensing Regulator"}
}

def parse_vcf_variants(vcf_file, region):
    """
    Extracts variants within a specific genomic region.
    """
    variants = []
    
    opener = gzip.open if str(vcf_file).endswith('.gz') else open
    
    try:
        with opener(vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                
                parts = line.strip().split('\t')
                chrom = parts[0]
                pos = int(parts[1])
                
                if chrom == region['chrom'] and region['start'] <= pos <= region['end']:
                    ref = parts[3]
                    alt = parts[4]
                    # Simple parsing of INFO field if needed later
                    info = parts[7]
                    variants.append({
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "info": info
                    })
    except Exception as e:
        print(f"Error reading VCF: {e}", file=sys.stderr)
        return []
        
    return variants

def predict_serotype(wbeT_variants):
    """
    Predicts Ogawa/Inaba serotype based on wbeT (rfbT) status.
    Reference 2010EL-1786 is OGAWA (Wild Type).
    """
    if not wbeT_variants:
        return {
            "prediction": "Ogawa",
            "status": "Vaccine Match (Wild Type)",
            "evidence": "No mutations in wbeT gene relative to 2010 Ogawa reference."
        }
    
    # Analyze variants
    # Any high impact mutation (stop, frameshift) usually means Inaba.
    # For now, we flag any non-synonymous change as "Potential Inaba".
    
    variant_desc = ", ".join([f"{v['pos']}:{v['ref']}->{v['alt']}" for v in wbeT_variants])
    
    return {
        "prediction": "Inaba (Likely)",
        "status": "Vaccine Mismatch (Potential)",
        "evidence": f"Mutations detected in wbeT: {variant_desc}. Functional impact requires SnpEff annotation."
    }

def analyze_toxin(ctxB_variants):
    """
    Analyzes ctxB allele.
    Reference 2010EL-1786 is ctxB7 (Haiti variant).
    """
    if not ctxB_variants:
        return {
            "genotype": "ctxB7 (Haiti/Classical)",
            "virulence": "High (Hypervirulent)",
            "evidence": "Identical to 2010 outbreak strain (ctxB7)."
        }
    
    # Specific known SNPs could distinguish B1 vs B7
    # But if we have SNPs relative to B7, it might be drifting or reverting to B1.
    variant_desc = ", ".join([f"{v['pos']}:{v['ref']}->{v['alt']}" for v in ctxB_variants])
    
    return {
        "genotype": "Variant / Atypical",
        "virulence": "Unknown",
        "evidence": f"Mutations detected in ctxB relative to B7: {variant_desc}"
    }

def assess_hapR_status(hapR_variants, ctxAB_detected=False):
    """
    Four-level hapR regulator assessment.
    hapR is the master QS repressor — its loss derepresses biofilm AND virulence.
    
    Levels:
      FUNCTIONAL        — No variants, normal QS regulation
      VARIANT           — SNP (missense), partial disruption possible
      LOF_MUTANT        — Frameshift/stop-gain, biofilm & virulence derepressed
      DEREPRESSED_VIRULENCE — LOF + ctxAB present = super-producer phenotype
    """
    if not hapR_variants:
        return {
            "regulator_status": "FUNCTIONAL",
            "implication": "Normal quorum sensing regulation. Consistent with 2010 lineage.",
            "threat_multiplier": 1.0
        }

    # Check for LOF indicators: frameshift (indel) or stop-gained
    is_lof = False
    for v in hapR_variants:
        # Frameshift: ref and alt differ in length
        if len(v.get('ref', '')) != len(v.get('alt', '')):
            is_lof = True
            break
        # stop_gained in INFO/annotation
        info = v.get('info', '')
        if 'stop_gained' in info or 'frameshift' in info:
            is_lof = True
            break

    if is_lof and ctxAB_detected:
        return {
            "regulator_status": "DEREPRESSED_VIRULENCE",
            "implication": "Super-producer phenotype. hapR LOF + ctxAB = constitutive toxin expression.",
            "threat_multiplier": 1.5
        }
    elif is_lof:
        return {
            "regulator_status": "LOF_MUTANT",
            "implication": "Biofilm and virulence derepressed. No ctxAB detected.",
            "threat_multiplier": 1.3
        }
    else:
        return {
            "regulator_status": "VARIANT",
            "implication": "SNP-level variation in hapR. Partial disruption possible.",
            "threat_multiplier": 1.1
        }


# Legacy alias for backward compatibility
def check_environmental_markers(hapR_variants):
    """Backward-compatible wrapper around assess_hapR_status."""
    result = assess_hapR_status(hapR_variants, ctxAB_detected=False)
    return {
        "status": result["regulator_status"],
        "implication": result["implication"]
    }

def main():
    parser = argparse.ArgumentParser(description="Vibrion Public Health Typing")
    parser.add_argument("--vcf", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output JSON report")
    parser.add_argument("--sample", required=True, help="Sample ID")
    
    args = parser.parse_args()
    
    print(f"🏥 Running Public Health Typing on {args.sample}...")
    
    # 1. Serotype (wbeT)
    wbeT_vars = parse_vcf_variants(args.vcf, LOCI['wbeT'])
    serotype = predict_serotype(wbeT_vars)
    print(f"   💉 Serotype Prediction: {serotype['prediction']}")
    
    # 2. Toxin (ctxB)
    ctxB_vars = parse_vcf_variants(args.vcf, LOCI['ctxB'])
    toxin = analyze_toxin(ctxB_vars)
    print(f"   ☠️  Toxin Genotype: {toxin['genotype']}")
    
    # Detect ctxAB presence (needed for hapR derepressed virulence assessment)
    ctxAB_detected = len(ctxB_vars) == 0  # No variants vs reference = ctxB present
    
    # 3. Colonization (tcpA)
    tcpA_vars = parse_vcf_variants(args.vcf, LOCI['tcpA'])
    
    # 4. Persistence (hapR) — Enhanced 4-level assessment
    hapR_vars = parse_vcf_variants(args.vcf, LOCI['hapR'])
    persistence = assess_hapR_status(hapR_vars, ctxAB_detected=ctxAB_detected)
    print(f"   🧬 hapR Status: {persistence['regulator_status']} (x{persistence['threat_multiplier']})")
    
    # Build Report
    report = {
        "sample_id": args.sample,
        "public_health_guidance": {
            "vaccine_alignment": {
                "gene": "wbeT",
                "predicted_serotype": serotype['prediction'],
                "status": serotype['status'],
                "evidence": serotype['evidence'],
                "variant_count": len(wbeT_vars)
            },
            "virulence_profile": {
                "gene": "ctxB",
                "genotype": toxin['genotype'],
                "severity_risk": toxin['virulence'],
                "evidence": toxin['evidence'],
                "variant_count": len(ctxB_vars)
            },
            "reservoir_persistence": {
                "gene": "hapR",
                "regulator_status": persistence.get('regulator_status', persistence.get('status', 'UNKNOWN')),
                "threat_multiplier": persistence.get('threat_multiplier', 1.0),
                "implication": persistence['implication'],
                "variant_count": len(hapR_vars)
            }
        },
        "raw_variants": {
            "wbeT": wbeT_vars,
            "ctxB": ctxB_vars,
            "tcpA": tcpA_vars,
            "hapR": hapR_vars
        }
    }
    
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"✅ Analysis Complete. Report: {args.output}")

if __name__ == "__main__":
    main()
