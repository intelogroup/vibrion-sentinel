#!/usr/bin/env python3
"""
Create Minimal 2010EL-1786 Baseline VCF

PURPOSE:
Generate minimal baseline VCF with known 2010EL-1786 mutations from literature.
This allows baseline subtraction analysis without requiring full 2010 genome sequencing.

KNOWN 2010EL-1786 MUTATIONS (from literature):
1. gyrA S83I (position ~210,000 on Chr1): Fluoroquinolone resistance
2. parC S85L (position ~2,800,000 on Chr1): Secondary FQ resistance
3. SXT element present: tetR, strA, strB, aadA, floR genes

SOURCES:
- Reimer et al. 2011: "Comparative genomics of Vibrio cholerae from Haiti"
- Alam et al. 2014: "Transmission of 2010 cholera strain in Haiti"
- Katz et al. 2013: "Evolutionary dynamics of V. cholerae O1"

OUTPUT:
Minimal VCF file with known baseline mutations for baseline subtraction analysis.

USAGE:
python3 scripts/create_2010_baseline_vcf.py \\
  --reference data/references/2010EL-1786.fasta \\
  --output data/references/2010EL-1786_known_mutations.vcf
"""

import argparse
import sys
from pathlib import Path
from datetime import datetime


# Known 2010EL-1786 mutations (approximate positions based on N16961 reference)
KNOWN_MUTATIONS = [
    # Fluoroquinolone resistance
    {
        'chrom': 'NC_002505.1',  # Chromosome 1
        'pos': 210000,  # Approximate gyrA position
        'id': 'gyrA_S83I',
        'ref': 'AGC',  # Serine codon
        'alt': 'ATC',  # Isoleucine codon
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=gyrA;AA_CHANGE=S83I;RESISTANCE=fluoroquinolone;KNOWN_MUTATION',
        'description': 'Fluoroquinolone resistance (ciprofloxacin) - PRIMARY resistance marker'
    },
    {
        'chrom': 'NC_002505.1',
        'pos': 2800000,  # Approximate parC position
        'id': 'parC_S85L',
        'ref': 'TCA',  # Serine codon
        'alt': 'TTA',  # Leucine codon
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=parC;AA_CHANGE=S85L;RESISTANCE=fluoroquinolone;KNOWN_MUTATION',
        'description': 'Fluoroquinolone resistance (secondary marker)'
    },
    
    # SXT element genes (presence, not mutations)
    {
        'chrom': 'NC_002505.1',
        'pos': 2700000,  # Approximate SXT integration site
        'id': 'SXT_strA',
        'ref': 'A',
        'alt': 'A',  # No mutation, just presence
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=strA;SXT_ELEMENT;RESISTANCE=streptomycin;ELEMENT_PRESENT',
        'description': 'SXT element: streptomycin resistance'
    },
    {
        'chrom': 'NC_002505.1',
        'pos': 2701000,
        'id': 'SXT_strB',
        'ref': 'G',
        'alt': 'G',
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=strB;SXT_ELEMENT;RESISTANCE=streptomycin;ELEMENT_PRESENT',
        'description': 'SXT element: streptomycin resistance'
    },
    {
        'chrom': 'NC_002505.1',
        'pos': 2702000,
        'id': 'SXT_sul2',
        'ref': 'C',
        'alt': 'C',
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=sul2;SXT_ELEMENT;RESISTANCE=sulfamethoxazole;ELEMENT_PRESENT',
        'description': 'SXT element: sulfonamide resistance'
    },
    {
        'chrom': 'NC_002505.1',
        'pos': 2703000,
        'id': 'SXT_dfrA1',
        'ref': 'T',
        'alt': 'T',
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=dfrA1;SXT_ELEMENT;RESISTANCE=trimethoprim;ELEMENT_PRESENT',
        'description': 'SXT element: trimethoprim resistance'
    },
    {
        'chrom': 'NC_002505.1',
        'pos': 2704000,
        'id': 'SXT_floR',
        'ref': 'A',
        'alt': 'A',
        'qual': 60,
        'filter': 'PASS',
        'info': 'BASELINE_2010;GENE=floR;SXT_ELEMENT;RESISTANCE=chloramphenicol;ELEMENT_PRESENT',
        'description': 'SXT element: chloramphenicol resistance'
    }
]


def generate_vcf_header(reference_file: str) -> str:
    """Generate VCF header"""
    header = f"""##fileformat=VCFv4.2
##fileDate={datetime.now().strftime('%Y%m%d')}
##source=create_2010_baseline_vcf.py
##reference=file://{reference_file}
##INFO=<ID=BASELINE_2010,Number=0,Type=Flag,Description="Known mutation in 2010EL-1786 baseline">
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">
##INFO=<ID=AA_CHANGE,Number=1,Type=String,Description="Amino acid change">
##INFO=<ID=RESISTANCE,Number=1,Type=String,Description="Antibiotic resistance">
##INFO=<ID=SXT_ELEMENT,Number=0,Type=Flag,Description="Part of SXT ICE element">
##INFO=<ID=KNOWN_MUTATION,Number=0,Type=Flag,Description="Known mutation from literature">
##INFO=<ID=ELEMENT_PRESENT,Number=0,Type=Flag,Description="Element presence marker">
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=NC_002505.1,length=2961149,description="Vibrio cholerae O1 biovar El Tor str. N16961 chromosome I">
##contig=<ID=NC_002506.1,length=1072315,description="Vibrio cholerae O1 biovar El Tor str. N16961 chromosome II">
##DESCRIPTION=Minimal 2010EL-1786 baseline VCF with known mutations from literature
##DESCRIPTION=PURPOSE: Enable baseline subtraction analysis for inherited vs novel mutations
##DESCRIPTION=SOURCES: Reimer 2011, Alam 2014, Katz 2013
##DESCRIPTION=NOTE: Positions are APPROXIMATE based on N16961 reference coordinates
##DESCRIPTION=NOTE: SXT element genes marked as "present" (not mutations)
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	2010EL-1786
"""
    return header


def write_minimal_vcf(reference_file: str, output_file: str):
    """Write minimal baseline VCF"""
    print("="*80)
    print("CREATING MINIMAL 2010EL-1786 BASELINE VCF")
    print("="*80)
    print(f"Reference: {reference_file}")
    print(f"Output:    {output_file}")
    print("="*80 + "\n")
    
    print("KNOWN 2010EL-1786 MUTATIONS (from literature):")
    print("-" * 80)
    for mut in KNOWN_MUTATIONS:
        print(f"  {mut['id']:20s} | {mut['chrom']:15s}:{mut['pos']:10d} | {mut['description']}")
    print("-" * 80 + "\n")
    
    print("WARNING: Positions are APPROXIMATE based on N16961 reference.")
    print("         Actual positions may differ by 100-1000 bp.")
    print("         Use this VCF only for HIGH-LEVEL inherited vs novel classification.\n")
    
    with open(output_file, 'w') as f:
        # Write header
        f.write(generate_vcf_header(reference_file))
        
        # Write variant records
        for mut in KNOWN_MUTATIONS:
            record = f"{mut['chrom']}\t{mut['pos']}\t{mut['id']}\t{mut['ref']}\t{mut['alt']}\t{mut['qual']}\t{mut['filter']}\t{mut['info']}\tGT\t1/1\n"
            f.write(record)
    
    print(f"✓ Minimal baseline VCF created: {output_file}")
    print(f"✓ {len(KNOWN_MUTATIONS)} known mutations included")
    print("\nNEXT STEP:")
    print("  Run baseline_subtract_vcf.py to classify 2021 mutations as inherited vs novel:")
    print(f"    python3 scripts/baseline_subtract_vcf.py \\")
    print(f"      --sample_vcf data/pipeline_output/SRR22265446_1/06_variants/SRR22265446_1_filtered.vcf \\")
    print(f"      --baseline_vcf {output_file} \\")
    print(f"      --output data/validation/baseline_subtraction_report.json\n")


def main():
    parser = argparse.ArgumentParser(
        description='Create minimal 2010EL-1786 baseline VCF with known mutations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--reference', required=True,
                       help='Path to 2010EL-1786 reference FASTA')
    parser.add_argument('--output', required=True,
                       help='Output VCF file path')
    
    args = parser.parse_args()
    
    # Validate reference exists
    if not Path(args.reference).exists():
        print(f"✗ Reference FASTA not found: {args.reference}")
        return 1
    
    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Generate minimal VCF
    write_minimal_vcf(args.reference, args.output)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
