#!/usr/bin/env python3
"""
Chromosome 1 vs Chromosome 2 SNP Distribution Analysis

Resolves the 67,848 SNP paradox by separating core genome (Chr1) 
from plastic super-integron (Chr2) variants.

Expected: ~80% SNPs on Chr2 if normal background variation
          >70% SNPs on Chr1 if true lineage replacement

Usage:
    python3 chr1_vs_chr2_analysis.py \
        --vcf data/pipeline_output/SRR22265446_1/06_variants/SRR22265446_1_filtered.vcf \
        --output data/validation/chr_snp_distribution.json
"""

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def parse_vcf(vcf_path):
    """Parse VCF and extract SNP positions by chromosome."""
    snps_by_chr = defaultdict(list)
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            qual = float(fields[5]) if fields[5] != '.' else 0
            filter_status = fields[6]
            
            # Classify variant type
            if len(ref) == 1 and len(alt) == 1:
                variant_type = 'SNP'
            elif len(ref) < len(alt):
                variant_type = 'INSERTION'
            elif len(ref) > len(alt):
                variant_type = 'DELETION'
            else:
                variant_type = 'COMPLEX'
            
            snps_by_chr[chrom].append({
                'position': pos,
                'ref': ref,
                'alt': alt,
                'quality': qual,
                'filter': filter_status,
                'type': variant_type
            })
    
    return snps_by_chr


def classify_chromosomes(snps_by_chr):
    """Classify SNPs by chromosome and calculate statistics."""
    
    # V. cholerae 2010EL-1786 reference genome
    chr_info = {
        'CP003069.1': {'name': 'Chromosome 1', 'type': 'core', 'size': 3031375},
        'CP003070.1': {'name': 'Chromosome 2', 'type': 'plastic', 'size': 1048055}
    }
    
    results = {
        'chr1': {'name': 'Chromosome 1 (Core Genome)', 'snps': 0, 'indels': 0, 'total': 0, 'variants': []},
        'chr2': {'name': 'Chromosome 2 (Plastic Super-integron)', 'snps': 0, 'indels': 0, 'total': 0, 'variants': []},
        'other': {'name': 'Other contigs', 'snps': 0, 'indels': 0, 'total': 0, 'variants': []}
    }
    
    total_variants = 0
    
    for chrom, variants in snps_by_chr.items():
        if chrom == 'CP003069.1':
            category = 'chr1'
        elif chrom == 'CP003070.1':
            category = 'chr2'
        else:
            category = 'other'
        
        for variant in variants:
            total_variants += 1
            results[category]['variants'].append(variant)
            results[category]['total'] += 1
            
            if variant['type'] == 'SNP':
                results[category]['snps'] += 1
            else:
                results[category]['indels'] += 1
    
    # Calculate percentages
    if total_variants > 0:
        for cat in ['chr1', 'chr2', 'other']:
            results[cat]['percentage'] = (results[cat]['total'] / total_variants) * 100
    
    # Calculate SNP density (SNPs per Mb)
    for cat, chr_id in [('chr1', 'CP003069.1'), ('chr2', 'CP003070.1')]:
        if chr_id in chr_info:
            size_mb = chr_info[chr_id]['size'] / 1_000_000
            results[cat]['snp_density'] = results[cat]['snps'] / size_mb if size_mb > 0 else 0
    
    return results, total_variants, chr_info


def interpret_distribution(results, total_variants):
    """Interpret the SNP distribution pattern."""
    
    chr1_pct = results['chr1']['percentage']
    chr2_pct = results['chr2']['percentage']
    
    interpretation = {
        'pattern': '',
        'likelihood': '',
        'explanation': '',
        'recommendation': ''
    }
    
    if chr2_pct > 80:
        interpretation['pattern'] = 'NORMAL_BACKGROUND_VARIATION'
        interpretation['likelihood'] = 'HIGH'
        interpretation['explanation'] = (
            f"Chr2 contains {chr2_pct:.1f}% of variants. This is consistent with "
            "normal background variation in the plastic super-integron region. "
            "High SNP count is EXPECTED due to Chr2's variable nature."
        )
        interpretation['recommendation'] = (
            "Δ=0.85 likely driven by functional mutations in surveillance loci, "
            "NOT by overall SNP count. Revise interpretation to focus on core genome "
            "mutations in virulence/resistance genes."
        )
    
    elif chr1_pct > 70:
        interpretation['pattern'] = 'LINEAGE_REPLACEMENT'
        interpretation['likelihood'] = 'HIGH'
        interpretation['explanation'] = (
            f"Chr1 contains {chr1_pct:.1f}% of variants. This concentration in the "
            "core genome suggests true phylogenetic divergence or lineage replacement. "
            "67,848 SNPs is NOT normal background—investigate origin."
        )
        interpretation['recommendation'] = (
            "CRITICAL: Align to 2016-2017 Haiti isolates immediately. "
            "If phylogenetic distance >50K SNPs from recent Haiti strains, "
            "this represents new introduction (Bangladesh/Pakistan origin?). "
            "Outbreak investigation warranted."
        )
    
    elif 50 <= chr1_pct <= 70:
        interpretation['pattern'] = 'MIXED_SIGNAL'
        interpretation['likelihood'] = 'MODERATE'
        interpretation['explanation'] = (
            f"Chr1 contains {chr1_pct:.1f}% of variants, Chr2 contains {chr2_pct:.1f}%. "
            "Mixed distribution suggests both background variation and "
            "potential core genome evolution. Requires further analysis."
        )
        interpretation['recommendation'] = (
            "Perform baseline subtraction (compare to 2010EL-1786 genotype) AND "
            "align to 2016-2017 isolates. Separate inherited baseline mutations "
            "from novel variants before making clinical claims."
        )
    
    else:
        interpretation['pattern'] = 'UNEXPECTED_DISTRIBUTION'
        interpretation['likelihood'] = 'UNKNOWN'
        interpretation['explanation'] = (
            f"Chr1 contains {chr1_pct:.1f}% of variants, Chr2 contains {chr2_pct:.1f}%. "
            "Distribution does not match expected patterns. Check for "
            "alignment artifacts or contamination."
        )
        interpretation['recommendation'] = (
            "Investigate data quality: (1) Check reference genome alignment quality, "
            "(2) Validate NT-500M rescue reads (BLAST), (3) Screen for non-Vibrio contamination."
        )
    
    return interpretation


def generate_report(results, total_variants, chr_info, interpretation, output_path):
    """Generate comprehensive JSON report."""
    
    report = {
        'summary': {
            'total_variants': total_variants,
            'chr1_variants': results['chr1']['total'],
            'chr2_variants': results['chr2']['total'],
            'chr1_percentage': round(results['chr1']['percentage'], 2),
            'chr2_percentage': round(results['chr2']['percentage'], 2),
            'chr1_snp_density': round(results['chr1']['snp_density'], 2),
            'chr2_snp_density': round(results['chr2']['snp_density'], 2)
        },
        'chr1_core_genome': {
            'description': 'Highly conserved core genome (~3.0 Mb)',
            'total_variants': results['chr1']['total'],
            'snps': results['chr1']['snps'],
            'indels': results['chr1']['indels'],
            'percentage_of_total': round(results['chr1']['percentage'], 2),
            'snp_density_per_mb': round(results['chr1']['snp_density'], 2),
            'expected_snp_count_endemic': '~100-300 SNPs from 2016-2017',
            'expected_snp_count_new_lineage': '>50,000 SNPs'
        },
        'chr2_plastic_region': {
            'description': 'Plastic super-integron with ~120 gene cassettes (~1.0 Mb)',
            'total_variants': results['chr2']['total'],
            'snps': results['chr2']['snps'],
            'indels': results['chr2']['indels'],
            'percentage_of_total': round(results['chr2']['percentage'], 2),
            'snp_density_per_mb': round(results['chr2']['snp_density'], 2),
            'expected_behavior': 'High background variation expected'
        },
        'interpretation': interpretation,
        'quality_flags': []
    }
    
    # Add quality flags
    if total_variants > 100000:
        report['quality_flags'].append({
            'flag': 'VERY_HIGH_SNP_COUNT',
            'severity': 'WARNING',
            'message': f'{total_variants} variants is unusually high. Check for contamination or misalignment.'
        })
    
    if results['chr1']['percentage'] > 70:
        report['quality_flags'].append({
            'flag': 'CORE_GENOME_ENRICHMENT',
            'severity': 'CRITICAL',
            'message': 'Majority of SNPs in core genome suggests lineage replacement.'
        })
    
    if results['chr2']['snp_density'] > results['chr1']['snp_density'] * 5:
        report['quality_flags'].append({
            'flag': 'EXTREME_CHR2_VARIATION',
            'severity': 'INFO',
            'message': 'Chr2 shows expected high variation (super-integron plasticity).'
        })
    
    # Write report
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report


def print_summary(report):
    """Print human-readable summary to console."""
    
    print("\n" + "="*80)
    print("CHROMOSOME 1 vs CHROMOSOME 2 SNP DISTRIBUTION ANALYSIS")
    print("="*80 + "\n")
    
    print(f"Total Variants: {report['summary']['total_variants']:,}")
    print(f"  - Chr1 (Core Genome):         {report['summary']['chr1_variants']:,} ({report['summary']['chr1_percentage']:.1f}%)")
    print(f"  - Chr2 (Plastic Super-integron): {report['summary']['chr2_variants']:,} ({report['summary']['chr2_percentage']:.1f}%)")
    print()
    
    print(f"SNP Density:")
    print(f"  - Chr1: {report['summary']['chr1_snp_density']:,.1f} SNPs/Mb")
    print(f"  - Chr2: {report['summary']['chr2_snp_density']:,.1f} SNPs/Mb")
    print()
    
    print("INTERPRETATION:")
    print(f"  Pattern:     {report['interpretation']['pattern']}")
    print(f"  Likelihood:  {report['interpretation']['likelihood']}")
    print()
    print(f"  Explanation: {report['interpretation']['explanation']}")
    print()
    print(f"  Recommendation: {report['interpretation']['recommendation']}")
    print()
    
    if report['quality_flags']:
        print("QUALITY FLAGS:")
        for flag in report['quality_flags']:
            print(f"  [{flag['severity']}] {flag['flag']}: {flag['message']}")
        print()
    
    print("="*80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Analyze SNP distribution across V. cholerae chromosomes'
    )
    parser.add_argument(
        '--vcf',
        required=True,
        help='Path to filtered VCF file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to output JSON report'
    )
    
    args = parser.parse_args()
    
    # Check input file exists
    if not Path(args.vcf).exists():
        print(f"ERROR: VCF file not found: {args.vcf}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Parsing VCF: {args.vcf}")
    snps_by_chr = parse_vcf(args.vcf)
    
    print("Classifying variants by chromosome...")
    results, total_variants, chr_info = classify_chromosomes(snps_by_chr)
    
    print("Interpreting distribution pattern...")
    interpretation = interpret_distribution(results, total_variants)
    
    print(f"Generating report: {args.output}")
    report = generate_report(results, total_variants, chr_info, interpretation, args.output)
    
    print_summary(report)
    
    print(f"✓ Analysis complete. Report saved to: {args.output}")


if __name__ == '__main__':
    main()
