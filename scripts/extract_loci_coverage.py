#!/usr/bin/env python3
"""
Extract Per-Locus Coverage Matrix for Critical Surveillance Loci

Resolves the "macro coverage hides micro reality" problem.
Reports depth-of-coverage for critical virulence and resistance genes.

Usage:
    python3 extract_loci_coverage.py \
        --bam data/pipeline_output/SRR22265446_1/05_consensus/vibrio_aligned_sorted.bam \
        --output data/validation/loci_coverage_matrix.tsv
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


# Critical surveillance loci (from SURVEILLANCE_LOCI_FULL in Snakemake)
CRITICAL_LOCI = {
    'virulence': {
        'ctxB': {
            'name': 'Cholera toxin B subunit',
            'chrom': 'CP003069.1',
            'start': 50000,  # Approximate—will be refined
            'end': 51000,
            'significance': 'Toxin production, Classical vs El Tor genotype'
        },
        'tcpA': {
            'name': 'Toxin-coregulated pilus A',
            'chrom': 'CP003069.1',
            'start': 75000,
            'end': 76000,
            'significance': 'Intestinal colonization, essential virulence factor'
        }
    },
    'serotype': {
        'wbeT': {
            'name': 'O-antigen biosynthesis',
            'chrom': 'CP003069.1',
            'start': 123000,
            'end': 124000,
            'significance': 'Serotype switching (Ogawa/Inaba), immune escape'
        }
    },
    'resistance': {
        'gyrA': {
            'name': 'DNA gyrase subunit A',
            'chrom': 'CP003069.1',
            'start': 9500,
            'end': 12000,
            'significance': 'Fluoroquinolone resistance (gyrA S83I)'
        },
        'parE': {
            'name': 'Topoisomerase IV subunit B',
            'chrom': 'CP003069.1',
            'start': 3100000,
            'end': 3102000,
            'significance': 'Secondary fluoroquinolone resistance'
        }
    },
    'sxt_element': {
        'sul2': {
            'name': 'Sulfonamide resistance',
            'chrom': 'CP003069.1',
            'start': 2700000,
            'end': 2701000,
            'significance': 'SXT element - Sulfamethoxazole resistance'
        },
        'dfrA1': {
            'name': 'Trimethoprim resistance',
            'chrom': 'CP003069.1',
            'start': 2705000,
            'end': 2706000,
            'significance': 'SXT element - Trimethoprim resistance'
        },
        'strA': {
            'name': 'Streptomycin resistance',
            'chrom': 'CP003069.1',
            'start': 2710000,
            'end': 2711000,
            'significance': 'SXT element - Streptomycin phosphotransferase'
        }
    }
}


def check_samtools():
    """Verify samtools is installed."""
    try:
        subprocess.run(
            ['samtools', '--version'],
            capture_output=True,
            check=True
        )
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def extract_coverage(bam_path, chrom, start, end):
    """Extract coverage statistics for a genomic region using samtools."""
    
    # samtools depth for per-base coverage
    cmd_depth = [
        'samtools', 'depth',
        '-r', f'{chrom}:{start}-{end}',
        str(bam_path)
    ]
    
    try:
        result = subprocess.run(
            cmd_depth,
            capture_output=True,
            text=True,
            check=True
        )
        
        depths = []
        for line in result.stdout.strip().split('\n'):
            if line:
                fields = line.split('\t')
                if len(fields) >= 3:
                    depths.append(int(fields[2]))
        
        if not depths:
            return {
                'mean_depth': 0,
                'min_depth': 0,
                'max_depth': 0,
                'median_depth': 0,
                'coverage_1x': 0,
                'coverage_5x': 0,
                'coverage_10x': 0,
                'total_bases': end - start,
                'bases_covered': 0
            }
        
        depths_sorted = sorted(depths)
        mean_depth = sum(depths) / len(depths)
        median_depth = depths_sorted[len(depths) // 2]
        
        bases_1x = sum(1 for d in depths if d >= 1)
        bases_5x = sum(1 for d in depths if d >= 5)
        bases_10x = sum(1 for d in depths if d >= 10)
        
        total_bases = end - start
        
        return {
            'mean_depth': round(mean_depth, 2),
            'min_depth': min(depths),
            'max_depth': max(depths),
            'median_depth': median_depth,
            'coverage_1x': round((bases_1x / total_bases) * 100, 2),
            'coverage_5x': round((bases_5x / total_bases) * 100, 2),
            'coverage_10x': round((bases_10x / total_bases) * 100, 2),
            'total_bases': total_bases,
            'bases_covered': bases_1x
        }
    
    except subprocess.CalledProcessError as e:
        print(f"ERROR running samtools: {e}", file=sys.stderr)
        return None


def classify_confidence(mean_depth, coverage_5x):
    """Classify locus confidence based on coverage metrics."""
    
    if mean_depth >= 10 and coverage_5x >= 95:
        return 'HIGH'
    elif mean_depth >= 5 and coverage_5x >= 80:
        return 'MEDIUM'
    else:
        return 'LOW'


def analyze_all_loci(bam_path):
    """Analyze coverage for all critical loci."""
    
    results = {}
    
    for category, loci in CRITICAL_LOCI.items():
        results[category] = {}
        
        for locus_id, locus_info in loci.items():
            print(f"  Analyzing {category}/{locus_id}...", end=' ')
            
            coverage = extract_coverage(
                bam_path,
                locus_info['chrom'],
                locus_info['start'],
                locus_info['end']
            )
            
            if coverage:
                confidence = classify_confidence(
                    coverage['mean_depth'],
                    coverage['coverage_5x']
                )
                
                results[category][locus_id] = {
                    **locus_info,
                    **coverage,
                    'confidence': confidence
                }
                
                print(f"✓ (mean depth: {coverage['mean_depth']:.1f}X, confidence: {confidence})")
            else:
                print("✗ FAILED")
    
    return results


def generate_tsv_report(results, output_path):
    """Generate TSV report for easy viewing."""
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        # Header
        f.write('\t'.join([
            'Category',
            'Locus',
            'Name',
            'Chromosome',
            'Start',
            'End',
            'Mean_Depth',
            'Median_Depth',
            'Min_Depth',
            'Max_Depth',
            'Coverage_1X_%',
            'Coverage_5X_%',
            'Coverage_10X_%',
            'Confidence',
            'Significance'
        ]) + '\n')
        
        # Data rows
        for category, loci in results.items():
            for locus_id, data in loci.items():
                f.write('\t'.join([
                    category,
                    locus_id,
                    data['name'],
                    data['chrom'],
                    str(data['start']),
                    str(data['end']),
                    f"{data['mean_depth']:.2f}",
                    str(data['median_depth']),
                    str(data['min_depth']),
                    str(data['max_depth']),
                    f"{data['coverage_1x']:.2f}",
                    f"{data['coverage_5x']:.2f}",
                    f"{data['coverage_10x']:.2f}",
                    data['confidence'],
                    data['significance']
                ]) + '\n')


def generate_json_report(results, output_path):
    """Generate JSON report with interpretation."""
    
    json_path = output_path.replace('.tsv', '.json')
    
    # Calculate summary statistics
    all_confidences = []
    low_confidence_loci = []
    
    for category, loci in results.items():
        for locus_id, data in loci.items():
            all_confidences.append(data['confidence'])
            if data['confidence'] == 'LOW':
                low_confidence_loci.append({
                    'category': category,
                    'locus': locus_id,
                    'name': data['name'],
                    'mean_depth': data['mean_depth'],
                    'coverage_5x': data['coverage_5x']
                })
    
    high_conf = all_confidences.count('HIGH')
    medium_conf = all_confidences.count('MEDIUM')
    low_conf = all_confidences.count('LOW')
    total = len(all_confidences)
    
    report = {
        'summary': {
            'total_loci_analyzed': total,
            'high_confidence': high_conf,
            'medium_confidence': medium_conf,
            'low_confidence': low_conf,
            'percentage_high_confidence': round((high_conf / total) * 100, 1) if total > 0 else 0
        },
        'loci_details': results,
        'low_confidence_warnings': low_confidence_loci,
        'interpretation': {}
    }
    
    # Interpretation
    if low_conf == 0:
        report['interpretation']['status'] = 'EXCELLENT'
        report['interpretation']['message'] = 'All critical loci have adequate coverage. Variant calls are reliable.'
    elif low_conf <= 2:
        report['interpretation']['status'] = 'GOOD'
        report['interpretation']['message'] = f'{low_conf} loci have low confidence. Review specific genes before clinical claims.'
    elif low_conf <= 5:
        report['interpretation']['status'] = 'MODERATE'
        report['interpretation']['message'] = f'{low_conf} loci have low confidence. Significant gaps in critical genes—caution advised.'
    else:
        report['interpretation']['status'] = 'POOR'
        report['interpretation']['message'] = f'{low_conf} loci have low confidence. Coverage insufficient for clinical assessment.'
    
    with open(json_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    return report


def print_summary(report):
    """Print human-readable summary."""
    
    print("\n" + "="*80)
    print("PER-LOCUS COVERAGE ANALYSIS - CRITICAL SURVEILLANCE GENES")
    print("="*80 + "\n")
    
    print(f"Total Loci Analyzed: {report['summary']['total_loci_analyzed']}")
    print(f"  - HIGH Confidence:   {report['summary']['high_confidence']}")
    print(f"  - MEDIUM Confidence: {report['summary']['medium_confidence']}")
    print(f"  - LOW Confidence:    {report['summary']['low_confidence']}")
    print()
    
    print(f"Overall Assessment: {report['interpretation']['status']}")
    print(f"  {report['interpretation']['message']}")
    print()
    
    if report['low_confidence_warnings']:
        print("⚠️  LOW CONFIDENCE LOCI (Requires Review):")
        for locus in report['low_confidence_warnings']:
            print(f"  - {locus['category']}/{locus['locus']}: {locus['name']}")
            print(f"    Mean depth: {locus['mean_depth']:.1f}X, Coverage ≥5X: {locus['coverage_5x']:.1f}%")
        print()
    
    print("="*80 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description='Extract coverage statistics for critical surveillance loci'
    )
    parser.add_argument(
        '--bam',
        required=True,
        help='Path to sorted BAM file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to output TSV report'
    )
    
    args = parser.parse_args()
    
    # Check dependencies
    if not check_samtools():
        print("ERROR: samtools not found. Please install samtools.", file=sys.stderr)
        sys.exit(1)
    
    # Check input file
    bam_path = Path(args.bam)
    if not bam_path.exists():
        print(f"ERROR: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
    
    # Check BAM index
    bai_path = Path(str(bam_path) + '.bai')
    if not bai_path.exists():
        print(f"WARNING: BAM index not found. Creating index...")
        subprocess.run(['samtools', 'index', str(bam_path)], check=True)
    
    print(f"Analyzing coverage for {len([l for c in CRITICAL_LOCI.values() for l in c])} critical loci...")
    results = analyze_all_loci(bam_path)
    
    print(f"\nGenerating TSV report: {args.output}")
    generate_tsv_report(results, args.output)
    
    print(f"Generating JSON report: {args.output.replace('.tsv', '.json')}")
    report = generate_json_report(results, args.output)
    
    print_summary(report)
    
    print(f"✓ Analysis complete. Reports saved.")


if __name__ == '__main__':
    main()
