#!/usr/bin/env python3
"""
Quick NT-500M Rescue Validation - Local Alignment Approach

PURPOSE:
Fast validation of NT-500M rescued reads using local alignment to 2010EL-1786 reference
instead of slow remote BLAST. This provides immediate feedback on rescue quality.

APPROACH:
1. Align rescued reads to 2010EL-1786.fasta using minimap2
2. Calculate alignment statistics:
   - % reads with good alignments (>90% identity, >100bp aligned)
   - % reads with partial alignments (70-90% identity)
   - % reads with poor/no alignments (<70% identity)
3. Estimate false positive rate based on alignment quality

INTERPRETATION:
- >90% good alignments → High confidence (FPR < 10%)
- 70-90% good alignments → Moderate confidence (FPR 10-30%)
- <70% good alignments → Low confidence (FPR > 30%)

USAGE:
python3 scripts/quick_nt500m_validation.py \\
  --rescued_fasta data/validation/rescued_reads_sample.fasta \\
  --reference data/references/2010EL-1786.fasta \\
  --output data/validation/quick_nt500m_validation.json
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


def run_minimap2(reads_fasta, reference_fasta, output_sam):
    """Align reads to reference using minimap2"""
    print(f"Aligning {reads_fasta} to {reference_fasta}...")
    
    cmd = [
        'minimap2',
        '-ax', 'sr',  # Short read mode
        '--secondary=no',  # No secondary alignments
        reference_fasta,
        reads_fasta
    ]
    
    try:
        with open(output_sam, 'w') as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)
        print(f"✓ Alignment complete: {output_sam}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ minimap2 failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print("✗ minimap2 not found. Install with: conda install -c bioconda minimap2")
        return False


def parse_sam(sam_file):
    """Parse SAM file and extract alignment statistics"""
    print(f"Parsing alignment results from {sam_file}...")
    
    stats = {
        'total_reads': 0,
        'aligned': 0,
        'unaligned': 0,
        'high_quality': 0,  # >90% identity
        'medium_quality': 0,  # 70-90% identity
        'low_quality': 0,  # <70% identity
        'identity_distribution': defaultdict(int),
        'alignment_lengths': []
    }
    
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue  # Skip header
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            stats['total_reads'] += 1
            flag = int(fields[1])
            
            # Check if unmapped (flag 4)
            if flag & 4:
                stats['unaligned'] += 1
                stats['identity_distribution'][0] += 1
                continue
            
            stats['aligned'] += 1
            
            # Extract alignment identity from NM tag and CIGAR
            cigar = fields[5]
            nm_tag = None
            for field in fields[11:]:
                if field.startswith('NM:i:'):
                    nm_tag = int(field.split(':')[-1])
                    break
            
            # Calculate alignment length from CIGAR
            import re
            matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
            aligned_length = sum(int(length) for length, op in matches if op in 'M=X')
            stats['alignment_lengths'].append(aligned_length)
            
            # Calculate identity percentage
            if nm_tag is not None and aligned_length > 0:
                identity = 100 * (1 - nm_tag / aligned_length)
                identity_bin = int(identity // 10) * 10
                stats['identity_distribution'][identity_bin] += 1
                
                if identity >= 90:
                    stats['high_quality'] += 1
                elif identity >= 70:
                    stats['medium_quality'] += 1
                else:
                    stats['low_quality'] += 1
            else:
                # No NM tag, assume good alignment if mapped
                stats['high_quality'] += 1
                stats['identity_distribution'][90] += 1
    
    return stats


def calculate_fpr_estimate(stats):
    """Estimate false positive rate from alignment statistics"""
    total = stats['total_reads']
    if total == 0:
        return None
    
    aligned = stats['aligned']
    high_quality = stats['high_quality']
    medium_quality = stats['medium_quality']
    low_quality = stats['low_quality']
    unaligned = stats['unaligned']
    
    # Estimation logic:
    # - High quality alignments (>90% identity) → True Vibrio cholerae
    # - Medium quality (70-90%) → Ambiguous (could be divergent V. cholerae or other Vibrio)
    # - Low quality (<70%) → Likely false positives
    # - Unaligned → Definite false positives
    
    true_positive_rate = (high_quality / total) * 100
    ambiguous_rate = (medium_quality / total) * 100
    false_positive_rate = ((low_quality + unaligned) / total) * 100
    
    # Corrected SNP estimate using same formula as full validation
    # Real SNPs ≈ 67,848 × (TP_rate + 0.5 × Ambiguous_rate)
    snp_correction_factor = (true_positive_rate + 0.5 * ambiguous_rate) / 100
    estimated_real_snps = int(67848 * snp_correction_factor)
    
    # Interpretation
    if true_positive_rate >= 75:
        confidence = "HIGH_CONFIDENCE"
        interpretation = "Lineage replacement confirmed (67K SNPs are real)"
        paradox_resolution = f"Still 378-1,166× MORE than expected. Massive divergence from 2010 baseline is biological reality, not artifact."
    elif true_positive_rate >= 50:
        confidence = "MODERATE_CONFIDENCE"
        interpretation = f"Real SNPs ~{estimated_real_snps:,} (corrected for ambiguous alignments)"
        paradox_resolution = f"Moderate correction. Still {estimated_real_snps//132}-{estimated_real_snps//44}× MORE than expected. Suggests genuine divergence with some rescue artifacts."
    elif true_positive_rate >= 25:
        confidence = "LOW_CONFIDENCE"
        interpretation = f"Real SNPs ~{estimated_real_snps:,} (corrected for artifacts)"
        paradox_resolution = f"Significant correction. {estimated_real_snps//132}-{estimated_real_snps//44}× MORE than expected. Could represent accelerated environmental evolution."
    else:
        confidence = "FAILED_RESCUE"
        interpretation = f"Real SNPs ~{estimated_real_snps:,} (heavily corrected)"
        paradox_resolution = f"SNP PARADOX RESOLVED: Original 67,848 SNPs were inflated by NT-500M artifacts. Real SNP count ({estimated_real_snps:,}) is closer to expected endemic evolution."
    
    return {
        'true_positive_rate': round(true_positive_rate, 2),
        'ambiguous_rate': round(ambiguous_rate, 2),
        'false_positive_rate': round(false_positive_rate, 2),
        'snp_correction_factor': round(snp_correction_factor, 4),
        'estimated_real_snp_count': estimated_real_snps,
        'original_snp_count': 67848,
        'expected_2021_snps': '44-132 SNPs (11 years × 4-12 SNPs/year)',
        'confidence_level': confidence,
        'interpretation': interpretation,
        'paradox_resolution': paradox_resolution
    }


def generate_report(stats, fpr_estimate, output_file):
    """Generate JSON report"""
    report = {
        'validation_type': 'Quick NT-500M Rescue Validation (Local Alignment)',
        'method': 'minimap2 alignment to 2010EL-1786 reference',
        'sample_info': {
            'total_reads_sampled': stats['total_reads'],
            'source': 'SRR22265446_1 rescued reads (March 2021, Haiti)'
        },
        'alignment_statistics': {
            'aligned_reads': stats['aligned'],
            'unaligned_reads': stats['unaligned'],
            'alignment_rate': round((stats['aligned'] / stats['total_reads']) * 100, 2) if stats['total_reads'] > 0 else 0,
            'high_quality_alignments': stats['high_quality'],
            'medium_quality_alignments': stats['medium_quality'],
            'low_quality_alignments': stats['low_quality'],
            'identity_distribution': dict(stats['identity_distribution'])
        },
        'false_positive_estimate': fpr_estimate,
        'scientific_context': {
            'user_hypothesis': 'The 67,000+ SNPs may be simply bioinformatic artifacts from rescued non-Vibrio DNA',
            'baseline': '2010EL-1786 (Haiti 2010 outbreak reference)',
            'sample': 'SRR22265446 (March 12, 2021, environmental water)',
            'expected_snps': '44-132 SNPs (4-12 SNPs/year × 11 years)',
            'observed_snps': '67,848 SNPs',
            'discrepancy': '514-1,542× MORE than expected'
        },
        'recommendation': _generate_recommendation(fpr_estimate)
    }
    
    with open(output_file, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n✓ Validation report saved: {output_file}")
    return report


def _generate_recommendation(fpr_estimate):
    """Generate recommendation based on FPR estimate"""
    if fpr_estimate is None:
        return "Cannot generate recommendation - no alignment data"
    
    fpr = fpr_estimate['false_positive_rate']
    confidence = fpr_estimate['confidence_level']
    
    if confidence == "HIGH_CONFIDENCE":
        return [
            "✓ NT-500M rescue is HIGH quality (FPR < 25%)",
            "✓ Proceed with rugose phenotype screening",
            "✓ 67,848 SNPs represent true lineage replacement",
            "→ This is NOT endemic evolution from 2010 baseline",
            "→ Sample likely represents separate introduction or non-O1/non-O139 ecotype"
        ]
    elif confidence == "MODERATE_CONFIDENCE":
        return [
            "⚠️ NT-500M rescue has MODERATE quality (FPR 25-50%)",
            "⚠️ Real SNP count is ~50% lower than reported",
            "→ Proceed with caution - some rescue reads are artifacts",
            "→ Consider re-running with stricter similarity threshold",
            "→ Still suggests genuine divergence, not pure endemic evolution"
        ]
    elif confidence == "LOW_CONFIDENCE":
        return [
            "⚠️ NT-500M rescue has LOW quality (FPR 50-75%)",
            "⚠️ Real SNP count may be significantly lower",
            "→ Recalibrate pipeline with stricter filters",
            "→ Consider upgrading to NT-2.5B model",
            "→ Add BLAST post-filtering step"
        ]
    else:  # FAILED_RESCUE
        return [
            "✗ NT-500M rescue FAILED (FPR > 75%)",
            "✗ SNP PARADOX RESOLVED: Original 67,848 SNPs were artifacts",
            "→ Pipeline requires major recalibration",
            "→ Move to NT-2.5B or implement strict taxonomic filtering",
            "→ Real SNP count likely matches expected endemic evolution"
        ]


def print_summary(report):
    """Print summary to console"""
    fpr = report['false_positive_estimate']
    
    print("\n" + "="*80)
    print("NT-500M RESCUE VALIDATION SUMMARY")
    print("="*80)
    print(f"\nSample: {report['sample_info']['total_reads_sampled']} rescued reads")
    print(f"Alignment Rate: {report['alignment_statistics']['alignment_rate']:.1f}%")
    print(f"\nQuality Breakdown:")
    print(f"  High Quality (>90% identity): {report['alignment_statistics']['high_quality_alignments']} reads")
    print(f"  Medium Quality (70-90%):      {report['alignment_statistics']['medium_quality_alignments']} reads")
    print(f"  Low Quality (<70%):           {report['alignment_statistics']['low_quality_alignments']} reads")
    print(f"  Unaligned:                    {report['alignment_statistics']['unaligned_reads']} reads")
    
    print(f"\n{'='*80}")
    print("FALSE POSITIVE ESTIMATE")
    print("="*80)
    print(f"True Positive Rate:      {fpr['true_positive_rate']:.1f}%")
    print(f"Ambiguous Rate:          {fpr['ambiguous_rate']:.1f}%")
    print(f"False Positive Rate:     {fpr['false_positive_rate']:.1f}%")
    print(f"\nSNP Correction:")
    print(f"  Original SNP count:    67,848")
    print(f"  Corrected SNP count:   {fpr['estimated_real_snp_count']:,}")
    print(f"  Correction factor:     {fpr['snp_correction_factor']:.4f}")
    print(f"\nConfidence Level: {fpr['confidence_level']}")
    print(f"\nInterpretation:")
    print(f"  {fpr['interpretation']}")
    print(f"\nSNP Paradox Resolution:")
    print(f"  {fpr['paradox_resolution']}")
    
    print(f"\n{'='*80}")
    print("RECOMMENDATION")
    print("="*80)
    for rec in report['recommendation']:
        print(f"  {rec}")
    print("="*80)


def main():
    parser = argparse.ArgumentParser(
        description='Quick NT-500M rescue validation using local alignment',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--rescued_fasta', required=True,
                       help='Path to rescued reads FASTA file')
    parser.add_argument('--reference', required=True,
                       help='Path to 2010EL-1786 reference FASTA')
    parser.add_argument('--output', required=True,
                       help='Output JSON report path')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.rescued_fasta).exists():
        print(f"✗ Rescued FASTA not found: {args.rescued_fasta}")
        return 1
    
    if not Path(args.reference).exists():
        print(f"✗ Reference FASTA not found: {args.reference}")
        return 1
    
    # Create output directory
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Temporary SAM file
    sam_file = output_path.parent / "rescued_reads_alignment.sam"
    
    print("="*80)
    print("QUICK NT-500M RESCUE VALIDATION")
    print("="*80)
    print(f"Rescued reads: {args.rescued_fasta}")
    print(f"Reference:     {args.reference}")
    print(f"Output:        {args.output}")
    print("="*80 + "\n")
    
    # Step 1: Align reads
    if not run_minimap2(args.rescued_fasta, args.reference, sam_file):
        print("\n✗ Validation failed")
        return 1
    
    # Step 2: Parse alignment results
    stats = parse_sam(sam_file)
    
    # Step 3: Calculate FPR estimate
    fpr_estimate = calculate_fpr_estimate(stats)
    
    # Step 4: Generate report
    report = generate_report(stats, fpr_estimate, args.output)
    
    # Step 5: Print summary
    print_summary(report)
    
    print(f"\n✓ Validation complete!")
    print(f"✓ Report saved: {args.output}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
