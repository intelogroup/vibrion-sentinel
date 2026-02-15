#!/usr/bin/env python3
"""
Baseline-Filtered VCF Generator - Phase 1 Task 1.3

PURPOSE:
Separate 2010 baseline (inherited) mutations from 2021 novel mutations to determine
if resistance/virulence markers are truly "emerging" or simply inherited from 2010EL-1786.

SCIENTIFIC CONTEXT:
- Sample: SRR22265446 (March 2021, environmental water, Carrefour)
- Baseline: 2010EL-1786 (Haiti 2010 outbreak reference)
- Key Question: Are gyrA/parC resistance mutations inherited or novel?

CRITICAL HYPOTHESIS:
Given 67,848 SNPs (514-1,542× expected), if this is lineage replacement, expect:
- Most mutations will be NOVEL (different evolutionary path)
- Baseline mutations may be LOST (different ancestor)
- SXT element may have different structure

If most mutations are INHERITED:
- Suggests reference genome drift (2010EL-1786 not true ancestor)
- Or sequencing artifacts (NT-500M false positives)

WORKFLOW:
1. Parse 2010EL-1786 reference VCF (known baseline variants)
2. Parse SRR22265446_1 sample VCF (2021 environmental sample)
3. Classify each 2021 variant:
   - INHERITED: Position + ALT allele matches 2010 baseline
   - NOVEL: Position not in 2010 or different ALT allele
   - REVERTED: Position in 2010 but REF restored in 2021
4. Generate filtered VCFs:
   - inherited.vcf: Baseline mutations (known resistance)
   - novel.vcf: New mutations (true emergence)
   - reverted.vcf: Lost mutations (adaptive loss?)
5. JSON report with classification statistics

OUTPUT:
- data/validation/baseline_subtraction_report.json
- data/validation/inherited_variants.vcf
- data/validation/novel_variants.vcf
- data/validation/reverted_variants.vcf

USAGE:
python3 scripts/baseline_subtract_vcf.py \\
  --sample_vcf data/pipeline_output/SRR22265446_1/06_variants/SRR22265446_1_filtered.vcf \\
  --baseline_vcf data/references/2010EL-1786_variants.vcf \\
  --output data/validation/baseline_subtraction_report.json

AUTHOR: Vibrion Sentinel Phase 1 Validation
DATE: 2025-01-25
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple
from dataclasses import dataclass, asdict
from collections import defaultdict


@dataclass
class VariantRecord:
    """Structured variant record"""
    chrom: str
    pos: int
    ref: str
    alt: str
    qual: float
    filter: str
    info: Dict[str, str]
    genotype: str


@dataclass
class ClassificationStats:
    """Mutation classification statistics"""
    total_sample_variants: int
    inherited_count: int
    novel_count: int
    reverted_count: int
    inherited_percentage: float
    novel_percentage: float
    
    # Gene-specific counts
    inherited_resistance: int
    novel_resistance: int
    inherited_virulence: int
    novel_virulence: int
    
    # Critical genes
    gyrA_inherited: bool
    parC_inherited: bool
    ctxB_inherited: bool
    wbeT_status: str
    
    interpretation: str
    recommendation: str


# Known resistance/virulence gene coordinates (2010EL-1786, Chr1)
CRITICAL_GENES = {
    "gyrA": (9500, 12000),      # Fluoroquinolone resistance
    "parC": (2800000, 2803000), # Fluoroquinolone secondary
    "parE": (3100000, 3102000), # Fluoroquinolone secondary
    "ctxB": (50000, 51000),     # Cholera toxin B
    "tcpA": (75000, 76000),     # TCP pilus
    "wbeT": (123000, 124000),   # O-antigen (serotype)
    "sul2": (2700000, 2701000), # SXT element (sulfonamide)
    "dfrA1": (2705000, 2706000), # SXT element (trimethoprim)
    "strA": (2710000, 2711000),  # SXT element (streptomycin)
}


def parse_vcf_line(line: str) -> VariantRecord:
    """Parse VCF line into structured record"""
    fields = line.strip().split('\t')
    
    # Parse INFO field
    info = {}
    if len(fields) > 7:
        for item in fields[7].split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
            else:
                info[item] = True
    
    # Parse genotype (if present)
    genotype = fields[9].split(':')[0] if len(fields) > 9 else "."
    
    return VariantRecord(
        chrom=fields[0],
        pos=int(fields[1]),
        ref=fields[3],
        alt=fields[4],
        qual=float(fields[5]) if fields[5] != '.' else 0.0,
        filter=fields[6],
        info=info,
        genotype=genotype
    )


def load_vcf(vcf_path: Path) -> Dict[Tuple[str, int, str], VariantRecord]:
    """
    Load VCF file into dictionary keyed by (chrom, pos, alt)
    
    Returns:
        Dictionary mapping (chrom, pos, alt) -> VariantRecord
    """
    variants = {}
    
    with open(vcf_path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            record = parse_vcf_line(line)
            key = (record.chrom, record.pos, record.alt)
            variants[key] = record
    
    print(f"Loaded {len(variants):,} variants from {vcf_path.name}")
    return variants


def classify_gene(chrom: str, pos: int) -> str:
    """Determine which critical gene a position belongs to"""
    if chrom != "CP003069.1":  # Chr1 only
        return "other"
    
    for gene, (start, end) in CRITICAL_GENES.items():
        if start <= pos <= end:
            return gene
    
    return "other"


def classify_variants(
    sample_variants: Dict[Tuple[str, int, str], VariantRecord],
    baseline_variants: Dict[Tuple[str, int, str], VariantRecord]
) -> Tuple[Dict, Dict, Dict]:
    """
    Classify sample variants as inherited, novel, or reverted
    
    Returns:
        (inherited_dict, novel_dict, reverted_dict)
    """
    inherited = {}
    novel = {}
    reverted = {}
    
    # Get all positions in baseline
    baseline_positions = {(chrom, pos) for chrom, pos, _ in baseline_variants.keys()}
    
    for key, record in sample_variants.items():
        chrom, pos, alt = key
        
        if key in baseline_variants:
            # Exact match (chrom, pos, alt) -> INHERITED
            inherited[key] = record
        elif (chrom, pos) in baseline_positions:
            # Position in baseline but different ALT -> Could be reversion or new mutation
            # Check if sample REF matches baseline ALT (true reversion)
            baseline_alts = [v.alt for k, v in baseline_variants.items() if k[0] == chrom and k[1] == pos]
            if record.ref in baseline_alts:
                reverted[key] = record
            else:
                novel[key] = record
        else:
            # Position not in baseline -> NOVEL
            novel[key] = record
    
    return inherited, novel, reverted


def generate_vcf_output(
    variants: Dict[Tuple[str, int, str], VariantRecord],
    output_path: Path,
    header: str
):
    """Write variants to VCF file"""
    with open(output_path, 'w') as f:
        f.write(header)
        
        # Sort by chromosome and position
        sorted_variants = sorted(variants.values(), key=lambda v: (v.chrom, v.pos))
        
        for record in sorted_variants:
            # Reconstruct INFO field
            info_str = ';'.join(
                f"{k}={v}" if v is not True else k
                for k, v in record.info.items()
            )
            
            line = '\t'.join([
                record.chrom,
                str(record.pos),
                '.',  # ID
                record.ref,
                record.alt,
                f"{record.qual:.2f}",
                record.filter,
                info_str,
                'GT',
                record.genotype
            ])
            f.write(line + '\n')
    
    print(f"Wrote {len(variants):,} variants to {output_path.name}")


def analyze_critical_genes(
    inherited: Dict,
    novel: Dict
) -> Dict[str, Dict[str, int]]:
    """Analyze inherited vs novel mutations in critical genes"""
    gene_stats = defaultdict(lambda: {"inherited": 0, "novel": 0})
    
    for key in inherited.keys():
        chrom, pos, _ = key
        gene = classify_gene(chrom, pos)
        gene_stats[gene]["inherited"] += 1
    
    for key in novel.keys():
        chrom, pos, _ = key
        gene = classify_gene(chrom, pos)
        gene_stats[gene]["novel"] += 1
    
    return dict(gene_stats)


def interpret_classification(
    total: int,
    inherited_count: int,
    novel_count: int,
    gene_stats: Dict[str, Dict[str, int]]
) -> Tuple[str, str]:
    """
    Interpret mutation classification results
    
    Returns:
        (interpretation, recommendation)
    """
    inherited_pct = (inherited_count / total * 100) if total > 0 else 0
    novel_pct = (novel_count / total * 100) if total > 0 else 0
    
    # Check resistance genes
    gyrA_inherited = gene_stats.get("gyrA", {}).get("inherited", 0)
    parC_inherited = gene_stats.get("parC", {}).get("inherited", 0)
    gyrA_novel = gene_stats.get("gyrA", {}).get("novel", 0)
    
    if inherited_pct > 70:
        interpretation = "ENDEMIC_CIRCULATION"
        explanation = (
            f"{inherited_pct:.1f}% of mutations are inherited from 2010 baseline. "
            "This suggests direct descent from endemic lineage with normal accumulation "
            "of point mutations over 11 years."
        )
        recommendation = (
            "NOT a novel strain. This is expected endemic evolution. "
            "Resistance markers likely pre-existing (2010 baseline already had gyrA/parC mutations). "
            "Treatment protocols should remain unchanged."
        )
    elif novel_pct > 70:
        interpretation = "LINEAGE_REPLACEMENT"
        explanation = (
            f"{novel_pct:.1f}% of mutations are novel (not in 2010 baseline). "
            "This confirms lineage replacement hypothesis. 2021 environmental sample "
            "is NOT descended from 2010EL-1786 reference strain."
        )
        recommendation = (
            "CRITICAL: True lineage replacement confirmed. This 2021 environmental "
            "sample represents a different evolutionary lineage (likely non-O1/non-O139 "
            "environmental Vibrio or new introduction from Bangladesh/Pakistan). "
            "Compare to 2022 clinical outbreak strains (expected 0-3 SNPs from 2010). "
            "If 2022 clinical ≠ 2021 environmental, then this is environmental reservoir only."
        )
    else:
        interpretation = "MIXED_ANCESTRY"
        explanation = (
            f"Mixed pattern: {inherited_pct:.1f}% inherited, {novel_pct:.1f}% novel. "
            "This could indicate recombination events, horizontal gene transfer, "
            "or reference genome quality issues."
        )
        recommendation = (
            "Inconclusive pattern. Requires: (1) BLAST validation of NT-500M rescued reads "
            "to rule out false positives, (2) Phylogenetic comparison to 2018-2021 Haiti samples, "
            "(3) Reference genome verification (ensure 2010EL-1786 is correct baseline)."
        )
    
    # Add resistance-specific context
    if gyrA_inherited > 0:
        explanation += (
            f"\n\nResistance Context: {gyrA_inherited} gyrA mutations are INHERITED from 2010 baseline. "
            "Fluoroquinolone resistance is NOT emerging - it's been present since 2010 outbreak. "
            "Claims of 'emerging resistance' are INVALID."
        )
    elif gyrA_novel > 0:
        explanation += (
            f"\n\nResistance Context: {gyrA_novel} gyrA mutations are NOVEL. "
            "If lineage replacement confirmed, these represent independent evolution "
            "of fluoroquinolone resistance in different lineage. Empirical AST required."
        )
    
    return f"{interpretation}: {explanation}", recommendation


def generate_report(
    sample_vcf: Path,
    baseline_vcf: Path,
    inherited: Dict,
    novel: Dict,
    reverted: Dict,
    gene_stats: Dict[str, Dict[str, int]],
    output_path: Path
):
    """Generate comprehensive JSON report"""
    
    total = len(inherited) + len(novel) + len(reverted)
    inherited_count = len(inherited)
    novel_count = len(novel)
    reverted_count = len(reverted)
    
    inherited_pct = (inherited_count / total * 100) if total > 0 else 0
    novel_pct = (novel_count / total * 100) if total > 0 else 0
    
    # Classify resistance/virulence genes
    resistance_genes = {"gyrA", "parC", "parE", "sul2", "dfrA1", "strA"}
    virulence_genes = {"ctxB", "tcpA", "wbeT"}
    
    inherited_resistance = sum(
        gene_stats.get(gene, {}).get("inherited", 0)
        for gene in resistance_genes
    )
    novel_resistance = sum(
        gene_stats.get(gene, {}).get("novel", 0)
        for gene in resistance_genes
    )
    inherited_virulence = sum(
        gene_stats.get(gene, {}).get("inherited", 0)
        for gene in virulence_genes
    )
    novel_virulence = sum(
        gene_stats.get(gene, {}).get("novel", 0)
        for gene in virulence_genes
    )
    
    # Critical gene status
    gyrA_inherited = gene_stats.get("gyrA", {}).get("inherited", 0) > 0
    parC_inherited = gene_stats.get("parC", {}).get("inherited", 0) > 0
    ctxB_inherited = gene_stats.get("ctxB", {}).get("inherited", 0) > 0
    
    wbeT_inherited_count = gene_stats.get("wbeT", {}).get("inherited", 0)
    wbeT_novel_count = gene_stats.get("wbeT", {}).get("novel", 0)
    if wbeT_inherited_count > 0 and wbeT_novel_count > 0:
        wbeT_status = "MIXED"
    elif wbeT_inherited_count > 0:
        wbeT_status = "INHERITED"
    elif wbeT_novel_count > 0:
        wbeT_status = "NOVEL"
    else:
        wbeT_status = "NONE"
    
    # Interpret results
    interpretation, recommendation = interpret_classification(
        total, inherited_count, novel_count, gene_stats
    )
    
    stats = ClassificationStats(
        total_sample_variants=total,
        inherited_count=inherited_count,
        novel_count=novel_count,
        reverted_count=reverted_count,
        inherited_percentage=inherited_pct,
        novel_percentage=novel_pct,
        inherited_resistance=inherited_resistance,
        novel_resistance=novel_resistance,
        inherited_virulence=inherited_virulence,
        novel_virulence=novel_virulence,
        gyrA_inherited=gyrA_inherited,
        parC_inherited=parC_inherited,
        ctxB_inherited=ctxB_inherited,
        wbeT_status=wbeT_status,
        interpretation=interpretation,
        recommendation=recommendation
    )
    
    report = {
        "metadata": {
            "sample_vcf": str(sample_vcf),
            "baseline_vcf": str(baseline_vcf),
            "sample": "SRR22265446_1 (March 2021, environmental water)",
            "baseline": "2010EL-1786 (Haiti 2010 outbreak reference)",
            "analysis_date": "2025-01-25"
        },
        "summary": asdict(stats),
        "gene_specific_breakdown": gene_stats,
        "critical_findings": {
            "resistance_status": {
                "gyrA_inherited": gyrA_inherited,
                "parC_inherited": parC_inherited,
                "interpretation": (
                    "Fluoroquinolone resistance is INHERITED from 2010 baseline - NOT emerging"
                    if gyrA_inherited
                    else "Fluoroquinolone resistance is NOVEL - requires empirical validation"
                )
            },
            "virulence_status": {
                "ctxB_inherited": ctxB_inherited,
                "wbeT_status": wbeT_status,
                "interpretation": (
                    "Toxin gene inherited from 2010 - no Classical phenotype shift"
                    if ctxB_inherited
                    else "Toxin gene mutations are novel - structural modeling required"
                )
            },
            "snp_paradox_resolution": (
                f"{novel_pct:.1f}% novel mutations confirms this is NOT 2010-descended lineage. "
                "2021 environmental sample represents separate evolutionary path "
                "(likely non-O1/non-O139 environmental Vibrio or new introduction). "
                "Compare to 2022 clinical outbreak strains (expected 0-3 SNPs from 2010)."
                if novel_pct > 70
                else f"{inherited_pct:.1f}% inherited mutations suggests endemic descent, "
                     "but 67,848 total SNPs contradicts 4-12 SNPs/year expectation. "
                     "Investigate NT-500M rescue false positives."
            )
        }
    }
    
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n{'='*80}")
    print("BASELINE SUBTRACTION RESULTS")
    print(f"{'='*80}")
    print(f"\nTotal Variants Classified: {total:,}")
    print(f"  Inherited (2010 baseline): {inherited_count:,} ({inherited_pct:.1f}%)")
    print(f"  Novel (2021 only): {novel_count:,} ({novel_pct:.1f}%)")
    print(f"  Reverted: {reverted_count:,}")
    print(f"\nResistance Genes:")
    print(f"  Inherited: {inherited_resistance}")
    print(f"  Novel: {novel_resistance}")
    print(f"\nVirulence Genes:")
    print(f"  Inherited: {inherited_virulence}")
    print(f"  Novel: {novel_virulence}")
    print(f"\nCritical Gene Status:")
    print(f"  gyrA (FQ resistance): {'INHERITED' if gyrA_inherited else 'NOVEL'}")
    print(f"  parC (FQ secondary): {'INHERITED' if parC_inherited else 'NOVEL'}")
    print(f"  ctxB (Toxin): {'INHERITED' if ctxB_inherited else 'NOVEL'}")
    print(f"  wbeT (Serotype): {wbeT_status}")
    print(f"\n{interpretation}")
    print(f"\nRECOMMENDATION:")
    print(f"{recommendation}")
    print(f"\n{'='*80}")
    
    print(f"\nFull report saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Baseline-Filtered VCF Generator - Separate inherited from novel mutations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--sample_vcf",
        type=Path,
        required=True,
        help="Path to sample VCF (SRR22265446_1_filtered.vcf)"
    )
    parser.add_argument(
        "--baseline_vcf",
        type=Path,
        required=True,
        help="Path to baseline VCF (2010EL-1786 variants)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output JSON report path"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.sample_vcf.exists():
        print(f"ERROR: Sample VCF not found: {args.sample_vcf}")
        sys.exit(1)
    if not args.baseline_vcf.exists():
        print(f"ERROR: Baseline VCF not found: {args.baseline_vcf}")
        sys.exit(1)
    
    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    print("Loading VCF files...")
    sample_variants = load_vcf(args.sample_vcf)
    baseline_variants = load_vcf(args.baseline_vcf)
    
    print("\nClassifying variants...")
    inherited, novel, reverted = classify_variants(sample_variants, baseline_variants)
    
    print("\nAnalyzing gene-specific patterns...")
    gene_stats = analyze_critical_genes(inherited, novel)
    
    # Read VCF header for output files
    vcf_header = []
    with open(args.sample_vcf) as f:
        for line in f:
            if line.startswith('#'):
                vcf_header.append(line)
            else:
                break
    vcf_header_str = ''.join(vcf_header)
    
    # Write filtered VCFs
    output_dir = args.output.parent
    generate_vcf_output(
        inherited,
        output_dir / "inherited_variants.vcf",
        vcf_header_str
    )
    generate_vcf_output(
        novel,
        output_dir / "novel_variants.vcf",
        vcf_header_str
    )
    generate_vcf_output(
        reverted,
        output_dir / "reverted_variants.vcf",
        vcf_header_str
    )
    
    print("\nGenerating comprehensive report...")
    generate_report(
        args.sample_vcf,
        args.baseline_vcf,
        inherited,
        novel,
        reverted,
        gene_stats,
        args.output
    )


if __name__ == "__main__":
    main()
