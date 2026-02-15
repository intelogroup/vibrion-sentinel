#!/usr/bin/env python3
"""
NT-500M Rescue Validation - Phase 1 Task 1.4 (CRITICAL)

PURPOSE:
BLAST-validate "rescued" reads to quantify false positive rate in NT-500M rescue step.
This directly addresses the "SNP Paradox" - if 25% of rescued reads are non-Vibrio DNA,
the 67,848 SNPs may be inflated by bioinformatic artifacts.

SCIENTIFIC CONTEXT - USER'S CRITICAL INSIGHT:
"The team must validate the NT-500M rescue to ensure that the 67,000+ SNPs are not
simply bioinformatic artifacts from 'rescued' non-Vibrio DNA."

If NT-500M false positive rate is 25%:
- Real coverage improvement: ~3-5% (not 10.1%)
- Real SNPs: ~50,000 (not 67,848)
- Interpretation: Still massive, but context changes

If NT-500M false positive rate is 50-75%:
- Real coverage improvement: ~2-3%
- Real SNPs: ~17,000-34,000
- Interpretation: Could match endemic evolution with some divergence

HYPOTHESIS TO TEST:
The 67,848 SNPs include:
1. True Vibrio cholerae O1 divergence (unknown %)
2. Non-O1/non-O139 V. cholerae (environmental variants)
3. Other Vibrio species (V. mimicus, V. parahaemolyticus, V. fluvialis)
4. Complete false positives (non-Vibrio bacteria)

WORKFLOW:
1. Extract 1000 random "rescued" reads from rescue output
2. BLAST each read against NCBI RefSeq bacterial database
3. Classify top hit by taxonomy:
   - V. cholerae O1/O139 (TRUE POSITIVE - target)
   - V. cholerae non-O1/non-O139 (AMBIGUOUS - environmental)
   - Other Vibrio spp. (FALSE POSITIVE - wrong genus)
   - Non-Vibrio bacteria (FALSE POSITIVE - contamination)
4. Calculate precision/recall statistics
5. Re-estimate "true" SNP count after FPR correction

OUTPUT:
- data/validation/nt500m_rescue_validation.json
- data/validation/rescued_reads_sample.fasta (1000 reads)
- data/validation/blast_results.tsv
- data/validation/taxonomy_breakdown.json

CRITICAL FOR:
- Resolving 67,848 vs 0-3 SNP paradox
- Validating +10.1% coverage claim
- Determining if this is true lineage replacement or artifacts

USAGE:
python3 scripts/validate_nt500m_rescue.py \\
  --rescued_fastq data/pipeline_output/SRR22265446_1/03_rescue/rescued_reads.fastq \\
  --sample_size 1000 \\
  --output data/validation/nt500m_rescue_validation.json

REQUIRES:
- NCBI BLAST+ suite (blastn)
- RefSeq bacterial database (or online BLAST API)
- Biopython (for FASTA/FASTQ parsing)

AUTHOR: Vibrion Sentinel Phase 1 Validation
DATE: 2025-01-25
"""

import argparse
import json
import random
import sys
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
from dataclasses import dataclass, asdict
from collections import Counter, defaultdict


@dataclass
class BlastHit:
    """BLAST search result"""
    query_id: str
    subject_id: str
    subject_title: str
    pct_identity: float
    alignment_length: int
    evalue: float
    bitscore: float
    species: str
    genus: str
    taxonomy_classification: str


@dataclass
class ValidationStats:
    """NT-500M rescue validation statistics"""
    total_rescued_reads: int
    sample_size: int
    
    # Taxonomy classification
    true_positive_vibrio_cholerae_O1: int
    ambiguous_vibrio_cholerae_non_O1: int
    false_positive_other_vibrio: int
    false_positive_non_vibrio: int
    no_significant_hit: int
    
    # Percentages
    true_positive_rate: float
    false_positive_rate: float
    ambiguous_rate: float
    
    # Coverage correction
    original_coverage_improvement: float  # +10.1%
    corrected_coverage_improvement: float  # After FPR correction
    
    # SNP correction
    original_snp_count: int  # 67,848
    estimated_real_snp_count: int  # After FPR correction
    snp_correction_factor: float
    
    interpretation: str
    snp_paradox_resolution: str


def extract_random_reads(
    fastq_path: Path,
    sample_size: int,
    output_fasta: Path
) -> List[Tuple[str, str]]:
    """
    Extract random subset of reads from FASTQ file
    
    Returns:
        List of (read_id, sequence) tuples
    """
    print(f"Counting reads in {fastq_path.name}...")
    
    # Count total reads
    total_reads = 0
    with open(fastq_path) as f:
        for line in f:
            if line.startswith('@'):
                total_reads += 1
    
    print(f"Total rescued reads: {total_reads:,}")
    
    if sample_size > total_reads:
        sample_size = total_reads
        print(f"⚠️  Sample size reduced to {sample_size} (all available reads)")
    
    # Generate random indices
    random_indices = sorted(random.sample(range(total_reads), sample_size))
    selected_reads = []
    
    print(f"Extracting {sample_size:,} random reads...")
    
    with open(fastq_path) as f:
        read_idx = 0
        next_target = random_indices.pop(0) if random_indices else None
        
        while next_target is not None:
            line = f.readline()
            if not line:
                break
            
            if line.startswith('@'):
                if read_idx == next_target:
                    read_id = line[1:].strip().split()[0]  # Remove @ and get ID only
                    sequence = f.readline().strip()
                    f.readline()  # Skip +
                    f.readline()  # Skip quality
                    
                    selected_reads.append((read_id, sequence))
                    next_target = random_indices.pop(0) if random_indices else None
                else:
                    # Skip this read
                    f.readline()  # sequence
                    f.readline()  # +
                    f.readline()  # quality
                
                read_idx += 1
    
    # Write to FASTA for BLAST
    with open(output_fasta, 'w') as f:
        for read_id, sequence in selected_reads:
            f.write(f">{read_id}\n{sequence}\n")
    
    print(f"Wrote {len(selected_reads):,} reads to {output_fasta.name}")
    return selected_reads, total_reads


def run_blast(
    fasta_path: Path,
    output_tsv: Path,
    database: str = "nt",
    max_target_seqs: int = 1,
    use_remote: bool = True
) -> bool:
    """
    Run BLAST search against RefSeq database
    
    Args:
        fasta_path: Input FASTA file
        output_tsv: Output TSV file
        database: BLAST database (nt, refseq_genomic, etc.)
        max_target_seqs: Number of top hits to return
        use_remote: Use NCBI remote BLAST (slower but no local DB required)
    
    Returns:
        True if successful, False otherwise
    """
    outfmt = "6 qseqid sseqid stitle pident length evalue bitscore staxids"
    
    cmd = [
        "blastn",
        "-query", str(fasta_path),
        "-out", str(output_tsv),
        "-outfmt", outfmt,
        "-max_target_seqs", str(max_target_seqs),
        "-evalue", "1e-5"
    ]
    
    if use_remote:
        cmd.extend(["-remote", "-db", database])
        print(f"\n⚠️  Using NCBI remote BLAST - this will be SLOW (~10-30 minutes for 1000 reads)")
        print(f"   Consider using local BLAST database for faster results")
    else:
        cmd.extend(["-db", database])
    
    print(f"\nRunning BLAST search...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200  # 2 hour timeout
        )
        
        if result.returncode != 0:
            print(f"ERROR: BLAST failed with return code {result.returncode}")
            print(f"STDERR: {result.stderr}")
            return False
        
        print(f"✓ BLAST search complete")
        return True
    
    except subprocess.TimeoutExpired:
        print(f"ERROR: BLAST timed out after 2 hours")
        return False
    except FileNotFoundError:
        print(f"ERROR: blastn not found. Install NCBI BLAST+ suite:")
        print(f"  conda install -c bioconda blast")
        print(f"  or: brew install blast")
        return False


def parse_blast_results(
    blast_tsv: Path
) -> List[BlastHit]:
    """Parse BLAST TSV output into structured hits"""
    hits = []
    
    with open(blast_tsv) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 7:
                continue
            
            query_id = fields[0]
            subject_id = fields[1]
            subject_title = fields[2]
            pct_identity = float(fields[3])
            alignment_length = int(fields[4])
            evalue = float(fields[5])
            bitscore = float(fields[6])
            taxid = fields[7] if len(fields) > 7 else "unknown"
            
            # Extract species from subject title
            # Format: "Vibrio cholerae strain X chromosome, complete genome"
            title_lower = subject_title.lower()
            
            if "vibrio cholerae" in title_lower:
                species = "Vibrio cholerae"
                genus = "Vibrio"
                # Determine O1/O139 status
                if any(x in title_lower for x in ["o1", "o139", "serogroup o1", "serogroup o139"]):
                    classification = "TRUE_POSITIVE_O1_O139"
                elif "non-o1" in title_lower or "non-o139" in title_lower:
                    classification = "AMBIGUOUS_NON_O1_O139"
                else:
                    classification = "AMBIGUOUS_UNKNOWN_SEROGROUP"
            elif "vibrio" in title_lower:
                species = subject_title.split()[0:2]  # Get genus + species
                species = ' '.join(species) if len(species) == 2 else subject_title
                genus = "Vibrio"
                classification = "FALSE_POSITIVE_OTHER_VIBRIO"
            else:
                species = subject_title.split()[0:2]
                species = ' '.join(species) if len(species) == 2 else subject_title
                genus = species.split()[0] if species else "unknown"
                classification = "FALSE_POSITIVE_NON_VIBRIO"
            
            hit = BlastHit(
                query_id=query_id,
                subject_id=subject_id,
                subject_title=subject_title,
                pct_identity=pct_identity,
                alignment_length=alignment_length,
                evalue=evalue,
                bitscore=bitscore,
                species=species,
                genus=genus,
                taxonomy_classification=classification
            )
            hits.append(hit)
    
    print(f"Parsed {len(hits):,} BLAST hits")
    return hits


def classify_hits(
    hits: List[BlastHit]
) -> Dict[str, int]:
    """Classify BLAST hits by taxonomy"""
    classification_counts = Counter(
        hit.taxonomy_classification for hit in hits
    )
    
    return dict(classification_counts)


def calculate_statistics(
    sample_size: int,
    total_rescued: int,
    classification_counts: Dict[str, int],
    original_snp_count: int = 67848,
    original_coverage_improvement: float = 10.1
) -> ValidationStats:
    """Calculate validation statistics and corrections"""
    
    # Count classifications
    true_positive = classification_counts.get("TRUE_POSITIVE_O1_O139", 0)
    ambiguous_non_o1 = classification_counts.get("AMBIGUOUS_NON_O1_O139", 0)
    ambiguous_unknown = classification_counts.get("AMBIGUOUS_UNKNOWN_SEROGROUP", 0)
    false_positive_vibrio = classification_counts.get("FALSE_POSITIVE_OTHER_VIBRIO", 0)
    false_positive_non_vibrio = classification_counts.get("FALSE_POSITIVE_NON_VIBRIO", 0)
    no_hit = sample_size - sum(classification_counts.values())
    
    ambiguous_total = ambiguous_non_o1 + ambiguous_unknown
    
    # Calculate rates
    true_positive_rate = (true_positive / sample_size * 100) if sample_size > 0 else 0
    false_positive_rate = ((false_positive_vibrio + false_positive_non_vibrio) / sample_size * 100) if sample_size > 0 else 0
    ambiguous_rate = (ambiguous_total / sample_size * 100) if sample_size > 0 else 0
    
    # Coverage correction
    # If X% are false positives, reduce coverage improvement proportionally
    corrected_coverage_improvement = original_coverage_improvement * (true_positive_rate / 100)
    
    # SNP correction
    # Assume false positive reads contribute random SNPs (high divergence)
    # Real SNP count ≈ Original SNPs × (TP rate + 0.5 × Ambiguous rate)
    # Rationale: Non-O1/O139 V. cholerae are still Vibrio but divergent
    snp_correction_factor = (true_positive_rate + 0.5 * ambiguous_rate) / 100
    estimated_real_snp_count = int(original_snp_count * snp_correction_factor)
    
    # Interpretation
    if true_positive_rate > 75:
        interpretation = "HIGH_CONFIDENCE_RESCUE"
        explanation = (
            f"{true_positive_rate:.1f}% of rescued reads are true V. cholerae O1/O139. "
            "NT-500M rescue is performing well. The 67,848 SNPs are REAL genetic divergence, "
            "not bioinformatic artifacts. Lineage replacement hypothesis remains valid."
        )
    elif true_positive_rate > 50:
        interpretation = "MODERATE_CONFIDENCE_RESCUE"
        explanation = (
            f"{true_positive_rate:.1f}% true positives, {false_positive_rate:.1f}% false positives. "
            f"After correction, estimated real SNP count: {estimated_real_snp_count:,}. "
            "Still represents massive divergence (514-1,542× expected), confirming lineage "
            "replacement, but magnitude is lower than initially reported."
        )
    elif true_positive_rate > 25:
        interpretation = "LOW_CONFIDENCE_RESCUE"
        explanation = (
            f"{true_positive_rate:.1f}% true positives, {false_positive_rate:.1f}% false positives. "
            f"After correction, estimated real SNP count: {estimated_real_snp_count:,}. "
            "NT-500M rescue has significant false positive rate. Recommend re-analysis with "
            "stricter similarity threshold (≥0.65 instead of 0.5493)."
        )
    else:
        interpretation = "FAILED_RESCUE"
        explanation = (
            f"Only {true_positive_rate:.1f}% true V. cholerae O1/O139. "
            f"After correction, estimated real SNP count: {estimated_real_snp_count:,}. "
            "NT-500M rescue is capturing primarily non-target DNA. The 67,848 SNPs "
            "may be INFLATED by false positives. SNP paradox may be RESOLVED: "
            f"{estimated_real_snp_count:,} SNPs is closer to expected endemic evolution."
        )
    
    # SNP paradox resolution
    if estimated_real_snp_count > 50000:
        paradox_resolution = (
            f"After FPR correction: {estimated_real_snp_count:,} SNPs. "
            "Still 378-1,166× MORE than expected (44-132 SNPs for 11 years). "
            "Lineage replacement hypothesis CONFIRMED. This is NOT endemic evolution."
        )
    elif estimated_real_snp_count > 10000:
        paradox_resolution = (
            f"After FPR correction: {estimated_real_snp_count:,} SNPs. "
            "75-227× MORE than expected. Suggests partial lineage replacement or "
            "recombination with divergent environmental Vibrio lineage."
        )
    elif estimated_real_snp_count > 1000:
        paradox_resolution = (
            f"After FPR correction: {estimated_real_snp_count:,} SNPs. "
            "7-22× MORE than expected. Could represent accelerated evolution in "
            "environmental reservoir (biofilm, chlorine stress) or misaligned Chr2 variants."
        )
    else:
        paradox_resolution = (
            f"After FPR correction: {estimated_real_snp_count:,} SNPs. "
            "Within 1-2× expected range for 11 years of endemic evolution. "
            "SNP PARADOX RESOLVED: Original 67,848 SNPs were NT-500M artifacts."
        )
    
    return ValidationStats(
        total_rescued_reads=total_rescued,
        sample_size=sample_size,
        true_positive_vibrio_cholerae_O1=true_positive,
        ambiguous_vibrio_cholerae_non_O1=ambiguous_total,
        false_positive_other_vibrio=false_positive_vibrio,
        false_positive_non_vibrio=false_positive_non_vibrio,
        no_significant_hit=no_hit,
        true_positive_rate=true_positive_rate,
        false_positive_rate=false_positive_rate,
        ambiguous_rate=ambiguous_rate,
        original_coverage_improvement=original_coverage_improvement,
        corrected_coverage_improvement=corrected_coverage_improvement,
        original_snp_count=original_snp_count,
        estimated_real_snp_count=estimated_real_snp_count,
        snp_correction_factor=snp_correction_factor,
        interpretation=f"{interpretation}: {explanation}",
        snp_paradox_resolution=paradox_resolution
    )


def generate_report(
    stats: ValidationStats,
    classification_counts: Dict[str, int],
    top_species: List[Tuple[str, int]],
    output_path: Path
):
    """Generate comprehensive validation report"""
    
    report = {
        "metadata": {
            "sample": "SRR22265446_1 (March 2021, environmental water)",
            "analysis": "NT-500M Rescue Validation (BLAST taxonomy)",
            "analysis_date": "2025-01-25",
            "critical_question": "Are 67,848 SNPs real or NT-500M artifacts?"
        },
        "summary": asdict(stats),
        "taxonomy_breakdown": classification_counts,
        "top_species": [
            {"species": species, "count": count, "percentage": f"{count/stats.sample_size*100:.1f}%"}
            for species, count in top_species
        ],
        "critical_findings": {
            "snp_paradox_status": stats.snp_paradox_resolution,
            "coverage_correction": {
                "original": f"+{stats.original_coverage_improvement}%",
                "corrected": f"+{stats.corrected_coverage_improvement:.1f}%",
                "interpretation": (
                    f"Real coverage improvement is {stats.corrected_coverage_improvement:.1f}% "
                    f"(not {stats.original_coverage_improvement}%)"
                )
            },
            "snp_correction": {
                "original": stats.original_snp_count,
                "corrected": stats.estimated_real_snp_count,
                "correction_factor": f"{stats.snp_correction_factor:.2f}",
                "interpretation": (
                    f"Estimated real SNPs: {stats.estimated_real_snp_count:,} "
                    f"({stats.snp_correction_factor*100:.1f}% of original)"
                )
            },
            "lineage_replacement_hypothesis": (
                "CONFIRMED - Real SNPs still vastly exceed endemic evolution"
                if stats.estimated_real_snp_count > 10000
                else "REFUTED - SNPs within expected range after FPR correction"
            )
        },
        "recommendation": (
            "Proceed with Phase 2 phylogenetic analysis using corrected SNP estimate"
            if stats.true_positive_rate > 50
            else "Re-run NT-500M rescue with stricter threshold (≥0.65) before phylogenetic analysis"
        )
    }
    
    with open(output_path, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\n{'='*80}")
    print("NT-500M RESCUE VALIDATION RESULTS")
    print(f"{'='*80}")
    print(f"\nSample Size: {stats.sample_size:,} / {stats.total_rescued_reads:,} rescued reads")
    print(f"\nTaxonomy Classification:")
    print(f"  ✓ V. cholerae O1/O139 (TRUE POSITIVE): {stats.true_positive_vibrio_cholerae_O1} ({stats.true_positive_rate:.1f}%)")
    print(f"  ? V. cholerae non-O1/O139 (AMBIGUOUS): {stats.ambiguous_vibrio_cholerae_non_O1} ({stats.ambiguous_rate:.1f}%)")
    print(f"  ✗ Other Vibrio spp. (FALSE POSITIVE): {stats.false_positive_other_vibrio}")
    print(f"  ✗ Non-Vibrio bacteria (FALSE POSITIVE): {stats.false_positive_non_vibrio}")
    print(f"  - No significant hit: {stats.no_significant_hit}")
    print(f"\nFalse Positive Rate: {stats.false_positive_rate:.1f}%")
    print(f"\nCoverage Correction:")
    print(f"  Original: +{stats.original_coverage_improvement}%")
    print(f"  Corrected: +{stats.corrected_coverage_improvement:.1f}%")
    print(f"\nSNP Correction:")
    print(f"  Original: {stats.original_snp_count:,} SNPs")
    print(f"  Estimated Real: {stats.estimated_real_snp_count:,} SNPs")
    print(f"  Correction Factor: {stats.snp_correction_factor:.2f}")
    print(f"\n{stats.interpretation}")
    print(f"\nSNP PARADOX RESOLUTION:")
    print(f"{stats.snp_paradox_resolution}")
    print(f"\n{'='*80}")
    
    print(f"\nFull report saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="NT-500M Rescue Validation - BLAST taxonomy verification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        "--rescued_fastq",
        type=Path,
        required=True,
        help="Path to rescued reads FASTQ file"
    )
    parser.add_argument(
        "--sample_size",
        type=int,
        default=1000,
        help="Number of random reads to validate (default: 1000)"
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output JSON report path"
    )
    parser.add_argument(
        "--blast_db",
        type=str,
        default="nt",
        help="BLAST database (default: nt)"
    )
    parser.add_argument(
        "--use_local_blast",
        action="store_true",
        help="Use local BLAST database (faster, requires setup)"
    )
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.rescued_fastq.exists():
        print(f"ERROR: Rescued FASTQ not found: {args.rescued_fastq}")
        sys.exit(1)
    
    # Create output directory
    args.output.parent.mkdir(parents=True, exist_ok=True)
    
    # Set random seed for reproducibility
    random.seed(42)
    
    # Extract random reads
    sample_fasta = args.output.parent / "rescued_reads_sample.fasta"
    selected_reads, total_rescued = extract_random_reads(
        args.rescued_fastq,
        args.sample_size,
        sample_fasta
    )
    
    # Run BLAST
    blast_output = args.output.parent / "blast_results.tsv"
    success = run_blast(
        sample_fasta,
        blast_output,
        database=args.blast_db,
        use_remote=not args.use_local_blast
    )
    
    if not success:
        print("\n⚠️  BLAST search failed. Cannot complete validation.")
        print("   To continue, either:")
        print("   1. Install local BLAST database and use --use_local_blast")
        print("   2. Wait for remote BLAST to complete (may take 30+ minutes)")
        sys.exit(1)
    
    # Parse BLAST results
    print("\nParsing BLAST results...")
    hits = parse_blast_results(blast_output)
    
    if not hits:
        print("\n⚠️  No BLAST hits found. Check BLAST output manually.")
        sys.exit(1)
    
    # Classify hits
    classification_counts = classify_hits(hits)
    
    # Get top species
    species_counts = Counter(hit.species for hit in hits)
    top_species = species_counts.most_common(10)
    
    # Calculate statistics
    print("\nCalculating validation statistics...")
    stats = calculate_statistics(
        sample_size=len(hits),
        total_rescued=total_rescued,
        classification_counts=classification_counts
    )
    
    # Generate report
    generate_report(
        stats,
        classification_counts,
        top_species,
        args.output
    )


if __name__ == "__main__":
    main()
