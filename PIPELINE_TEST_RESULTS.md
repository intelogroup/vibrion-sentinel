# Vibrion Sentinel Pipeline Test Results

## Test Date: 2026-02-06

## Test Configuration
- **Sample**: SRR22265446 (test_sample.fastq.gz - 20,000 reads)
- **Config**: `workflow/pipeline_test_single.yaml`
- **Snakemake Version**: 9.16.3
- **Test Type**: Dry-run (workflow validation)

## Pipeline Validation Results

### ✅ Dry-Run Successful
The pipeline DAG (Directed Acyclic Graph) was successfully built with **38 jobs** identified.

### Pipeline Stages Identified

| Stage | Job Count | Description |
|-------|-----------|-------------|
| **0. Quality Control** | 1 | fastp_qc - Read quality filtering (Q20, 4bp window) |
| **1. Decontamination** | 1 | hostile_clean - Remove human/non-target DNA |
| **2. Classification** | 2 | kraken2_classify, kraken2_broad_audit - Taxonomic classification |
| **2b. Serogroup** | 1 | detect_serogroup - Serogroup identification |
| **2c. Phage** | 1 | phage_sentinel - Phage detection |
| **3. Read Extraction** | 2 | extract_vibrio, extract_unclassified - Separate Vibrio reads |
| **3b. Rescue** | 2 | rescue_vibrio_reads, merge_vibrio_reads - MMseqs2 read rescue |
| **4. Alignment** | 2 | align_to_reference, align_sv_supplementary - Map reads to reference |
| **5. Structural Validation** | 2 | verify_coverage_integrity, validate_structural_variants |
| **6. Consensus** | 2 | generate_consensus, polish_consensus - Generate polished genome |
| **6b. Platform Detection** | 1 | detect_platform - Detect sequencing platform |
| **7. Variant Calling** | 2 | call_variants, snpeff_annotate - SNP/indel detection |
| **8. Loci Extraction** | 2 | extract_surveillance_loci, extract_reference_loci - Extract key regions |
| **9. AMR Detection** | 2 | detect_amr, detect_amr_rgi - Antimicrobial resistance genes |
| **10. Phenotype Prediction** | 1 | predict_phenotype - Predict bacterial phenotype |
| **11. CTX Integration** | 1 | detect_ctx_integration - Cholera toxin integration |
| **12. SXT Assembly** | 1 | assemble_sxt - SXT element assembly |
| **13. Phylogeny** | 1 | cholera_seq_tree - Phylogenetic tree construction |
| **14. Triage (Multi-tier)** | 5 | sketch_sample_triage, sourmash_divergence, strain_triage, global_triage, global_reference_triage |
| **15. Local AI Triage** | 1 | local_triage - HyenaDNA local analysis |
| **16. Cloud Triage** | 1 | evo2_analyze - Evo2 cloud analysis (conditional) |
| **17. Quality Validation** | 2 | validate_checksum, forensic_confirmation |
| **18. Final Report** | 1 | generate_comprehensive_report - Surveillance report generation |

**Total Jobs**: 38

### Expected Output Files

For sample `SRR22265446`, the pipeline will generate:

```
data/pipeline_output/SRR22265446/
├── 00_fastp/
│   ├── SRR22265446_1.fastq.gz (cleaned reads)
│   ├── fastp.json
│   └── fastp.html
├── 01_decontamination/
│   └── SRR22265446_clean.fastq.gz
├── 02_kraken2/
│   ├── output.txt
│   ├── report.txt
│   └── broad_report.txt
├── 02_serogroup/
│   └── serogroup_report.json
├── 02_phage/
│   └── phage_report.json
├── 03_vibrio/
│   ├── vibrio.fastq.gz
│   ├── stats.json
│   └── unclassified.fastq.gz
├── 03_rescue/
│   ├── rescued.fastq.gz
│   └── rescue_stats.json
├── 04_alignment/
│   ├── aligned.bam
│   ├── aligned.bam.bai
│   ├── sv_supplementary.bam
│   └── coverage_integrity.json
├── 05_variants/
│   ├── SRR22265446_variants.vcf.gz
│   └── snp_report.json
├── 05_loci/
│   ├── surveillance_loci.fasta
│   └── reference_loci.fasta
├── 05_structural/
│   └── structural_validation.json
├── 06_amr/
│   ├── amr_report.json
│   └── rgi_report.json
├── 06_phenotype/
│   └── phenotype_report.json
├── 07_triage/
│   ├── tier0_sourmash.json
│   ├── triage_decision.json
│   ├── local_triage.json
│   └── global_match.json
├── 07_evo2/
│   └── evo2_analysis.json
├── 08_comprehensive_report/
│   ├── surveillance_report.md (Main forensic report)
│   └── surveillance_report.json
├── 09_consensus/
│   ├── SRR22265446_polished.fasta (Final genome)
│   ├── platform_detection.json
│   ├── qc_checksum.json
│   ├── ctx_integration.json
│   ├── sxt_assembly.json
│   └── forensic_biotype.json
└── 10_phylogeny/
    ├── tree.nwk
    └── tree.png
```

## Issues Fixed

1. **Missing Configuration Parameters**: Added `mmseqs_min_seq_id`, `kraken_broad_db`, `snpeff_db`, and `evo2` section to test config
2. **Snakefile Bug**: Fixed line 1162 in `workflow/Snakefile` where `sentinel_loci` was trying to format a config string with `{sample}` wildcard at parse time
3. **Missing Directories**: Created `data/global_references/` directory
4. **Sample Naming**: Created symlink from `test_sample.fastq.gz` to `SRR22265446_1.fastq.gz`

## To Run Full Pipeline

To execute the pipeline (not just dry-run), you would need:

### 1. Install Conda Environment
```bash
cd /home/runner/work/vibrion-sentinel/vibrion-sentinel
mamba env create -f environment.yml
conda activate vibrion-sentinel
```

### 2. Download Databases (~10GB)
```bash
bash scripts/setup_databases.sh
```

Required databases:
- Kraken2 standard database (~8GB)
- Kraken2 serogroup database
- MMseqs2 SwissProt database
- Reference genomes (Haiti 2010, 2022, etc.)

### 3. Run Pipeline
```bash
snakemake --snakefile workflow/Snakefile \
  --configfile workflow/pipeline_test_single.yaml \
  --cores 4 \
  --use-conda
```

## System Requirements

- **CPU**: 4+ cores recommended
- **RAM**: 8GB minimum, 16GB recommended
- **Storage**: ~15GB for databases + outputs
- **Time**: ~30-60 minutes for small test sample (depends on database loading)

## Conclusion

The pipeline structure is **valid and ready for execution**. The dry-run successfully identified all 38 jobs and their dependencies. The workflow follows a logical progression from raw reads to a comprehensive surveillance report with multiple quality control and triage checkpoints.

The pipeline implements a sophisticated multi-tier triage system:
- **Tier 0**: Sourmash k-mer sketching (fast pre-filter)
- **Tier 1**: HyenaDNA local AI analysis
- **Tier 2**: Evo2 cloud analysis (conditional escalation)

This design allows for cost-effective surveillance with intelligent escalation to more expensive analyses only when necessary.
