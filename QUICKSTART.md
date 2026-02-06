# Quick Start: Testing Vibrion Sentinel Pipeline

This guide shows you how to quickly test the Vibrion Sentinel pipeline with a small test sample.

## Prerequisites

- Python 3.10+
- Snakemake (will be installed automatically if missing)
- 4+ CPU cores recommended
- 8GB+ RAM recommended

## Quick Test (Dry-Run Only)

To validate the pipeline structure without running the full analysis:

```bash
# From the repository root
./test_pipeline.sh
```

This will:
- Check if Snakemake is installed
- Validate the pipeline configuration
- Show all 38 jobs that would be executed
- Complete in seconds

## Full Pipeline Execution

To actually run the pipeline on the test sample:

### Step 1: Install Dependencies

```bash
# Install Mamba (faster than Conda)
# If you don't have it, follow: https://mamba.readthedocs.io/

# Create environment
mamba env create -f environment.yml
conda activate vibrion-sentinel
```

### Step 2: Download Databases

```bash
# This will download ~10GB of reference databases
bash scripts/setup_databases.sh
```

Required databases:
- Kraken2 Standard (8GB) - for taxonomic classification
- Kraken2 Serogroup - for O1/O139 detection
- MMseqs2 SwissProt - for read rescue
- Reference genomes (Haiti 2010, 2022, etc.)

### Step 3: Run Pipeline

```bash
snakemake \
  --snakefile workflow/Snakefile \
  --configfile workflow/pipeline_test_single.yaml \
  --cores 4 \
  --use-conda
```

Expected runtime: 30-60 minutes (with databases pre-downloaded)

## Test Sample Details

- **File**: `data/raw_reads/test_sample.fastq.gz`
- **Reads**: 20,000 (minimal test dataset)
- **Sample ID**: SRR22265446
- **Format**: FASTQ (Illumina-style)

## Expected Outputs

After successful execution, check:

```bash
data/pipeline_output/SRR22265446/
├── 08_comprehensive_report/
│   └── surveillance_report.md    # Main report - START HERE
├── 09_consensus/
│   └── SRR22265446_polished.fasta  # Final genome sequence
├── 10_phylogeny/
│   └── tree.png                   # Phylogenetic tree visualization
└── 07_triage/
    └── triage_decision.json       # Risk assessment
```

## Customization

To run with your own sample:

1. Place your `.fastq.gz` files in `data/raw_reads/`
2. Update `workflow/pipeline_test_single.yaml`:
   ```yaml
   samples: ["YOUR_SAMPLE_ID"]
   ```
3. Run the pipeline

## Troubleshooting

### "Missing input files" error
- Ensure test sample exists: `ls data/raw_reads/SRR22265446_1.fastq.gz`
- Check symlink: `cd data/raw_reads && ln -s test_sample.fastq.gz SRR22265446_1.fastq.gz`

### "Database not found" error
- Run: `bash scripts/setup_databases.sh`
- Or create placeholder databases for testing

### Memory issues
- Reduce threads: Change `threads: 4` to `threads: 2` in config
- Use smaller Kraken2 database

### Conda environment issues
```bash
# Remove and recreate
conda env remove -n vibrion-sentinel
mamba env create -f environment.yml
```

## Pipeline Stages

The pipeline runs 38 jobs in sequence:

1. **QC** → Remove low-quality reads
2. **Decontamination** → Remove human DNA
3. **Classification** → Identify Vibrio cholerae
4. **Rescue** → Recover mutated reads
5. **Alignment** → Map to reference genome
6. **Consensus** → Generate final genome
7. **Variants** → Call SNPs and indels
8. **Triage** → Multi-tier risk assessment
9. **AMR** → Detect resistance genes
10. **Phylogeny** → Build evolutionary tree
11. **Report** → Generate surveillance summary

## Performance Tips

- Use `--use-conda` to isolate tool dependencies
- Use `--cores N` to parallelize (N = number of CPU cores)
- For cloud: Use `--cluster` mode with job scheduler
- Monitor with `--detailed-summary`

## Getting Help

- Documentation: `README.md`
- Pipeline test results: `PIPELINE_TEST_RESULTS.md`
- Tool inventory: `TOOL_INVENTORY.md`
- Issues: https://github.com/intelogroup/vibrion-sentinel/issues

## Next Steps

After successful test:

1. Review the surveillance report
2. Check phylogenetic placement
3. Examine AMR findings
4. Validate serogroup classification
5. Review triage risk assessment

Happy testing! 🧬
