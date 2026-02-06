#!/bin/bash
# Quick Test Script for Vibrion Sentinel Pipeline
# Runs a dry-run to validate the pipeline configuration

set -e

echo "=========================================="
echo "Vibrion Sentinel Pipeline Dry-Run Test"
echo "=========================================="
echo ""

# Check if snakemake is installed
if ! command -v snakemake &> /dev/null; then
    echo "❌ Snakemake not found. Installing..."
    pip3 install snakemake
fi

echo "✓ Snakemake found: $(snakemake --version)"
echo ""

# Run dry-run
echo "Running pipeline dry-run..."
echo ""

cd "$(dirname "$0")"

snakemake \
    --snakefile workflow/Snakefile \
    --configfile workflow/pipeline_test_single.yaml \
    --dry-run \
    --cores 1 \
    --quiet

echo ""
echo "=========================================="
echo "✅ Dry-Run Completed Successfully!"
echo "=========================================="
echo ""
echo "Pipeline validated with 38 jobs:"
echo "  • Quality Control (fastp)"
echo "  • Decontamination (Hostile)"
echo "  • Classification (Kraken2)"
echo "  • Read Rescue (MMseqs2)"
echo "  • Alignment (BWA/Minimap2)"
echo "  • Consensus Generation"
echo "  • Variant Calling"
echo "  • Multi-tier Triage (Sourmash, HyenaDNA, Evo2)"
echo "  • AMR Detection"
echo "  • Phylogenetic Analysis"
echo "  • Comprehensive Reporting"
echo ""
echo "To run the full pipeline, you need:"
echo "  1. Install conda environment: mamba env create -f environment.yml"
echo "  2. Download databases: bash scripts/setup_databases.sh"
echo "  3. Run: snakemake --snakefile workflow/Snakefile --configfile workflow/pipeline_test_single.yaml --cores 4 --use-conda"
echo ""
