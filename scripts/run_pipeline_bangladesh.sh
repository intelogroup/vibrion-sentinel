#!/bin/bash
# Run Full Pipeline for Bangladesh Blind Test (Hostile -> Kraken -> Triage -> Evo2)

# Ensure Environment Variables
if [ -z "$NVIDIA_API_KEY" ]; then
    echo "❌ Error: NVIDIA_API_KEY is not set. Please export it."
    exit 1
fi

echo "🧬 Starting Full Vibrion Sentinel Pipeline Run"
echo "   Target Sample: Bangladesh_O139_Challenge"
echo "   Steps: Hostile (Human Filter) -> Kraken2 (ID) -> Assembly -> AI Triage"

# Run Snakemake
# --use-conda: Use isolated environments for tools (including hostile)
# --cores 4: Parallelize
# -R: Force re-run for this sample
snakemake --use-conda --cores 4 -R evo2_analyze \
    --config nvidia_api_key=$NVIDIA_API_KEY \
    data/pipeline_output/Bangladesh_O139_Challenge/07_evo2/evo2_analysis.json

echo "✅ Pipeline run command sent. Check logs for progress."
