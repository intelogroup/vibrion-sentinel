#!/bin/bash
# Test script for enhanced tree visualization

echo "🌳 Testing Enhanced Phylogenetic Tree Visualization"
echo "=================================================="

TREE_FILE="data/pipeline_output/haiti_golden10k/10_phylogeny/tree.nwk"
METADATA_FILE="data/metadata/haiti_phylogeny_metadata.json"
OUTPUT_FILE="data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced.png"

# Check if tree exists
if [ ! -f "$TREE_FILE" ]; then
    echo "❌ Error: Tree file not found: $TREE_FILE"
    echo "Run the pipeline first: snakemake --configfile workflow/test_config.yaml --cores 4"
    exit 1
fi

# Check if metadata exists
if [ ! -f "$METADATA_FILE" ]; then
    echo "⚠️  Warning: Metadata file not found: $METADATA_FILE"
    echo "Will use auto-generated metadata from tip labels"
fi

# Activate conda environment
echo "📦 Activating phylogeny environment..."
if conda env list | grep -q "phylogeny-env"; then
    eval "$(conda shell.bash hook)"
    conda activate phylogeny-env
else
    echo "⚠️  Creating phylogeny-env (first time setup)..."
    conda env create -f workflow/envs/phylogeny.yaml
    conda activate phylogeny-env
fi

# Run enhanced visualization
echo "🎨 Generating enhanced tree visualization..."
Rscript workflow/scripts/enhanced_tree_viz.R \
    "$TREE_FILE" \
    "$METADATA_FILE" \
    "$OUTPUT_FILE" \
    14 \
    10

# Check if successful
if [ -f "$OUTPUT_FILE" ]; then
    echo "✅ Enhanced tree visualization complete!"
    echo "📁 Output files:"
    ls -lh data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced*.png
    echo ""
    echo "🖼️  Open with: open $OUTPUT_FILE"
else
    echo "❌ Error: Visualization failed"
    exit 1
fi
