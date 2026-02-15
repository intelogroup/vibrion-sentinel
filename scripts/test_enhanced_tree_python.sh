#!/bin/bash
# Quick test for Python-based enhanced tree visualization

set -e

echo "🌳 Testing Enhanced Phylogenetic Tree Visualization (Python)"
echo "=========================================================="

TREE_FILE="data/pipeline_output/haiti_golden10k/10_phylogeny/tree.nwk"
METADATA_FILE="data/metadata/haiti_phylogeny_metadata.json"
OUTPUT_FILE="data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced_python.png"

# Check if tree exists
if [ ! -f "$TREE_FILE" ]; then
    echo "❌ Error: Tree file not found: $TREE_FILE"
    exit 1
fi

echo "✓ Found tree file"
echo "✓ Using metadata: $METADATA_FILE"

# Run Python visualization
echo ""
echo "🎨 Generating enhanced tree visualization..."
python3 workflow/scripts/enhanced_tree_viz_python.py \
    "$TREE_FILE" \
    "$METADATA_FILE" \
    "$OUTPUT_FILE" \
    --layout circular

# Check if successful
if [ -f "$OUTPUT_FILE" ]; then
    echo ""
    echo "✅ SUCCESS! Enhanced tree visualization complete!"
    echo "📁 Output files:"
    ls -lh data/pipeline_output/haiti_golden10k/10_phylogeny/tree_enhanced_python*.png 2>/dev/null || echo "  $OUTPUT_FILE"
    echo ""
    echo "🖼️  Open with:"
    echo "   open $OUTPUT_FILE"
else
    echo "❌ Error: Visualization failed"
    exit 1
fi
