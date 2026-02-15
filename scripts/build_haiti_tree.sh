#!/bin/bash
# Build and visualize Haiti-specific phylogeny (2010-2022)

set -e

echo "🇭🇹 Building Haiti Cholera Phylogeny (2010-2022)"
echo "================================================"

SAMPLE_DIR="data/pipeline_output/haiti_golden10k/10_phylogeny"
HAITI_TREE="$SAMPLE_DIR/haiti_2010_2022_tree.nwk"
HAITI_METADATA="$SAMPLE_DIR/haiti_metadata.json"
HAITI_VIZ="$SAMPLE_DIR/haiti_2010_2022_tree.png"

echo ""
echo "Step 1: Building Haiti phylogenetic tree..."
python3 workflow/scripts/build_haiti_phylogeny.py \
    --sample-dir "$SAMPLE_DIR" \
    --output "$HAITI_TREE" \
    --method auto

if [ ! -f "$HAITI_TREE" ]; then
    echo "❌ Error: Tree building failed"
    exit 1
fi

echo ""
echo "Step 2: Creating enhanced visualization..."
python3 workflow/scripts/enhanced_tree_viz_python.py \
    "$HAITI_TREE" \
    "$HAITI_METADATA" \
    "$HAITI_VIZ" \
    --layout circular \
    --title "Haiti Cholera Outbreak Evolution (2010-2022)"

if [ ! -f "$HAITI_VIZ" ]; then
    echo "❌ Error: Visualization failed"
    exit 1
fi

echo ""
echo "✅ SUCCESS! Haiti phylogeny complete"
echo ""
echo "📁 Output files:"
ls -lh "$SAMPLE_DIR"/haiti_*.{nwk,json,png,fasta} 2>/dev/null | awk '{print "  "$9" ("$5")"}'

echo ""
echo "🖼️  View tree:"
echo "   open $HAITI_VIZ"
echo ""
echo "📊 Tree structure:"
echo "   cat $HAITI_TREE"
