#!/bin/bash

echo "📊 NT v2-500M Download Status"
echo "=============================="
echo ""

MODEL_DIR="backend/models/nucleotide-transformer-v2-500m"
MODEL_FILE="$MODEL_DIR/model.safetensors"

if [ ! -d "$MODEL_DIR" ]; then
    echo "❌ Model directory not found"
    echo "Download may not have started yet"
    exit 1
fi

echo "📂 Files:"
ls -lh "$MODEL_DIR" | tail -n +2 | while read -r line; do
    echo "   $line"
done

echo ""

if [ -f "$MODEL_FILE" ]; then
    CURRENT_SIZE=$(stat -f%z "$MODEL_FILE" 2>/dev/null || stat -c%s "$MODEL_FILE" 2>/dev/null)
    CURRENT_MB=$((CURRENT_SIZE / 1024 / 1024))
    TARGET_MB=1920  # ~2GB expected
    
    PERCENT=$((CURRENT_MB * 100 / TARGET_MB))
    
    echo "📥 model.safetensors: ${CURRENT_MB}MB / ${TARGET_MB}MB (${PERCENT}%)"
    
    if [ $CURRENT_MB -ge $TARGET_MB ]; then
        echo "✅ Download appears complete!"
        echo ""
        echo "Test with: .venv312/bin/python backend/scripts/test_nt_v2_model.py"
    else
        echo "⏳ Still downloading..."
        echo ""
        echo "Monitor: tail -f /tmp/nt500m_v2_download.log"
    fi
else
    echo "⏳ Waiting for model.safetensors to start downloading..."
fi

echo ""
echo "Download process:"
ps aux | grep "download_nt500m_v2.sh" | grep -v grep || echo "   (not running - may have completed or failed)"
