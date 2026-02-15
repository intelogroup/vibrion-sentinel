#!/bin/bash
# Rebuild Kraken2 database the proper way
# Uses NCBI taxonomy but custom genomes

set -euo pipefail

DB_PATH="/Users/kalinovdameus/Developer/Vibrion/data/kraken2_haiti_custom"
LIBRARY_PATH="/Users/kalinovdameus/Developer/Vibrion/data/kraken2_library"

echo "=========================================="
echo "Rebuilding Kraken2 Haiti Custom Database"
echo "=========================================="
echo ""

# Step 1: Download NCBI taxonomy
echo "[1/3] Downloading NCBI taxonomy..."
kraken2-build --download-taxonomy --db "${DB_PATH}" --use-ftp 2>&1 | tail -10
echo "✓ Taxonomy downloaded"
echo ""

# Step 2: Add genomes to library
echo "[2/3] Adding genomes to library..."
total=0

for category in haiti_2010_lineage novc_environmental v_mimicus other_vibrio; do
    echo "  Processing ${category}..."
    for fasta in "${LIBRARY_PATH}/${category}"/*.fasta; do
        if [ -f "${fasta}" ]; then
            genome=$(basename "${fasta}" .fasta)
            echo "    Adding ${genome}..."
            kraken2-build --add-to-library "${fasta}" --db "${DB_PATH}" 2>&1 | grep -E "Added|Skipping" || true
            ((total++))
        fi
    done
done

echo ""
echo "✓ Added ${total} genomes"
echo ""

# Step 3: Build the database
echo "[3/3] Building database (10-20 minutes)..."
echo "  Parameters: k=35, m=31, s=7, threads=8"
echo ""

kraken2-build --build \
    --db "${DB_PATH}" \
    --threads 8 \
    --kmer-len 35 \
    --minimizer-len 31 \
    --minimizer-spaces 7 \
    --max-db-size 4000000000 \
    2>&1 | tail -30

echo ""
echo "=========================================="
echo "Database Build Complete!"
echo "=========================================="
echo ""

du -sh "${DB_PATH}"
kraken2-inspect --db "${DB_PATH}" 2>&1 | head -20

echo ""
echo "Database ready at: ${DB_PATH}"
