#!/bin/bash
set -e

echo "🛡️  Vibrion Sentinel - Database Setup"
echo "======================================"

DB_DIR="data"
mkdir -p $DB_DIR/kraken2_db
mkdir -p $DB_DIR/mmseqs_db

# Check if Conda env is active
if ! command -v kraken2 &> /dev/null; then
    echo "❌ Error: Kraken2 not found. Please activate the conda environment:"
    echo "   conda activate vibrion-sentinel"
    exit 1
fi

# 1. Kraken2 Database (Standard 8GB)
echo "---------------------------------------------------"
echo "⬇️  Downloading Kraken2 Standard Database (8GB)..."
echo "   Source: BenLangmead AWS Index (2024)"
echo "---------------------------------------------------"

K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240605.tar.gz"
K2_TAR="$DB_DIR/k2_standard_08gb.tar.gz"

if [ ! -f "$DB_DIR/kraken2_db/hash.k2d" ]; then
    # Download
    curl -L -o "$K2_TAR" "$K2_URL"
    
    # Extract
    echo "📦 Extracting Kraken2 Database..."
    tar -xzf "$K2_TAR" -C "$DB_DIR/kraken2_db"
    
    # Cleanup
    rm "$K2_TAR"
    echo "✅ Kraken2 Database Ready."
else
    echo "✅ Kraken2 Database already exists. Skipping."
fi

# 2. MMseqs2 Database (SwissProt)
echo "---------------------------------------------------"
echo "⬇️  Downloading SwissProt for MMseqs2 Rescue..."
echo "---------------------------------------------------"

if [ ! -f "$DB_DIR/mmseqs_db/swissprot.dbtype" ]; then
    mmseqs databases swissprot "$DB_DIR/mmseqs_db/swissprot" tmp_mmseqs --remove-tmp-files 1
    echo "✅ MMseqs2 Database Ready."
else
    echo "✅ MMseqs2 Database already exists. Skipping."
fi

echo "======================================"
echo "🎉 Setup Complete. You can now run the pipeline."