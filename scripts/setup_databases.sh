#!/bin/bash
set -e

echo "🛡️  Vibrion Sentinel - Database Setup"
echo "======================================"

DB_DIR="data"
mkdir -p $DB_DIR/kraken2_db
mkdir -p $DB_DIR/mmseqs_db

# 1. Kraken2 Database (Standard 8GB)
echo "⬇️  Downloading Kraken2 Standard Database (8GB)..."
# (Placeholder for actual URL - users typically need to download this manually or use a specific link)
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz
# tar -xzf k2_standard_08gb_20230605.tar.gz -C $DB_DIR/kraken2_db
echo "⚠️  NOTE: Please download the Kraken2 Standard 8GB database to data/kraken2_db/"

# 2. MMseqs2 Database (SwissProt)
echo "⬇️  Downloading SwissProt for MMseqs2 Rescue..."
mmseqs databases swissprot $DB_DIR/mmseqs_db/swissprot tmp_mmseqs --remove-tmp-files 1

echo "✅ Database setup complete."
