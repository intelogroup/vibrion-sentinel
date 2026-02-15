#!/bin/bash
set -e

DB_DIR="/Users/kalinovdameus/Developer/Vibrion/data/kraken2_serogroup"
THREADS=8

echo "Building Kraken2 Serogroup Database..."

# Add files to library
# Note: Kraken2 requires the files to be added one by one or via a loop if there are many
# Actually it can take multiple files but sometimes shell limits apply.
echo "Adding sequences to library..."
for f in $DB_DIR/library/added/*.fasta; do
    kraken2-build --db $DB_DIR --add-to-library "$f"
done

echo "Building the database (k=35, l=31)..."
kraken2-build --db $DB_DIR --build --threads $THREADS --kmer-len 35 --minimizer-len 31

echo "Database build complete!"
echo "To use: kraken2 --db $DB_DIR reads.fastq"
