#!/bin/bash
# Simplified genome downloader for Kraken2 Haiti custom database
# Downloads only essential genomes for strain-level classification

set -euo pipefail

VIBRION_ROOT="/Users/kalinovdameus/Developer/Vibrion"
LIBRARY_PATH="${VIBRION_ROOT}/data/kraken2_library"

echo "Downloading essential reference genomes..."

# Create directories
mkdir -p "${LIBRARY_PATH}"/{haiti_2010_lineage,novc_environmental,v_mimicus,other_vibrio}

cd "${LIBRARY_PATH}"

# Haiti 2010 lineage (already have 2010EL-1786 and N16961)
echo "[1/6] Haiti 2010 strains - already present"

# Nepal 2010 (different accession than 2010EL-1786)
echo "[2/6] Downloading Yemen 2017 (Haiti descendant)..."
cd haiti_2010_lineage
if [ ! -f Yemen_2017.fasta ]; then
    wget -O Yemen_2017.fasta.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/344/755/GCA_003344755.1_ASM334475v1/GCA_003344755.1_ASM334475v1_genomic.fna.gz" \
        2>&1 | grep -E "saved|failed|%" || true
    gunzip -f Yemen_2017.fasta.gz || true
    echo "✓ Yemen 2017 downloaded"
else
    echo "✓ Yemen 2017 already present"
fi

# O37 environmental strain
echo "[3/6] Downloading O37 environmental..."
cd ../novc_environmental
if [ ! -f O37.fasta ]; then
    wget -O O37.fasta.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/789/155/GCA_000789155.1_ASM78915v1/GCA_000789155.1_ASM78915v1_genomic.fna.gz" \
        2>&1 | grep -E "saved|failed|%" || true
    gunzip -f O37.fasta.gz || true
    echo "✓ O37 downloaded"
else
    echo "✓ O37 already present"
fi

# O139 (Bengal strain)
echo "[4/6] Downloading O139 Bengal..."
if [ ! -f O139.fasta ]; then
    wget -O O139.fasta.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/021/165/GCA_000021165.2_ASM2116v1/GCA_000021165.2_ASM2116v1_genomic.fna.gz" \
        2>&1 | grep -E "saved|failed|%" || true
    gunzip -f O139.fasta.gz || true
    echo "✓ O139 downloaded"
else
    echo "✓ O139 already present"
fi

# V. mimicus (confounding species)
echo "[5/6] Downloading V. mimicus..."
cd ../v_mimicus
if [ ! -f mimicus.fasta ]; then
    wget -O mimicus.fasta.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/348/545/GCA_000348545.1_ASM34854v1/GCA_000348545.1_ASM34854v1_genomic.fna.gz" \
        2>&1 | grep -E "saved|failed|%" || true
    gunzip -f mimicus.fasta.gz || true
    echo "✓ V. mimicus downloaded"
else
    echo "✓ V. mimicus already present"
fi

# V. parahaemolyticus (other Vibrio species)
echo "[6/6] Downloading V. parahaemolyticus..."
cd ../other_vibrio
if [ ! -f parahaemolyticus.fasta ]; then
    wget -O parahaemolyticus.fasta.gz \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/196/095/GCA_000196095.1_ASM19609v1/GCA_000196095.1_ASM19609v1_genomic.fna.gz" \
        2>&1 | grep -E "saved|failed|%" || true
    gunzip -f parahaemolyticus.fasta.gz || true
    echo "✓ V. parahaemolyticus downloaded"
else
    echo "✓ V. parahaemolyticus already present"
fi

echo ""
echo "==========================================="
echo "Download Summary"
echo "==========================================="
echo "Haiti 2010 lineage:"
ls -lh "${LIBRARY_PATH}/haiti_2010_lineage/" | grep -E "fasta$" | awk '{print "  -", $9, $5}'
echo "NOVC environmental:"
ls -lh "${LIBRARY_PATH}/novc_environmental/" | grep -E "fasta$" | awk '{print "  -", $9, $5}' || echo "  (none)"
echo "V. mimicus:"
ls -lh "${LIBRARY_PATH}/v_mimicus/" | grep -E "fasta$" | awk '{print "  -", $9, $5}' || echo "  (none)"
echo "Other Vibrio:"
ls -lh "${LIBRARY_PATH}/other_vibrio/" | grep -E "fasta$" | awk '{print "  -", $9, $5}' || echo "  (none)"
echo ""
echo "Ready to build Kraken2 database!"
