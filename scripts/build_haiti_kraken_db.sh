#!/usr/bin/env bash
# Vibrion Sentinel — Haiti custom Kraken2 database builder
#
# Builds the small specialist database referenced by the pipeline as
# data/kraken2_haiti_custom ("Specialist: Haiti Vibrio + Decoys").
# v1 seed: three V. cholerae 7PET/El Tor representatives + one E. coli decoy.
# Expand GENOMES below as the reference set grows.
#
# Requirements: kraken2-build, python3. Output ~50-200MB.
# Usage: bash scripts/build_haiti_kraken_db.sh [output_dir] [threads]

set -euo pipefail

OUT_DIR="${1:-data/kraken2_haiti_custom}"
THREADS="${2:-2}"
WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT

mkdir -p "$OUT_DIR"

# accession -> label : efetch fasta URLs (complete genomes / chromosomes)
declare -A GENOMES=(
  # 2010EL-1786 — Haiti 2010 outbreak, 7PET lineage (chr1 + chr2)
  ["CP003069.1"]="haiti_2010EL1786_chr1"
  ["CP003070.1"]="haiti_2010EL1786_chr2"
  # N16961 — El Tor reference (chr1 + chr2)
  ["AE003852.1"]="eltor_N16961_chr1"
  ["AE003853.1"]="eltor_N16961_chr2"
  # E. coli K-12 MG1655 — decoy (common contaminant / near-neighbor outgroup)
  ["U00096.3"]="decoy_ecoli_K12"
)

echo "▶ Fetching ${#GENOMES[@]} reference sequences..."
mkdir -p "$OUT_DIR/library/added"
for acc in "${!GENOMES[@]}"; do
  label="${GENOMES[$acc]}"
  echo "  $acc ($label)"
  curl -sL "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${acc}&rettype=fasta&retmode=text" \
    -o "$WORKDIR/${label}.fasta"
  # sanity: must look like fasta
  head -c 1 "$WORKDIR/${label}.fasta" | grep -q ">" || { echo "  ✗ failed to fetch $acc"; exit 1; }
  kraken2-build --add-to-library "$WORKDIR/${label}.fasta" --db "$OUT_DIR" >/dev/null
done

echo "▶ Downloading taxonomy..."
kraken2-build --download-taxonomy --db "$OUT_DIR" >/dev/null

echo "▶ Building database ($THREADS threads)..."
kraken2-build --build --db "$OUT_DIR" --threads "$THREADS"

echo "✅ Haiti custom Kraken2 DB ready: $OUT_DIR"
du -sh "$OUT_DIR"
echo ""
echo "Manifest (expand as the reference set grows):"
for acc in "${!GENOMES[@]}"; do echo "  - $acc  ${GENOMES[$acc]}"; done
