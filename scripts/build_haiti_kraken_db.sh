#!/usr/bin/env bash
# Vibrion Sentinel — Haiti custom Kraken2 database builder
#
# Builds the small specialist database referenced by the pipeline as
# data/kraken2_haiti_custom ("Specialist: Haiti Vibrio + near-neighbors +
# decoys"). Scoped for cultured bacterial isolate WGS (no human genome):
# the realistic off-target space here is other Vibrio/Aeromonas species
# documented as biochemically misidentified as V. cholerae in the lab
# (https://pmc.ncbi.nlm.nih.gov/articles/PMC5857081/), plus common WGS
# lab contaminants. Without them in the DB, a non-cholerae read has
# nothing to match and falls into "unclassified" -- which extract_vibrio
# deliberately keeps (to rescue divergent cholerae strains), so a wrong-
# species sample was sailing through undetected until it failed later,
# expensively, at the mapping stage.
# Expand GENOMES below as the reference set grows.
#
# Requirements: kraken2-build, python3. Output ~650-700MB as of the 25
# sequences below (verified: `du -sh` after a real build) -- re-check
# this comment if the reference set grows meaningfully.
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
  # Near-neighbor Vibrio/Aeromonas species documented as biochemically
  # misidentified as V. cholerae in the lab (see script header comment).
  # V. parahaemolyticus RIMD 2210633 (chr1+chr2)
  ["BA000031.2"]="vparahaemolyticus_chr1"
  ["BA000032.2"]="vparahaemolyticus_chr2"
  # V. vulnificus CMCP6 (chr1+chr2)
  ["AE016795.3"]="vvulnificus_chr1"
  ["AE016796.2"]="vvulnificus_chr2"
  # V. mimicus SCCF01 (chr1+chr2)
  ["CP016383.1"]="vmimicus_chr1"
  ["CP016384.1"]="vmimicus_chr2"
  # V. fluvialis F8658 (chr1+chr2+plasmid)
  ["CP046855.1"]="vfluvialis_chr1"
  ["CP046856.1"]="vfluvialis_chr2"
  ["CP046854.1"]="vfluvialis_plasmid"
  # V. metschnikovii CIP 69.14 — draft (no complete genome exists); largest
  # 9 of 11 scaffolds, ~98.5% of the ~3.8 Mb assembly (GCA_000176155.1)
  ["ACZO01000001.1"]="vmetschnikovii_s1"
  ["ACZO01000002.1"]="vmetschnikovii_s2"
  ["ACZO01000003.1"]="vmetschnikovii_s3"
  ["ACZO01000004.1"]="vmetschnikovii_s4"
  ["ACZO01000005.1"]="vmetschnikovii_s5"
  ["ACZO01000006.1"]="vmetschnikovii_s6"
  ["ACZO01000007.1"]="vmetschnikovii_s7"
  ["ACZO01000008.1"]="vmetschnikovii_s8"
  ["ACZO01000009.1"]="vmetschnikovii_s9"
  # A. hydrophila ATCC 7966 (single chromosome)
  ["CP000462.1"]="ahydrophila"
  # A. caviae NUITM-VA2 (single chromosome)
  ["AP025280.1"]="acaviae"
)

# accession -> NCBI taxid (embedded in fasta headers as |kraken:taxid| so that
# kraken2-build needs no accession2taxid map files). Verified directly
# against NCBI's taxonomy API, not a search-engine summary -- one such
# summary incorrectly gave A. caviae the taxid 562, which is E. coli's.
declare -A TAXIDS=(
  ["CP003069.1"]=666 ["CP003070.1"]=666
  ["AE003852.1"]=666 ["AE003853.1"]=666
  ["U00096.3"]=562
  ["BA000031.2"]=670 ["BA000032.2"]=670
  ["AE016795.3"]=672 ["AE016796.2"]=672
  ["CP016383.1"]=674 ["CP016384.1"]=674
  ["CP046855.1"]=676 ["CP046856.1"]=676 ["CP046854.1"]=676
  ["ACZO01000001.1"]=28172 ["ACZO01000002.1"]=28172 ["ACZO01000003.1"]=28172
  ["ACZO01000004.1"]=28172 ["ACZO01000005.1"]=28172 ["ACZO01000006.1"]=28172
  ["ACZO01000007.1"]=28172 ["ACZO01000008.1"]=28172 ["ACZO01000009.1"]=28172
  ["CP000462.1"]=644
  ["AP025280.1"]=648
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
  # embed taxid: lets kraken2-build skip the multi-GB accession2taxid maps
  sed -i -E "1s/^(>\\S+)/\\1|kraken:taxid|${TAXIDS[$acc]}/" "$WORKDIR/${label}.fasta"
  kraken2-build --add-to-library "$WORKDIR/${label}.fasta" --db "$OUT_DIR" >/dev/null
done

echo "▶ Downloading taxonomy (HTTPS: rsync is blocked on CI runners)..."
mkdir -p "$OUT_DIR/taxonomy"
curl -sL "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" -o "$WORKDIR/taxdump.tar.gz"
tar -xzf "$WORKDIR/taxdump.tar.gz" -C "$OUT_DIR/taxonomy"

echo "▶ Building database ($THREADS threads)..."
kraken2-build --build --db "$OUT_DIR" --threads "$THREADS"

echo "✅ Haiti custom Kraken2 DB ready: $OUT_DIR"
du -sh "$OUT_DIR"
echo ""
echo "Manifest (expand as the reference set grows):"
for acc in "${!GENOMES[@]}"; do echo "  - $acc  ${GENOMES[$acc]}"; done
