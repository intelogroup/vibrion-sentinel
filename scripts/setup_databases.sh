#!/usr/bin/env bash
# Vibrion Sentinel — Database & Reference Setup
#
# Downloads all required databases and reference genomes.
# Run once after cloning the repository.
#
# Requirements: conda activate vibrion (or mamba), wget, ~12GB free disk space
#
# Usage: bash scripts/setup_databases.sh [--skip-kraken] [--skip-mmseqs]

set -euo pipefail
cd "$(dirname "$0")/.."

SKIP_KRAKEN=false
SKIP_MMSEQS=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-kraken) SKIP_KRAKEN=true; shift ;;
        --skip-mmseqs) SKIP_MMSEQS=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

mkdir -p data/references data/global_references data/kraken2_standard_8gb \
         data/kraken2_serogroup data/mmseqs_db data/core_alignment \
         data/raw_reads data/pipeline_output

echo "============================================================"
echo " Vibrion Sentinel — Database Setup"
echo "============================================================"

# ──────────────────────────────────────────────────────────────
# 1. Core reference genomes (small, always downloaded)
# ──────────────────────────────────────────────────────────────
echo ""
echo "▶ [1/4] Downloading core reference genomes (~100MB)..."

REF_BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all"

# 2010EL-1786 — Haiti sentinel reference (7PET pandemic lineage)
if [[ ! -f data/references/2010EL-1786.fasta ]]; then
    echo "  Fetching 2010EL-1786 (CP003069.1 + CP003070.1)..."
    python3 -c "
from urllib.request import urlretrieve
import os, gzip, shutil

accessions = {
    'CP003069.1': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=CP003069.1&rettype=fasta&retmode=text',
    'CP003070.1': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=CP003070.1&rettype=fasta&retmode=text',
}
parts = []
for acc, url in accessions.items():
    print(f'    Fetching {acc}...')
    from urllib.request import urlopen
    parts.append(urlopen(url).read().decode())

with open('data/references/2010EL-1786.fasta', 'w') as f:
    f.write(''.join(parts))
print('  ✅ 2010EL-1786.fasta written')
"
else
    echo "  ✅ 2010EL-1786.fasta already present"
fi

# Haiti 2022 Resurgence reference
if [[ ! -f data/references/Haiti_2022_Resurgence.fasta ]]; then
    echo "  Fetching Haiti 2022 Resurgence (OQ860440.1)..."
    python3 -c "
from urllib.request import urlopen
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=OQ860440.1&rettype=fasta&retmode=text'
data = urlopen(url).read().decode()
open('data/references/Haiti_2022_Resurgence.fasta','w').write(data)
print('  ✅ Haiti_2022_Resurgence.fasta written')
"
else
    echo "  ✅ Haiti_2022_Resurgence.fasta already present"
fi

# ──────────────────────────────────────────────────────────────
# 2. Kraken2 standard 8GB database
# ──────────────────────────────────────────────────────────────
if [[ "$SKIP_KRAKEN" == "false" ]]; then
    echo ""
    echo "▶ [2/4] Downloading Kraken2 standard 8GB database (~8GB)..."
    echo "  Source: https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz"

    if [[ ! -f data/kraken2_standard_8gb/hash.k2d ]]; then
        wget -q --show-progress \
            "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz" \
            -O /tmp/k2_standard_8gb.tar.gz
        echo "  Extracting..."
        tar -xzf /tmp/k2_standard_8gb.tar.gz -C data/kraken2_standard_8gb/
        rm /tmp/k2_standard_8gb.tar.gz
        echo "  ✅ Kraken2 standard 8GB database ready"
    else
        echo "  ✅ Kraken2 standard 8GB already present"
    fi
else
    echo ""
    echo "⏭  [2/4] Skipping Kraken2 download (--skip-kraken)"
fi

# ──────────────────────────────────────────────────────────────
# 3. Kraken2 serogroup database (custom, lightweight ~20MB)
# ──────────────────────────────────────────────────────────────
echo ""
echo "▶ [3/4] Building Kraken2 serogroup database (~20MB)..."

if [[ ! -f data/kraken2_serogroup/hash.k2d ]]; then
    echo "  Building from O1/O139 serogroup references..."
    mkdir -p data/kraken2_serogroup/library/added
    kraken2-build --download-taxonomy --db data/kraken2_serogroup 2>/dev/null || true
    kraken2-build --add-to-library data/references/2010EL-1786.fasta \
        --db data/kraken2_serogroup 2>/dev/null || true
    kraken2-build --build --db data/kraken2_serogroup --threads 4 2>/dev/null || \
        echo "  ⚠  kraken2-build not found — activate conda env first: conda activate vibrion"
    echo "  ✅ Serogroup DB built"
else
    echo "  ✅ Kraken2 serogroup DB already present"
fi

# ──────────────────────────────────────────────────────────────
# 4. MMseqs2 SwissProt database (for unclassified read rescue)
# ──────────────────────────────────────────────────────────────
if [[ "$SKIP_MMSEQS" == "false" ]]; then
    echo ""
    echo "▶ [4/4] Downloading MMseqs2 SwissProt database (~1GB)..."

    if [[ ! -f data/mmseqs_db/swissprot ]]; then
        mmseqs databases UniProtKB/Swiss-Prot data/mmseqs_db/swissprot \
            /tmp/mmseqs_tmp --threads 4 2>/dev/null || \
            echo "  ⚠  mmseqs not found — activate conda env first: conda activate vibrion"
        echo "  ✅ MMseqs2 SwissProt ready"
    else
        echo "  ✅ MMseqs2 SwissProt already present"
    fi
else
    echo ""
    echo "⏭  [4/4] Skipping MMseqs2 download (--skip-mmseqs)"
fi

# ──────────────────────────────────────────────────────────────
# 5. Index the primary reference with BWA
# ──────────────────────────────────────────────────────────────
echo ""
echo "▶ Indexing references with BWA..."
for ref in data/references/2010EL-1786.fasta data/references/Haiti_2022_Resurgence.fasta; do
    if [[ -f "$ref" && ! -f "${ref}.bwt" ]]; then
        echo "  Indexing $(basename $ref)..."
        bwa index "$ref" 2>/dev/null || echo "  ⚠  bwa not found — activate conda env"
    fi
done
echo "  ✅ BWA index done"

echo ""
echo "============================================================"
echo " Setup complete. Disk usage:"
echo "============================================================"
du -sh data/references/ data/global_references/ data/kraken2_standard_8gb/ \
       data/mmseqs_db/ data/core_alignment/ 2>/dev/null || true
echo ""
echo "Next steps:"
echo "  1. Export your NVIDIA API key (optional, for Evo2 cloud escalation):"
echo "       export NVIDIA_API_KEY='nvapi-...'"
echo "  2. Place your FASTQ reads in data/raw_reads/<SAMPLE_ID>.fastq.gz"
echo "  3. Run the pipeline:"
echo "       bash scripts/run_pipeline.sh --config workflow/test_config.yaml --cores 8"
