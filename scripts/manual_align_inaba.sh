#!/bin/bash
set -e

# Config
CONDA_ENV="/Users/kalinovdameus/Developer/Vibrion/.snakemake/conda/ff8213fe4b44fd04c84ecf895c15e02a_"
export PATH="$CONDA_ENV/bin:$PATH"

SAMPLE="SRR23509888"
INPUT="data/pipeline_output/validation_run/$SAMPLE/03_vibrio/${SAMPLE}_vibrio_complete.fastq.gz"
REF="data/references/2010EL-1786.fasta"
OUTDIR="data/pipeline_output/validation_run/$SAMPLE/04_alignment"
BAM="$OUTDIR/${SAMPLE}_aligned.sorted.bam"

mkdir -p "$OUTDIR"

echo "Running Minimap2..."
minimap2 -ax sr -t 8 "$REF" "$INPUT" | samtools sort -o "$BAM" -

echo "Indexing BAM..."
samtools index "$BAM"

echo "Done. BAM created at $BAM"
