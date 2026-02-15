#!/bin/bash
set -e

# Config
CONDA_ENV="/Users/kalinovdameus/Developer/Vibrion/.snakemake/conda/ff8213fe4b44fd04c84ecf895c15e02a_"
export PATH="$CONDA_ENV/bin:$PATH"
REF="data/references/2010EL-1786.fasta"

for SAMPLE in SRR23509888 SRR23509871; do
    echo "---------------------------------------------------"
    echo "Processing $SAMPLE..."
    
    BASE_DIR="data/pipeline_output/validation_run/$SAMPLE"
    ALIGN_DIR="$BASE_DIR/04_alignment"
    STRUCT_DIR="$BASE_DIR/05_structural"
    INPUT="$BASE_DIR/03_vibrio/${SAMPLE}_vibrio_complete.fastq.gz"
    BAM="$ALIGN_DIR/${SAMPLE}_aligned.sorted.bam"
    REPORT="$STRUCT_DIR/structural_validation.json"
    
    mkdir -p "$ALIGN_DIR"
    mkdir -p "$STRUCT_DIR"
    
    # 1. Align if needed
    if [ ! -f "$BAM" ]; then
        echo "  -> Aligning..."
        if [ ! -f "$INPUT" ]; then
             echo "  ERROR: Input $INPUT missing!"
             continue
        fi
        minimap2 -ax sr -t 8 "$REF" "$INPUT" 2>/dev/null | samtools sort -o "$BAM" -
        samtools index "$BAM"
    else
        echo "  -> BAM exists."
    fi
    
    # 2. Validate
    echo "  -> Validating Structural Variants..."
    python3 workflow/scripts/validate_structural_variants.py \
        --bam "$BAM" \
        --vcf "data/references/2010EL-1786_known_mutations.vcf" \
        --reference "$REF" \
        --output "$REPORT"
        
    echo "  -> Report: $REPORT"
done
