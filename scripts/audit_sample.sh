#!/bin/bash
# audit_sample.sh - Comprehensive audit of a single sample
# Usage: ./scripts/audit_sample.sh SRR22265437

set -euo pipefail

SAMPLE=$1
OUTPUT_DIR="data/pipeline_output/${SAMPLE}"

if [ ! -d "$OUTPUT_DIR" ]; then
    echo "❌ Sample directory not found: ${OUTPUT_DIR}"
    exit 1
fi

echo ""
echo "🔍 Auditing sample: ${SAMPLE}"
echo "========================================"

# 1. Check alignment integrity
echo ""
echo "1️⃣  Alignment Integrity"
alignment_manifest="${OUTPUT_DIR}/04_alignment/alignment_manifest.json"
if [ -f "$alignment_manifest" ]; then
    bam_file="${OUTPUT_DIR}/04_alignment/${SAMPLE}_aligned.sorted.bam"
    if [ -f "$bam_file" ]; then
        bam_hash=$(jq -r '.hashes.bam_sha256' "$alignment_manifest")
        actual_hash=$(shasum -a 256 "$bam_file" | cut -d' ' -f1)
        if [ "$bam_hash" = "$actual_hash" ]; then
            echo "   ✅ BAM integrity verified"
        else
            echo "   ❌ BAM file modified! Expected: $bam_hash, Got: $actual_hash"
        fi
        
        aligner=$(jq -r '.tools.aligner' "$alignment_manifest")
        tool_version=$(jq -r ".tools.${aligner}" "$alignment_manifest")
        echo "   Aligner: ${aligner} (${tool_version})"
    else
        echo "   ⚠️  BAM file not found"
    fi
else
    echo "   ⚠️  Alignment manifest not found"
fi

# 2. Check consensus quality
echo ""
echo "2️⃣  Consensus Quality"
consensus_manifest="${OUTPUT_DIR}/09_consensus/consensus_manifest.json"
if [ -f "$consensus_manifest" ]; then
    coverage=$(jq -r '.qc_metrics.coverage_percentage' "$consensus_manifest")
    depth=$(jq -r '.qc_metrics.mean_depth' "$consensus_manifest")
    hetero=$(jq -r '.qc_metrics.heterogeneous_sites' "$consensus_manifest")
    confident=$(jq -r '.qc_metrics.confident_positions' "$consensus_manifest")
    
    echo "   Coverage: ${coverage}%"
    echo "   Mean depth: ${depth}x"
    echo "   Heterogeneous sites: ${hetero}"
    echo "   Confident positions: ${confident}"
    
    # Quality flags
    if (( $(echo "$coverage >= 95.0" | bc -l) )); then
        echo "   ✅ Coverage meets threshold (≥95%)"
    else
        echo "   ⚠️  Low coverage (<95%)"
    fi
    
    if (( $(echo "$depth >= 30.0" | bc -l) )); then
        echo "   ✅ Depth meets threshold (≥30x)"
    else
        echo "   ⚠️  Low depth (<30x)"
    fi
else
    echo "   ⚠️  Consensus manifest not found"
fi

# 3. Check polishing applied
echo ""
echo "3️⃣  Polishing Status"
polishing_manifest="${OUTPUT_DIR}/09_consensus/polishing_manifest.json"
if [ -f "$polishing_manifest" ]; then
    applied=$(jq -r '.qc_metrics.polishing_applied' "$polishing_manifest")
    polisher=$(jq -r '.qc_metrics.polisher' "$polishing_manifest")
    platform=$(jq -r '.qc_metrics.platform' "$polishing_manifest")
    global_depth=$(jq -r '.qc_metrics.global_depth' "$polishing_manifest")
    
    echo "   Platform: ${platform}"
    echo "   Polisher: ${polisher}"
    echo "   Polishing applied: ${applied}"
    echo "   Global depth: ${global_depth}x"
    
    if [ "$applied" = "true" ]; then
        echo "   ✅ Polishing was applied"
    else
        echo "   ⚠️  Polishing skipped (likely low coverage)"
    fi
    
    # Check polished genome integrity
    polished_file="${OUTPUT_DIR}/09_consensus/${SAMPLE}_polished.fasta"
    if [ -f "$polished_file" ]; then
        polished_hash=$(jq -r '.hashes.polished_sha256' "$polishing_manifest")
        actual_hash=$(shasum -a 256 "$polished_file" | cut -d' ' -f1)
        if [ "$polished_hash" = "$actual_hash" ]; then
            echo "   ✅ Polished genome integrity verified"
        else
            echo "   ❌ Polished genome modified!"
        fi
    fi
else
    echo "   ⚠️  Polishing manifest not found"
fi

# 4. Check rescue stats
echo ""
echo "4️⃣  Rescue Performance"
rescue_stats="${OUTPUT_DIR}/04_alignment/${SAMPLE}_unmapped_rescue_stats.json"
if [ -f "$rescue_stats" ]; then
    if [ -s "$rescue_stats" ] && [ "$(cat "$rescue_stats")" != "{}" ]; then
        total=$(jq -r '.total_unmapped // 0' "$rescue_stats")
        rescued=$(jq -r '.rescued_by_evo2 // 0' "$rescue_stats")
        rescue_rate=$(jq -r '.rescue_rate // 0' "$rescue_stats")
        
        echo "   Total unmapped: ${total}"
        echo "   Rescued by Evo2: ${rescued}"
        echo "   Rescue rate: ${rescue_rate}"
        
        if (( $(echo "$rescue_rate > 0.0" | bc -l) )); then
            echo "   ✅ Evo2 rescue active and effective"
        fi
    else
        echo "   ℹ️  Mapping rescue not enabled or no unmapped reads"
    fi
else
    echo "   ⚠️  Rescue stats not found"
fi

# 5. Chain of custody verification
echo ""
echo "5️⃣  Chain of Custody"
if [ -f "$alignment_manifest" ] && [ -f "$consensus_manifest" ]; then
    alignment_bam_hash=$(jq -r '.hashes.bam_sha256' "$alignment_manifest")
    consensus_bam_hash=$(jq -r '.hashes.bam_sha256' "$consensus_manifest")
    
    if [ "$alignment_bam_hash" = "$consensus_bam_hash" ]; then
        echo "   ✅ Alignment → Consensus: Chain verified"
    else
        echo "   ❌ Chain broken between alignment and consensus!"
    fi
fi

if [ -f "$consensus_manifest" ] && [ -f "$polishing_manifest" ]; then
    consensus_hash=$(jq -r '.hashes.consensus_sha256' "$consensus_manifest")
    polishing_draft_hash=$(jq -r '.hashes.draft_sha256' "$polishing_manifest")
    
    if [ "$consensus_hash" = "$polishing_draft_hash" ]; then
        echo "   ✅ Consensus → Polishing: Chain verified"
    else
        echo "   ❌ Chain broken between consensus and polishing!"
    fi
fi

# 6. Timestamp tracking
echo ""
echo "6️⃣  Processing Timeline"
if [ -f "$alignment_manifest" ]; then
    alignment_time=$(jq -r '.timestamp' "$alignment_manifest")
    echo "   Alignment: ${alignment_time}"
fi
if [ -f "$consensus_manifest" ]; then
    consensus_time=$(jq -r '.timestamp' "$consensus_manifest")
    echo "   Consensus: ${consensus_time}"
fi
if [ -f "$polishing_manifest" ]; then
    polishing_time=$(jq -r '.timestamp' "$polishing_manifest")
    echo "   Polishing: ${polishing_time}"
fi

echo ""
echo "========================================"
echo "✅ Audit complete for ${SAMPLE}"
echo ""
