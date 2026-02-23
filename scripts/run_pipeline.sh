#!/usr/bin/env bash
# Vibrion Pipeline Runner
#
# Usage:
#   bash scripts/run_pipeline.sh --config workflow/my_config.yaml [--cores N] [--dry-run]
#
# Requirements:
#   - NVIDIA_API_KEY env var (optional; evo2_analyze runs in local mode if unset)
#   - Conda environments pre-built: snakemake --use-conda --conda-create-envs-only
#
# Lessons applied:
#   - nohup + tee: pipeline survives terminal closure
#   - --rerun-incomplete: recovers cleanly from mid-run interruptions
#   - --keep-going: partial rule failures don't abort the whole run
#   - nvidia_api_key resolved here, never put ${VAR} literals in YAML
#     (Snakemake expands them as wildcards in params: blocks → WildcardError)

set -euo pipefail
cd "$(dirname "$0")/.."

# ---------- Defaults ----------
CONFIG=""
CORES=4
DRY_RUN=false

# ---------- Argument parsing ----------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)    CONFIG="$2"; shift 2 ;;
        --cores)     CORES="$2";  shift 2 ;;
        --dry-run)   DRY_RUN=true; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

[[ -z "$CONFIG" ]] && { echo "ERROR: --config <path> is required"; exit 1; }
[[ -f "$CONFIG" ]] || { echo "ERROR: config not found: $CONFIG"; exit 1; }

# ---------- API key resolution ----------
# Never embed ${NVIDIA_API_KEY} in YAML — pass it here via --config override
RESOLVED_API_KEY="${NVIDIA_API_KEY:-DUMMY_KEY_NOT_SET}"
if [[ "$RESOLVED_API_KEY" == "DUMMY_KEY_NOT_SET" ]]; then
    echo "⚠  NVIDIA_API_KEY not set — evo2_analyze will run in local mode"
fi

# ---------- Derive sample name for log ----------
SAMPLE=$(grep "^samples:" -A1 "$CONFIG" | tail -1 | tr -d ' -[]"' | cut -d'#' -f1)
[[ -z "$SAMPLE" ]] && SAMPLE="pipeline"
mkdir -p logs
LOG="logs/${SAMPLE}_run.log"

# ---------- Dry-run mode ----------
if [[ "$DRY_RUN" == "true" ]]; then
    echo "🔍 Dry-run for: $CONFIG"
    snakemake \
        --snakefile workflow/Snakefile \
        --configfile "$CONFIG" \
        --config nvidia_api_key="${RESOLVED_API_KEY}" \
        --dry-run \
        --cores "$CORES"
    exit 0
fi

# ---------- Full run ----------
echo "🚀 Launching pipeline"
echo "   Config:  $CONFIG"
echo "   Cores:   $CORES"
echo "   Log:     $LOG"

# nohup ensures the process survives terminal closure
nohup snakemake \
    --snakefile workflow/Snakefile \
    --configfile "$CONFIG" \
    --config nvidia_api_key="${RESOLVED_API_KEY}" \
    --use-conda \
    --cores "$CORES" \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    2>&1 | tee "$LOG" &

echo "✅ Snakemake PID: $!  (tail -f $LOG to follow)"
