#!/usr/bin/env bash
# run_bench.sh — Run both variants in separate processes with cooldown
# Usage: bash run_bench.sh [REPS]
#   REPS defaults to 100

set -euo pipefail

REPS="${1:-100}"
SCRIPT="zsolqa_loopnest.py"
CSV_3D="bench_3d.csv"
CSV_SPLIT="bench_split.csv"
SLEEP_SEC=10

echo "========================================"
echo "  Condensation kernel benchmark"
echo "  REPS=${REPS}  COOLDOWN=${SLEEP_SEC}s"
echo "========================================"
echo ""

# --- Phase 1: 3D variant ---
echo "[$(date '+%H:%M:%S')] Running 3D variant (${REPS} reps)..."
python3 "${SCRIPT}" --variant=3d --reps="${REPS}" --csv="${CSV_3D}"

echo "[$(date '+%H:%M:%S')] Cooling down for ${SLEEP_SEC}s..."
sleep "${SLEEP_SEC}"

# --- Phase 2: Split variant ---
echo "[$(date '+%H:%M:%S')] Running Split variant (${REPS} reps)..."
python3 "${SCRIPT}" --variant=split --reps="${REPS}" --csv="${CSV_SPLIT}"

echo ""
echo "[$(date '+%H:%M:%S')] Both runs complete."
echo ""

# --- Phase 3: Generate violin plot ---
echo "[$(date '+%H:%M:%S')] Generating violin plot..."
python3 plot_violin.py --csv "${CSV_3D}" "${CSV_SPLIT}" --output condense_violin.png
echo "[$(date '+%H:%M:%S')] Done. Output: condense_violin.png"