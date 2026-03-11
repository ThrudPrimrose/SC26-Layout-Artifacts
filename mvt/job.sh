#!/bin/bash
#SBATCH --job-name=mvt_bench
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --output=job_output_mvt.txt
#SBATCH --error=job_error_mvt.txt

# sweep_mvt.sh — sweep block/tile sizes for mvt modes 1-6
# Usage: ./sweep_mvt.sh [N=8192] [CSV=results_mvt.csv] [NRUNS=50] [NWARM=5]
#
# Modes:
#   1 — A layout variants   (row-major, col-major)
#   2 — loop orderings      (ij_ji, ji_ji, ij_ij, ji_ij)
#   3 — blocked A SZ_B×SZ_B (RR/RC/CR/CC inner×outer), swept over BLOCK_SIZES
#   4 — tiled loops SZ_T×SZ_T (row-major A), swept over BLOCK_SIZES
#   5 — OpenBLAS 2× dgemv
#   6 — fused single-pass:
#         fused_tiled              (row-major A, SZ_T×SZ_T, heap-private y2)
#         fused_blk_{RR,RC,CR,CC}  (blocked A,  heap-private y2)
#       Both swept over BLOCK_SIZES (SZ_B=SZ_T per binary).

set -euo pipefail

# ── script location — SRC is always next to this script ─────────────────
SCRIPT_DIR="."
SRC="mvt_bench.cpp"

# ── environment ──────────────────────────────────────────────────────────
spack load gcc@14
spack load openblas
spack load cuda   # for hwloc / numactl if needed

# ── parameters (override via positional args) ────────────────────────────
N=${1:-16384}
CSV=${2:-results_mvt.csv}
NRUNS=${3:-50}
NWARM=${4:-5}

OPENBLAS_ROOT=$(spack location -i openblas)
CXX_FLAGS="-O3 -fopenmp -std=c++17 -march=native -mtune=native -ffast-math"
INC="-I${OPENBLAS_ROOT}/include"
LIB="-L${OPENBLAS_ROOT}/lib"
CXX="g++ ${CXX_FLAGS} ${INC} ${LIB}"
LINK="-lopenblas"

BLOCK_SIZES=(16 32 64 128 256 512)

echo "============================================================"
echo "  MVT benchmark sweep"
echo "  N=${N}  NRUNS=${NRUNS}  NWARM=${NWARM}  -> ${CSV}"
echo "  Block/tile sizes: ${BLOCK_SIZES[*]}"
echo "  SRC=${SRC}"
echo "============================================================"

# ── modes 1, 2, 5 — layout/ordering/openblas, compile once ──────────────
FIXED_BIN="${SCRIPT_DIR}/.mvt_fixed"
echo ""
echo "[compile] fixed binary (modes 1, 2, 5) ..."
$CXX -DSZ_N=${N} -DSZ_B=32 -DSZ_T=32 \
     -DNRUNS=${NRUNS} -DNWARM=${NWARM} \
     -o "${FIXED_BIN}" "${SRC}" ${LINK} || {
    echo "ERROR: fixed compilation failed"; exit 1
}

echo "[run] mode 5 — OpenBLAS 2x dgemv ..."
"${FIXED_BIN}" 5 "${CSV}"

echo "[run] mode 1 — A layout variants (row-major, col-major) ..."
"${FIXED_BIN}" 1 "${CSV}"

echo "[run] mode 2 — loop orderings ..."
"${FIXED_BIN}" 2 "${CSV}"

rm -f "${FIXED_BIN}"

# ── modes 3, 4, 6 — sweep over BLOCK_SIZES, one binary per size ─────────
# SZ_B and SZ_T are set equal per binary.
# mode 3 uses SZ_B, mode 4 uses SZ_T, mode 6 fused_tiled uses SZ_T and
# fused_blk_* use SZ_B.
echo ""
echo "[sweep] modes 3, 4, 6 over block/tile sizes ..."
for BS in "${BLOCK_SIZES[@]}"; do
    BIN="${SCRIPT_DIR}/.mvt_bs_${BS}"
    echo "  [compile] SZ_B=${BS} SZ_T=${BS} ..."
    $CXX -DSZ_N=${N} -DSZ_B=${BS} -DSZ_T=${BS} \
         -DNRUNS=${NRUNS} -DNWARM=${NWARM} \
         -o "${BIN}" "${SRC}" ${LINK} 2>/dev/null || {
        echo "  SKIP SZ_B/SZ_T=${BS} (compilation failed)"; continue
    }
    echo "  [run] mode 3 — blocked A (RR/RC/CR/CC) SZ_B=${BS} ..."
    "${BIN}" 3 "${CSV}"
    echo "  [run] mode 4 — tiled loops SZ_T=${BS} ..."
    "${BIN}" 4 "${CSV}"
    rm -f "${BIN}"
done

echo ""
echo "============================================================"
echo "  DONE — $(wc -l < "${CSV}") rows in ${CSV}"
echo "============================================================"