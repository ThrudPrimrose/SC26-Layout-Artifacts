#!/bin/bash
#SBATCH --job-name=mvt_bench
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --time=06:00:00
#SBATCH --output=job_output_mvt.txt
#SBATCH --error=job_error_mvt.txt

# sweep_mvt.sh — sweep block/tile sizes for mvt modes 1-6
# Usage: ./sweep_mvt.sh [N=8192] [CSV=results_mvt.csv] [NRUNS=50] [NWARM=5]
#
# Modes:
#   1 — A layout variants   (row-major, col-major)
#   2 — loop orderings      (ij_ji, ji_ji, ij_ij, ji_ij)
#   3 — blocked A B×B       (RR/RC/CR/CC inner×outer), swept over BLOCK_SIZES
#   4 — tiled loops T×T     (row-major A), swept over BLOCK_SIZES
#   5 — OpenBLAS 2× dgemv
#   6 — fused single-pass:
#         fused_tiled            (row-major A, T×T, heap-private y2)
#         fused_blk_{RR,RC,CR,CC}(blocked A, heap-private y2)
#       Both swept over BLOCK_SIZES (B=T per binary).

set -euo pipefail

# ── environment ─────────────────────────────────────────────────────────
spack load gcc@14
spack load openblas
spack load cuda   # for hwloc / numactl if needed

# ── parameters (override via positional args) ────────────────────────────
N=${1:-8192}
CSV=${2:-results_mvt.csv}
NRUNS=${3:-50}
NWARM=${4:-5}
SRC=mvt_bench.cpp

OPENBLAS_ROOT=$(spack location -i openblas)
CXX_FLAGS="-O3 -fopenmp -std=c++17 -march=native -mtune=native -ffast-math"
INC="-I${OPENBLAS_ROOT}/include"
LIB="-L${OPENBLAS_ROOT}/lib -lopenblas"
CXX="g++ ${CXX_FLAGS} ${INC} ${LIB}"

BLOCK_SIZES=(32 64 128 256)

echo "============================================================"
echo "  MVT benchmark sweep"
echo "  N=${N}  NRUNS=${NRUNS}  NWARM=${NWARM}  -> ${CSV}"
echo "  Block/tile sizes: ${BLOCK_SIZES[*]}"
echo "============================================================"

# ── modes 1, 2, 5 — layout/ordering/openblas, compile once ──────────────
FIXED_BIN=".mvt_fixed"
echo ""
echo "[compile] fixed binary (modes 1, 2, 5) ..."
$CXX -DSZ_N=${N} -DB=32 -DT=32 \
     -DNRUNS=${NRUNS} -DNWARM=${NWARM} \
     -o "${FIXED_BIN}" "${SRC}" || {
    echo "ERROR: fixed compilation failed"; exit 1
}

echo "[run] mode 5 — OpenBLAS 2x dgemv ..."
./"${FIXED_BIN}" 5 "${CSV}"

echo "[run] mode 1 — A layout variants (row-major, col-major) ..."
./"${FIXED_BIN}" 1 "${CSV}"

echo "[run] mode 2 — loop orderings ..."
./"${FIXED_BIN}" 2 "${CSV}"

rm -f "${FIXED_BIN}"

# ── modes 3, 4, 6 — sweep over BLOCK_SIZES, one binary per size ─────────
# B and T are set equal per binary; mode 3 uses B, mode 4 uses T,
# mode 6 fused_tiled uses T and fused_blk_* use B.
echo ""
echo "[sweep] modes 3, 4, 6 over block/tile sizes ..."
for BS in "${BLOCK_SIZES[@]}"; do
    BIN=".mvt_bs_${BS}"
    echo "  [compile] B=${BS} T=${BS} ..."
    $CXX -DSZ_N=${N} -DB=${BS} -DT=${BS} \
         -DNRUNS=${NRUNS} -DNWARM=${NWARM} \
         -o "${BIN}" "${SRC}" 2>/dev/null || {
        echo "  SKIP B/T=${BS} (compilation failed)"; continue
    }
    echo "  [run] mode 3 — blocked A (RR/RC/CR/CC) B=${BS} ..."
    ./"${BIN}" 3 "${CSV}"
    echo "  [run] mode 4 — tiled loops T=${BS} ..."
    ./"${BIN}" 4 "${CSV}"
    echo "  [run] mode 6 — fused tiled T=${BS} + fused_blk_* B=${BS} ..."
    ./"${BIN}" 6 "${CSV}"
    rm -f "${BIN}"
done

echo ""
echo "============================================================"
echo "  DONE — $(wc -l < "${CSV}") rows in ${CSV}"
echo "============================================================"