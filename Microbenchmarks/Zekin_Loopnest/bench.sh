#!/usr/bin/env bash
# run_benchmark.sh — compile and run ICON ekinh CUDA kernel benchmarks
set -euo pipefail

SM="sm_90a"

echo "Target architecture: ${SM}"

# ── compile ───────────────────────────────────────────────────────────────
NVCC=${NVCC:-nvcc}
BIN=./benchmark
export LD_LIBRARY_PATH="$(spack location -i cuda@12.9)/lib64:${LD_LIBRARY_PATH:-}"

"${NVCC}" \
    -O3 \
    -arch="${SM}" \
    --use_fast_math \
    -lineinfo \
    -Xcompiler "-O3,-march=native,-mtune=native" \
    benchmark.cu \
    -lcurand \
    -o "${BIN}"

echo "Compiled → ${BIN}"
echo ""

# ── run ───────────────────────────────────────────────────────────────────
"${BIN}" 2>&1 | tee benchmark_results.txt

echo ""
echo "Results saved to benchmark_results.txt"