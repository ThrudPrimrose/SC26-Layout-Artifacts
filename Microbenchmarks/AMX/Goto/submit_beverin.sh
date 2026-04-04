#!/bin/bash
#SBATCH --job-name=goto_numa_sweep
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=02:00:00
#SBATCH --output=goto_numa_%j.out
#SBATCH --error=goto_numa_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-bind=local
#SBATCH --exclusive

# ═══════════════════════════════════════════════════════════════════════════
# GOTO GEMM sweep on MI300A — Intel compiler toolchain throughout
#
# Using icpx (Intel LLVM C++) ensures consistent OpenMP runtime (libiomp5)
# across the compiler, MKL, and the application.  No GCC/libgomp conflict.
# ═══════════════════════════════════════════════════════════════════════════

# ── Toolchain setup ──────────────────────────────────────────────────────
spack load intel-oneapi-compilers
spack load intel-oneapi-mkl

# Intel compiler + MKL paths
INTEL_ROOT=$(spack location -i intel-oneapi-compilers)/compiler/latest
export MKLROOT=$(spack location -i intel-oneapi-mkl)/mkl/latest
export PATH=${INTEL_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${INTEL_ROOT}/lib:${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH:-}

# ── OpenMP configuration ─────────────────────────────────────────────────
T=96
export OMP_NUM_THREADS=$T
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_DISPLAY_ENV=TRUE
export OMP_STACKSIZE=2M
export MKL_DYNAMIC=FALSE

echo "Running on $(hostname) — $T threads"
echo "Compiler: $(icpx --version 2>&1 | head -1)"
echo "MKLROOT:  $MKLROOT"
echo ""
numactl --hardware 2>/dev/null || true
echo ""

OUTDIR="results"
mkdir -p "$OUTDIR"

# ── Build with Intel compiler ────────────────────────────────────────────
#
#  icpx flags:
#    -qopenmp         Intel OpenMP (libiomp5), consistent with MKL intel_thread
#    -march=native    auto-detect Zen4 AVX-512
#    -O3 -ffast-math  aggressive optimization
#    -fp-model fast   relaxed floating point (allows FMA fusion)
#
#  MKL flags (Intel threading layer — matches -qopenmp / libiomp5):
#    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#
CFLAGS="-O3 -qopenmp -march=native -ffast-math -std=c++17"
MKLFLAGS="-I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm"

echo "═══ build ═══"
echo "icpx $CFLAGS $MKLFLAGS -o bench_goto_vec bench_goto_vec.cpp"
icpx $CFLAGS $MKLFLAGS -o bench_goto_vec bench_goto_vec.cpp
echo "Build OK"
echo ""

echo "Linked libraries:"
ldd ./bench_goto_vec | grep -iE "omp|mkl|iomp"
echo ""

# ── Sanity test ──────────────────────────────────────────────────────────
echo "═══ sanity test (256³, 4 threads) ═══"
OMP_NUM_THREADS=4 OMP_DISPLAY_ENV=FALSE ./bench_goto_vec 256 256 256 \
    -t 4 -p 1x1 -w 1 -i 2 -o /dev/null 2>&1
if [ $? -ne 0 ]; then
    echo "SANITY TEST FAILED — aborting sweep"
    exit 1
fi
echo "Sanity test passed"
echo ""

M=8192; N=8192; K=8192
W=3; I=20

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 1: NUMA grid configs
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 1: NUMA grids ($T threads) ═══"

for GRID in 1x1 2x1 1x2 2x2 4x1 1x4 4x2 2x4 4x4; do
    PX=${GRID%x*}; PY=${GRID#*x}; ND=$((PX * PY))
    [ $((T % ND)) -ne 0 ] && continue
    [ $((M % (PX * 32))) -ne 0 ] && continue
    [ $((N % (PY * 32))) -ne 0 ] && continue

    OUT="$OUTDIR/grid_${GRID}.csv"
    echo "  ${GRID} (${ND} domains × $((T/ND)) cores/dom)"
    ./bench_goto_vec $M $N $K -t $T -p $GRID -w $W -i $I \
        -mc 1024 -kc 256 -nc 4096 -o "$OUT"
    echo ""
done

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 2: Blocking MC × KC (2x2 grid)
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 2: MC × KC (2x2 grid, $T threads) ═══"

for MC in 256 512 1024 2048; do
    for KC in 64 128 256 512; do
        OUT="$OUTDIR/block_mc${MC}_kc${KC}.csv"
        echo "  MC=$MC KC=$KC"
        ./bench_goto_vec $M $N $K -t $T -p 2x2 -w $W -i $I \
            -mc $MC -kc $KC -nc 4096 -o "$OUT"
        echo ""
    done
done

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 3: NUMA distribution (2x2 and 4x4)
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 3: NUMA distribution ═══"

for GRID in 2x2 4x4; do
    PX=${GRID%x*}; PY=${GRID#*x}; ND=$((PX * PY))
    [ $((T % ND)) -ne 0 ] && continue
    [ $((M % (PX * 32))) -ne 0 ] && continue
    [ $((N % (PY * 32))) -ne 0 ] && continue

    OUT="$OUTDIR/numa_${GRID}.csv"
    echo "  ${GRID}"
    ./bench_goto_vec $M $N $K -t $T -p $GRID -w $W -i $I \
        -mc 1024 -kc 256 -nc 4096 -o "$OUT"
    echo ""
done

# ═══════════════════════════════════════════════════════════════════════════
# Merge + summary
# ═══════════════════════════════════════════════════════════════════════════
MASTER="$OUTDIR/all.csv"
head -1 "$(ls "$OUTDIR"/*.csv | head -1)" > "$MASTER"
for f in "$OUTDIR"/*.csv; do [ "$f" = "$MASTER" ] && continue; tail -n +2 "$f" >> "$MASTER"; done

echo "═══ summary ═══"
for f in "$OUTDIR"/*.csv; do
    [ "$f" = "$MASTER" ] && continue
    echo "[$(basename "$f" .csv)]"
    awk -F, 'NR>1 && $7>0{k=$1; s[k]+=$8; n[k]++}
        END{for(k in s) printf "  %-40s %8.2f GF/s\n",k,s[k]/n[k]}' "$f" | sort -k2 -rn
    echo ""
done

echo "Results: $OUTDIR/"