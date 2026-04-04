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
# GOTO GEMM sweep on MI300A — all 96 cores, always
#
# Sweeps:
#   1. NUMA grid configs:     1x1, 2x1, 1x2, 2x2, 4x1, 1x4, 4x2, 2x4, 4x4
#   2. Blocking parameters:   MC × KC (L2/L1 cache mapping)
#   3. NUMA distribution:     Static vs SUMMA(KC) vs SUMMA(4KC) vs Cannon
# ═══════════════════════════════════════════════════════════════════════════


T=96
export OMP_NUM_THREADS=$T
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_DISPLAY_ENV=TRUE
export OMP_STACKSIZE=2M
export MKL_DYNAMIC=FALSE

ulimit -s unlimited

spack load gcc/ktd4slj

echo "Running on $(hostname) — $T threads"
numactl --hardware 2>/dev/null || true
echo ""

OUTDIR="results"
mkdir -p "$OUTDIR"

source $(spack location -i intel-oneapi-mkl)/setvars.sh
export MKLROOT=$(spack location -i intel-oneapi-mkl)/mkl/latest

# ── Build ────────────────────────────────────────────────────────────────
CFLAGS="-O3 -fopenmp -march=native -mtune=native -ffast-math -fno-vect-cost-model -std=c++17"
MKLFLAGS="-I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm"
 
echo "═══ build ═══"
g++ $CFLAGS $MKLFLAGS -o bench_goto_vec bench_goto_vec.cpp
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
echo "Plot:    python3 plot_sweep.py $OUTDIR/"