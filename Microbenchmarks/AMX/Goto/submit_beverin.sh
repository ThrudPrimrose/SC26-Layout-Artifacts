#!/bin/bash
#SBATCH --job-name=goto_sweep
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=02:00:00
#SBATCH --output=goto_%j.out
#SBATCH --error=goto_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-bind=local
#SBATCH --exclusive

# ═══════════════════════════════════════════════════════════════════════════
# Focused sweep: KC is the only parameter that matters.
#
# MC has no effect in the no-pack implementation (IC partitions at BM=32).
# NC has no effect (JC loop runs once when N_local ≤ NC).
#
# The data shows goto_flat scales linearly with KC:
#   KC=64  → 3000 GF/s (128 C-writeback passes)
#   KC=512 → 5200 GF/s (16 C-writeback passes)
# Extrapolating: KC=4096 → ~2 passes → should approach MKL's 5800.
#
# L2 constraint: each thread's working set in the k-loop is
#   BM × KC × 2B (A row strip) + KC × BN × 2B (B column strip being streamed)
#   = 32 × KC × 2 + KC × 32 × 2 = 128 × KC bytes
# For KC=4096: 512 KB.  Fits in Zen4's 1 MB L2 with room to spare.
# For KC=8192 (= K): 1 MB.  Exactly fills L2 — only one PC step needed.
# ═══════════════════════════════════════════════════════════════════════════

spack load intel-oneapi-compilers
spack load intel-oneapi-mkl

INTEL_ROOT=$(spack location -i intel-oneapi-compilers)/compiler/latest
export MKLROOT=$(spack location -i intel-oneapi-mkl)/mkl/latest
export PATH=${INTEL_ROOT}/bin:${PATH}
export LD_LIBRARY_PATH=${INTEL_ROOT}/lib:${MKLROOT}/lib/intel64:${LD_LIBRARY_PATH:-}

T=96
export OMP_NUM_THREADS=$T
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_DISPLAY_ENV=TRUE
export OMP_STACKSIZE=2M
export MKL_DYNAMIC=FALSE

echo "Running on $(hostname) — $T threads"
numactl --hardware 2>/dev/null | head -20 || true
echo ""

OUTDIR="results"
mkdir -p "$OUTDIR"

CFLAGS="-O3 -g -qopenmp -march=native -ffast-math -std=c++17"
MKLFLAGS="-I${MKLROOT}/include -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm"

echo "═══ build ═══"
icpx $CFLAGS $MKLFLAGS -o bench_goto_vec bench_goto_vec.cpp
echo "Build OK"
ldd ./bench_goto_vec | grep -iE "omp|mkl"
echo ""

# Quick sanity
echo "═══ sanity (256³) ═══"
OMP_NUM_THREADS=4 OMP_DISPLAY_ENV=FALSE \
    ./bench_goto_vec 256 256 256 -t 4 -w 1 -i 2 -o /dev/null 2>&1
[ $? -ne 0 ] && echo "FAILED" && exit 1
echo "OK"
echo ""

M=8192; N=8192; K=8192
W=3; I=20

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 1: KC scan (the only parameter that matters)
#
# For each KC value, both goto_flat and cannon_2x2 are benchmarked.
# goto_flat does K/KC writeback passes over C (256 MB).
# cannon_2x2 has 2 CANNON steps × K/(2*KC) passes per step.
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 1: KC scan ($T threads, 8192³) ═══"

for KC in 32 64 128 256 512 1024 2048 4096 8192; do
    # KC must divide K and divide BK (=32)
    [ $((K % KC)) -ne 0 ] && continue
    [ $((KC % 32)) -ne 0 ] && continue

    OUT="$OUTDIR/kc_${KC}.csv"
    echo "  KC=$KC  (PC steps: $((K/KC)))"
    ./bench_goto_vec $M $N $K -t $T -kc $KC -w $W -i $I -o "$OUT"
    echo ""
done

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 2: Problem size scaling at best KC
#
# Test whether the KC=large advantage holds at different matrix sizes.
# Uses KC=2048 (expected sweet spot from sweep 1).
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 2: Problem size scaling (KC=2048) ═══"

for N_SZ in 1024 2048 4096 8192 16384; do
    [ $((N_SZ % 64)) -ne 0 ] && continue

    OUT="$OUTDIR/size_${N_SZ}_kc2048.csv"
    echo "  ${N_SZ}³  KC=2048"
    ./bench_goto_vec $N_SZ $N_SZ $N_SZ -t $T -kc 2048 -w $W -i $I -o "$OUT"
    echo ""
done

# ═══════════════════════════════════════════════════════════════════════════
# SWEEP 3: Non-square problems (tall-skinny, short-wide)
#
# Tests whether CANNON's 2D decomposition helps when M ≠ N.
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Sweep 3: Non-square (KC=2048) ═══"

for SHAPE in "16384 8192 8192" "8192 16384 8192" "8192 8192 16384" "4096 4096 16384"; do
    set -- $SHAPE
    SM=$1; SN=$2; SK=$3
    [ $((SM % 64)) -ne 0 ] && continue
    [ $((SN % 64)) -ne 0 ] && continue
    [ $((SK % 64)) -ne 0 ] && continue

    OUT="$OUTDIR/shape_${SM}x${SN}x${SK}_kc2048.csv"
    echo "  ${SM}×${SN}×${SK}  KC=2048"
    ./bench_goto_vec $SM $SN $SK -t $T -kc 2048 -w $W -i $I -o "$OUT"
    echo ""
done

# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════
echo "═══ Summary ═══"
echo ""

# KC sweep summary (most important)
echo "KC scan (8192³):"
echo "  KC    goto_flat     cannon_2x2    MKL"
echo "  ----  ----------    ----------    ----------"
for KC in 32 64 128 256 512 1024 2048 4096 8192; do
    F="$OUTDIR/kc_${KC}.csv"
    [ -f "$F" ] || continue
    awk -F, 'NR>1{s[$1]+=$8; n[$1]++}
        END{
            printf "  %-5s", "'$KC'";
            for (k in s) {
                if (k == "goto_flat") gf = s[k]/n[k];
                else if (k == "cannon_2x2") cn = s[k]/n[k];
                else if (k == "mkl_sgemm") mk = s[k]/n[k];
            }
            printf " %8.0f      %8.0f      %8.0f\n", gf, cn, mk;
        }' "$F"
done
echo ""

# All results
echo "All results:"
for f in "$OUTDIR"/*.csv; do
    [ "$(basename "$f")" = "all.csv" ] && continue
    echo "  [$(basename "$f" .csv)]"
    awk -F, 'NR>1 && $7>0{k=$1; s[k]+=$8; n[k]++}
        END{for(k in s) printf "    %-20s %8.0f GF/s\n",k,s[k]/n[k]}' "$f" | sort -k2 -rn
done
echo ""
echo "Results: $OUTDIR/"