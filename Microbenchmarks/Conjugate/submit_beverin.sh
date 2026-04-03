#!/bin/bash
#SBATCH --job-name=conj_bench
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=00:50:00
#SBATCH --output=conj_%j.out
#SBATCH --error=conj_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-bind=local
# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=96
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"
set -e

CFLAGS="-O3 -fopenmp -mtune=native -ftree-vectorize -fno-vect-cost-model -march=native -ffast-math -std=c++17"
HIPFLAGS="--offload-arch=gfx942 -O3 -ffast-math -std=c++17 ${CFLAGS}"

echo "═══ build ═══"
echo "[1/4] conj_cpu_inplace.cpp    (CPU in-place)"
g++ $CFLAGS -o conj_cpu_inplace conj_cpu_inplace.cpp -lnuma

echo "[2/4] conj_cpu_oop.cpp        (CPU out-of-place)"
g++ $CFLAGS -o conj_cpu_oop conj_cpu_oop.cpp -lnuma

echo "[3/4] conj_gpu_inplace.cpp    (GPU in-place)"
hipcc $HIPFLAGS -o conj_gpu_inplace conj_gpu_inplace.cpp

echo "[4/4] conj_gpu_oop.cpp        (GPU out-of-place)"
hipcc $HIPFLAGS -o conj_gpu_oop conj_gpu_oop.cpp

echo ""
echo "═══ run ═══"

echo "--- CPU in-place ---"
./conj_cpu_inplace
echo ""

echo "--- CPU out-of-place ---"
./conj_cpu_oop
echo ""

echo "--- GPU in-place ---"
./conj_gpu_inplace
echo ""

echo "--- GPU out-of-place ---"
./conj_gpu_oop
echo ""

echo "═══ CSV files ═══"
echo "  results_cpu_inplace.csv"
echo "  results_cpu_oop.csv"
echo "  results_gpu_inplace.csv"
echo "  results_gpu_oop.csv"
echo ""

echo "═══ summary (averages) ═══"
for f in results_cpu_inplace.csv results_cpu_oop.csv results_gpu_inplace.csv results_gpu_oop.csv; do
    [ -f "$f" ] || continue
    echo "[$f]"
    awk -F, 'NR>1{s[$1","$2]+=$4; g[$1","$2]+=$5; n[$1","$2]++}
            END{for(k in s) printf "  K=%-2s %-14s avg=%8.4f ms  %7.1f GB/s\n",
                substr(k,1,index(k,",")-1), substr(k,index(k,",")+1),
                s[k]/n[k], g[k]/n[k]}' "$f" | sort -t= -k1,1n -k2
    echo ""
done

mkdir -p results/beverin
mv results_cpu_inplace.csv results_cpu_oop.csv \
   results_gpu_inplace.csv results_gpu_oop.csv \
   results/beverin/