#!/bin/bash
#SBATCH --job-name=conjugate_bench
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=00:50:00
#SBATCH --output=conjugate_%j.out
#SBATCH --error=conjugate_%j.err
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
HIPFLAGS="-D__HIP_PLATFORM_AMD__ -fvect-cost-model=cheap -march=native -ffast-math --offload-arch=gfx942 -O3 -ffast-math -std=c++17 -fopenmp=libgomp"

echo "═══ build ═══"
echo "[1/4] conjugate_inplace.cpp    (CPU in-place)"
g++ $CFLAGS -o conjugate_cpu_inplace conjugate_inplace.cpp -lnuma

echo "[2/4] conjugate.cpp        (CPU out-of-place)"
g++ $CFLAGS -o conjugate_cpu_oop conjugate.cpp -lnuma

echo "[3/4] conjugate_inplace_hip.cpp    (GPU in-place)"
hipcc $HIPFLAGS -o conjugate_gpu_inplace conjugate_inplace_hip.cpp

echo "[4/4] conjugate_hip.cpp        (GPU out-of-place)"
hipcc $HIPFLAGS -o conjugate_gpu_oop conjugate_hip.cpp

echo ""
echo "═══ run ═══"

echo "--- CPU in-place ---"
./conjugate_cpu_inplace
echo ""

echo "--- CPU out-of-place ---"
./conjugate_cpu_oop
echo ""

echo "--- GPU in-place ---"
./conjugate_gpu_inplace
echo ""

echo "--- GPU out-of-place ---"
./conjugate_gpu_oop
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