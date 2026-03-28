#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=nbody_daint_%j.out
#SBATCH --error=nbody_daint_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288
#SBATCH --mem-bind=local

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

spack load gcc/76jw6nu # 14.3
spack loadd cuda@12.9
# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
set -e

CFLAGS="-O3 -fopenmp -ftree-vectorize -fvect-cost-model=cheap -march=native -ffast-math -std=c++17"
NVFLAGS="-O3 -use_fast_math -arch=native --expt-relaxed-constexpr -std=c++17"

echo "═══ build ═══"

echo "[1/4] conjugate_inplace.cpp  (CPU in-place)"
g++ $CFLAGS -o conj_ip_cpu conjugate_inplace.cpp -lnuma

echo "[2/4] conjugate.cpp           (CPU out-of-place)"
g++ $CFLAGS -o conj_oop_cpu conjugate.cpp -lnuma

echo "[3/4] conjugate_inplace.cu   (GPU in-place)"
nvcc $NVFLAGS -Xcompiler "$CFLAGS" -o conj_ip_gpu conjugate_inplace.cu

echo "[4/4] conjugate.cu           (GPU out-of-place)"
nvcc $NVFLAGS -Xcompiler "$CFLAGS" -o conj_oop_gpu conjugate.cu


echo ""
echo "═══ run ═══"

echo "--- CPU in-place ---"
./conj_ip_cpu
echo ""

echo "--- CPU out-of-place ---"
./conj_oop_cpu
echo ""

echo "--- GPU in-place ---"
./conj_ip_gpu
echo ""

echo "--- GPU out-of-place ---"
./conj_oop_gpu
echo ""

echo "═══ CSV files ═══"
echo "  results_cpu_ip.csv   (CPU in-place,  per-run)"
echo "  results_cpu.csv      (CPU out-of-place, per-run)"
echo "  results_gpu_ip.csv   (GPU in-place,  per-run)"
echo "  results_gpu_oop.csv   (GPU out-of-place,  per-run)"

echo ""
echo "═══ summary (averages) ═══"
for f in results_cpu_ip.csv results_cpu.csv results_gpu_ip.csv results_gpu_oop.csv; do
    [ -f "$f" ] || continue
    echo "[$f]"
    awk -F, 'NR>1{s[$1","$2]+=$4; g[$1","$2]+=$5; n[$1","$2]++}
        END{for(k in s) printf "  %-20s avg=%8.4f ms  %7.1f GB/s\n",
            k, s[k]/n[k], g[k]/n[k]}' "$f" | sort
    echo ""
done

mkdir -p results/daint
mv results_cpu_ip.csv results_cpu.csv results_gpu_ip.csv results_gpu_oop.csv results/daint/