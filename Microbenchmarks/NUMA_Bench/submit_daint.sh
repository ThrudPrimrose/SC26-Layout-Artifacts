#!/bin/bash
#SBATCH --job-name=numa_bench
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=numa_bench_daint_%j.out
#SBATCH --error=numa_bench_daint_%j.err
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
spack load cuda@12.9

CFLAGS="-O3 -fopenmp -mtune=native -ftree-vectorize -fno-vect-cost-model -march=native -ffast-math -std=c++17"
HIPFLAGS="--offload-arch=gfx942 -O3 -ffast-math -std=c++17"

g++ $CFLAGS -o bench_numa numa_bench.cpp -lnuma

./bench_numa