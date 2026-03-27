#!/bin/bash
#SBATCH --job-name=numa_bench
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=numa_bench_%j.out
#SBATCH --error=numa_bench_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288

set -e
PLATFORM="daint"

# === OpenMP (Grace: 4 NUMA × 72 Neoverse V2 cores) ===
export OMP_NUM_THREADS=288
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PROC_BIND=close
export OMP_SCHEDULE=static
export OMP_DISPLAY_ENV=TRUE
export SLURM_CPU_BIND=cores

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

# === Toolchain ===
spack load gcc/76jw6nu
spack load cutensor
spack load cuda@12.9

source ${SCRATCH}/yakup-dev-env/bin/activate

export CUDA_HOME=$(spack location -i cuda@12.9)
export CUTENSOR_HOME=$(spack location -i cutensor)
export C_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$CUTENSOR_HOME/bin:$PATH

# === Build ===
echo "Building..."
gcc -O3 -std=c++17 -march=native -fopenmp -ffast-math -mtune=native \
    -ftree-vectorize -o numa_calibrate numa_calibrate.c -lm
gcc -O3 -march=native -fopenmp -o numa_triad numa_triad.c -lnuma

echo ""
echo "=== Daint (Grace Neoverse V2, 4×72 cores) ==="
echo "  Expected: β=64 (4KB pages / 64B CL)"
echo ""

# === Topology ===
echo "--- System topology ---"
numactl --hardware 2>/dev/null || true
echo ""
lscpu | grep -E "^(Thread|Core|Socket|NUMA|CPU\(s\)|Model name)" || true
echo ""
getconf PAGESIZE && echo "Page size: $(getconf PAGESIZE) bytes" || true
echo ""

# === Run all benchmarks ===
for bench in stride_bw stride_lat numa_bw numa_lat numa_matrix; do
    echo "================================================================"
    echo "  ${bench}"
    echo "================================================================"
    ./numa_calibrate $bench | tee ${bench}_${PLATFORM}.txt
done

echo "================================================================"
echo "  numa_triad"
echo "================================================================"
./numa_triad | tee numa_triad_${PLATFORM}.txt

echo ""
echo "=== Done. Results: *_${PLATFORM}.txt ==="