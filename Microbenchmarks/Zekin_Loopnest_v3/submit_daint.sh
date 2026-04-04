#!/bin/bash
#SBATCH --job-name=zekin_d
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=02:30:00
#SBATCH --output=zekin_d_gpu_%j.out
#SBATCH --error=zekin_d_gpu_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288
#SBATCH --account=g177-1
#SBATCH --exclusive
# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE


echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

spack load gcc/76jw6nu # 14.3
spack load cutensor
spack load cuda@12.9

# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
set -e
export CUTENSOR_HOME=$(spack location -i cutensor)
source ${SCRATCH}/yakup-dev-env/bin/activate
export CUDA_HOME=$(spack location -i cuda@12.9)
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH

export C_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$CUTENSOR_HOME/bin:$PATH

export ICON_DATA_PATH=/capstor/scratch/cscs/ybudanaz/icon-artifacts/velocity/data_r02b05

nvcc -O3 -std=c++17 \
    -arch=sm_90 \
    -Xcompiler=-fopenmp  \
    -Xcompiler="-march=native -ffast-math -mtune=native" \
    --use_fast_math \
    -o bench_gpu_d bench_gpu.cu

./bench_gpu_d