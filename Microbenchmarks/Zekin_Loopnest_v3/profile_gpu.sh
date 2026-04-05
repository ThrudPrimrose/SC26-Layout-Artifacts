#!/bin/bash
#SBATCH --job-name=profile_gpu
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=02:30:00
#SBATCH --output=profile_gpu_%j.out
#SBATCH --error=profile_gpu_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288
#SBATCH --account=g177-1
#SBATCH --exclusive
# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE
echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"
spack load gcc/76jw6nu # 14.3
spack load cutensor
spack load cuda@12.9
# -------------------------------
# Environment
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
# -------------------------------
# Compile with debug info for source import
# -------------------------------
echo "=== Compiling profile_gpu.cu ==="
nvcc -O3 -std=c++17 \
    -arch=sm_90 \
    -lineinfo \
    -Xcompiler=-fopenmp \
    -Xcompiler="-march=native -ffast-math -mtune=native" \
    --use_fast_math \
    -o profile_gpu profile_gpu.cu
echo "=== Compilation done ==="
# -------------------------------
# Dry run (verify correctness, no profiling)
# -------------------------------
echo "=== Dry run ==="
./profile_gpu
echo "=== Dry run done ==="
# -------------------------------
# Profile with nsight compute
# -------------------------------
echo "=== Profiling with ncu ==="
ncu --set full \
    --import-source=yes \
    -f -o ncu_profile_6kernels \
    ./profile_gpu
echo "=== Profiling done ==="
echo "Report: ncu_profile_6kernels.ncu-rep"
echo "View:   ncu-ui ncu_profile_6kernels.ncu-rep"
