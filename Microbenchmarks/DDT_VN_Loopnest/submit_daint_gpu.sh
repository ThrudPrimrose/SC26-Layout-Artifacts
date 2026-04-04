#!/bin/bash
#SBATCH --job-name=ddtvn_dg
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --output=ddtvn_dg_%j.out
#SBATCH --error=ddtvn_dg_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PROC_BIND=close
export OMP_SCHEDULE=static
export OMP_DISPLAY_ENV=TRUE


echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

# -------------------------------
# Toolchain
# -------------------------------
spack load gcc/76jw6nu  # 14.3
spack load cutensor
spack load cuda@12.9

set -e

export CUDA_HOME=$(spack location -i cuda@12.9)
export CUTENSOR_HOME=$(spack location -i cutensor)
source ${SCRATCH}/yakup-dev-env/bin/activate

export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$CUTENSOR_HOME/bin:$PATH
export ICON_DATA_PATH=/capstor/scratch/cscs/ybudanaz/beverin/icon-artifacts/velocity/data_r02b05

# -------------------------------
# Build
# -------------------------------
echo "Building bench_ddt_vn_gpu_sweep for sm_90 (GH200)"

nvcc -O3 -std=c++17 \
    -arch=sm_90 \
    --use_fast_math \
    -Xcompiler="-fopenmp,-march=native,-mtune=native,-ffast-math" \
    -o bench_ddt_vn_gpu_sweep ddt_vn_gpu.cu

echo "Build succeeded"

# -------------------------------
# Run (default nlev=90)
# -------------------------------
./bench_ddt_vn_gpu_sweep 90