#!/bin/bash
#SBATCH --job-name=ddtvn_bc
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=08:00:00
#SBATCH --output=ddtvn_bc_%j.out
#SBATCH --error=ddtvn_bc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=96
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_PROC_BIND=close
export OMP_SCHEDULE=static
export OMP_DISPLAY_ENV=TRUE
export SLURM_CPU_BIND=cores

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

export BEVERIN=1

# -------------------------------
# Toolchain
# -------------------------------
export ROCM_HOME=/opt/rocm
export LLVM_HOME=/opt/rocm/llvm/
export PATH=$ROCM_HOME/bin:$LLVM_HOME/bin:$SCRATCH/bin:$PATH
export LD_LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LLVM_HOME/lib:$LLVM_HOME/lib64:$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export CPATH=$ROCM_HOME/include:$LLVM_HOME/include:$CPATH
export LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LLVM_HOME/lib:$LLVM_HOME/lib64:$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH

spack load python@3.13.8
source ${SCRATCH}/yakup-dev-env/bin/activate

# -------------------------------
# Build
# -------------------------------
echo "Building bench_ddt_vn_cpu_sweep"

set -e

g++ -O3 -std=c++17 \
    -march=native \
    -mtune=native \
    -fopenmp \
    -ffast-math \
    -fno-vect-cost-model \
    -ftree-vectorize \
    -o bench_ddt_vn_cpu_sweep ddt_vn_cpu.cpp

echo "Build succeeded"

# -------------------------------
# Run (default nlev sweep: 65, 90)
# -------------------------------
./bench_ddt_vn_cpu_sweep