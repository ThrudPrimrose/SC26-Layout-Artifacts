#!/bin/bash
#SBATCH --job-name=madd_daint
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=madd_daint_%j.out
#SBATCH --error=madd_daint_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PROC_BIND="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PLACES=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE


echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

spack load gcc/76jw6nu # 14.3
spack load cuda@12.9
spack load cutensor

# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
set -e

source ${SCRATCH}/yakup-dev-env/bin/activate
export CUDA_HOME=$(spack location -i cuda@12.9)
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH

export CUTENSOR_HOME=$(spack location -i cutensor)

export C_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$CUTENSOR_HOME/bin:$PATH

export OPENBLAS_HOME=$(spack location -i openblas@0.3.29)
export C_INCLUDE_PATH=$SCRATCH/include:$OPENBLAS_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$OPENBLAS_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$OPENBLAS_HOME/lib:$OPENBLAS_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$OPENBLAS_HOME/lib:$OPENBLAS_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$OPENBLAS_HOME/bin:$PATH

export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=0

g++ -O3 -march=native -mtune=native -fopenmp -ffast-math -fno-vect-cost-model  -fprefetch-loop-arrays -funroll-loops -ftree-loop-distribution -falign-loops=64 -std=c++17 -o bench_cpu bench_cpu.cpp -lnuma

./bench_cpu