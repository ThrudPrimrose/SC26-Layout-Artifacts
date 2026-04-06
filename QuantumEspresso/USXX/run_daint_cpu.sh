#!/bin/bash
#SBATCH --job-name=addusxx_cpu_daint
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=00:30:00
#SBATCH --output=addusxx_cpu_daint_%j.out
#SBATCH --error=addusxx_cpu_daint_%j.err
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --account=g177-1
#SBATCH --cpus-per-task=288
#SBATCH --exclusive

export OMP_NUM_THREADS=288
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PROC_BIND=close
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

set -e

spack load gcc/76jw6nu

export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=0

echo "Building CPU-only"

# CPU source files: main_cpu.cpp, usxx_kernels_cpu.cpp
g++ -O3 -std=c++17 \
    -march=native \
    -mtune=native \
    -fopenmp \
    -ffast-math \
    -fno-vect-cost-model \
    -ftree-vectorize \
    -ffast-math -fno-math-errno -fno-trapping-math -ffinite-math-only \
    -o addusxx_cpu_daint \
    main_cpu.cpp \
    -lgomp

./addusxx_cpu_daint