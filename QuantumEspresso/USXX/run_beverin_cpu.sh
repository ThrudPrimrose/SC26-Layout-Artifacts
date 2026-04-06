#!/bin/bash
#SBATCH --job-name=addusxx_cpu_beverin
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=04:00:00
#SBATCH --output=addusxx_cpu_beverin_%j.out
#SBATCH --error=addusxx_cpu_beverin_%j.err
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive

export OMP_NUM_THREADS=96
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_SCHEDULE=static
export OMP_PROC_BIND=close
export SLURM_CPU_BIND=cores
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=1

spack load gcc/ktd4slj
export GCCHOME=$(spack location -i gcc/ktd4slj)
export C_INCLUDE_PATH=$GCCHOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$GCCHOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$GCCHOME/lib:$GCCHOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$GCCHOME/lib:$GCCHOME/lib64:$LD_LIBRARY_PATH
export PATH=$GCCHOME/bin:$PATH

echo "Building CPU-only"

# CPU source files: main_cpu.cpp, usxx_kernels_cpu.cpp
g++ -O3 -std=c++17 \
    -march=native \
    -mtune=native \
    -fopenmp \
    -ffast-math \
    -fno-vect-cost-model \
    -ftree-vectorize \
    -o addusxx_cpu_beverin \
    main_cpu.cpp \
    -lgomp -lnuma

./addusxx_cpu_beverin