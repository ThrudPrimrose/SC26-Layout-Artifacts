#!/bin/bash
#SBATCH --job-name=addusxx_gpu_daint
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=addusxx_gpu_daint_%j.out
#SBATCH --error=addusxx_gpu_daint_%j.err
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
spack load cuda@12.9

export CUDA_HOME=$(spack location -i cuda@12.9)
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH

export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=0

echo "Building GPU (nvcc)"

nvcc -O3 -std=c++17 \
    -arch=native \
    -Xcompiler -fopenmp \
    -Xcompiler -march=native \
    -Xcompiler "-ffast-math -fno-math-errno -fno-trapping-math -ffinite-math-only" \
    --use_fast_math \
    -o addusxx_gpu_daint \
    main.cu 
    
./addusxx_gpu_daint