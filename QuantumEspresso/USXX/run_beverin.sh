#!/bin/bash
#SBATCH --job-name=addusxx_gpu_beverin
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=04:00:00
#SBATCH --output=addusxx_gpu_beverin_%j.out
#SBATCH --error=addusxx_gpu_beverin_%j.err
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

export __HIP_PLATFORM_AMD__=1
export HIP_PLATFORM_AMD=1

export ROCM_HOME=/opt/rocm
export HIP_PATH=$ROCM_HOME
export HIPCC=$ROCM_HOME/bin/hipcc
export PATH=$ROCM_HOME/bin:$PATH
export LD_LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LD_LIBRARY_PATH
export CPATH=$ROCM_HOME/include:$CPATH
export LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LIBRARY_PATH
export HCC_AMDGPU_TARGET=gfx942

export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=1

export ARCH=${HIP_ARCH:-gfx942}

hipcc -O3 -std=c++17 \
    --offload-arch=$ARCH \
    -march=native \
    -mtune=native \
    -ffast-math \
    -fno-math-errno \
    -ffinite-math-only \
    -fno-trapping-math \
    -freciprocal-math \
    -ffp-contract=fast \
    -munsafe-fp-atomics \
    -mllvm -amdgpu-early-inline-all=true \
    -mllvm -amdgpu-function-calls=false \
    -fgpu-flush-denormals-to-zero \
    -funroll-loops \
    -D__HIP_PLATFORM_AMD__=1 -DHIP_PLATFORM_AMD=1 \
    -fgpu-rdc \
    -o addusxx_gpu_beverin \
    main_hip.cpp

./addusxx_gpu_beverin