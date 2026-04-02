#!/bin/bash
#SBATCH --job-name=b_gpu_madd
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=04:00:00
#SBATCH --output=beverin_gpu_madd_%j.out
#SBATCH --error=beverin_gpu_madd_%j.err
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive

export OMP_NUM_THREADS=96
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_SCHEDULE=static
export OMP_PROC_BIND=close

# Strong pinning via Slurm
export SLURM_CPU_BIND=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

export _DACE_NO_SYNC=1
export __HIP_PLATFORM_AMD__=1
export HIP_PLATFORM_AMD=1

# -------------------------------
# Workload parameters (BIG!)
# -------------------------------


export ROCM_HOME=/opt/rocm
export HIP_PATH=$ROCM_HOME
export HIPCC=$ROCM_HOME/bin/hipcc
export PATH=$ROCM_HOME/bin:$PATH
export LD_LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LD_LIBRARY_PATH
export CPATH=$ROCM_HOME/include:$CPATH
export LIBRARY_PATH=$ROCM_HOME/lib:$ROCM_HOME/lib64:$LIBRARY_PATH
export CFLAGS="-I$ROCM_HOME/include"
export LDFLAGS="-L$ROCM_HOME/lib -L$ROCM_HOME/lib64"
export CUPY_INSTALL_USE_HIP=1
export HCC_AMDGPU_TARGET=gfx942

spack load python@3.13.8
source ${SCRATCH}/yakup-dev-env/bin/activate

export CFLAGS="-I$(python3.13 -c "import sysconfig; print(sysconfig.get_path('include'))") ${CFLAGS}"
export C_INCLUDE_PATH="$(python3.13 -c "import sysconfig; print(sysconfig.get_path('include'))"):${C_INCLUDE_PATH}"

export HCC_AMDGPU_TARGET=gfx942
export CUPY_HIPCC_GENERATE_CODE=--offload-arch=gfx942


export C_INCLUDE_PATH=$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$PATH
export BEVERIN=1

spack load gcc/ktd4slj

hipcc -O3 -std=c++17 \
    --offload-arch=$ARCH \
    -fopenmp \
    -march=native \
    -mtune=native \
    -ffast-math \
    -munsafe-fp-atomics \
    -mllvm -amdgpu-early-inline-all=true \
    -mllvm -amdgpu-function-calls=false \
    -fgpu-flush-denormals-to-zero \
    -D__HIP_PLATFORM_AMD__=1 -DHIP_PLATFORM_AMD=1 \
    -ffast-math --offload-arch=$ARCH -fopenmp \
    -o bench_gpu bench_gpu_hip.cpp

./bench_gpu
