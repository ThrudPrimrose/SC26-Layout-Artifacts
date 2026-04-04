#!/bin/bash
#SBATCH --job-name=ddtvn_bg
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=08:00:00
#SBATCH --output=ddtvn_bg_%j.out
#SBATCH --error=ddtvn_bg_%j.err
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

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

export _DACE_NO_SYNC=1
export __HIP_PLATFORM_AMD__=1
export HIP_PLATFORM_AMD=1

# -------------------------------
# ROCm / HIP environment
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
export LLVM_HOME=/opt/rocm/llvm/
export PATH=$LLVM_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_HOME/lib:$LLVM_HOME/lib64:$LD_LIBRARY_PATH
export CPATH=$LLVM_HOME/include:$CPATH
export LIBRARY_PATH=$LLVM_HOME/lib:$LLVM_HOME/lib64:$LIBRARY_PATH
export ICON_DATA_PATH=/capstor/scratch/cscs/ybudanaz/beverin/icon-artifacts/velocity/data_r02b05

# -------------------------------
# Build
# -------------------------------
ARCH=${HIP_ARCH:-gfx942}
echo "Building bench_ddt_vn_gpu_sweep for $ARCH"

set -e

hipcc -O3 -std=c++17 \
    --offload-arch=$ARCH \
    -fopenmp \
    -march=native \
    -mtune=native \
    -ffast-math \
    -munsafe-fp-atomics \
    -fgpu-flush-denormals-to-zero \
    -D__HIP_PLATFORM_AMD__=1 -DHIP_PLATFORM_AMD=1 \
    -o bench_ddt_vn_gpu_sweep ddt_vn_gpu.cpp

echo "Build succeeded"

# -------------------------------
# Run (default nlev=90)
# -------------------------------
./bench_ddt_vn_gpu_sweep 90