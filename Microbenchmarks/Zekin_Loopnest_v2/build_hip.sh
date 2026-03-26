#!/bin/bash
set -e
ARCH=${HIP_ARCH:-gfx942}
echo "Building for $ARCH"
hipcc -O3 -std=c++17 \
    --offload-arch=$ARCH \
    -fopenmp \
    -march=native \
    -ffast-math \
    -funroll-loops \
    -munsafe-fp-atomics \
    -mllvm -amdgpu-early-inline-all=true \
    -D__HIP_PLATFORM_AMD__=1 -DHIP_PLATFORM_AMD=1 \
    -ffast-math --offload-arch=$ARCH -fopenmp \
    -o bench_gpu_hip bench_gpu_hip.cpp
echo "Done: ./bench_gpu_hip"

export LLVM_HOME=/opt/rocm/llvm/
export PATH=$LLVM_HOME/bin:$PATH
export LD_LIBRARY_PATH=$LLVM_HOME/lib:$LLVM_HOME/lib64:$LD_LIBRARY_PATH
export CPATH=$LLVM_HOME/include:$CPATH
export LIBRARY_PATH=$LLVM_HOME/lib:$LLVM_HOME/lib64:$LIBRARY_PATH

#!/bin/bash
set -e

export OMP_NUM_THREADS=96
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

echo "Building bench_cpu (OMP_NUM_THREADS=$OMP_NUM_THREADS)"
g++ -O3 -std=c++17 \
    -march=native \
    -fopenmp \
    -ffast-math \
    -funroll-loops \
    -ftree-vectorize \
    -o bench_cpu bench_cpu.cpp
echo "Done: OMP_NUM_THREADS=96 OMP_PROC_BIND=spread OMP_PLACES=cores ./bench_cpu"