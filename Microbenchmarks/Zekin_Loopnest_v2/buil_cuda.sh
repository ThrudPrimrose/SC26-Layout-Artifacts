#!/bin/bash
set -e
ARCH=${HIP_ARCH:-gfx942}
echo "Building for $ARCH"
nvcc -O3 -std=c++17 \
    -arch=sm_90 \
    -Xcompiler=-fopenmp \
    -o bench_gpu bench_gpu.cu
echo "Done: ./bench_gpu"


#!/bin/bash
set -e

export OMP_NUM_THREADS=288
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
echo "Done: OMP_NUM_THREADS=$OMP_NUM_THREADS OMP_PROC_BIND=$OMP_PROC_BIND OMP_PLACES=$OMP_PLACES ./bench_cpu"