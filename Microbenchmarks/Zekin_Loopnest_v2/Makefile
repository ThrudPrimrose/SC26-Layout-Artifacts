# Adjust ARCH for your GPU:
#   sm_70 = V100,  sm_80 = A100,  sm_90 = H100
ARCH ?= sm_80

CXX       = g++
CXXFLAGS  = -O3 -fopenmp -march=native -std=c++17

NVCC      = nvcc
NVCCFLAGS = -O3 -arch=$(ARCH) -std=c++17 -Xcompiler "-fopenmp -O3 -march=native"

all: bench_cpu bench_gpu

bench_cpu: bench_cpu.cpp bench_common.h
	$(CXX) $(CXXFLAGS) -o $@ bench_cpu.cpp

bench_gpu: bench_gpu.cu bench_common.h
	$(NVCC) $(NVCCFLAGS) -o $@ bench_gpu.cu -lgomp

clean:
	rm -f bench_cpu bench_gpu z_v_grad_w_cpu.csv z_v_grad_w_gpu.csv

run_cpu: bench_cpu
	./bench_cpu

run_gpu: bench_gpu
	./bench_gpu

run: run_cpu run_gpu

.PHONY: all clean run run_cpu run_gpu
