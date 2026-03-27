#!/bin/bash
#!/bin/bash
#SBATCH --job-name=numa_bench
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=01:30:00
#SBATCH --output=numa_bench_%j.out
#SBATCH --error=numa_bench_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --exclusive

PLATFORM="beverin"
BINARY="./numa_calibrate"


echo "Building..."
gcc -O3 -std=c++17 \
    -march=native \
    -fopenmp \
    -ffast-math \
    -mtune=native \
    -fno-vect-cost-model \
    -ftree-vectorize \
     -o numa_calibrate numa_calibrate.c -lm


gcc -O3 -std=c++17 \
    -march=native \
    -fopenmp \
    -ffast-math \
    -mtune=native \
    -fno-vect-cost-model \
    -ftree-vectorize \
     -o numa_triad numa_triad.c -lm

# ================================================================
# Platform-specific configuration
# ================================================================

echo "=== Beverin (MI300A Zen 4 cores) ==="
echo "  Expected: β=64 (4KB pages / 64B CL)"
echo "  Expected: α≈?, γ≈? (cross-CCD)"
echo ""


# ================================================================
# Print topology
# ================================================================

spack load python@3.13.8
source ${SCRATCH}/yakup-dev-env/bin/activate

# === ROCm / HIP ===
export ROCM_HOME=/opt/rocm
export HIP_PATH=$ROCM_HOME
export HIPCC=$ROCM_HOME/bin/hipcc
export __HIP_PLATFORM_AMD__=1
export HIP_PLATFORM_AMD=1
export ARCH=gfx942
export HCC_AMDGPU_TARGET=$ARCH
export CUPY_HIPCC_GENERATE_CODE=--offload-arch=$ARCH
export CUPY_INSTALL_USE_HIP=1

# === LLVM (ROCm bundled) ===
export LLVM_HOME=$ROCM_HOME/llvm

# === GPU runtime ===
export HSA_ENABLE_SDMA=0
export HSA_XNACK=1
export GPU_MAX_HEAP_SIZE=100

# === DaCe ===
export _DACE_NO_SYNC=1

# === OpenMP ===
export OMP_NUM_THREADS=96
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_SCHEDULE=static
export OMP_PROC_BIND=close
export SLURM_CPU_BIND=cores

# === Python includes ===
_PYINC=$(python3.13 -c "import sysconfig; print(sysconfig.get_path('include'))")

# === Paths (ROCm + LLVM + $SCRATCH + Python) ===
export PATH=$SCRATCH/bin:$ROCM_HOME/bin:$LLVM_HOME/bin:$PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$ROCM_HOME/lib:$ROCM_HOME/lib64:$LLVM_HOME/lib:$LLVM_HOME/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$ROCM_HOME/lib:$ROCM_HOME/lib64:$LLVM_HOME/lib:$LLVM_HOME/lib64:$LIBRARY_PATH
export CPATH=$_PYINC:$SCRATCH/include:$ROCM_HOME/include:$LLVM_HOME/include:$CPATH
export C_INCLUDE_PATH=$_PYINC:$SCRATCH/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CPLUS_INCLUDE_PATH
export CFLAGS="-I$ROCM_HOME/include -I$_PYINC"
export LDFLAGS="-L$ROCM_HOME/lib -L$ROCM_HOME/lib64"

export BEVERIN=1


echo "--- System topology ---"
numactl --hardware 2>/dev/null || echo "(numactl not available)"
echo ""
lscpu | grep -E "^(Thread|Core|Socket|NUMA|CPU\(s\)|Model name)" || true
echo ""
getconf PAGESIZE 2>/dev/null && echo "Page size: $(getconf PAGESIZE) bytes" || true
echo ""

# ================================================================
# Run each benchmark
# ================================================================

echo "================================================================"
echo "  Benchmark 1: Stride-sweep bandwidth (β detection)"
echo "================================================================"
$BINARY stride_bw | tee stride_bw_${PLATFORM}.txt

echo "================================================================"
echo "  Benchmark 2: Pointer-chase latency sweep (α calibration)"
echo "================================================================"
$BINARY stride_lat | tee stride_lat_${PLATFORM}.txt

echo "================================================================"
echo "  Benchmark 3: Local vs Remote NUMA bandwidth (γ_bw)"
echo "================================================================"
$BINARY numa_bw | tee numa_bw_${PLATFORM}.txt

echo "================================================================"
echo "  Benchmark 4: Local vs Remote NUMA latency (α, γ_lat)"
echo "================================================================"
$BINARY numa_lat | tee numa_lat_${PLATFORM}.txt

echo "================================================================"
echo "  Benchmark 5: NUMA BW matrix (thread → owner)"
echo "================================================================"
$BINARY numa_matrix | tee numa_matrix_${PLATFORM}.txt


echo "================================================================"
echo "  Benchmark 5: NUMA Triad copy"
echo "================================================================"
$BINARY numa_triad | tee numa_triad_${PLATFORM}.txt

echo ""
echo "=== Summary ==="
echo "Results saved to: stride_bw_${PLATFORM}.txt, stride_lat_${PLATFORM}.txt,"
echo "  numa_bw_${PLATFORM}.txt, numa_lat_${PLATFORM}.txt, numa_matrix_${PLATFORM}.txt"
echo ""
echo "To extract parameters:"
echo "  β  = stride (in CLs) where stride_lat shows a latency jump"
echo "  α  = lat(stride < β) / lat(stride = β)"
echo "  γ  = lat_remote / lat_local  (from numa_lat)"

# MI300A Zen 4:  β = 4,  α = 0.38,  γ = 2.0
# Grace Neoverse V2: β = 1024,  α = 0.52,  γ = 2.7