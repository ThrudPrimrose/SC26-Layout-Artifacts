#!/bin/bash
#SBATCH --job-name=conjugate_bench
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=03:30:00
#SBATCH --output=conjugate_%j.out
#SBATCH --error=conjugate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288
#SBATCH --exclusive
#SBATCH --account=g177-1
#SBATCH --gpus-per-task=1

echo "Running on $(hostname)"
spack load gcc/76jw6nu
spack load cuda@12.9

# ═══════════════════════════════════════════════════════════
#  System configuration checks
# ═══════════════════════════════════════════════════════════
echo "═══ system checks ═══"
echo "Kernel cmdline (relevant):"
cat /proc/cmdline | tr ' ' '\n' | grep -E "numa_balancing|init_on_alloc|preempt"
echo ""
echo "numa_balancing = $(cat /proc/sys/kernel/numa_balancing 2>/dev/null || echo 'N/A')"
echo "page_size      = $(getconf PAGESIZE)"
echo "cores          = $(nproc)"
echo "NUMA nodes:"
numactl --hardware 2>/dev/null | grep -E "available|node .* cpus" || true
echo ""

# ═══════════════════════════════════════════════════════════
#  Build
# ═══════════════════════════════════════════════════════════
CFLAGS="-O3 -fopenmp -mtune=native -ftree-vectorize -fno-vect-cost-model -march=native -ffast-math -std=c++17"
NVFLAGS="-O3 -arch=sm_90 -std=c++17"

echo "═══ build ═══"
echo "[1/5] conjugate_inplace.cpp    (CPU in-place)"
g++ $CFLAGS -o conjugate_cpu_inplace conjugate_inplace.cpp -lnuma

echo "[2/5] conjugate.cpp            (CPU out-of-place)"
g++ $CFLAGS -o conjugate_cpu_oop conjugate.cpp -lnuma

echo "[3/5] conjugate_inplace.cu     (GPU in-place)"
nvcc $NVFLAGS -o conjugate_gpu_inplace conjugate_inplace.cu

echo "[4/5] conjugate.cu             (GPU out-of-place)"
nvcc $NVFLAGS -o conjugate_gpu_oop conjugate.cu

echo "[5/5] conj_prof_arm.cpp        (CPU profiling)"
g++ $CFLAGS -o conj_prof_arm conj_prof_arm.cpp -lnuma

echo ""

# ═══════════════════════════════════════════════════════════
#  NUMA variance experiment: 1 node vs 4 nodes
#  Same binary, same node, only NUMA scope changes.
#  If 1-NUMA is clean and 4-NUMA is noisy, the cross-chip
#  NVLink coherency hypothesis is confirmed.
# ═══════════════════════════════════════════════════════════
echo "═══ NUMA variance experiment ═══"

echo "--- 1 NUMA node (72 threads, node 0 only) ---"
export OMP_NUM_THREADS=72
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:72:1"
numactl --cpunodebind=0 --membind=0 ./conjugate_cpu_inplace
mv results_cpu_inplace.csv results_cpu_inplace_1numa.csv
echo ""

echo "--- 4 NUMA nodes (288 threads, all nodes) ---"
export OMP_NUM_THREADS=288
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
numactl --cpunodebind=0-3 --membind=0-3 ./conjugate_cpu_inplace
mv results_cpu_inplace.csv results_cpu_inplace_4numa.csv
echo ""

echo "── variance comparison ──"
for f in results_cpu_inplace_1numa.csv results_cpu_inplace_4numa.csv; do
    echo "[$f]"
    awk -F, 'NR>1 {
        k=$1","$2
        n[k]++
        vals[k,n[k]] = $5
    } END {
        for (k in n) {
            cnt = n[k]
            for (i = 1; i <= cnt; i++)
                for (j = i+1; j <= cnt; j++)
                    if (vals[k,i]+0 > vals[k,j]+0) {
                        t = vals[k,i]; vals[k,i] = vals[k,j]; vals[k,j] = t
                    }
            p5  = vals[k, int(cnt*0.05)+1]
            med = vals[k, int(cnt*0.5)+1]
            p95 = vals[k, int(cnt*0.95)+1]
            spread = (med > 0) ? (p95 - p5) / med * 100 : 0
            printf "  P=%-2s %-14s  med=%7.1f GB/s  P5=%7.1f  P95=%7.1f  spread=%5.1f%%\n",
                   substr(k,1,index(k,",")-1),
                   substr(k,index(k,",")+1), med, p5, p95, spread
        }
    }' "$f" | sort
    echo ""
done

# ═══════════════════════════════════════════════════════════
#  Main benchmarks (4 NUMA nodes, 288 threads)
# ═══════════════════════════════════════════════════════════
echo "═══ main benchmarks (4 NUMA, 288 threads) ═══"
export OMP_NUM_THREADS=288
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_DISPLAY_ENV=TRUE

echo "--- ARM PMC profiling ---"
./conj_prof_arm
echo ""

echo "--- CPU in-place ---"
./conjugate_cpu_inplace
echo ""

echo "--- CPU out-of-place ---"
./conjugate_cpu_oop
echo ""

echo "--- GPU in-place ---"
./conjugate_gpu_inplace
echo ""

echo "--- GPU out-of-place ---"
./conjugate_gpu_oop
echo ""

# ═══════════════════════════════════════════════════════════
#  Summary
# ═══════════════════════════════════════════════════════════
echo "═══ summary (averages) ═══"
for f in results_cpu_inplace.csv results_cpu_oop.csv results_gpu_inplace.csv results_gpu_oop.csv; do
    [ -f "$f" ] || continue
    echo "[$f]"
    awk -F, 'NR>1{s[$1","$2]+=$4; g[$1","$2]+=$5; n[$1","$2]++}
            END{for(k in s) printf "  P=%-2s %-14s avg=%8.4f ms  %7.1f GB/s\n",
                substr(k,1,index(k,",")-1), substr(k,index(k,",")+1),
                s[k]/n[k], g[k]/n[k]}' "$f" | sort -t= -k1,1n -k2
    echo ""
done

# ═══════════════════════════════════════════════════════════
#  Collect results
# ═══════════════════════════════════════════════════════════
mkdir -p results/daint
mv results_cpu_inplace.csv results_cpu_oop.csv \
   results_gpu_inplace.csv results_gpu_oop.csv \
   results_cpu_inplace_1numa.csv results_cpu_inplace_4numa.csv \
   results_cpu_inplace_prof.csv results_cpu_inplace_pmc.csv \
   results/daint/ 2>/dev/null

echo ""
echo "═══ done ═══"
echo "Results in results/daint/"
ls -la results/daint/