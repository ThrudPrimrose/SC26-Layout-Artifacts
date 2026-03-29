#!/bin/bash
#SBATCH --job-name=conj_bench
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=00:50:00
#SBATCH --output=conj_%j.out
#SBATCH --error=conj_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192
#SBATCH --mem-bind=local

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=96
export OMP_PROC_BIND=close
export OMP_PLACES="{0}:24:1,{24}:24:1,{48}:24:1,{72}:24:1"
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

set -e

CFLAGS="-O3 -fopenmp -mtune=native -ftree-vectorize -fno-vect-cost-model -march=native -ffast-math -std=c++17"
HIPFLAGS="--offload-arch=gfx942 -O3 -ffast-math -std=c++17"

echo "═══ build ═══"
echo "[0/5] bench_triad.cpp       (CPU in-place)"
g++ $CFLAGS -o bench_triad bench_triad.cpp -lnuma
g++ $CFLAGS -o bench_triad_local bench_triad_local.cpp -lnuma

echo "[1/5] conjugate_inplace.cpp       (CPU in-place)"
g++ $CFLAGS -o conj_ip_cpu conjugate_inplace.cpp -lnuma

echo "[2/5] conjugate.cpp               (CPU out-of-place)"
g++ $CFLAGS -o conj_oop_cpu conjugate.cpp -lnuma

echo "[3/5] conjugate_inplace_hip.cpp   (GPU in-place)"
hipcc $HIPFLAGS -o conj_ip_gpu conjugate_inplace_hip.cpp

echo "[4/5] conjugate_hip.cpp           (GPU out-of-place)"
hipcc $HIPFLAGS -o conj_oop_gpu conjugate_hip.cpp

echo "[5/5] bench_streams.cpp           (CPU stream concurrency)"
g++ $CFLAGS -o bench_streams bench_streams.cpp

echo ""
echo "═══ run ═══"

echo "--- Triad ---"
./bench_triad_local
./bench_triad
echo ""

echo "--- CPU in-place ---"
./conj_ip_cpu
echo ""

echo "--- CPU out-of-place ---"
./conj_oop_cpu
echo ""

echo "--- GPU in-place ---"
./conj_ip_gpu
echo ""

echo "--- GPU out-of-place ---"
./conj_oop_gpu
echo ""

echo "--- CPU stream concurrency ---"
./bench_streams
echo ""

echo "═══ CSV files ═══"
echo "  results_cpu_ip.csv      (CPU in-place,        per-run)"
echo "  results_cpu.csv         (CPU out-of-place,    per-run)"
echo "  results_gpu_ip.csv      (GPU in-place,        per-run)"
echo "  results_gpu_oop.csv     (GPU out-of-place,    per-run)"
echo "  results_streams.csv     (CPU stream conc.,    per-run)"
echo ""

echo "═══ summary (averages) ═══"
for f in results_cpu_ip.csv results_cpu.csv results_gpu_ip.csv results_gpu_oop.csv results_streams.csv; do
    [ -f "$f" ] || continue
    echo "[$f]"
    awk -F, 'NR>1{s[$1","$2]+=$4; g[$1","$2]+=$5; n[$1","$2]++}
        END{for(k in s) printf "  %-20s avg=%8.4f ms  %7.1f GB/s\n",
            k, s[k]/n[k], g[k]/n[k]}' "$f" | sort
    echo ""
done

mkdir -p results/beverin
mv results_cpu_ip.csv results/beverin/
mv results_cpu.csv results/beverin/
mv results_gpu_ip.csv results/beverin/
mv results_gpu_oop.csv results/beverin/
mv results_streams.csv results/beverin/