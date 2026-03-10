#!/bin/bash
#SBATCH --job-name=syr2k
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --time=04:00:00
#SBATCH --output=job_output_syr2k.txt
#SBATCH --error=job_error_syr2k.txt


# sweep.sh — sweep block/tile sizes for syr2k modes 1-5
# Usage: ./sweep.sh [N=1024] [M=1024] [CSV=results.csv] [NRUNS=100]
spack load gcc@14
spack load openblas
spack load cuda

N=${1:-8192}; M=${2:-8192}; CSV=${3:-results.csv}; NRUNS=${4:-50}
SRC=syr2k_bench.cpp
CXX="g++ -O3 -fopenmp -std=c++17 -march=native -mtune=native -ffast-math -I$(spack location -i openblas)/include -L$(spack location -i openblas)/lib "
SIZES=(32 64 128 256)

echo "N=$N M=$M NRUNS=$NRUNS -> $CSV"

# Fixed modes (1,2,5) – compile once
FIXED_BIN=".syr2k_fixed"
echo "Compiling fixed binary for modes 1,2,5..."
$CXX -DSZ_N=$N -DSZ_M=$M -DTI=32 -DTJ=32 -DTK=32 -DNRUNS=$NRUNS -o "$FIXED_BIN" "$SRC" -lopenblas  || {
    echo "ERROR: compilation failed for fixed binary"
    exit 1
}
echo "Running mode 5 (OpenBLAS)..."
./"$FIXED_BIN" 5 "$CSV"
echo "Running mode 1 (layout permutations)..."
./"$FIXED_BIN" 1 "$CSV"
echo "Running mode 2 (loop orderings)..."
./"$FIXED_BIN" 2 "$CSV"
rm -f "$FIXED_BIN"

# Sweep modes 3 & 4 over block/tile sizes
for BI in "${SIZES[@]}"; do
 for BJ in "${SIZES[@]}"; do
  for BK in "${SIZES[@]}"; do
    BIN=".syr2k_${BI}_${BJ}_${BK}"
    echo "Compiling for BI=$BI BJ=$BJ BK=$BK ..."
    $CXX -DSZ_N=$N -DSZ_M=$M -DBI=$BI -DBJ=$BJ -DBK=$BK \
          -DTI=$BI -DTJ=$BJ -DTK=$BK -DNRUNS=$NRUNS \
          -o "$BIN" "$SRC" -lopenblas  2>/dev/null || {
              echo "SKIP $BI $BJ $BK (compilation failed)"
              continue
          }
    echo "  Running mode 3 (blocked arrays + loops)..."
    ./"$BIN" 3 "$CSV"
    echo "  Running mode 4 (tiled loops)..."
    ./"$BIN" 4 "$CSV"
    rm -f "$BIN"
  done
 done
done

echo "done. $(wc -l < "$CSV") rows in $CSV"