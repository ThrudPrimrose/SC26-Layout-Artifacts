#!/bin/bash
# sweep.sh — sweep block/tile sizes for syr2k modes 3 & 4
# Usage: ./sweep.sh [N=1024] [M=1024] [CSV=results.csv] [NRUNS=100]

N=${1:-1024}; M=${2:-1024}; CSV=${3:-results.csv}; NRUNS=${4:-100}
SRC=syr2k_bench.cpp
CXX="g++ -O3 -fopenmp -std=c++17 -march=native -mtune=native -ffast-math"
SIZES=(32 64 128)

echo "N=$N M=$M NRUNS=$NRUNS -> $CSV"

for BI in "${SIZES[@]}"; do
 for BJ in "${SIZES[@]}"; do
  for BK in "${SIZES[@]}"; do
    BIN=".syr2k_${BI}_${BJ}_${BK}"
    $CXX -DSZ_N=$N -DSZ_M=$M -DBI=$BI -DBJ=$BJ -DBK=$BK \
          -DTI=$BI -DTJ=$BJ -DTK=$BK -DNRUNS=$NRUNS \
          -o "$BIN" "$SRC" 2>/dev/null || { echo "SKIP $BI $BJ $BK"; continue; }
    echo "BI=$BI BJ=$BJ BK=$BK"
    ./"$BIN" 3 "$CSV"
    ./"$BIN" 4 "$CSV"
    rm -f "$BIN"
  done
 done
done

echo "done. $(wc -l < "$CSV") rows in $CSV"