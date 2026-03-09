# compile with desired sizes
g++ -O3 -fopenmp -march=native -mtune=native -ffast-math \
    -std=c++17 -DSZ_N=1024 -DSZ_M=1024 \
    -DBI=64 -DBJ=64 -DBK=64 -DTI=32 -DTJ=32 -DTK=32 \
    -o syr2k_bench syr2k_bench.cpp

# run a mode (results append to CSV)
./syr2k_bench 1 results.csv   # 8 layout permutations (RRR..CCC)
./syr2k_bench 2 results.csv   # 6 loop orderings (ikj,ijk,kij,kji,jik,jki)
./syr2k_bench 3 results.csv   # blocked arrays+loops
./syr2k_bench 4 results.csv   # tiled loops, row-major arrays