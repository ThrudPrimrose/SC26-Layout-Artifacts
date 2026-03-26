export OMP_NUM_THREADS=96
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

g++ -O3 -fopenmp -march=native -std=c++17 profile_cpu.cpp -o profile_cpu_simple

perf stat -e cycles,instructions,cache-misses,cache-references,\
L1-dcache-load-misses,L1-dcache-loads,\
LLC-load-misses,LLC-loads,\
dTLB-load-misses,dTLB-loads \
./profile_cpu_simple