export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_NUM_THREADS=2

g++ -O3 -fopenmp -march=native \
    -mamx-tile -mamx-bf16 -mamx-int8 -mavx512f -mavx512bf16 -D_NTHREADS=$OMP_NUM_THREADS \
    -I. -L/usr/lib/x86_64-linux-gnu/openblas-pthread/ \
    -o bench1 amx_gemm_base.cpp -lopenblas

g++ -O3 -fopenmp -march=native -mamx-tile -mamx-bf16 -mamx-int8 -D_NTHREADS=$OMP_NUM_THREADS \
    -I. -o bench2 amx_gemm_layouts.cpp \
    -L/usr/lib/x86_64-linux-gnu/openblas-pthread/ -lopenblas

./bench1 256 256 256
./bench2 256 256 256 -t $OMP_NUM_THREADS -i 100