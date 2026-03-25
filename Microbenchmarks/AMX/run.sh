#source /opt/intel/oneapi/setvars.sh

for NTHR in 1 2; do
    export OMP_NUM_THREADS=$NTHR
    export MKL_NUM_THREADS=$NTHR
    export MKL_DYNAMIC=FALSE
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    echo "Running with $NTHR threads"

    g++-14 -O3 -fopenmp -march=graniterapids -mamx-tile -mamx-bf16 -mamx-int8 \
        -I. -I${MKLROOT}/include -D_NTHREADS=$OMP_NUM_THREADS \
        -o bench1 amx_gemm_base.cpp \
        -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
        -lgomp -lpthread -lm

    g++-14 -O3 -fopenmp -march=graniterapids -mamx-tile -mamx-bf16 -mamx-int8 \
        -I. -I${MKLROOT}/include -D_NTHREADS=$OMP_NUM_THREADS \
        -o bench2 amx_gemm_layouts.cpp \
        -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core \
        -lgomp -lpthread -lm

    for S in 1024 2048 4096; do
        echo "Running bench1 with size $S"
        ./bench1 $S $S $S
        ./bench2 $S $S $S -t 2 -i 100

        mv bench_amx_gemm.csv bench_amx_gemm_${NTHR}_threads_${S}_MNK.csv
        mv bench_amx_layouts.csv bench_amx_layouts_${NTHR}_threads_${S}_MNK.csv
    done
done



