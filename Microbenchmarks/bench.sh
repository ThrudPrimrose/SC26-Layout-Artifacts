# compile with desired sizes
nvcc -O3 -Xcompiler=-fopenmp -march=native -mtune=native -ffast-math \
    -std=c++17 \
    -o segmented_reduction segmented_reduction.cu

# run a mode (results append to CSV)
./segmented_reduction