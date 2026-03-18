#Compile with desired sizes
nvcc -arch=sm_90a -O3 -Xcompiler=-fopenmp  \
    -std=c++17 \
    -o segmented_reduction segmented_reduction.cu

# run a mode (results append to CSV)
./segmented_reduction
