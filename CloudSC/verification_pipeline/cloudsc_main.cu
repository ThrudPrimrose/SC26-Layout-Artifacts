#include <iostream>
#include <vector>
#include <cuda_runtime.h>
#include <chrono>

// DaCe generated header
#include "cloudsc_py.h"

int main() {
    // Symbols (SDFG Constants/Symbols)
    int klon = 128;
    int klev = 64;
    int nclv = 5;
    int ncldqi = 2;
    int ncldql = 1;
    int ncldqr = 3;
    int ncldqs = 4;
    int ncldqv = 5;

    // Allocate host memory for the mandatory arrays (based on SDFG)
    // For simplicity, we allocate everything with double precision as per the SDFG inspection
    std::vector<double> pt(klon * klev, 280.0);
    std::vector<double> pq(klon * klev, 0.01);
    std::vector<double> tendency_loc_t(klon * klev, 0.0);
    
    // Scalars (passed by pointer/array in DaCe generated init)
    int kidia = 1;
    int kfdia = klon;
    double ptsphy = 3600.0;

    // Allocate GPU memory
    double *d_pt, *d_pq, *d_tendency_loc_t;
    cudaMalloc(&d_pt, klon * klev * sizeof(double));
    cudaMalloc(&d_pq, klon * klev * sizeof(double));
    cudaMalloc(&d_tendency_loc_t, klon * klev * sizeof(double));

    cudaMemcpy(d_pt, pt.data(), klon * klev * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_pq, pq.data(), klon * klev * sizeof(double), cudaMemcpyHostToDevice);

    std::cout << "Initializing DaCe CloudSC program..." << std::endl;
    // We initialize the DaCe runtime. 
    // The handle depends on the specific DaCe generated entry point.
    auto *handle = __dace_init_cloudsc_py(
        &kidia, &kfdia, &ptsphy, 
        d_pt, d_pq, 
        d_tendency_loc_t, // and many other arrays...
        // Note: The actual call needs ALL the non-transient arrays from the SDFG.
        // For a full standalone, we would need to provide all of them here.
        klon, klev, nclv, ncldqi, ncldql, ncldqr, ncldqs, ncldqv
    );

    std::cout << "Launching CloudSC kernel..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    
    __program_cloudsc_py(handle, 
        &kidia, &kfdia, &ptsphy, 
        d_pt, d_pq, 
        d_tendency_loc_t, // ... additional arrays
        klon, klev, nclv, ncldqi, ncldql, ncldqr, ncldqs, ncldqv
    );
    
    cudaDeviceSynchronize();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Execution time: " << diff.count() << " seconds" << std::endl;

    __dace_exit_cloudsc_py(handle);

    // Free resources
    cudaFree(d_pt);
    cudaFree(d_pq);
    cudaFree(d_tendency_loc_t);

    std::cout << "Done." << std::endl;
    return 0;
}
