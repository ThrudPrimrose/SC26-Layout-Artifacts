
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#define CK(x) do { hipError_t _e=(x); if(_e){ \
    fprintf(stderr,"HIP error %s:%d: %s\n",__FILE__,__LINE__,hipGetErrorString(_e)); \
    exit(1); } } while(0)
__global__ void kc(const double* s, double* d, long long n) {
    long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i < n) d[i] = s[i];
}
int main() {
    const long long SZ = 128LL * 1024 * 1024;
    const size_t bytes = SZ * sizeof(double);
    double *a, *b;
    CK(hipMalloc(&a, bytes)); CK(hipMalloc(&b, bytes));
    hipEvent_t s, e; CK(hipEventCreate(&s)); CK(hipEventCreate(&e));
    const int threads = 256;
    const int blocks  = (int)((SZ + threads - 1) / threads);
    kc<<<blocks, threads>>>(a, b, SZ);
    CK(hipGetLastError()); CK(hipDeviceSynchronize());
    CK(hipEventRecord(s));
    for (int i = 0; i < 20; i++) { kc<<<blocks, threads>>>(a, b, SZ); CK(hipGetLastError()); }
    CK(hipEventRecord(e)); CK(hipEventSynchronize(e));
    float ms; CK(hipEventElapsedTime(&ms, s, e));
    if (ms < 1e-3f) { fprintf(stderr, "timing error: ms=%.6f\n", ms); return 1; }
    printf("%.1f\n", 20.0 * 2 * SZ * sizeof(double) / (ms / 1000.0) / 1e9);
    CK(hipFree(a)); CK(hipFree(b));
    return 0;
}