#pragma once

// ============================================================
// Thin compat layer: cuda* API -> hip* API on AMD
// On NVIDIA this header is a no-op.
// ============================================================

#if defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)

#include <hip/hip_runtime.h>

// Memory
#define cudaMalloc      hipMalloc
#define cudaFree        hipFree
#define cudaMemcpy      hipMemcpy
#define cudaMemset      hipMemset
#define cudaMemcpyHostToDevice   hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost   hipMemcpyDeviceToHost

// Events
#define cudaEvent_t             hipEvent_t
#define cudaEventCreate         hipEventCreate
#define cudaEventRecord         hipEventRecord
#define cudaEventSynchronize    hipEventSynchronize
#define cudaEventElapsedTime    hipEventElapsedTime
#define cudaEventDestroy        hipEventDestroy

// Sync / errors
#define cudaDeviceSynchronize   hipDeviceSynchronize
#define cudaGetLastError        hipGetLastError
#define cudaError_t             hipError_t
#define cudaSuccess             hipSuccess

#elif defined(__CUDACC__)

#include <cuda_runtime.h>

#endif