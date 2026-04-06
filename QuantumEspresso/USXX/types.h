#pragma once
#include <complex>
#include <cmath>

// ============================================================
// Platform detection and includes
// ============================================================
#if defined(__CUDACC__)
    // NVIDIA CUDA
    #include <cuda_runtime.h>
    #include <cuComplex.h>
    using Complex_DP = cuDoubleComplex;
    #define GPU_HOST_DEVICE __host__ __device__
#elif defined(__HIP_PLATFORM_AMD__) || defined(__HIPCC__)
    // AMD HIP
    #include <hip/hip_runtime.h>
    #include <hip/hip_complex.h>
    using Complex_DP = hipDoubleComplex;
    #define GPU_HOST_DEVICE __host__ __device__
#else
    // Plain host compiler (e.g. g++ compiling a .cpp that includes this)
    struct Complex_DP { double x, y; };
    #define GPU_HOST_DEVICE
#endif

using DP = double;
using HostComplex = std::complex<double>;

// ============================================================
// Complex helpers — use .x/.y fields directly so we don't
// depend on cuCmul/hipCmul which vary across SDK versions.
// Both cuDoubleComplex and hipDoubleComplex are { double x, y; }.
// ============================================================
GPU_HOST_DEVICE inline Complex_DP make_cmplx(double r, double i) {
    Complex_DP c;
    c.x = r; c.y = i;
    return c;
}

GPU_HOST_DEVICE inline double creal_val(Complex_DP a) { return a.x; }
GPU_HOST_DEVICE inline double cimag_val(Complex_DP a) { return a.y; }

GPU_HOST_DEVICE inline Complex_DP cmul(Complex_DP a, Complex_DP b) {
    Complex_DP c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

GPU_HOST_DEVICE inline Complex_DP cadd(Complex_DP a, Complex_DP b) {
    Complex_DP c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

GPU_HOST_DEVICE inline Complex_DP csub(Complex_DP a, Complex_DP b) {
    Complex_DP c;
    c.x = a.x - b.x;
    c.y = a.y - b.y;
    return c;
}

GPU_HOST_DEVICE inline Complex_DP cconj(Complex_DP a) {
    Complex_DP c;
    c.x = a.x;
    c.y = -a.y;
    return c;
}

GPU_HOST_DEVICE inline double cabs_val(Complex_DP a) {
    return sqrt(a.x * a.x + a.y * a.y);
}

static constexpr double PI  = 3.14159265358979323846;
static constexpr double TPI = 2.0 * PI;