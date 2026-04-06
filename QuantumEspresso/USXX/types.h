#pragma once
#include <complex>
#include <cuda_runtime.h>
#include <cuComplex.h>

using DP = double;
using Complex_DP = cuDoubleComplex;

// Host-side complex
using HostComplex = std::complex<double>;

__host__ __device__ inline Complex_DP make_cmplx(double r, double i) {
    return make_cuDoubleComplex(r, i);
}

__host__ __device__ inline Complex_DP cmul(Complex_DP a, Complex_DP b) {
    return cuCmul(a, b);
}

__host__ __device__ inline Complex_DP cadd(Complex_DP a, Complex_DP b) {
    return cuCadd(a, b);
}

__host__ __device__ inline Complex_DP cconj(Complex_DP a) {
    return cuConj(a);
}

__host__ __device__ inline double cabs_val(Complex_DP a) {
    double r = cuCreal(a);
    double i = cuCimag(a);
    return sqrt(r * r + i * i);
}

__host__ __device__ inline Complex_DP csub(Complex_DP a, Complex_DP b) {
    return make_cuDoubleComplex(cuCreal(a) - cuCreal(b), cuCimag(a) - cuCimag(b));
}

static constexpr double PI  = 3.14159265358979323846;
static constexpr double TPI = 2.0 * PI;
