// fp16_operators.h — Mixed-type arithmetic operators for dace::float16 / double.
// Resolves ambiguous operator errors when fp16 and fp64 values appear in the
// same expression (e.g. inside DaCe-generated tasklets with mixed-precision inputs).
// All mixed ops promote to float (not double) so GPU ALU stays in fp32.
#pragma once

#include <dace/types.h>
#include <type_traits>

#ifdef __CUDACC__

__device__ __host__ inline dace::float16 abs(dace::float16 a) {
    return __habs(a);
}
__device__ __host__ inline bool operator<(dace::float16 a, double b) {
   return a < static_cast<dace::float16>(b);
} 
__device__ __host__ inline bool operator/(float a, dace::float16 b) {
    return static_cast<dace::float16>(a) / b;
}

#endif // __CUDACC__
