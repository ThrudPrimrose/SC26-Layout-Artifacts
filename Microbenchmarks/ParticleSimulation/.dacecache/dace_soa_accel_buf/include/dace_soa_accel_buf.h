#include <dace/dace.h>
typedef void * dace_soa_accel_bufHandle_t;
extern "C" dace_soa_accel_bufHandle_t __dace_init_dace_soa_accel_buf(int64_t n);
extern "C" int __dace_exit_dace_soa_accel_buf(dace_soa_accel_bufHandle_t handle);
extern "C" void __program_dace_soa_accel_buf(dace_soa_accel_bufHandle_t handle, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ fx, double * __restrict__ fy, double * __restrict__ fz, double * __restrict__ m, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, int64_t n, double rsoft2);
