#include <dace/dace.h>
typedef void * dace_aos_accel_bufHandle_t;
extern "C" dace_aos_accel_bufHandle_t __dace_init_dace_aos_accel_buf(int64_t n);
extern "C" int __dace_exit_dace_aos_accel_buf(dace_aos_accel_bufHandle_t handle);
extern "C" void __program_dace_aos_accel_buf(dace_aos_accel_bufHandle_t handle, double * __restrict__ a, double * __restrict__ fx, double * __restrict__ fy, double * __restrict__ fz, double * __restrict__ m, double * __restrict__ r, int64_t n, double rsoft2);
