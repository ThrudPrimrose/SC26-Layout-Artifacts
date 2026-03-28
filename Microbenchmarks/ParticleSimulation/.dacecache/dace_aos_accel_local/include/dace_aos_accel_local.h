#include <dace/dace.h>
typedef void * dace_aos_accel_localHandle_t;
extern "C" dace_aos_accel_localHandle_t __dace_init_dace_aos_accel_local(int64_t n);
extern "C" int __dace_exit_dace_aos_accel_local(dace_aos_accel_localHandle_t handle);
extern "C" void __program_dace_aos_accel_local(dace_aos_accel_localHandle_t handle, double * __restrict__ a, double * __restrict__ m, double * __restrict__ r, int64_t n, double rsoft2);
