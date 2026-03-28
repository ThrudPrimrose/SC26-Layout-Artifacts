#include <dace/dace.h>
typedef void * dace_aosoa_accel_localHandle_t;
extern "C" dace_aosoa_accel_localHandle_t __dace_init_dace_aosoa_accel_local(int64_t nb, int64_t vl);
extern "C" int __dace_exit_dace_aosoa_accel_local(dace_aosoa_accel_localHandle_t handle);
extern "C" void __program_dace_aosoa_accel_local(dace_aosoa_accel_localHandle_t handle, double * __restrict__ a, double * __restrict__ m, double * __restrict__ r, int64_t N_val, int64_t nb, double rsoft2, int64_t vl);
