#include <dace/dace.h>
typedef void * dace_aosoa_kickHandle_t;
extern "C" dace_aosoa_kickHandle_t __dace_init_dace_aosoa_kick(int64_t nb, int64_t vl);
extern "C" int __dace_exit_dace_aosoa_kick(dace_aosoa_kickHandle_t handle);
extern "C" void __program_dace_aosoa_kick(dace_aosoa_kickHandle_t handle, double * __restrict__ a, double * __restrict__ v, double dt, int64_t nb, int64_t vl);
