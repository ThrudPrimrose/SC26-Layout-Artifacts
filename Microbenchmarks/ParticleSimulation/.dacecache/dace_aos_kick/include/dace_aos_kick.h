#include <dace/dace.h>
typedef void * dace_aos_kickHandle_t;
extern "C" dace_aos_kickHandle_t __dace_init_dace_aos_kick(int64_t n);
extern "C" int __dace_exit_dace_aos_kick(dace_aos_kickHandle_t handle);
extern "C" void __program_dace_aos_kick(dace_aos_kickHandle_t handle, double * __restrict__ a, double * __restrict__ v, double dt, int64_t n);
