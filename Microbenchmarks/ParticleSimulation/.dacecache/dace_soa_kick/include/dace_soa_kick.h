#include <dace/dace.h>
typedef void * dace_soa_kickHandle_t;
extern "C" dace_soa_kickHandle_t __dace_init_dace_soa_kick(int64_t n);
extern "C" int __dace_exit_dace_soa_kick(dace_soa_kickHandle_t handle);
extern "C" void __program_dace_soa_kick(dace_soa_kickHandle_t handle, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n);
