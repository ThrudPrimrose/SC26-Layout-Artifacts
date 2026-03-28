#include <dace/dace.h>
typedef void * dace_soa_driftHandle_t;
extern "C" dace_soa_driftHandle_t __dace_init_dace_soa_drift(int64_t n);
extern "C" int __dace_exit_dace_soa_drift(dace_soa_driftHandle_t handle);
extern "C" void __program_dace_soa_drift(dace_soa_driftHandle_t handle, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n);
