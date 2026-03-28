#include <dace/dace.h>
typedef void * dace_aosoa_driftHandle_t;
extern "C" dace_aosoa_driftHandle_t __dace_init_dace_aosoa_drift(int64_t nb, int64_t vl);
extern "C" int __dace_exit_dace_aosoa_drift(dace_aosoa_driftHandle_t handle);
extern "C" void __program_dace_aosoa_drift(dace_aosoa_driftHandle_t handle, double * __restrict__ r, double * __restrict__ v, double dt, int64_t nb, int64_t vl);
