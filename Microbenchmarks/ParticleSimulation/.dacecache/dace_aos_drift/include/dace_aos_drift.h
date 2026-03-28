#include <dace/dace.h>
typedef void * dace_aos_driftHandle_t;
extern "C" dace_aos_driftHandle_t __dace_init_dace_aos_drift(int64_t n);
extern "C" int __dace_exit_dace_aos_drift(dace_aos_driftHandle_t handle);
extern "C" void __program_dace_aos_drift(dace_aos_driftHandle_t handle, double * __restrict__ r, double * __restrict__ v, double dt, int64_t n);
