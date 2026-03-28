#include <dace/dace.h>
typedef void * condense_3d_cpuHandle_t;
extern "C" condense_3d_cpuHandle_t __dace_init_condense_3d_cpu(int klev, int klon);
extern "C" int __dace_exit_condense_3d_cpu(condense_3d_cpuHandle_t handle);
extern "C" void __program_condense_3d_cpu(condense_3d_cpuHandle_t handle, double * __restrict__ za, double * __restrict__ zdqs, double * __restrict__ zqsmix, double * __restrict__ zqv, double * __restrict__ zqxfg, double * __restrict__ zsolqa, double * __restrict__ ztp1, int klev, int klon, double r4ies, double r4les, double r5alscp, double r5alvcp, double retv, double rlmin, double rthomo, double rtice, double rtwat, double rtwat_rtice_r);
