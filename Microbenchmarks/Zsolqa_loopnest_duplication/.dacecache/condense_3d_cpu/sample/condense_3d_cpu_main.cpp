#include <cstdlib>
#include "../include/condense_3d_cpu.h"

int main(int argc, char **argv) {
    condense_3d_cpuHandle_t handle;
    int err;
    int klev = 42;
    int klon = 42;
    double r4ies = 42;
    double r4les = 42;
    double r5alscp = 42;
    double r5alvcp = 42;
    double retv = 42;
    double rlmin = 42;
    double rthomo = 42;
    double rtice = 42;
    double rtwat = 42;
    double rtwat_rtice_r = 42;
    double * __restrict__ za = (double*) calloc((klev * klon), sizeof(double));
    double * __restrict__ zdqs = (double*) calloc(klon, sizeof(double));
    double * __restrict__ zqsmix = (double*) calloc((klev * klon), sizeof(double));
    double * __restrict__ zqv = (double*) calloc((klev * klon), sizeof(double));
    double * __restrict__ zqxfg = (double*) calloc((5 * klon), sizeof(double));
    double * __restrict__ zsolqa = (double*) calloc((25 * klon), sizeof(double));
    double * __restrict__ ztp1 = (double*) calloc((klev * klon), sizeof(double));


    handle = __dace_init_condense_3d_cpu(klev, klon);
    __program_condense_3d_cpu(handle, za, zdqs, zqsmix, zqv, zqxfg, zsolqa, ztp1, klev, klon, r4ies, r4les, r5alscp, r5alvcp, retv, rlmin, rthomo, rtice, rtwat, rtwat_rtice_r);
    err = __dace_exit_condense_3d_cpu(handle);

    free(za);
    free(zdqs);
    free(zqsmix);
    free(zqv);
    free(zqxfg);
    free(zsolqa);
    free(ztp1);


    return err;
}
