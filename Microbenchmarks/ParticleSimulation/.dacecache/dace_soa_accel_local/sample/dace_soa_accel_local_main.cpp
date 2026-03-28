#include <cstdlib>
#include "../include/dace_soa_accel_local.h"

int main(int argc, char **argv) {
    dace_soa_accel_localHandle_t handle;
    int err;
    int64_t n = 42;
    double rsoft2 = 42;
    double * __restrict__ ax = (double*) calloc(n, sizeof(double));
    double * __restrict__ ay = (double*) calloc(n, sizeof(double));
    double * __restrict__ az = (double*) calloc(n, sizeof(double));
    double * __restrict__ m = (double*) calloc(n, sizeof(double));
    double * __restrict__ rx = (double*) calloc(n, sizeof(double));
    double * __restrict__ ry = (double*) calloc(n, sizeof(double));
    double * __restrict__ rz = (double*) calloc(n, sizeof(double));


    handle = __dace_init_dace_soa_accel_local(n);
    __program_dace_soa_accel_local(handle, ax, ay, az, m, rx, ry, rz, n, rsoft2);
    err = __dace_exit_dace_soa_accel_local(handle);

    free(ax);
    free(ay);
    free(az);
    free(m);
    free(rx);
    free(ry);
    free(rz);


    return err;
}
