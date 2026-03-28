#include <cstdlib>
#include "../include/dace_soa_accel_buf.h"

int main(int argc, char **argv) {
    dace_soa_accel_bufHandle_t handle;
    int err;
    int64_t n = 42;
    double rsoft2 = 42;
    double * __restrict__ ax = (double*) calloc(n, sizeof(double));
    double * __restrict__ ay = (double*) calloc(n, sizeof(double));
    double * __restrict__ az = (double*) calloc(n, sizeof(double));
    double * __restrict__ fx = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fy = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fz = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ m = (double*) calloc(n, sizeof(double));
    double * __restrict__ rx = (double*) calloc(n, sizeof(double));
    double * __restrict__ ry = (double*) calloc(n, sizeof(double));
    double * __restrict__ rz = (double*) calloc(n, sizeof(double));


    handle = __dace_init_dace_soa_accel_buf(n);
    __program_dace_soa_accel_buf(handle, ax, ay, az, fx, fy, fz, m, rx, ry, rz, n, rsoft2);
    err = __dace_exit_dace_soa_accel_buf(handle);

    free(ax);
    free(ay);
    free(az);
    free(fx);
    free(fy);
    free(fz);
    free(m);
    free(rx);
    free(ry);
    free(rz);


    return err;
}
