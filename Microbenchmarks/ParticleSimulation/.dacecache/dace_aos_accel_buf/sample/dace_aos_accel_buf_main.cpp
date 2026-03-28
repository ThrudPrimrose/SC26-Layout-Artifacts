#include <cstdlib>
#include "../include/dace_aos_accel_buf.h"

int main(int argc, char **argv) {
    dace_aos_accel_bufHandle_t handle;
    int err;
    int64_t n = 42;
    double rsoft2 = 42;
    double * __restrict__ a = (double*) calloc((3 * n), sizeof(double));
    double * __restrict__ fx = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fy = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fz = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ m = (double*) calloc(n, sizeof(double));
    double * __restrict__ r = (double*) calloc((3 * n), sizeof(double));


    handle = __dace_init_dace_aos_accel_buf(n);
    __program_dace_aos_accel_buf(handle, a, fx, fy, fz, m, r, n, rsoft2);
    err = __dace_exit_dace_aos_accel_buf(handle);

    free(a);
    free(fx);
    free(fy);
    free(fz);
    free(m);
    free(r);


    return err;
}
