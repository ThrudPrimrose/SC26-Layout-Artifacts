#include <cstdlib>
#include "../include/dace_aosoa_accel_buf.h"

int main(int argc, char **argv) {
    dace_aosoa_accel_bufHandle_t handle;
    int err;
    int64_t N_val = 42;
    int64_t n = 42;
    int64_t nb = 42;
    double rsoft2 = 42;
    int64_t vl = 42;
    double * __restrict__ a = (double*) calloc(((3 * nb) * vl), sizeof(double));
    double * __restrict__ fx = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fy = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ fz = (double*) calloc((n * n), sizeof(double));
    double * __restrict__ m = (double*) calloc((nb * vl), sizeof(double));
    double * __restrict__ r = (double*) calloc(((3 * nb) * vl), sizeof(double));


    handle = __dace_init_dace_aosoa_accel_buf(n, nb, vl);
    __program_dace_aosoa_accel_buf(handle, a, fx, fy, fz, m, r, N_val, n, nb, rsoft2, vl);
    err = __dace_exit_dace_aosoa_accel_buf(handle);

    free(a);
    free(fx);
    free(fy);
    free(fz);
    free(m);
    free(r);


    return err;
}
