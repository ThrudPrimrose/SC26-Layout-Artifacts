#include <cstdlib>
#include "../include/dace_aos_accel_local.h"

int main(int argc, char **argv) {
    dace_aos_accel_localHandle_t handle;
    int err;
    int64_t n = 42;
    double rsoft2 = 42;
    double * __restrict__ a = (double*) calloc((3 * n), sizeof(double));
    double * __restrict__ m = (double*) calloc(n, sizeof(double));
    double * __restrict__ r = (double*) calloc((3 * n), sizeof(double));


    handle = __dace_init_dace_aos_accel_local(n);
    __program_dace_aos_accel_local(handle, a, m, r, n, rsoft2);
    err = __dace_exit_dace_aos_accel_local(handle);

    free(a);
    free(m);
    free(r);


    return err;
}
