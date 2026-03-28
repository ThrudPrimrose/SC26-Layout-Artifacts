#include <cstdlib>
#include "../include/dace_aos_drift.h"

int main(int argc, char **argv) {
    dace_aos_driftHandle_t handle;
    int err;
    double dt = 42;
    int64_t n = 42;
    double * __restrict__ r = (double*) calloc((3 * n), sizeof(double));
    double * __restrict__ v = (double*) calloc((3 * n), sizeof(double));


    handle = __dace_init_dace_aos_drift(n);
    __program_dace_aos_drift(handle, r, v, dt, n);
    err = __dace_exit_dace_aos_drift(handle);

    free(r);
    free(v);


    return err;
}
