#include <cstdlib>
#include "../include/dace_aos_kick.h"

int main(int argc, char **argv) {
    dace_aos_kickHandle_t handle;
    int err;
    double dt = 42;
    int64_t n = 42;
    double * __restrict__ a = (double*) calloc((3 * n), sizeof(double));
    double * __restrict__ v = (double*) calloc((3 * n), sizeof(double));


    handle = __dace_init_dace_aos_kick(n);
    __program_dace_aos_kick(handle, a, v, dt, n);
    err = __dace_exit_dace_aos_kick(handle);

    free(a);
    free(v);


    return err;
}
