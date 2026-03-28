#include <cstdlib>
#include "../include/dace_aosoa_kick.h"

int main(int argc, char **argv) {
    dace_aosoa_kickHandle_t handle;
    int err;
    double dt = 42;
    int64_t nb = 42;
    int64_t vl = 42;
    double * __restrict__ a = (double*) calloc(((3 * nb) * vl), sizeof(double));
    double * __restrict__ v = (double*) calloc(((3 * nb) * vl), sizeof(double));


    handle = __dace_init_dace_aosoa_kick(nb, vl);
    __program_dace_aosoa_kick(handle, a, v, dt, nb, vl);
    err = __dace_exit_dace_aosoa_kick(handle);

    free(a);
    free(v);


    return err;
}
