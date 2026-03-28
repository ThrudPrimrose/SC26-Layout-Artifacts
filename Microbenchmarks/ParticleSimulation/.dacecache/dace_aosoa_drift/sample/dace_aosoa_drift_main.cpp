#include <cstdlib>
#include "../include/dace_aosoa_drift.h"

int main(int argc, char **argv) {
    dace_aosoa_driftHandle_t handle;
    int err;
    double dt = 42;
    int64_t nb = 42;
    int64_t vl = 42;
    double * __restrict__ r = (double*) calloc(((3 * nb) * vl), sizeof(double));
    double * __restrict__ v = (double*) calloc(((3 * nb) * vl), sizeof(double));


    handle = __dace_init_dace_aosoa_drift(nb, vl);
    __program_dace_aosoa_drift(handle, r, v, dt, nb, vl);
    err = __dace_exit_dace_aosoa_drift(handle);

    free(r);
    free(v);


    return err;
}
