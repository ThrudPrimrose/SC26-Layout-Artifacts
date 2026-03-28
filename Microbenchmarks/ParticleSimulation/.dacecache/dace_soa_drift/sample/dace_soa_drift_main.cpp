#include <cstdlib>
#include "../include/dace_soa_drift.h"

int main(int argc, char **argv) {
    dace_soa_driftHandle_t handle;
    int err;
    double dt = 42;
    int64_t n = 42;
    double * __restrict__ rx = (double*) calloc(n, sizeof(double));
    double * __restrict__ ry = (double*) calloc(n, sizeof(double));
    double * __restrict__ rz = (double*) calloc(n, sizeof(double));
    double * __restrict__ vx = (double*) calloc(n, sizeof(double));
    double * __restrict__ vy = (double*) calloc(n, sizeof(double));
    double * __restrict__ vz = (double*) calloc(n, sizeof(double));


    handle = __dace_init_dace_soa_drift(n);
    __program_dace_soa_drift(handle, rx, ry, rz, vx, vy, vz, dt, n);
    err = __dace_exit_dace_soa_drift(handle);

    free(rx);
    free(ry);
    free(rz);
    free(vx);
    free(vy);
    free(vz);


    return err;
}
