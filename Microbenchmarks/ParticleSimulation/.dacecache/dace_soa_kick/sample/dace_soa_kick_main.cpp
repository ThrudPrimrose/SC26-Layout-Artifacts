#include <cstdlib>
#include "../include/dace_soa_kick.h"

int main(int argc, char **argv) {
    dace_soa_kickHandle_t handle;
    int err;
    double dt = 42;
    int64_t n = 42;
    double * __restrict__ ax = (double*) calloc(n, sizeof(double));
    double * __restrict__ ay = (double*) calloc(n, sizeof(double));
    double * __restrict__ az = (double*) calloc(n, sizeof(double));
    double * __restrict__ vx = (double*) calloc(n, sizeof(double));
    double * __restrict__ vy = (double*) calloc(n, sizeof(double));
    double * __restrict__ vz = (double*) calloc(n, sizeof(double));


    handle = __dace_init_dace_soa_kick(n);
    __program_dace_soa_kick(handle, ax, ay, az, vx, vy, vz, dt, n);
    err = __dace_exit_dace_soa_kick(handle);

    free(ax);
    free(ay);
    free(az);
    free(vx);
    free(vy);
    free(vz);


    return err;
}
