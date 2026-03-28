#include <cstdlib>
#include "../include/dace_aosoa_accel_local.h"

int main(int argc, char **argv) {
    dace_aosoa_accel_localHandle_t handle;
    int err;
    int64_t N_val = 42;
    int64_t nb = 42;
    double rsoft2 = 42;
    int64_t vl = 42;
    double * __restrict__ a = (double*) calloc(((3 * nb) * vl), sizeof(double));
    double * __restrict__ m = (double*) calloc((nb * vl), sizeof(double));
    double * __restrict__ r = (double*) calloc(((3 * nb) * vl), sizeof(double));


    handle = __dace_init_dace_aosoa_accel_local(nb, vl);
    __program_dace_aosoa_accel_local(handle, a, m, r, N_val, nb, rsoft2, vl);
    err = __dace_exit_dace_aosoa_accel_local(handle);

    free(a);
    free(m);
    free(r);


    return err;
}
