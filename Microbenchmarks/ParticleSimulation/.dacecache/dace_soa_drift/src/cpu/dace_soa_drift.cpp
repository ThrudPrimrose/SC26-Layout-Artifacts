/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_soa_drift_state_t {
    dace::perf::Report report;
};


#include <chrono>
#include <chrono>
void __program_dace_soa_drift_internal(dace_soa_drift_state_t*__state, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                double vx_slice_times_dt;
                double rx_slice_plus_vx_slice_dt;
                double vy_slice_times_dt;
                double ry_slice_plus_vy_slice_dt;
                double vz_slice_times_dt;
                double rz_slice_plus_vz_slice_dt;
                {
                    double __in1 = vx[i];
                    double __in2 = dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    vx_slice_times_dt = __out;
                }
                {
                    double __in2 = vx_slice_times_dt;
                    double __in1 = rx[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    rx_slice_plus_vx_slice_dt = __out;
                }
                {
                    double __inp = rx_slice_plus_vx_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_386_8)
                    __out = __inp;
                    ///////////////////

                    rx[i] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = vy[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    vy_slice_times_dt = __out;
                }
                {
                    double __in2 = vy_slice_times_dt;
                    double __in1 = ry[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    ry_slice_plus_vy_slice_dt = __out;
                }
                {
                    double __inp = ry_slice_plus_vy_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_387_8)
                    __out = __inp;
                    ///////////////////

                    ry[i] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = vz[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    vz_slice_times_dt = __out;
                }
                {
                    double __in2 = vz_slice_times_dt;
                    double __in1 = rz[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    rz_slice_plus_vz_slice_dt = __out;
                }
                {
                    double __inp = rz_slice_plus_vz_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_388_8)
                    __out = __inp;
                    ///////////////////

                    rz[i] = __out;
                }
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_soa_drift", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_soa_drift/perf", __HASH_dace_soa_drift);
}

DACE_EXPORTED void __program_dace_soa_drift(dace_soa_drift_state_t *__state, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n)
{
    __program_dace_soa_drift_internal(__state, rx, ry, rz, vx, vy, vz, dt, n);
}

DACE_EXPORTED dace_soa_drift_state_t *__dace_init_dace_soa_drift(int64_t n)
{

    int __result = 0;
    dace_soa_drift_state_t *__state = new dace_soa_drift_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_soa_drift(dace_soa_drift_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
