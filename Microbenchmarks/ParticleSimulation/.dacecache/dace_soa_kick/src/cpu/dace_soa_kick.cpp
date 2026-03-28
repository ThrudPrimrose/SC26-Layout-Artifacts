/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_soa_kick_state_t {
    dace::perf::Report report;
};


#include <chrono>
#include <chrono>
void __program_dace_soa_kick_internal(dace_soa_kick_state_t*__state, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                double ax_slice_times_dt;
                double vx_slice_plus_ax_slice_dt;
                double ay_slice_times_dt;
                double vy_slice_plus_ay_slice_dt;
                double az_slice_times_dt;
                double vz_slice_plus_az_slice_dt;
                {
                    double __in1 = ax[i];
                    double __in2 = dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ax_slice_times_dt = __out;
                }
                {
                    double __in2 = ax_slice_times_dt;
                    double __in1 = vx[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    vx_slice_plus_ax_slice_dt = __out;
                }
                {
                    double __inp = vx_slice_plus_ax_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_371_8)
                    __out = __inp;
                    ///////////////////

                    vx[i] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = ay[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    ay_slice_times_dt = __out;
                }
                {
                    double __in2 = ay_slice_times_dt;
                    double __in1 = vy[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    vy_slice_plus_ay_slice_dt = __out;
                }
                {
                    double __inp = vy_slice_plus_ay_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_372_8)
                    __out = __inp;
                    ///////////////////

                    vy[i] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = az[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    az_slice_times_dt = __out;
                }
                {
                    double __in2 = az_slice_times_dt;
                    double __in1 = vz[i];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    vz_slice_plus_az_slice_dt = __out;
                }
                {
                    double __inp = vz_slice_plus_az_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_373_8)
                    __out = __inp;
                    ///////////////////

                    vz[i] = __out;
                }
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_soa_kick", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_soa_kick/perf", __HASH_dace_soa_kick);
}

DACE_EXPORTED void __program_dace_soa_kick(dace_soa_kick_state_t *__state, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ vx, double * __restrict__ vy, double * __restrict__ vz, double dt, int64_t n)
{
    __program_dace_soa_kick_internal(__state, ax, ay, az, vx, vy, vz, dt, n);
}

DACE_EXPORTED dace_soa_kick_state_t *__dace_init_dace_soa_kick(int64_t n)
{

    int __result = 0;
    dace_soa_kick_state_t *__state = new dace_soa_kick_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_soa_kick(dace_soa_kick_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
