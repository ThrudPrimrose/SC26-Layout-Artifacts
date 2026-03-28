/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_aos_drift_state_t {
    dace::perf::Report report;
};


#include <chrono>
#include <chrono>
void __program_dace_aos_drift_internal(dace_aos_drift_state_t*__state, double * __restrict__ r, double * __restrict__ v, double dt, int64_t n)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                double v_slice_times_dt;
                double r_slice_plus_v_slice_dt;
                double v_slice_times_dt_0;
                double r_slice_plus_v_slice_dt_0;
                double v_slice_times_dt_1;
                double r_slice_plus_v_slice_dt_1;
                {
                    double __in1 = v[(3 * i)];
                    double __in2 = dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    v_slice_times_dt = __out;
                }
                {
                    double __in2 = v_slice_times_dt;
                    double __in1 = r[(3 * i)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    r_slice_plus_v_slice_dt = __out;
                }
                {
                    double __inp = r_slice_plus_v_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_410_8)
                    __out = __inp;
                    ///////////////////

                    r[(3 * i)] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = v[((3 * i) + 1)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    v_slice_times_dt_0 = __out;
                }
                {
                    double __in2 = v_slice_times_dt_0;
                    double __in1 = r[((3 * i) + 1)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    r_slice_plus_v_slice_dt_0 = __out;
                }
                {
                    double __inp = r_slice_plus_v_slice_dt_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_411_8)
                    __out = __inp;
                    ///////////////////

                    r[((3 * i) + 1)] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = v[((3 * i) + 2)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    v_slice_times_dt_1 = __out;
                }
                {
                    double __in2 = v_slice_times_dt_1;
                    double __in1 = r[((3 * i) + 2)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    r_slice_plus_v_slice_dt_1 = __out;
                }
                {
                    double __inp = r_slice_plus_v_slice_dt_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_412_8)
                    __out = __inp;
                    ///////////////////

                    r[((3 * i) + 2)] = __out;
                }
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_aos_drift", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_aos_drift/perf", __HASH_dace_aos_drift);
}

DACE_EXPORTED void __program_dace_aos_drift(dace_aos_drift_state_t *__state, double * __restrict__ r, double * __restrict__ v, double dt, int64_t n)
{
    __program_dace_aos_drift_internal(__state, r, v, dt, n);
}

DACE_EXPORTED dace_aos_drift_state_t *__dace_init_dace_aos_drift(int64_t n)
{

    int __result = 0;
    dace_aos_drift_state_t *__state = new dace_aos_drift_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_aos_drift(dace_aos_drift_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
