/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_aos_kick_state_t {
    dace::perf::Report report;
};


#include <chrono>
#include <chrono>
void __program_dace_aos_kick_internal(dace_aos_kick_state_t*__state, double * __restrict__ a, double * __restrict__ v, double dt, int64_t n)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                double a_slice_times_dt;
                double v_slice_plus_a_slice_dt;
                double a_slice_times_dt_0;
                double v_slice_plus_a_slice_dt_0;
                double a_slice_times_dt_1;
                double v_slice_plus_a_slice_dt_1;
                {
                    double __in1 = a[(3 * i)];
                    double __in2 = dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    a_slice_times_dt = __out;
                }
                {
                    double __in2 = a_slice_times_dt;
                    double __in1 = v[(3 * i)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    v_slice_plus_a_slice_dt = __out;
                }
                {
                    double __inp = v_slice_plus_a_slice_dt;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_399_8)
                    __out = __inp;
                    ///////////////////

                    v[(3 * i)] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = a[((3 * i) + 1)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    a_slice_times_dt_0 = __out;
                }
                {
                    double __in2 = a_slice_times_dt_0;
                    double __in1 = v[((3 * i) + 1)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    v_slice_plus_a_slice_dt_0 = __out;
                }
                {
                    double __inp = v_slice_plus_a_slice_dt_0;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_400_8)
                    __out = __inp;
                    ///////////////////

                    v[((3 * i) + 1)] = __out;
                }
                {
                    double __in2 = dt;
                    double __in1 = a[((3 * i) + 2)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    a_slice_times_dt_1 = __out;
                }
                {
                    double __in2 = a_slice_times_dt_1;
                    double __in1 = v[((3 * i) + 2)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    v_slice_plus_a_slice_dt_1 = __out;
                }
                {
                    double __inp = v_slice_plus_a_slice_dt_1;
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_401_8)
                    __out = __inp;
                    ///////////////////

                    v[((3 * i) + 2)] = __out;
                }
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_aos_kick", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_aos_kick/perf", __HASH_dace_aos_kick);
}

DACE_EXPORTED void __program_dace_aos_kick(dace_aos_kick_state_t *__state, double * __restrict__ a, double * __restrict__ v, double dt, int64_t n)
{
    __program_dace_aos_kick_internal(__state, a, v, dt, n);
}

DACE_EXPORTED dace_aos_kick_state_t *__dace_init_dace_aos_kick(int64_t n)
{

    int __result = 0;
    dace_aos_kick_state_t *__state = new dace_aos_kick_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_aos_kick(dace_aos_kick_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
