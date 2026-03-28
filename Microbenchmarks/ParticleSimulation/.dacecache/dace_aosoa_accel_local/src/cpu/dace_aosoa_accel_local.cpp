/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_aosoa_accel_local_state_t {
    dace::perf::Report report;
};


#include <chrono>

#include <chrono>
#include <chrono>
#include <chrono>
inline void dace_aosoa_accel_local_584_4_0_0_2(dace_aosoa_accel_local_state_t *__state, const double&  __tmp_588_14_r, const double&  __tmp_589_14_r, const double&  __tmp_590_14_r, double*  __tmp_593_27_r, double*  __tmp_594_27_r, double*  __tmp_595_27_r, const double&  __tmp_596_46_r, double* __restrict__ __tmp_599_28_r, double&  __tmp_602_8_w, double&  __tmp_603_8_w, double&  __tmp_604_8_w, int64_t nb, int64_t vl) {
    double axi;
    double ayi;
    double azi;
    double rxi;
    double ryi;
    double rzi;
    int64_t bi_j;
    int64_t li_j;

    {

        {
            double __out;

            ///////////////////
            // Tasklet code (assign_585_8)
            __out = 0.0;
            ///////////////////

            axi = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_586_8)
            __out = 0.0;
            ///////////////////

            ayi = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_587_8)
            __out = 0.0;
            ///////////////////

            azi = __out;
        }
        {
            double __inp = __tmp_588_14_r;
            double __out;

            ///////////////////
            // Tasklet code (assign_588_8)
            __out = __inp;
            ///////////////////

            rxi = __out;
        }
        {
            double __inp = __tmp_589_14_r;
            double __out;

            ///////////////////
            // Tasklet code (assign_589_8)
            __out = __inp;
            ///////////////////

            ryi = __out;
        }
        {
            double __inp = __tmp_590_14_r;
            double __out;

            ///////////////////
            // Tasklet code (assign_590_8)
            __out = __inp;
            ///////////////////

            rzi = __out;
        }

    }
    for (bi_j = 0; (bi_j < nb); bi_j = (bi_j + 1)) {
        for (li_j = 0; (li_j < vl); li_j = (li_j + 1)) {
            {
                double __tmp_593_27_r_slice[1]  DACE_ALIGN(64);
                double dx[1]  DACE_ALIGN(64);
                double __tmp_594_27_r_slice[1]  DACE_ALIGN(64);
                double dy[1]  DACE_ALIGN(64);
                double __tmp_595_27_r_slice[1]  DACE_ALIGN(64);
                double dz[1]  DACE_ALIGN(64);
                double dx_times_dx[1]  DACE_ALIGN(64);
                double dy_times_dy[1]  DACE_ALIGN(64);
                double dx_dx_plus_dy_dy[1]  DACE_ALIGN(64);
                double dz_times_dz[1]  DACE_ALIGN(64);
                double dx_dx_dy_dy_plus_dz_dz[1]  DACE_ALIGN(64);
                double dsq[1]  DACE_ALIGN(64);
                double sqrt_dsq[1]  DACE_ALIGN(64);
                double inv_dist[1]  DACE_ALIGN(64);
                double inv_dist_times_inv_dist[1]  DACE_ALIGN(64);
                double inv_dist3[1]  DACE_ALIGN(64);
                double __tmp_599_28_r_slice[1]  DACE_ALIGN(64);
                double m_slice_times_dx[1]  DACE_ALIGN(64);
                double m_slice_dx_times_inv_dist3[1]  DACE_ALIGN(64);
                double axi_minus_m_slice_dx_inv_dist3[1]  DACE_ALIGN(64);
                double __tmp_599_28_r_slice_0[1]  DACE_ALIGN(64);
                double m_slice_times_dy[1]  DACE_ALIGN(64);
                double m_slice_dy_times_inv_dist3[1]  DACE_ALIGN(64);
                double ayi_minus_m_slice_dy_inv_dist3[1]  DACE_ALIGN(64);
                double __tmp_599_28_r_slice_1[1]  DACE_ALIGN(64);
                double m_slice_times_dz[1]  DACE_ALIGN(64);
                double m_slice_dz_times_inv_dist3[1]  DACE_ALIGN(64);
                double azi_minus_m_slice_dz_inv_dist3[1]  DACE_ALIGN(64);


                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_593_27_r + (((3 * bi_j) * vl) + li_j), __tmp_593_27_r_slice, 1);
                {
                    double __in1 = rxi;
                    double __in2 = __tmp_593_27_r_slice[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    dx[0] = __out;
                }
                {
                    double __in1 = dx[0];
                    double __in2 = dx[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    dx_times_dx[0] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_594_27_r + (((3 * bi_j) * vl) + li_j), __tmp_594_27_r_slice, 1);
                {
                    double __in1 = ryi;
                    double __in2 = __tmp_594_27_r_slice[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    dy[0] = __out;
                }
                {
                    double __in1 = dy[0];
                    double __in2 = dy[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    dy_times_dy[0] = __out;
                }
                {
                    double __in2 = dy_times_dy[0];
                    double __in1 = dx_times_dx[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    dx_dx_plus_dy_dy[0] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_595_27_r + (((3 * bi_j) * vl) + li_j), __tmp_595_27_r_slice, 1);
                {
                    double __in1 = rzi;
                    double __in2 = __tmp_595_27_r_slice[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    dz[0] = __out;
                }
                {
                    double __in1 = dz[0];
                    double __in2 = dz[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    dz_times_dz[0] = __out;
                }
                {
                    double __in2 = dz_times_dz[0];
                    double __in1 = dx_dx_plus_dy_dy[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    dx_dx_dy_dy_plus_dz_dz[0] = __out;
                }
                {
                    double __in2 = __tmp_596_46_r;
                    double __in1 = dx_dx_dy_dy_plus_dz_dz[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    dsq[0] = __out;
                }
                {
                    double __in1 = dsq[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_numpy_sqrt_)
                    __out = sqrt(__in1);
                    ///////////////////

                    sqrt_dsq[0] = __out;
                }
                {
                    double __in2 = sqrt_dsq[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    inv_dist[0] = __out;
                }
                {
                    double __in1 = inv_dist[0];
                    double __in2 = inv_dist[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    inv_dist_times_inv_dist[0] = __out;
                }
                {
                    double __in1 = inv_dist_times_inv_dist[0];
                    double __in2 = inv_dist[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    inv_dist3[0] = __out;
                }

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_599_28_r + ((bi_j * vl) + li_j), __tmp_599_28_r_slice, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_599_28_r + ((bi_j * vl) + li_j), __tmp_599_28_r_slice_0, 1);

                dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
                __tmp_599_28_r + ((bi_j * vl) + li_j), __tmp_599_28_r_slice_1, 1);
                {
                    double __in1 = __tmp_599_28_r_slice[0];
                    double __in2 = dx[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_times_dx[0] = __out;
                }
                {
                    double __in1 = m_slice_times_dx[0];
                    double __in2 = inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_dx_times_inv_dist3[0] = __out;
                }
                {
                    double __in1 = __tmp_599_28_r_slice_0[0];
                    double __in2 = dy[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_times_dy[0] = __out;
                }
                {
                    double __in1 = m_slice_times_dy[0];
                    double __in2 = inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_dy_times_inv_dist3[0] = __out;
                }
                {
                    double __in1 = __tmp_599_28_r_slice_1[0];
                    double __in2 = dz[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_times_dz[0] = __out;
                }
                {
                    double __in1 = m_slice_times_dz[0];
                    double __in2 = inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    m_slice_dz_times_inv_dist3[0] = __out;
                }
                {
                    double __in1 = axi;
                    double __in2 = m_slice_dx_times_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    axi_minus_m_slice_dx_inv_dist3[0] = __out;
                }
                {
                    double __inp = axi_minus_m_slice_dx_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_599_16)
                    __out = __inp;
                    ///////////////////

                    axi = __out;
                }
                {
                    double __in1 = ayi;
                    double __in2 = m_slice_dy_times_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ayi_minus_m_slice_dy_inv_dist3[0] = __out;
                }
                {
                    double __inp = ayi_minus_m_slice_dy_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_600_16)
                    __out = __inp;
                    ///////////////////

                    ayi = __out;
                }
                {
                    double __in1 = azi;
                    double __in2 = m_slice_dz_times_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    azi_minus_m_slice_dz_inv_dist3[0] = __out;
                }
                {
                    double __inp = azi_minus_m_slice_dz_inv_dist3[0];
                    double __out;

                    ///////////////////
                    // Tasklet code (assign_601_16)
                    __out = __inp;
                    ///////////////////

                    azi = __out;
                }

            }

        }

    }
    {

        {
            double __inp = axi;
            double __out;

            ///////////////////
            // Tasklet code (assign_602_8)
            __out = __inp;
            ///////////////////

            __tmp_602_8_w = __out;
        }
        {
            double __inp = ayi;
            double __out;

            ///////////////////
            // Tasklet code (assign_603_8)
            __out = __inp;
            ///////////////////

            __tmp_603_8_w = __out;
        }
        {
            double __inp = azi;
            double __out;

            ///////////////////
            // Tasklet code (assign_604_8)
            __out = __inp;
            ///////////////////

            __tmp_604_8_w = __out;
        }

    }
}

void __program_dace_aosoa_accel_local_internal(dace_aosoa_accel_local_state_t*__state, double * __restrict__ a, double * __restrict__ m, double * __restrict__ r, int64_t N_val, int64_t nb, double rsoft2, int64_t vl)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto bi_i = 0; bi_i < nb; bi_i += 1) {
                for (auto li_i = 0; li_i < vl; li_i += 1) {
                    dace_aosoa_accel_local_584_4_0_0_2(__state, r[(((3 * bi_i) * vl) + li_i)], r[((((3 * bi_i) * vl) + li_i) + vl)], r[((((3 * bi_i) * vl) + li_i) + (2 * vl))], &r[0], &r[vl], &r[(2 * vl)], rsoft2, &m[0], a[(((3 * bi_i) * vl) + li_i)], a[((((3 * bi_i) * vl) + li_i) + vl)], a[((((3 * bi_i) * vl) + li_i) + (2 * vl))], nb, vl);
                }
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_aosoa_accel_local", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_aosoa_accel_local/perf", __HASH_dace_aosoa_accel_local);
}

DACE_EXPORTED void __program_dace_aosoa_accel_local(dace_aosoa_accel_local_state_t *__state, double * __restrict__ a, double * __restrict__ m, double * __restrict__ r, int64_t N_val, int64_t nb, double rsoft2, int64_t vl)
{
    __program_dace_aosoa_accel_local_internal(__state, a, m, r, N_val, nb, rsoft2, vl);
}

DACE_EXPORTED dace_aosoa_accel_local_state_t *__dace_init_dace_aosoa_accel_local(int64_t nb, int64_t vl)
{

    int __result = 0;
    dace_aosoa_accel_local_state_t *__state = new dace_aosoa_accel_local_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_aosoa_accel_local(dace_aosoa_accel_local_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
