/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct dace_soa_accel_buf_state_t {
    dace::perf::Report report;
};


#include <chrono>

#include <chrono>
#include <chrono>
#include <chrono>
inline void dace_soa_accel_buf_498_4_0_0_12(dace_soa_accel_buf_state_t *__state, double* __restrict__ __tmp_503_24_r, double* __restrict__ __tmp_504_24_r, double* __restrict__ __tmp_505_24_r, double&  __tmp_506_8_w, double&  __tmp_507_8_w, double&  __tmp_508_8_w, int64_t n) {
    double axi;
    double ayi;
    double azi;
    int64_t j;

    {

        {
            double __out;

            ///////////////////
            // Tasklet code (assign_499_8)
            __out = 0.0;
            ///////////////////

            axi = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_500_8)
            __out = 0.0;
            ///////////////////

            ayi = __out;
        }
        {
            double __out;

            ///////////////////
            // Tasklet code (assign_501_8)
            __out = 0.0;
            ///////////////////

            azi = __out;
        }

    }
    for (j = 0; (j < n); j = (j + 1)) {
        {
            double __tmp_503_24_r_slice[1]  DACE_ALIGN(64);
            double axi_plus_fx_slice[1]  DACE_ALIGN(64);
            double __tmp_504_24_r_slice[1]  DACE_ALIGN(64);
            double ayi_plus_fy_slice[1]  DACE_ALIGN(64);
            double __tmp_505_24_r_slice[1]  DACE_ALIGN(64);
            double azi_plus_fz_slice[1]  DACE_ALIGN(64);


            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
            __tmp_503_24_r + j, __tmp_503_24_r_slice, 1);
            {
                double __in1 = axi;
                double __in2 = __tmp_503_24_r_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (_Add_)
                __out = (__in1 + __in2);
                ///////////////////

                axi_plus_fx_slice[0] = __out;
            }
            {
                double __inp = axi_plus_fx_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (assign_503_12)
                __out = __inp;
                ///////////////////

                axi = __out;
            }

            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
            __tmp_504_24_r + j, __tmp_504_24_r_slice, 1);
            {
                double __in1 = ayi;
                double __in2 = __tmp_504_24_r_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (_Add_)
                __out = (__in1 + __in2);
                ///////////////////

                ayi_plus_fy_slice[0] = __out;
            }
            {
                double __inp = ayi_plus_fy_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (assign_504_12)
                __out = __inp;
                ///////////////////

                ayi = __out;
            }

            dace::CopyND<double, 1, false, 1>::template ConstDst<1>::Copy(
            __tmp_505_24_r + j, __tmp_505_24_r_slice, 1);
            {
                double __in1 = azi;
                double __in2 = __tmp_505_24_r_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (_Add_)
                __out = (__in1 + __in2);
                ///////////////////

                azi_plus_fz_slice[0] = __out;
            }
            {
                double __inp = azi_plus_fz_slice[0];
                double __out;

                ///////////////////
                // Tasklet code (assign_505_12)
                __out = __inp;
                ///////////////////

                azi = __out;
            }

        }

    }
    {

        {
            double __inp = axi;
            double __out;

            ///////////////////
            // Tasklet code (assign_506_8)
            __out = __inp;
            ///////////////////

            __tmp_506_8_w = __out;
        }
        {
            double __inp = ayi;
            double __out;

            ///////////////////
            // Tasklet code (assign_507_8)
            __out = __inp;
            ///////////////////

            __tmp_507_8_w = __out;
        }
        {
            double __inp = azi;
            double __out;

            ///////////////////
            // Tasklet code (assign_508_8)
            __out = __inp;
            ///////////////////

            __tmp_508_8_w = __out;
        }

    }
}

void __program_dace_soa_accel_buf_internal(dace_soa_accel_buf_state_t*__state, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ fx, double * __restrict__ fy, double * __restrict__ fz, double * __restrict__ m, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, int64_t n, double rsoft2)
{
    __state->report.reset();
    auto __dace_tbegin_0 = std::chrono::high_resolution_clock::now();

    {

        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                for (auto j = 0; j < n; j += 1) {
                    double dx;
                    double dy;
                    double dz;
                    double dx_times_dx;
                    double dy_times_dy;
                    double dx_dx_plus_dy_dy;
                    double dz_times_dz;
                    double dx_dx_dy_dy_plus_dz_dz;
                    double dsq;
                    double sqrt_dsq;
                    double inv_dist;
                    double inv_dist_times_inv_dist;
                    double inv_dist3;
                    double neg_m_slice;
                    double expr_times_dx;
                    double fx_slice;
                    double neg_m_slice_0;
                    double expr_times_dy;
                    double fy_slice;
                    double neg_m_slice_1;
                    double expr_times_dz;
                    double fz_slice;
                    {
                        double __in1 = rx[i];
                        double __in2 = rx[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        dx = __out;
                    }
                    {
                        double __in1 = dx;
                        double __in2 = dx;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        dx_times_dx = __out;
                    }
                    {
                        double __in1 = ry[i];
                        double __in2 = ry[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        dy = __out;
                    }
                    {
                        double __in1 = dy;
                        double __in2 = dy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        dy_times_dy = __out;
                    }
                    {
                        double __in2 = dy_times_dy;
                        double __in1 = dx_times_dx;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        dx_dx_plus_dy_dy = __out;
                    }
                    {
                        double __in1 = rz[i];
                        double __in2 = rz[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_Sub_)
                        __out = (__in1 - __in2);
                        ///////////////////

                        dz = __out;
                    }
                    {
                        double __in1 = dz;
                        double __in2 = dz;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        dz_times_dz = __out;
                    }
                    {
                        double __in2 = dz_times_dz;
                        double __in1 = dx_dx_plus_dy_dy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        dx_dx_dy_dy_plus_dz_dz = __out;
                    }
                    {
                        double __in1 = dx_dx_dy_dy_plus_dz_dz;
                        double __in2 = rsoft2;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Add_)
                        __out = (__in1 + __in2);
                        ///////////////////

                        dsq = __out;
                    }
                    {
                        double __in1 = dsq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_numpy_sqrt_)
                        __out = sqrt(__in1);
                        ///////////////////

                        sqrt_dsq = __out;
                    }
                    {
                        double __in2 = sqrt_dsq;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Div_)
                        __out = (1.0 / __in2);
                        ///////////////////

                        inv_dist = __out;
                    }
                    {
                        double __in1 = inv_dist;
                        double __in2 = inv_dist;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        inv_dist_times_inv_dist = __out;
                    }
                    {
                        double __in1 = inv_dist_times_inv_dist;
                        double __in2 = inv_dist;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        inv_dist3 = __out;
                    }
                    {
                        double __in = m[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        neg_m_slice = __out;
                    }
                    {
                        double __in1 = neg_m_slice;
                        double __in2 = dx;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        expr_times_dx = __out;
                    }
                    {
                        double __in1 = expr_times_dx;
                        double __in2 = inv_dist3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        fx_slice = __out;
                    }
                    {
                        double __inp = fx_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_495_8)
                        __out = __inp;
                        ///////////////////

                        fx[((i * n) + j)] = __out;
                    }
                    {
                        double __in = m[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        neg_m_slice_0 = __out;
                    }
                    {
                        double __in1 = neg_m_slice_0;
                        double __in2 = dy;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        expr_times_dy = __out;
                    }
                    {
                        double __in1 = expr_times_dy;
                        double __in2 = inv_dist3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        fy_slice = __out;
                    }
                    {
                        double __inp = fy_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_496_8)
                        __out = __inp;
                        ///////////////////

                        fy[((i * n) + j)] = __out;
                    }
                    {
                        double __in = m[j];
                        double __out;

                        ///////////////////
                        // Tasklet code (_USub_)
                        __out = (- __in);
                        ///////////////////

                        neg_m_slice_1 = __out;
                    }
                    {
                        double __in1 = neg_m_slice_1;
                        double __in2 = dz;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        expr_times_dz = __out;
                    }
                    {
                        double __in1 = expr_times_dz;
                        double __in2 = inv_dist3;
                        double __out;

                        ///////////////////
                        // Tasklet code (_Mult_)
                        __out = (__in1 * __in2);
                        ///////////////////

                        fz_slice = __out;
                    }
                    {
                        double __inp = fz_slice;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_497_8)
                        __out = __inp;
                        ///////////////////

                        fz[((i * n) + j)] = __out;
                    }
                }
            }
        }
        {
            #pragma omp parallel for
            for (auto i = 0; i < n; i += 1) {
                dace_soa_accel_buf_498_4_0_0_12(__state, &fx[(i * n)], &fy[(i * n)], &fz[(i * n)], ax[i], ay[i], az[i], n);
            }
        }

    }
    auto __dace_tend_0 = std::chrono::high_resolution_clock::now();
    unsigned long int __dace_ts_start_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_0.time_since_epoch()).count();
    unsigned long int __dace_ts_end_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_0.time_since_epoch()).count();
    __state->report.add_completion("SDFG dace_soa_accel_buf", "Timer", __dace_ts_start_0, __dace_ts_end_0, 0, -1, -1);
    __state->report.save(".dacecache/dace_soa_accel_buf/perf", __HASH_dace_soa_accel_buf);
}

DACE_EXPORTED void __program_dace_soa_accel_buf(dace_soa_accel_buf_state_t *__state, double * __restrict__ ax, double * __restrict__ ay, double * __restrict__ az, double * __restrict__ fx, double * __restrict__ fy, double * __restrict__ fz, double * __restrict__ m, double * __restrict__ rx, double * __restrict__ ry, double * __restrict__ rz, int64_t n, double rsoft2)
{
    __program_dace_soa_accel_buf_internal(__state, ax, ay, az, fx, fy, fz, m, rx, ry, rz, n, rsoft2);
}

DACE_EXPORTED dace_soa_accel_buf_state_t *__dace_init_dace_soa_accel_buf(int64_t n)
{

    int __result = 0;
    dace_soa_accel_buf_state_t *__state = new dace_soa_accel_buf_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_dace_soa_accel_buf(dace_soa_accel_buf_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
