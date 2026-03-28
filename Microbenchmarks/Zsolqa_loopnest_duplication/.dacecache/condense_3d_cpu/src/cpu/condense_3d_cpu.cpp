/* DaCe AUTO-GENERATED FILE. DO NOT MODIFY */
#include <dace/dace.h>
#include "../../include/hash.h"

struct condense_3d_cpu_state_t {
    dace::perf::Report report;
};


#include <chrono>

#include <chrono>
#include <chrono>
#include <chrono>
inline void loop_body_1_0_0(condense_3d_cpu_state_t *__state, const double&  r4ies, const double&  r4les, const double&  r5alscp, const double&  r5alvcp, const double&  retv, const double&  rlmin, const double&  rthomo, const double&  rtice, const double&  rtwat, const double&  rtwat_rtice_r, double* __restrict__ za, double* __restrict__ zdqs, double* __restrict__ zqsmix, double* __restrict__ zqv, double* __restrict__ ztp1, double* __restrict__ zqxfg, double* __restrict__ zsolqa, int64_t jk, int64_t jl, int klon) {
    double cdm;
    double cdm_part;
    double za_slice_times_expr;
    double za_index;
    double zdqs_index;
    double neg_rlmin;
    double neg_zdqs_slice;
    double lc;
    double za_index_2;
    double ztp1_index_2;
    bool __tmp8;


    za_index = za[((jk * klon) + jl)];
    if ((za_index > 1e-14)) {

        zdqs_index = zdqs[jl];
        neg_rlmin = (- rlmin);
        if ((zdqs_index <= neg_rlmin)) {

            neg_zdqs_slice = (- zdqs[jl]);

            lc = max(neg_zdqs_slice, 0.0);
            za_index_2 = za[((jk * klon) + jl)];
            {
                double __tmp5;
                double af;
                double __tmp2;
                double expr_rtice_times_rtwat_rtice_r;
                double expr_minus_rtice;
                double ztp1_slice_r4ies_pow_2;
                double za_slice_times_zqsmix_slice;
                double __tmp6;
                double retv_times_zqsmix_slice;
                double af_r5alvcp_ztp1_slice_r4les_2_plus_1_0_af_r5alscp_ztp1_slice_r4ies_2;
                double ztp1_slice_minus_r4les;
                double zcor;
                double zqv_slice_minus_za_slice_zqsmix_slice;
                double af_times_r5alvcp;
                double zqv_slice_minus_zqsmix_slice;
                double af_r5alvcp_div_ztp1_slice_r4les_2;
                double ztp1_slice_minus_r4ies;
                double zcor_times_zqsmix_slice;
                double ztp1_slice_r4les_pow_2;
                double zcor_zqsmix_slice_times_af_r5alvcp_ztp1_slice_r4les_2_1_0_af_r5alscp_ztp1_slice_r4ies_2;
                double max_rtice_expr;
                double min_rtwat_ztp1_slice;
                double expr_rtice_rtwat_rtice_r_pow_2;
                double __tmp4;
                double __tmp3;

                {
                    double __in_a = rtwat;
                    double __in_b = ztp1[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(__in_a, __in_b);
                    ///////////////////

                    min_rtwat_ztp1_slice = __out;
                }
                {
                    double __in_a = rtice;
                    double __in_b = min_rtwat_ztp1_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, __in_b);
                    ///////////////////

                    max_rtice_expr = __out;
                }
                {
                    double __in1 = max_rtice_expr;
                    double __in2 = rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    expr_minus_rtice = __out;
                }
                {
                    double __in2 = rtwat_rtice_r;
                    double __in1 = expr_minus_rtice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    expr_rtice_times_rtwat_rtice_r = __out;
                }
                {
                    double __in1 = expr_rtice_times_rtwat_rtice_r;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    expr_rtice_rtwat_rtice_r_pow_2 = __out;
                }
                {
                    double __in_b = expr_rtice_rtwat_rtice_r_pow_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(1, __in_b);
                    ///////////////////

                    af = __out;
                }
                {
                    double __in2 = af;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp3 = __out;
                }
                {
                    double __in1 = retv;
                    double __in2 = zqsmix[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    retv_times_zqsmix_slice = __out;
                }
                {
                    double __in2 = retv_times_zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (1.0 - __in2);
                    ///////////////////

                    __tmp2 = __out;
                }
                {
                    double __in2 = __tmp2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (1.0 / __in2);
                    ///////////////////

                    zcor = __out;
                }
                {
                    double __in1 = zcor;
                    double __in2 = zqsmix[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zcor_times_zqsmix_slice = __out;
                }
                {
                    double __in2 = zqsmix[((jk * klon) + jl)];
                    double __in1 = zqv[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqv_slice_minus_zqsmix_slice = __out;
                }
                {
                    double __in2 = r5alvcp;
                    double __in1 = af;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    af_times_r5alvcp = __out;
                }
                {
                    double __in2 = r4les;
                    double __in1 = ztp1[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_r4les = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_r4les;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_r4les_pow_2 = __out;
                }
                {
                    double __in2 = ztp1_slice_r4les_pow_2;
                    double __in1 = af_times_r5alvcp;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    af_r5alvcp_div_ztp1_slice_r4les_2 = __out;
                }
                {
                    double __in2 = r5alscp;
                    double __in1 = __tmp3;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    __tmp4 = __out;
                }
                {
                    double __in2 = r4ies;
                    double __in1 = ztp1[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    ztp1_slice_minus_r4ies = __out;
                }
                {
                    double __in1 = ztp1_slice_minus_r4ies;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Pow_)
                    __out = (dace::math::ipow(__in1, 2));
                    ///////////////////

                    ztp1_slice_r4ies_pow_2 = __out;
                }
                {
                    double __in2 = ztp1_slice_r4ies_pow_2;
                    double __in1 = __tmp4;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    __tmp5 = __out;
                }
                {
                    double __in2 = __tmp5;
                    double __in1 = af_r5alvcp_div_ztp1_slice_r4les_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (__in1 + __in2);
                    ///////////////////

                    af_r5alvcp_ztp1_slice_r4les_2_plus_1_0_af_r5alscp_ztp1_slice_r4ies_2 = __out;
                }
                {
                    double __in2 = af_r5alvcp_ztp1_slice_r4les_2_plus_1_0_af_r5alscp_ztp1_slice_r4ies_2;
                    double __in1 = zcor_times_zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    zcor_zqsmix_slice_times_af_r5alvcp_ztp1_slice_r4les_2_1_0_af_r5alscp_ztp1_slice_r4ies_2 = __out;
                }
                {
                    double __in2 = zcor_zqsmix_slice_times_af_r5alvcp_ztp1_slice_r4les_2_1_0_af_r5alscp_ztp1_slice_r4ies_2;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Add_)
                    __out = (1.0 + __in2);
                    ///////////////////

                    __tmp6 = __out;
                }
                {
                    double __in2 = __tmp6;
                    double __in1 = zqv_slice_minus_zqsmix_slice;
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    cdm = __out;
                }
                {
                    double __in2 = zqsmix[((jk * klon) + jl)];
                    double __in1 = za[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    za_slice_times_zqsmix_slice = __out;
                }
                {
                    double __in2 = za_slice_times_zqsmix_slice;
                    double __in1 = zqv[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Sub_)
                    __out = (__in1 - __in2);
                    ///////////////////

                    zqv_slice_minus_za_slice_zqsmix_slice = __out;
                }
                {
                    double __in1 = zqv_slice_minus_za_slice_zqsmix_slice;
                    double __in2 = za[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Div_)
                    __out = (__in1 / __in2);
                    ///////////////////

                    cdm_part = __out;
                }

            }
            if ((! (za_index_2 > 0.99))) {
                {

                    {
                        double __inp = cdm_part;
                        double __out;

                        ///////////////////
                        // Tasklet code (assign_85_24)
                        __out = __inp;
                        ///////////////////

                        cdm = __out;
                    }

                }
            }
            {
                double max_expr_0_0;
                double min_lc_cdm;

                {
                    double __in_b = cdm;
                    double __out;

                    ///////////////////
                    // Tasklet code (__min2)
                    __out = min(max(neg_zdqs_slice, 0.0), __in_b);
                    ///////////////////

                    min_lc_cdm = __out;
                }
                {
                    double __in_a = min_lc_cdm;
                    double __out;

                    ///////////////////
                    // Tasklet code (__max2)
                    __out = max(__in_a, 0.0);
                    ///////////////////

                    max_expr_0_0 = __out;
                }
                {
                    double __in2 = max_expr_0_0;
                    double __in1 = za[((jk * klon) + jl)];
                    double __out;

                    ///////////////////
                    // Tasklet code (_Mult_)
                    __out = (__in1 * __in2);
                    ///////////////////

                    za_slice_times_expr = __out;
                }

            }
            lc = za_slice_times_expr;
            if ((lc >= rlmin)) {

                ztp1_index_2 = ztp1[((jk * klon) + jl)];

                __tmp8 = (ztp1_index_2 > rthomo);
                if (__tmp8) {
                    {
                        double zsolqa_slice_plus_lc;
                        double zsolqa_slice_minus_lc;
                        double zqxfg_slice_plus_lc;

                        {
                            double __in1 = zsolqa[(jl + (20 * klon))];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + lc);
                            ///////////////////

                            zsolqa_slice_plus_lc = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_lc;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_89_28)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl + (20 * klon))] = __out;
                        }
                        {
                            double __in1 = zsolqa[(jl + (4 * klon))];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - lc);
                            ///////////////////

                            zsolqa_slice_minus_lc = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_lc;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_90_28)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl + (4 * klon))] = __out;
                        }
                        {
                            double __in1 = zqxfg[jl];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + lc);
                            ///////////////////

                            zqxfg_slice_plus_lc = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_lc;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_91_28)
                            __out = __inp;
                            ///////////////////

                            zqxfg[jl] = __out;
                        }

                    }
                } else {
                    {
                        double zsolqa_slice_minus_lc_0;
                        double zsolqa_slice_plus_lc_0;
                        double zqxfg_slice_plus_lc_0;

                        {
                            double __in1 = zsolqa[(jl + (21 * klon))];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + lc);
                            ///////////////////

                            zsolqa_slice_plus_lc_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_plus_lc_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_93_28)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl + (21 * klon))] = __out;
                        }
                        {
                            double __in1 = zsolqa[(jl + (9 * klon))];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Sub_)
                            __out = (__in1 - lc);
                            ///////////////////

                            zsolqa_slice_minus_lc_0 = __out;
                        }
                        {
                            double __inp = zsolqa_slice_minus_lc_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_94_28)
                            __out = __inp;
                            ///////////////////

                            zsolqa[(jl + (9 * klon))] = __out;
                        }
                        {
                            double __in1 = zqxfg[(jl + klon)];
                            double __out;

                            ///////////////////
                            // Tasklet code (_Add_)
                            __out = (__in1 + lc);
                            ///////////////////

                            zqxfg_slice_plus_lc_0 = __out;
                        }
                        {
                            double __inp = zqxfg_slice_plus_lc_0;
                            double __out;

                            ///////////////////
                            // Tasklet code (assign_95_28)
                            __out = __inp;
                            ///////////////////

                            zqxfg[(jl + klon)] = __out;
                        }

                    }
                }
            }
        }
    }
}

void __program_condense_3d_cpu_internal(condense_3d_cpu_state_t*__state, double * __restrict__ za, double * __restrict__ zdqs, double * __restrict__ zqsmix, double * __restrict__ zqv, double * __restrict__ zqxfg, double * __restrict__ zsolqa, double * __restrict__ ztp1, int klev, int klon, double r4ies, double r4les, double r5alscp, double r5alvcp, double retv, double rlmin, double rthomo, double rtice, double rtwat, double rtwat_rtice_r)
{
    __state->report.reset();
    int64_t jk;

    for (jk = 0; (jk < klev); jk = (jk + 1)) {
        {

            auto __dace_tbegin_1_0 = std::chrono::high_resolution_clock::now();
            {
                #pragma omp parallel for
                for (auto jl = 0; jl < klon; jl += 1) {
                    loop_body_1_0_0(__state, r4ies, r4les, r5alscp, r5alvcp, retv, rlmin, rthomo, rtice, rtwat, rtwat_rtice_r, &za[0], &zdqs[0], &zqsmix[0], &zqv[0], &ztp1[0], &zqxfg[0], &zsolqa[0], jk, jl, klon);
                }
            }
            auto __dace_tend_1_0 = std::chrono::high_resolution_clock::now();
            unsigned long int __dace_ts_start_1_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tbegin_1_0.time_since_epoch()).count();
            unsigned long int __dace_ts_end_1_0 = std::chrono::duration_cast<std::chrono::microseconds>(__dace_tend_1_0.time_since_epoch()).count();
            __state->report.add_completion("State single_state_body", "Timer", __dace_ts_start_1_0, __dace_ts_end_1_0, 1, 0, -1);

        }

    }
    __state->report.save(".dacecache/condense_3d_cpu/perf", __HASH_condense_3d_cpu);
}

DACE_EXPORTED void __program_condense_3d_cpu(condense_3d_cpu_state_t *__state, double * __restrict__ za, double * __restrict__ zdqs, double * __restrict__ zqsmix, double * __restrict__ zqv, double * __restrict__ zqxfg, double * __restrict__ zsolqa, double * __restrict__ ztp1, int klev, int klon, double r4ies, double r4les, double r5alscp, double r5alvcp, double retv, double rlmin, double rthomo, double rtice, double rtwat, double rtwat_rtice_r)
{
    __program_condense_3d_cpu_internal(__state, za, zdqs, zqsmix, zqv, zqxfg, zsolqa, ztp1, klev, klon, r4ies, r4les, r5alscp, r5alvcp, retv, rlmin, rthomo, rtice, rtwat, rtwat_rtice_r);
}

DACE_EXPORTED condense_3d_cpu_state_t *__dace_init_condense_3d_cpu(int klev, int klon)
{

    int __result = 0;
    condense_3d_cpu_state_t *__state = new condense_3d_cpu_state_t;

    if (__result) {
        delete __state;
        return nullptr;
    }

    return __state;
}

DACE_EXPORTED int __dace_exit_condense_3d_cpu(condense_3d_cpu_state_t *__state)
{

    int __err = 0;
    delete __state;
    return __err;
}
