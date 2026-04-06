#include "types.h"
#include "data_loading.h"
#include "usxx_kernels_aos_gpu.h"
#include "usxx_kernels_soa_gpu.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

// Array dimensions
static constexpr int RHOC_SIZE   = 120000;
static constexpr int BECPC_SIZE  = 122;
static constexpr int QGM_ROWS   = 55191;
static constexpr int QGM_COLS   = 397;
static constexpr int EIGTS1_DIM1 = 217;
static constexpr int EIGTS1_DIM2 = 10;
static constexpr int EIGTS2_DIM1 = 109;
static constexpr int EIGTS2_DIM2 = 10;
static constexpr int NAT_MAX    = 10;
static constexpr int NTYP_MAX   = 3;
static constexpr int MILL_DIM1  = 3;
static constexpr int MILL_DIM2  = 156121;
static constexpr int NL_SIZE    = 55191;
static constexpr int IJTOH_N1   = 19;
static constexpr int IJTOH_N2   = 19;
static constexpr int IJTOH_N3   = 3;

static constexpr int NUM_ITERS  = 100;
static constexpr int NUM_WARMUP = 5;
static constexpr int TOTAL_ITERS = NUM_ITERS + NUM_WARMUP;

static const int TBLOCK_SIZES[] = { 32, 64, 128, 256, 512, 1024 };
static const int COARSEN_FACTORS[] = { 1, 2, 4, 8, 16 };
static constexpr int N_TBLOCK = sizeof(TBLOCK_SIZES) / sizeof(TBLOCK_SIZES[0]);
static constexpr int N_COARSEN = sizeof(COARSEN_FACTORS) / sizeof(COARSEN_FACTORS[0]);

// ---- Host arrays ----
static int nkb, ngms, nat, ntyp;
static Complex_DP rhoc[RHOC_SIZE], becphi_c[BECPC_SIZE], becpsi_c[BECPC_SIZE];
static Complex_DP qgm[QGM_COLS * QGM_ROWS];
static Complex_DP eigts1[EIGTS1_DIM2 * EIGTS1_DIM1];
static Complex_DP eigts2[EIGTS2_DIM2 * EIGTS2_DIM1];
static Complex_DP eigts3[EIGTS2_DIM2 * EIGTS2_DIM1];
static Complex_DP rhoc_out[RHOC_SIZE], rhoc_out_sim[RHOC_SIZE];
static DP xkq[3], xk[3], tau[NAT_MAX * 3];
static int nij_type[NTYP_MAX], ityp[NAT_MAX], ofsbeta[NAT_MAX], nh[NTYP_MAX];
static int ijtoh[IJTOH_N3 * IJTOH_N2 * IJTOH_N1];
static int mill[MILL_DIM2 * MILL_DIM1];
static int dfftt__nl[NL_SIZE], upf_tvanp[NTYP_MAX];

static Complex_DP qgm_T[QGM_ROWS * QGM_COLS];
static Complex_DP eigts1_T[EIGTS1_DIM1 * EIGTS1_DIM2];
static Complex_DP eigts2_T[EIGTS2_DIM1 * EIGTS2_DIM2];
static Complex_DP eigts3_T[EIGTS2_DIM1 * EIGTS2_DIM2];
static int dfftt__nl_sorted[NL_SIZE], dfftt__nl_ix[NL_SIZE];

static int file_id = 12;

static void transpose_2d(const Complex_DP* A, Complex_DP* B, int m, int n) {
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            B[i * n + j] = A[j * m + i];
}

static void init_addusxx_data() {
    char filename[256];
    #define LOAD_INT(var, name) \
        snprintf(filename, sizeof(filename), "./bin/addusxx_g__" name "_%d.bin", file_id); \
        data_load_int(var, filename);
    #define LOAD_CMPLX(var, cnt, name) \
        snprintf(filename, sizeof(filename), "./bin/addusxx_g__" name "_%d.bin", file_id); \
        data_load_cmplx_array(var, cnt, filename);
    #define LOAD_REAL(var, cnt, name) \
        snprintf(filename, sizeof(filename), "./bin/addusxx_g__" name "_%d.bin", file_id); \
        data_load_real_array(var, cnt, filename);
    #define LOAD_INTA(var, cnt, name) \
        snprintf(filename, sizeof(filename), "./bin/addusxx_g__" name "_%d.bin", file_id); \
        data_load_int_array(var, cnt, filename);

    printf("Loading data...\n");
    LOAD_INT(nkb, "nkb"); LOAD_INT(ngms, "ngms"); LOAD_INT(nat, "nat"); LOAD_INT(ntyp, "ntyp");
    LOAD_CMPLX(rhoc, RHOC_SIZE, "rhoc");
    LOAD_CMPLX(becphi_c, BECPC_SIZE, "becphi_c");
    LOAD_CMPLX(becpsi_c, BECPC_SIZE, "becpsi_c");
    LOAD_CMPLX(qgm, QGM_ROWS * QGM_COLS, "qgm");
    LOAD_CMPLX(eigts1, EIGTS1_DIM1 * EIGTS1_DIM2, "eigts1");
    LOAD_CMPLX(eigts2, EIGTS2_DIM1 * EIGTS2_DIM2, "eigts2");
    LOAD_CMPLX(eigts3, EIGTS2_DIM1 * EIGTS2_DIM2, "eigts3");
    LOAD_CMPLX(rhoc_out, RHOC_SIZE, "rhoc_out");
    LOAD_REAL(xkq, 3, "xkq"); LOAD_REAL(xk, 3, "xk"); LOAD_REAL(tau, 3 * NAT_MAX, "tau");
    LOAD_INTA(nij_type, NTYP_MAX, "nij_type");
    LOAD_INTA(ityp, NAT_MAX, "ityp");
    LOAD_INTA(ofsbeta, NAT_MAX, "ofsbeta");
    LOAD_INTA(nh, NTYP_MAX, "nh");
    LOAD_INTA(ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, "ijtoh");
    LOAD_INTA(mill, MILL_DIM1 * MILL_DIM2, "mill");
    LOAD_INTA(dfftt__nl, NL_SIZE, "dfftt__nl");
    LOAD_INTA(upf_tvanp, NTYP_MAX, "upf_tvanp");

    transpose_2d(qgm, qgm_T, QGM_ROWS, QGM_COLS);
    transpose_2d(eigts1, eigts1_T, EIGTS1_DIM1, EIGTS1_DIM2);
    transpose_2d(eigts2, eigts2_T, EIGTS2_DIM1, EIGTS2_DIM2);
    transpose_2d(eigts3, eigts3_T, EIGTS2_DIM1, EIGTS2_DIM2);

    struct NLEntry { int val; int orig_idx; };
    NLEntry* entries = new NLEntry[ngms];
    for (int i = 0; i < ngms; i++) { entries[i].val = dfftt__nl[i]; entries[i].orig_idx = i + 1; }
    std::sort(entries, entries + ngms, [](const NLEntry& a, const NLEntry& b) { return a.val < b.val; });
    for (int i = 0; i < ngms; i++) { dfftt__nl_sorted[i] = entries[i].val; dfftt__nl_ix[i] = entries[i].orig_idx; }
    delete[] entries;
}

// ============================================================
// Device arrays — AoS
// ============================================================
struct DeviceArrays {
    Complex_DP *d_rhoc, *d_becphi_c, *d_becpsi_c;
    Complex_DP *d_qgm, *d_eigts1, *d_eigts2, *d_eigts3;
    Complex_DP *d_qgm_T, *d_eigts1_T, *d_eigts2_T, *d_eigts3_T;
    DP *d_xkq, *d_xk, *d_tau;
    int *d_ityp, *d_ofsbeta, *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
    Complex_DP *d_eigqts;
};

#define CUDA_AC(d_ptr, h_ptr, count, type) \
    hipMalloc(&da.d_ptr, (count) * sizeof(type)); \
    hipMemcpy(da.d_ptr, h_ptr, (count) * sizeof(type), hipMemcpyHostToDevice);

static DeviceArrays allocate_device_arrays() {
    DeviceArrays da;
    CUDA_AC(d_rhoc, rhoc, RHOC_SIZE, Complex_DP);
    CUDA_AC(d_becphi_c, becphi_c, BECPC_SIZE, Complex_DP);
    CUDA_AC(d_becpsi_c, becpsi_c, BECPC_SIZE, Complex_DP);
    CUDA_AC(d_qgm, qgm, QGM_ROWS * QGM_COLS, Complex_DP);
    CUDA_AC(d_eigts1, eigts1, EIGTS1_DIM1 * EIGTS1_DIM2, Complex_DP);
    CUDA_AC(d_eigts2, eigts2, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_AC(d_eigts3, eigts3, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_AC(d_qgm_T, qgm_T, QGM_ROWS * QGM_COLS, Complex_DP);
    CUDA_AC(d_eigts1_T, eigts1_T, EIGTS1_DIM1 * EIGTS1_DIM2, Complex_DP);
    CUDA_AC(d_eigts2_T, eigts2_T, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_AC(d_eigts3_T, eigts3_T, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_AC(d_xkq, xkq, 3, DP); CUDA_AC(d_xk, xk, 3, DP);
    CUDA_AC(d_tau, tau, 3 * NAT_MAX, DP);
    CUDA_AC(d_ityp, ityp, NAT_MAX, int);
    CUDA_AC(d_ofsbeta, ofsbeta, NAT_MAX, int);
    CUDA_AC(d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, int);
    CUDA_AC(d_mill, mill, MILL_DIM1 * MILL_DIM2, int);
    CUDA_AC(d_dfftt__nl, dfftt__nl, NL_SIZE, int);
    CUDA_AC(d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE, int);
    CUDA_AC(d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE, int);
    hipMalloc(&da.d_eigqts, nat * sizeof(Complex_DP));
    return da;
}

static void free_device_arrays(DeviceArrays& da) {
    hipFree(da.d_rhoc); hipFree(da.d_becphi_c); hipFree(da.d_becpsi_c);
    hipFree(da.d_qgm); hipFree(da.d_eigts1); hipFree(da.d_eigts2); hipFree(da.d_eigts3);
    hipFree(da.d_qgm_T); hipFree(da.d_eigts1_T); hipFree(da.d_eigts2_T); hipFree(da.d_eigts3_T);
    hipFree(da.d_xkq); hipFree(da.d_xk); hipFree(da.d_tau);
    hipFree(da.d_ityp); hipFree(da.d_ofsbeta); hipFree(da.d_ijtoh);
    hipFree(da.d_mill); hipFree(da.d_dfftt__nl);
    hipFree(da.d_dfftt__nl_sorted); hipFree(da.d_dfftt__nl_ix);
    hipFree(da.d_eigqts);
}

static void reset_rhoc_device(DeviceArrays& da) {
    hipMemcpy(da.d_rhoc, rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyHostToDevice);
}

// ============================================================
// Device arrays — SoA
// ============================================================
struct DeviceArraysSoA {
    double *d_rhoc_re, *d_rhoc_im;
    double *d_becphi_re, *d_becphi_im, *d_becpsi_re, *d_becpsi_im;
    double *d_qgm_re, *d_qgm_im, *d_eigts1_re, *d_eigts1_im;
    double *d_eigts2_re, *d_eigts2_im, *d_eigts3_re, *d_eigts3_im;
    double *d_qgm_T_re, *d_qgm_T_im, *d_eigts1_T_re, *d_eigts1_T_im;
    double *d_eigts2_T_re, *d_eigts2_T_im, *d_eigts3_T_re, *d_eigts3_T_im;
    DP *d_xkq, *d_xk, *d_tau;
    int *d_ityp, *d_ofsbeta, *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
    double *d_eigqts_re, *d_eigqts_im;
};

static void alloc_soa_pair(double*& d_re, double*& d_im, const Complex_DP* h_aos, int n) {
    double* h_re = new double[n]; double* h_im = new double[n];
    aos_to_soa(h_aos, h_re, h_im, n);
    hipMalloc(&d_re, n * sizeof(double)); hipMalloc(&d_im, n * sizeof(double));
    hipMemcpy(d_re, h_re, n * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(d_im, h_im, n * sizeof(double), hipMemcpyHostToDevice);
    delete[] h_re; delete[] h_im;
}

#define SOA_INT_AC(d_ptr, h_ptr, count) \
    hipMalloc(&ds.d_ptr, (count) * sizeof(int)); \
    hipMemcpy(ds.d_ptr, h_ptr, (count) * sizeof(int), hipMemcpyHostToDevice);
#define SOA_DP_AC(d_ptr, h_ptr, count) \
    hipMalloc(&ds.d_ptr, (count) * sizeof(DP)); \
    hipMemcpy(ds.d_ptr, h_ptr, (count) * sizeof(DP), hipMemcpyHostToDevice);

static DeviceArraysSoA allocate_device_arrays_soa() {
    DeviceArraysSoA ds;
    alloc_soa_pair(ds.d_rhoc_re, ds.d_rhoc_im, rhoc, RHOC_SIZE);
    alloc_soa_pair(ds.d_becphi_re, ds.d_becphi_im, becphi_c, BECPC_SIZE);
    alloc_soa_pair(ds.d_becpsi_re, ds.d_becpsi_im, becpsi_c, BECPC_SIZE);
    alloc_soa_pair(ds.d_qgm_re, ds.d_qgm_im, qgm, QGM_ROWS * QGM_COLS);
    alloc_soa_pair(ds.d_eigts1_re, ds.d_eigts1_im, eigts1, EIGTS1_DIM1 * EIGTS1_DIM2);
    alloc_soa_pair(ds.d_eigts2_re, ds.d_eigts2_im, eigts2, EIGTS2_DIM1 * EIGTS2_DIM2);
    alloc_soa_pair(ds.d_eigts3_re, ds.d_eigts3_im, eigts3, EIGTS2_DIM1 * EIGTS2_DIM2);
    alloc_soa_pair(ds.d_qgm_T_re, ds.d_qgm_T_im, qgm_T, QGM_ROWS * QGM_COLS);
    alloc_soa_pair(ds.d_eigts1_T_re, ds.d_eigts1_T_im, eigts1_T, EIGTS1_DIM1 * EIGTS1_DIM2);
    alloc_soa_pair(ds.d_eigts2_T_re, ds.d_eigts2_T_im, eigts2_T, EIGTS2_DIM1 * EIGTS2_DIM2);
    alloc_soa_pair(ds.d_eigts3_T_re, ds.d_eigts3_T_im, eigts3_T, EIGTS2_DIM1 * EIGTS2_DIM2);
    SOA_DP_AC(d_xkq, xkq, 3); SOA_DP_AC(d_xk, xk, 3); SOA_DP_AC(d_tau, tau, 3 * NAT_MAX);
    SOA_INT_AC(d_ityp, ityp, NAT_MAX);
    SOA_INT_AC(d_ofsbeta, ofsbeta, NAT_MAX);
    SOA_INT_AC(d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3);
    SOA_INT_AC(d_mill, mill, MILL_DIM1 * MILL_DIM2);
    SOA_INT_AC(d_dfftt__nl, dfftt__nl, NL_SIZE);
    SOA_INT_AC(d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE);
    SOA_INT_AC(d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE);
    hipMalloc(&ds.d_eigqts_re, nat * sizeof(double));
    hipMalloc(&ds.d_eigqts_im, nat * sizeof(double));
    return ds;
}

static void free_device_arrays_soa(DeviceArraysSoA& ds) {
    hipFree(ds.d_rhoc_re); hipFree(ds.d_rhoc_im);
    hipFree(ds.d_becphi_re); hipFree(ds.d_becphi_im);
    hipFree(ds.d_becpsi_re); hipFree(ds.d_becpsi_im);
    hipFree(ds.d_qgm_re); hipFree(ds.d_qgm_im);
    hipFree(ds.d_eigts1_re); hipFree(ds.d_eigts1_im);
    hipFree(ds.d_eigts2_re); hipFree(ds.d_eigts2_im);
    hipFree(ds.d_eigts3_re); hipFree(ds.d_eigts3_im);
    hipFree(ds.d_qgm_T_re); hipFree(ds.d_qgm_T_im);
    hipFree(ds.d_eigts1_T_re); hipFree(ds.d_eigts1_T_im);
    hipFree(ds.d_eigts2_T_re); hipFree(ds.d_eigts2_T_im);
    hipFree(ds.d_eigts3_T_re); hipFree(ds.d_eigts3_T_im);
    hipFree(ds.d_xkq); hipFree(ds.d_xk); hipFree(ds.d_tau);
    hipFree(ds.d_ityp); hipFree(ds.d_ofsbeta); hipFree(ds.d_ijtoh);
    hipFree(ds.d_mill); hipFree(ds.d_dfftt__nl);
    hipFree(ds.d_dfftt__nl_sorted); hipFree(ds.d_dfftt__nl_ix);
    hipFree(ds.d_eigqts_re); hipFree(ds.d_eigqts_im);
}

static void reset_rhoc_soa(DeviceArraysSoA& ds) {
    double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
    aos_to_soa(rhoc, h_re, h_im, RHOC_SIZE);
    hipMemcpy(ds.d_rhoc_re, h_re, RHOC_SIZE * sizeof(double), hipMemcpyHostToDevice);
    hipMemcpy(ds.d_rhoc_im, h_im, RHOC_SIZE * sizeof(double), hipMemcpyHostToDevice);
    delete[] h_re; delete[] h_im;
}

// ============================================================
// Correctness
// ============================================================
static bool check_aos(DeviceArrays& da) {
    hipMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyDeviceToHost);
    for (int i = 0; i < RHOC_SIZE; i++)
        if (cabs_val(csub(rhoc_out_sim[i], rhoc_out[i])) >= 1e-8) return false;
    return true;
}

static bool check_soa(DeviceArraysSoA& ds) {
    double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
    hipMemcpy(h_re, ds.d_rhoc_re, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
    hipMemcpy(h_im, ds.d_rhoc_im, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
    bool ok = true;
    for (int i = 0; i < RHOC_SIZE && ok; i++) {
        double dr = h_re[i] - creal_val(rhoc_out[i]);
        double di = h_im[i] - cimag_val(rhoc_out[i]);
        if (sqrt(dr*dr+di*di) >= 1e-8) ok = false;
    }
    delete[] h_re; delete[] h_im;
    return ok;
}

// ============================================================
// Profiling — returns all iteration times
// ============================================================
static std::vector<float> profile_kernel(
    std::function<void()> reset_fn,
    std::function<void()> kernel_fn)
{
    hipEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventCreate(&start_ev[i]); hipEventCreate(&stop_ev[i]);
    }
    reset_fn();
    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventRecord(start_ev[i], 0);
        kernel_fn();
        hipEventRecord(stop_ev[i], 0);
        hipEventSynchronize(stop_ev[i]);
    }
    std::vector<float> times;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        hipEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        times.push_back(ms[i]);
    }
    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventDestroy(start_ev[i]); hipEventDestroy(stop_ev[i]);
    }
    return times;
}

static void write_csv(FILE* csv, const char* vname, int tbs, int cf,
                      const std::vector<float>& times) {
    for (int r = 0; r < (int)times.size(); r++)
        fprintf(csv, "%s,%d,%d,%d,%.6f\n", vname, tbs, cf, r, times[r]);
}

static float median(const std::vector<float>& v) {
    std::vector<float> s = v;
    std::sort(s.begin(), s.end());
    return s[s.size()/2];
}

// ============================================================
// Sweep helper macro — reduces boilerplate per variant
// ============================================================
#define RUN_AOS(vbase, launch_call) do { \
    for (bool unroll : {false, true}) { \
        char vn[64]; snprintf(vn, 64, "%s%s", vbase, unroll ? "_u" : ""); \
        reset_rhoc_device(da); \
        launch_call; \
        bool ok = check_aos(da); \
        auto reset = [&](){ reset_rhoc_device(da); }; \
        auto kern  = [&](){ launch_call; }; \
        auto times = profile_kernel(reset, kern); \
        write_csv(csv, vn, tbs, cf, times); \
        printf("%-34s %6d %4d %10.4f %s\n", vn, tbs, cf, median(times), ok?"PASS":"FAIL"); \
    } \
} while(0)

#define RUN_SOA(vbase, launch_call) do { \
    for (bool unroll : {false, true}) { \
        char vn[64]; snprintf(vn, 64, "%s%s", vbase, unroll ? "_u" : ""); \
        reset_rhoc_soa(ds); \
        launch_call; \
        bool ok = check_soa(ds); \
        auto reset = [&](){ reset_rhoc_soa(ds); }; \
        auto kern  = [&](){ launch_call; }; \
        auto times = profile_kernel(reset, kern); \
        write_csv(csv, vn, tbs, cf, times); \
        printf("%-34s %6d %4d %10.4f %s\n", vn, tbs, cf, median(times), ok?"PASS":"FAIL"); \
    } \
} while(0)

// ============================================================
// GPU sweep
// ============================================================
static void sweep_gpu(FILE* csv, DeviceArrays& da, DeviceArraysSoA& ds) {
    printf("\n=== GPU SWEEP: variant x tblock x coarsen x unroll ===\n");
    printf("%-34s %6s %4s %10s %s\n", "variant", "tblock", "cf", "median_ms", "ok");
    printf("------------------------------------------------------------\n");

    for (int ti = 0; ti < N_TBLOCK; ti++) {
        int tbs = TBLOCK_SIZES[ti];
        for (int ci = 0; ci < N_COARSEN; ci++) {
            int cf = COARSEN_FACTORS[ci];

            // AoS baseline
            RUN_AOS("gpu_baseline_aos",
                addusxx_g_gpu(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_ityp, da.d_ofsbeta, da.d_ijtoh,
                    da.d_qgm, da.d_eigts1, da.d_eigts2, da.d_eigts3,
                    da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll, da.d_eigqts));

            // AoS eigts-transposed
            RUN_AOS("gpu_eigts_t_aos",
                addusxx_g_gpu_eigts_t(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_ityp, da.d_ofsbeta, da.d_ijtoh,
                    da.d_qgm, da.d_eigts1_T, da.d_eigts2_T, da.d_eigts3_T,
                    da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll, da.d_eigqts));

            // AoS shared-bec
            RUN_AOS("gpu_shared_bec_aos",
                addusxx_g_gpu_shared_bec(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_ityp, da.d_ofsbeta, da.d_ijtoh,
                    da.d_qgm, da.d_eigts1_T, da.d_eigts2_T, da.d_eigts3_T,
                    da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll, da.d_eigqts));

            // AoS optimized (sorted fused)
            RUN_AOS("gpu_optimized_aos",
                addusxx_g_gpu_optimized(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_ityp, da.d_ofsbeta, da.d_ijtoh,
                    da.d_qgm_T, da.d_eigts1_T, da.d_eigts2_T, da.d_eigts3_T,
                    da.d_mill, da.d_dfftt__nl_sorted, da.d_dfftt__nl_ix,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll, da.d_eigqts));

            // SoA baseline
            RUN_SOA("gpu_baseline_soa",
                addusxx_g_gpu_soa(ds.d_rhoc_re, ds.d_rhoc_im,
                    ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im,
                    ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_ityp, ds.d_ofsbeta, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_re, ds.d_eigts1_im,
                    ds.d_eigts2_re, ds.d_eigts2_im,
                    ds.d_eigts3_re, ds.d_eigts3_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll,
                    ds.d_eigqts_re, ds.d_eigqts_im));

            // SoA eigts-transposed
            RUN_SOA("gpu_eigts_t_soa",
                addusxx_g_gpu_eigts_t_soa(ds.d_rhoc_re, ds.d_rhoc_im,
                    ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im,
                    ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_ityp, ds.d_ofsbeta, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll,
                    ds.d_eigqts_re, ds.d_eigqts_im));

            // SoA shared-bec
            RUN_SOA("gpu_shared_bec_soa",
                addusxx_g_gpu_shared_bec_soa(ds.d_rhoc_re, ds.d_rhoc_im,
                    ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im,
                    ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_ityp, ds.d_ofsbeta, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll,
                    ds.d_eigqts_re, ds.d_eigqts_im));

            // SoA optimized (sorted fused)
            RUN_SOA("gpu_optimized_soa",
                addusxx_g_gpu_optimized_soa(ds.d_rhoc_re, ds.d_rhoc_im,
                    ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im,
                    ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_ityp, ds.d_ofsbeta, ds.d_ijtoh,
                    ds.d_qgm_T_re, ds.d_qgm_T_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl_sorted, ds.d_dfftt__nl_ix,
                    upf_tvanp, nij_type, nh, tbs, cf, unroll,
                    ds.d_eigqts_re, ds.d_eigqts_im));

            fflush(csv);
        }
    }
}

int main() {
    init_addusxx_data();

    DeviceArrays da = allocate_device_arrays();
    DeviceArraysSoA ds = allocate_device_arrays_soa();

    FILE* csv = fopen("addusxx_gpu_sweep.csv", "w");
    fprintf(csv, "variant,tblock,coarsen,rep,time_ms\n");

    sweep_gpu(csv, da, ds);

    fclose(csv);
    printf("\nDone. CSV: addusxx_gpu_sweep.csv\n");

    free_device_arrays(da);
    free_device_arrays_soa(ds);
    return 0;
}
