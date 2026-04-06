#include "types.h"
#include "data_loading.h"
#include "usxx_kernels.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

// Array dimensions (matching Fortran ALLOCATEs)
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
static const int COARSEN_FACTORS[] = { 1, 2, 4, 8 };
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

    for (int i = 0; i < ngms; i++) {
        if (dfftt__nl[dfftt__nl_ix[i] - 1] != dfftt__nl_sorted[i]) {
            printf("Error in dfftt__nl_sorted/ix.\n"); break;
        }
    }
}

// ============================================================
// Device arrays — AoS (includes scratch buffers)
// ============================================================
struct DeviceArrays {
    Complex_DP *d_rhoc, *d_becphi_c, *d_becpsi_c;
    Complex_DP *d_qgm, *d_eigts1, *d_eigts2, *d_eigts3;
    Complex_DP *d_qgm_T, *d_eigts1_T, *d_eigts2_T, *d_eigts3_T;
    DP *d_xkq, *d_xk, *d_tau;
    int *d_upf_tvanp, *d_nij_type, *d_ityp, *d_ofsbeta, *d_nh;
    int *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
    // Scratch buffers (pre-allocated)
    Complex_DP *d_eigqts;        // size: nat
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
    CUDA_AC(d_xkq, xkq, 3, DP);
    CUDA_AC(d_xk, xk, 3, DP);
    CUDA_AC(d_tau, tau, 3 * NAT_MAX, DP);
    CUDA_AC(d_upf_tvanp, upf_tvanp, NTYP_MAX, int);
    CUDA_AC(d_nij_type, nij_type, NTYP_MAX, int);
    CUDA_AC(d_ityp, ityp, NAT_MAX, int);
    CUDA_AC(d_ofsbeta, ofsbeta, NAT_MAX, int);
    CUDA_AC(d_nh, nh, NTYP_MAX, int);
    CUDA_AC(d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, int);
    CUDA_AC(d_mill, mill, MILL_DIM1 * MILL_DIM2, int);
    CUDA_AC(d_dfftt__nl, dfftt__nl, NL_SIZE, int);
    CUDA_AC(d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE, int);
    CUDA_AC(d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE, int);
    // Scratch
    hipMalloc(&da.d_eigqts, nat * sizeof(Complex_DP));
    return da;
}

static void free_device_arrays(DeviceArrays& da) {
    hipFree(da.d_rhoc); hipFree(da.d_becphi_c); hipFree(da.d_becpsi_c);
    hipFree(da.d_qgm); hipFree(da.d_eigts1); hipFree(da.d_eigts2); hipFree(da.d_eigts3);
    hipFree(da.d_qgm_T); hipFree(da.d_eigts1_T); hipFree(da.d_eigts2_T); hipFree(da.d_eigts3_T);
    hipFree(da.d_xkq); hipFree(da.d_xk); hipFree(da.d_tau);
    hipFree(da.d_upf_tvanp); hipFree(da.d_nij_type); hipFree(da.d_ityp);
    hipFree(da.d_ofsbeta); hipFree(da.d_nh); hipFree(da.d_ijtoh);
    hipFree(da.d_mill); hipFree(da.d_dfftt__nl);
    hipFree(da.d_dfftt__nl_sorted); hipFree(da.d_dfftt__nl_ix);
    hipFree(da.d_eigqts);
}

static void reset_rhoc_device(DeviceArrays& da) {
    hipMemcpy(da.d_rhoc, rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyHostToDevice);
}

// ============================================================
// Device arrays — SoA (includes scratch buffers)
// ============================================================
struct DeviceArraysSoA {
    double *d_rhoc_re, *d_rhoc_im;
    double *d_becphi_re, *d_becphi_im, *d_becpsi_re, *d_becpsi_im;
    double *d_qgm_re, *d_qgm_im, *d_eigts1_re, *d_eigts1_im;
    double *d_eigts2_re, *d_eigts2_im, *d_eigts3_re, *d_eigts3_im;
    double *d_qgm_T_re, *d_qgm_T_im, *d_eigts1_T_re, *d_eigts1_T_im;
    double *d_eigts2_T_re, *d_eigts2_T_im, *d_eigts3_T_re, *d_eigts3_T_im;
    DP *d_xkq, *d_xk, *d_tau;
    int *d_upf_tvanp, *d_nij_type, *d_ityp, *d_ofsbeta, *d_nh;
    int *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
    // Scratch buffers (pre-allocated)
    double *d_eigqts_re, *d_eigqts_im;    // size: nat each
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
    SOA_INT_AC(d_upf_tvanp, upf_tvanp, NTYP_MAX);
    SOA_INT_AC(d_nij_type, nij_type, NTYP_MAX);
    SOA_INT_AC(d_ityp, ityp, NAT_MAX);
    SOA_INT_AC(d_ofsbeta, ofsbeta, NAT_MAX);
    SOA_INT_AC(d_nh, nh, NTYP_MAX);
    SOA_INT_AC(d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3);
    SOA_INT_AC(d_mill, mill, MILL_DIM1 * MILL_DIM2);
    SOA_INT_AC(d_dfftt__nl, dfftt__nl, NL_SIZE);
    SOA_INT_AC(d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE);
    SOA_INT_AC(d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE);
    // Scratch
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
    hipFree(ds.d_upf_tvanp); hipFree(ds.d_nij_type); hipFree(ds.d_ityp);
    hipFree(ds.d_ofsbeta); hipFree(ds.d_nh); hipFree(ds.d_ijtoh);
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
// Correctness checks
// ============================================================
static bool check_correctness(const Complex_DP* sim, const Complex_DP* ref, int n, const char* label) {
    double maxdiff = 0.0;
    for (int i = 0; i < n; i++) {
        double d = cabs_val(csub(sim[i], ref[i]));
        if (d > maxdiff) maxdiff = d;
    }
    if (maxdiff >= 1.0e-8) {
        printf("  %s -- FAILED. Max diff: %12.4e\n", label, maxdiff);
        return false;
    }
    printf("  %s -- PASSED\n", label);
    return true;
}

static bool check_correctness_soa(DeviceArraysSoA& ds, const Complex_DP* ref, int n, const char* label) {
    double* h_re = new double[n]; double* h_im = new double[n];
    hipMemcpy(h_re, ds.d_rhoc_re, n * sizeof(double), hipMemcpyDeviceToHost);
    hipMemcpy(h_im, ds.d_rhoc_im, n * sizeof(double), hipMemcpyDeviceToHost);
    double maxdiff = 0.0;
    for (int i = 0; i < n; i++) {
        double dr = h_re[i] - creal_val(ref[i]);
        double di = h_im[i] - cimag_val(ref[i]);
        double d = sqrt(dr * dr + di * di);
        if (d > maxdiff) maxdiff = d;
    }
    delete[] h_re; delete[] h_im;
    if (maxdiff >= 1.0e-8) {
        printf("  %s -- FAILED. Max diff: %12.4e\n", label, maxdiff);
        return false;
    }
    printf("  %s -- PASSED\n", label);
    return true;
}

// ============================================================
// Generic GPU profiling
// ============================================================
static float profile_kernel(
    std::function<void()> reset_fn,
    std::function<void()> kernel_fn)
{
    hipEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventCreate(&start_ev[i]);
        hipEventCreate(&stop_ev[i]);
    }

    reset_fn();
    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventRecord(start_ev[i], 0);
        kernel_fn();
        hipEventRecord(stop_ev[i], 0);
        hipEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        hipEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }

    for (int i = 0; i < TOTAL_ITERS; i++) {
        hipEventDestroy(start_ev[i]);
        hipEventDestroy(stop_ev[i]);
    }
    return total / NUM_ITERS;
}

// ============================================================
// CPU baseline
// ============================================================
static void profile_cpu_original() {
    printf("\n=== (1) CPU baseline ===\n");
    memcpy(rhoc_out_sim, rhoc, RHOC_SIZE * sizeof(Complex_DP));
    addusxx_g_cpu(rhoc_out_sim, xkq, xk, tau, becphi_c, becpsi_c,
                  nkb, ngms, nat, ntyp, upf_tvanp, nij_type,
                  ityp, ofsbeta, nh, ijtoh, qgm, eigts1, eigts2, eigts3,
                  mill, dfftt__nl);
    check_correctness(rhoc_out_sim, rhoc_out, RHOC_SIZE, "CPU");

    auto reset = [&]() { memcpy(rhoc_out_sim, rhoc, RHOC_SIZE * sizeof(Complex_DP)); };
    auto kernel = [&]() {
        addusxx_g_cpu(rhoc_out_sim, xkq, xk, tau, becphi_c, becpsi_c,
                      nkb, ngms, nat, ntyp, upf_tvanp, nij_type,
                      ityp, ofsbeta, nh, ijtoh, qgm, eigts1, eigts2, eigts3,
                      mill, dfftt__nl);
    };
    float avg = profile_kernel(reset, kernel);
    printf("  time: %13.6f ms\n", avg);
}

// ============================================================
// GPU sweep — all 4 variants x tblock_sizes x coarsen_factors
// ============================================================
static void sweep_gpu(DeviceArrays& da, DeviceArraysSoA& ds) {
    printf("\n=== GPU SWEEP: variant x tblock_size x coarsen ===\n");
    printf("%-30s %8s %8s %13s %8s\n", "variant", "tblock", "coarsen", "time_ms", "correct");
    printf("----------------------------------------------------------------------\n");

    FILE* csv = fopen("addusxx_gpu_sweep.csv", "w");
    fprintf(csv, "variant,tblock,coarsen,time_ms,correct\n");

    for (int ti = 0; ti < N_TBLOCK; ti++) {
        int tbs = TBLOCK_SIZES[ti];
        for (int ci = 0; ci < N_COARSEN; ci++) {
            int cf = COARSEN_FACTORS[ci];

            // ---- (2) GPU baseline AoS ----
            {
                reset_rhoc_device(da);
                addusxx_g_gpu(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                    da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1, da.d_eigts2,
                    da.d_eigts3, da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    da.d_eigqts);
                hipMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++)
                    if (cabs_val(csub(rhoc_out_sim[i], rhoc_out[i])) >= 1.0e-8) ok = false;

                auto reset = [&]() { reset_rhoc_device(da); };
                auto kernel = [&]() {
                    addusxx_g_gpu(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                        da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                        da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                        da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1, da.d_eigts2,
                        da.d_eigts3, da.d_mill, da.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        da.d_eigqts);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_baseline_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_baseline_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (2b) GPU eigts-transposed AoS ----
            {
                reset_rhoc_device(da);
                addusxx_g_gpu_eigts_transposed(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                    da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1_T,
                    da.d_eigts2_T, da.d_eigts3_T, da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    da.d_eigqts);
                hipMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++)
                    if (cabs_val(csub(rhoc_out_sim[i], rhoc_out[i])) >= 1.0e-8) ok = false;

                auto reset = [&]() { reset_rhoc_device(da); };
                auto kernel = [&]() {
                    addusxx_g_gpu_eigts_transposed(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                        da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                        da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                        da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1_T,
                        da.d_eigts2_T, da.d_eigts3_T, da.d_mill, da.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        da.d_eigqts);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_eigts_t_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_eigts_t_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (2c) GPU shared-bec AoS ----
            {
                reset_rhoc_device(da);
                addusxx_g_gpu_shared_bec(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                    da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1_T,
                    da.d_eigts2_T, da.d_eigts3_T, da.d_mill, da.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    da.d_eigqts);
                hipMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++)
                    if (cabs_val(csub(rhoc_out_sim[i], rhoc_out[i])) >= 1.0e-8) ok = false;

                auto reset = [&]() { reset_rhoc_device(da); };
                auto kernel = [&]() {
                    addusxx_g_gpu_shared_bec(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                        da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                        da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                        da.d_nh, da.d_ijtoh, da.d_qgm, da.d_eigts1_T,
                        da.d_eigts2_T, da.d_eigts3_T, da.d_mill, da.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        da.d_eigqts);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_shared_bec_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_shared_bec_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (3) GPU optimized AoS ----
            {
                reset_rhoc_device(da);
                addusxx_g_gpu_optimized(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                    da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                    da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                    da.d_nh, da.d_ijtoh, da.d_qgm_T, da.d_eigts1_T,
                    da.d_eigts2_T, da.d_eigts3_T, da.d_mill,
                    da.d_dfftt__nl_sorted, da.d_dfftt__nl_ix,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    da.d_eigqts);
                hipMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++)
                    if (cabs_val(csub(rhoc_out_sim[i], rhoc_out[i])) >= 1.0e-8) ok = false;

                auto reset = [&]() { reset_rhoc_device(da); };
                auto kernel = [&]() {
                    addusxx_g_gpu_optimized(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                        da.d_becphi_c, da.d_becpsi_c, nkb, ngms, nat, ntyp,
                        da.d_upf_tvanp, da.d_nij_type, da.d_ityp, da.d_ofsbeta,
                        da.d_nh, da.d_ijtoh, da.d_qgm_T, da.d_eigts1_T,
                        da.d_eigts2_T, da.d_eigts3_T, da.d_mill,
                        da.d_dfftt__nl_sorted, da.d_dfftt__nl_ix,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        da.d_eigqts);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_optimized_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_optimized_aos", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (4) GPU baseline SoA ----
            {
                reset_rhoc_soa(ds);
                addusxx_g_gpu_soa(
                    ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                    ds.d_nh, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_re, ds.d_eigts1_im,
                    ds.d_eigts2_re, ds.d_eigts2_im,
                    ds.d_eigts3_re, ds.d_eigts3_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    ds.d_eigqts_re, ds.d_eigqts_im);
                double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
                hipMemcpy(h_re, ds.d_rhoc_re, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                hipMemcpy(h_im, ds.d_rhoc_im, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++) {
                    double dr = h_re[i] - creal_val(rhoc_out[i]);
                    double di = h_im[i] - cimag_val(rhoc_out[i]);
                    if (sqrt(dr*dr + di*di) >= 1.0e-8) ok = false;
                }
                delete[] h_re; delete[] h_im;

                auto reset = [&]() { reset_rhoc_soa(ds); };
                auto kernel = [&]() {
                    addusxx_g_gpu_soa(
                        ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                        ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                        nkb, ngms, nat, ntyp,
                        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                        ds.d_nh, ds.d_ijtoh,
                        ds.d_qgm_re, ds.d_qgm_im,
                        ds.d_eigts1_re, ds.d_eigts1_im,
                        ds.d_eigts2_re, ds.d_eigts2_im,
                        ds.d_eigts3_re, ds.d_eigts3_im,
                        ds.d_mill, ds.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        ds.d_eigqts_re, ds.d_eigqts_im);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_baseline_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_baseline_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (4b) GPU eigts-transposed SoA ----
            {
                reset_rhoc_soa(ds);
                addusxx_g_gpu_eigts_transposed_soa(
                    ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                    ds.d_nh, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    ds.d_eigqts_re, ds.d_eigqts_im);
                double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
                hipMemcpy(h_re, ds.d_rhoc_re, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                hipMemcpy(h_im, ds.d_rhoc_im, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++) {
                    double dr = h_re[i] - creal_val(rhoc_out[i]);
                    double di = h_im[i] - cimag_val(rhoc_out[i]);
                    if (sqrt(dr*dr + di*di) >= 1.0e-8) ok = false;
                }
                delete[] h_re; delete[] h_im;

                auto reset = [&]() { reset_rhoc_soa(ds); };
                auto kernel = [&]() {
                    addusxx_g_gpu_eigts_transposed_soa(
                        ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                        ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                        nkb, ngms, nat, ntyp,
                        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                        ds.d_nh, ds.d_ijtoh,
                        ds.d_qgm_re, ds.d_qgm_im,
                        ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                        ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                        ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                        ds.d_mill, ds.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        ds.d_eigqts_re, ds.d_eigqts_im);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_eigts_t_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_eigts_t_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (4c) GPU shared-bec SoA ----
            {
                reset_rhoc_soa(ds);
                addusxx_g_gpu_shared_bec_soa(
                    ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                    ds.d_nh, ds.d_ijtoh,
                    ds.d_qgm_re, ds.d_qgm_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    ds.d_eigqts_re, ds.d_eigqts_im);
                double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
                hipMemcpy(h_re, ds.d_rhoc_re, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                hipMemcpy(h_im, ds.d_rhoc_im, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++) {
                    double dr = h_re[i] - creal_val(rhoc_out[i]);
                    double di = h_im[i] - cimag_val(rhoc_out[i]);
                    if (sqrt(dr*dr + di*di) >= 1.0e-8) ok = false;
                }
                delete[] h_re; delete[] h_im;

                auto reset = [&]() { reset_rhoc_soa(ds); };
                auto kernel = [&]() {
                    addusxx_g_gpu_shared_bec_soa(
                        ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                        ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                        nkb, ngms, nat, ntyp,
                        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                        ds.d_nh, ds.d_ijtoh,
                        ds.d_qgm_re, ds.d_qgm_im,
                        ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                        ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                        ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                        ds.d_mill, ds.d_dfftt__nl,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        ds.d_eigqts_re, ds.d_eigqts_im);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_shared_bec_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_shared_bec_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }

            // ---- (5) GPU optimized SoA ----
            {
                reset_rhoc_soa(ds);
                addusxx_g_gpu_optimized_soa(
                    ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                    ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                    nkb, ngms, nat, ntyp,
                    ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                    ds.d_nh, ds.d_ijtoh,
                    ds.d_qgm_T_re, ds.d_qgm_T_im,
                    ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                    ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                    ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                    ds.d_mill, ds.d_dfftt__nl_sorted, ds.d_dfftt__nl_ix,
                    upf_tvanp, nij_type, nh, tbs, cf,
                    ds.d_eigqts_re, ds.d_eigqts_im);
                double* h_re = new double[RHOC_SIZE]; double* h_im = new double[RHOC_SIZE];
                hipMemcpy(h_re, ds.d_rhoc_re, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                hipMemcpy(h_im, ds.d_rhoc_im, RHOC_SIZE * sizeof(double), hipMemcpyDeviceToHost);
                bool ok = true;
                for (int i = 0; i < RHOC_SIZE && ok; i++) {
                    double dr = h_re[i] - creal_val(rhoc_out[i]);
                    double di = h_im[i] - cimag_val(rhoc_out[i]);
                    if (sqrt(dr*dr + di*di) >= 1.0e-8) ok = false;
                }
                delete[] h_re; delete[] h_im;

                auto reset = [&]() { reset_rhoc_soa(ds); };
                auto kernel = [&]() {
                    addusxx_g_gpu_optimized_soa(
                        ds.d_rhoc_re, ds.d_rhoc_im, ds.d_xkq, ds.d_xk, ds.d_tau,
                        ds.d_becphi_re, ds.d_becphi_im, ds.d_becpsi_re, ds.d_becpsi_im,
                        nkb, ngms, nat, ntyp,
                        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp, ds.d_ofsbeta,
                        ds.d_nh, ds.d_ijtoh,
                        ds.d_qgm_T_re, ds.d_qgm_T_im,
                        ds.d_eigts1_T_re, ds.d_eigts1_T_im,
                        ds.d_eigts2_T_re, ds.d_eigts2_T_im,
                        ds.d_eigts3_T_re, ds.d_eigts3_T_im,
                        ds.d_mill, ds.d_dfftt__nl_sorted, ds.d_dfftt__nl_ix,
                        upf_tvanp, nij_type, nh, tbs, cf,
                        ds.d_eigqts_re, ds.d_eigqts_im);
                };
                float avg = profile_kernel(reset, kernel);
                printf("%-30s %8d %8d %13.6f %8s\n", "gpu_optimized_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
                fprintf(csv, "%s,%d,%d,%.6f,%s\n", "gpu_optimized_soa", tbs, cf, avg, ok ? "PASS" : "FAIL");
            }
        }
    }

    fclose(csv);
    printf("\nCSV: addusxx_gpu_sweep.csv\n");
}

int main() {
    init_addusxx_data();

    DeviceArrays da = allocate_device_arrays();
    DeviceArraysSoA ds = allocate_device_arrays_soa();

    profile_cpu_original();
    sweep_gpu(da, ds);

    free_device_arrays(da);
    free_device_arrays_soa(ds);
    return 0;
}
