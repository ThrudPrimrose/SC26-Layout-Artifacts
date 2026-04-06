#include "types.h"
#include "data_loading.h"
#include "usxx_kernels.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>

// Array dimensions (matching Fortran ALLOCATEs)
static constexpr int RHOC_SIZE   = 120000;
static constexpr int BECPC_SIZE  = 122;
static constexpr int QGM_ROWS   = 55191;
static constexpr int QGM_COLS   = 397;
static constexpr int EIGTS1_DIM1 = 217;  // -108:108
static constexpr int EIGTS1_DIM2 = 10;
static constexpr int EIGTS2_DIM1 = 109;  // -54:54
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

// ---- Host arrays ----
static int nkb, ngms, nat, ntyp;

static Complex_DP rhoc[RHOC_SIZE], becphi_c[BECPC_SIZE], becpsi_c[BECPC_SIZE];
static Complex_DP qgm[QGM_COLS * QGM_ROWS]; // col-major
static Complex_DP eigts1[EIGTS1_DIM2 * EIGTS1_DIM1];
static Complex_DP eigts2[EIGTS2_DIM2 * EIGTS2_DIM1];
static Complex_DP eigts3[EIGTS2_DIM2 * EIGTS2_DIM1];
static Complex_DP rhoc_out[RHOC_SIZE], rhoc_out_sim[RHOC_SIZE];

static DP xkq[3], xk[3], tau[NAT_MAX * 3];

static int nij_type[NTYP_MAX], ityp[NAT_MAX], ofsbeta[NAT_MAX], nh[NTYP_MAX];
static int ijtoh[IJTOH_N3 * IJTOH_N2 * IJTOH_N1];
static int mill[MILL_DIM2 * MILL_DIM1];
static int dfftt__nl[NL_SIZE], upf_tvanp[NTYP_MAX];

// Transposed arrays
static Complex_DP qgm_T[QGM_ROWS * QGM_COLS]; // (397, 55191) col-major
static Complex_DP eigts1_T[EIGTS1_DIM1 * EIGTS1_DIM2]; // (10, 217) col-major
static Complex_DP eigts2_T[EIGTS2_DIM1 * EIGTS2_DIM2];
static Complex_DP eigts3_T[EIGTS2_DIM1 * EIGTS2_DIM2];

// Sorted dfftt__nl
static int dfftt__nl_sorted[NL_SIZE], dfftt__nl_ix[NL_SIZE];

static int file_id = 12;

// ---- Helper: Fortran-style TRANSPOSE for col-major 2D ----
// Fortran TRANSPOSE(A(m,n)) -> B(n,m)
// col-major A: A[j*m + i], col-major B: B[i*n + j]
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

    LOAD_INT(nkb, "nkb");
    LOAD_INT(ngms, "ngms");
    LOAD_INT(nat, "nat");
    LOAD_INT(ntyp, "ntyp");

    LOAD_CMPLX(rhoc, RHOC_SIZE, "rhoc");
    LOAD_CMPLX(becphi_c, BECPC_SIZE, "becphi_c");
    LOAD_CMPLX(becpsi_c, BECPC_SIZE, "becpsi_c");
    LOAD_CMPLX(qgm, QGM_ROWS * QGM_COLS, "qgm");
    LOAD_CMPLX(eigts1, EIGTS1_DIM1 * EIGTS1_DIM2, "eigts1");
    LOAD_CMPLX(eigts2, EIGTS2_DIM1 * EIGTS2_DIM2, "eigts2");
    LOAD_CMPLX(eigts3, EIGTS2_DIM1 * EIGTS2_DIM2, "eigts3");
    LOAD_CMPLX(rhoc_out, RHOC_SIZE, "rhoc_out");

    LOAD_REAL(xkq, 3, "xkq");
    LOAD_REAL(xk, 3, "xk");
    LOAD_REAL(tau, 3 * NAT_MAX, "tau");

    LOAD_INTA(nij_type, NTYP_MAX, "nij_type");
    LOAD_INTA(ityp, NAT_MAX, "ityp");
    LOAD_INTA(ofsbeta, NAT_MAX, "ofsbeta");
    LOAD_INTA(nh, NTYP_MAX, "nh");
    LOAD_INTA(ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, "ijtoh");
    LOAD_INTA(mill, MILL_DIM1 * MILL_DIM2, "mill");
    LOAD_INTA(dfftt__nl, NL_SIZE, "dfftt__nl");
    LOAD_INTA(upf_tvanp, NTYP_MAX, "upf_tvanp");

    // Populate transposed arrays
    transpose_2d(qgm, qgm_T, QGM_ROWS, QGM_COLS);
    transpose_2d(eigts1, eigts1_T, EIGTS1_DIM1, EIGTS1_DIM2);
    transpose_2d(eigts2, eigts2_T, EIGTS2_DIM1, EIGTS2_DIM2);
    transpose_2d(eigts3, eigts3_T, EIGTS2_DIM1, EIGTS2_DIM2);

    // Sort dfftt__nl (selection-sort matching Fortran logic)
    // Produces a sorted array and index permutation (1-based, like Fortran)
    struct NLEntry { int val; int orig_idx; };
    NLEntry* entries = new NLEntry[ngms];
    for (int i = 0; i < ngms; i++) {
        entries[i].val = dfftt__nl[i];
        entries[i].orig_idx = i + 1; // 1-based
    }
    std::sort(entries, entries + ngms, [](const NLEntry& a, const NLEntry& b) {
        return a.val < b.val;
    });
    for (int i = 0; i < ngms; i++) {
        dfftt__nl_sorted[i] = entries[i].val;
        dfftt__nl_ix[i] = entries[i].orig_idx;
    }
    delete[] entries;

    // Verify
    for (int i = 0; i < ngms; i++) {
        if (dfftt__nl[dfftt__nl_ix[i] - 1] != dfftt__nl_sorted[i]) {
            printf("Error in setting dfftt__nl_sorted and dfftt__nl_ix.\n");
            break;
        }
    }
}

static bool check_correctness(const Complex_DP* sim, const Complex_DP* ref, int n, const char* label) {
    double maxdiff = 0.0;
    for (int i = 0; i < n; i++) {
        double d = cabs_val(csub(sim[i], ref[i]));
        if (d > maxdiff) maxdiff = d;
    }
    if (maxdiff >= 1.0e-8) {
        printf("%s -- FAILED correctness check. Max diff: %12.4e\n", label, maxdiff);
        return false;
    }
    printf("%s -- PASSED correctness test.\n", label);
    return true;
}

// ============================================================
// Helper: copy all input arrays to device (returns struct of device ptrs)
// ============================================================
struct DeviceArrays {
    Complex_DP *d_rhoc, *d_becphi_c, *d_becpsi_c;
    Complex_DP *d_qgm, *d_eigts1, *d_eigts2, *d_eigts3;
    Complex_DP *d_qgm_T, *d_eigts1_T, *d_eigts2_T, *d_eigts3_T;
    DP *d_xkq, *d_xk, *d_tau;
    int *d_upf_tvanp, *d_nij_type, *d_ityp, *d_ofsbeta, *d_nh;
    int *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
};

// SoA device arrays: each complex array split into _re and _im
struct DeviceArraysSoA {
    double *d_rhoc_re, *d_rhoc_im;
    double *d_becphi_re, *d_becphi_im;
    double *d_becpsi_re, *d_becpsi_im;
    // Original layout
    double *d_qgm_re, *d_qgm_im;
    double *d_eigts1_re, *d_eigts1_im;
    double *d_eigts2_re, *d_eigts2_im;
    double *d_eigts3_re, *d_eigts3_im;
    // Transposed layout
    double *d_qgm_T_re, *d_qgm_T_im;
    double *d_eigts1_T_re, *d_eigts1_T_im;
    double *d_eigts2_T_re, *d_eigts2_T_im;
    double *d_eigts3_T_re, *d_eigts3_T_im;
    // Non-complex arrays shared with DeviceArrays
    DP *d_xkq, *d_xk, *d_tau;
    int *d_upf_tvanp, *d_nij_type, *d_ityp, *d_ofsbeta, *d_nh;
    int *d_ijtoh, *d_mill, *d_dfftt__nl;
    int *d_dfftt__nl_sorted, *d_dfftt__nl_ix;
};

#define CUDA_ALLOC_COPY(d_ptr, h_ptr, count, type) \
    cudaMalloc(&da.d_ptr, (count) * sizeof(type)); \
    cudaMemcpy(da.d_ptr, h_ptr, (count) * sizeof(type), cudaMemcpyHostToDevice);

#define CUDA_ALLOC_COPY_RAW(d_ptr, h_ptr, count, type) \
    cudaMalloc(&d_ptr, (count) * sizeof(type)); \
    cudaMemcpy(d_ptr, h_ptr, (count) * sizeof(type), cudaMemcpyHostToDevice);

static DeviceArrays allocate_device_arrays() {
    DeviceArrays da;
    CUDA_ALLOC_COPY(d_rhoc, rhoc_out_sim, RHOC_SIZE, Complex_DP);
    CUDA_ALLOC_COPY(d_becphi_c, becphi_c, BECPC_SIZE, Complex_DP);
    CUDA_ALLOC_COPY(d_becpsi_c, becpsi_c, BECPC_SIZE, Complex_DP);
    CUDA_ALLOC_COPY(d_qgm, qgm, QGM_ROWS * QGM_COLS, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts1, eigts1, EIGTS1_DIM1 * EIGTS1_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts2, eigts2, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts3, eigts3, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_qgm_T, qgm_T, QGM_ROWS * QGM_COLS, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts1_T, eigts1_T, EIGTS1_DIM1 * EIGTS1_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts2_T, eigts2_T, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_eigts3_T, eigts3_T, EIGTS2_DIM1 * EIGTS2_DIM2, Complex_DP);
    CUDA_ALLOC_COPY(d_xkq, xkq, 3, DP);
    CUDA_ALLOC_COPY(d_xk, xk, 3, DP);
    CUDA_ALLOC_COPY(d_tau, tau, 3 * NAT_MAX, DP);
    CUDA_ALLOC_COPY(d_upf_tvanp, upf_tvanp, NTYP_MAX, int);
    CUDA_ALLOC_COPY(d_nij_type, nij_type, NTYP_MAX, int);
    CUDA_ALLOC_COPY(d_ityp, ityp, NAT_MAX, int);
    CUDA_ALLOC_COPY(d_ofsbeta, ofsbeta, NAT_MAX, int);
    CUDA_ALLOC_COPY(d_nh, nh, NTYP_MAX, int);
    CUDA_ALLOC_COPY(d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, int);
    CUDA_ALLOC_COPY(d_mill, mill, MILL_DIM1 * MILL_DIM2, int);
    CUDA_ALLOC_COPY(d_dfftt__nl, dfftt__nl, NL_SIZE, int);
    CUDA_ALLOC_COPY(d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE, int);
    CUDA_ALLOC_COPY(d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE, int);
    return da;
}

static void free_device_arrays(DeviceArrays& da) {
    cudaFree(da.d_rhoc); cudaFree(da.d_becphi_c); cudaFree(da.d_becpsi_c);
    cudaFree(da.d_qgm); cudaFree(da.d_eigts1); cudaFree(da.d_eigts2); cudaFree(da.d_eigts3);
    cudaFree(da.d_qgm_T); cudaFree(da.d_eigts1_T); cudaFree(da.d_eigts2_T); cudaFree(da.d_eigts3_T);
    cudaFree(da.d_xkq); cudaFree(da.d_xk); cudaFree(da.d_tau);
    cudaFree(da.d_upf_tvanp); cudaFree(da.d_nij_type); cudaFree(da.d_ityp);
    cudaFree(da.d_ofsbeta); cudaFree(da.d_nh); cudaFree(da.d_ijtoh);
    cudaFree(da.d_mill); cudaFree(da.d_dfftt__nl);
    cudaFree(da.d_dfftt__nl_sorted); cudaFree(da.d_dfftt__nl_ix);
}

// Reset d_rhoc to rhoc (initial charge density) on device
static void reset_rhoc_device(DeviceArrays& da) {
    cudaMemcpy(da.d_rhoc, rhoc, RHOC_SIZE * sizeof(Complex_DP), cudaMemcpyHostToDevice);
}

// ============================================================
// SoA device array helpers
// ============================================================
// Allocate SoA double pair on device from AoS Complex_DP host array
static void alloc_soa_pair(double*& d_re, double*& d_im, const Complex_DP* h_aos, int n) {
    double* h_re = new double[n];
    double* h_im = new double[n];
    aos_to_soa(h_aos, h_re, h_im, n);
    cudaMalloc(&d_re, n * sizeof(double));
    cudaMalloc(&d_im, n * sizeof(double));
    cudaMemcpy(d_re, h_re, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_im, h_im, n * sizeof(double), cudaMemcpyHostToDevice);
    delete[] h_re;
    delete[] h_im;
}

static DeviceArraysSoA allocate_device_arrays_soa() {
    DeviceArraysSoA ds;
    alloc_soa_pair(ds.d_rhoc_re, ds.d_rhoc_im, rhoc, RHOC_SIZE); // will be reset before use
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
    // Non-complex arrays
    CUDA_ALLOC_COPY_RAW(ds.d_xkq, xkq, 3, DP);
    CUDA_ALLOC_COPY_RAW(ds.d_xk, xk, 3, DP);
    CUDA_ALLOC_COPY_RAW(ds.d_tau, tau, 3 * NAT_MAX, DP);
    CUDA_ALLOC_COPY_RAW(ds.d_upf_tvanp, upf_tvanp, NTYP_MAX, int);
    CUDA_ALLOC_COPY_RAW(ds.d_nij_type, nij_type, NTYP_MAX, int);
    CUDA_ALLOC_COPY_RAW(ds.d_ityp, ityp, NAT_MAX, int);
    CUDA_ALLOC_COPY_RAW(ds.d_ofsbeta, ofsbeta, NAT_MAX, int);
    CUDA_ALLOC_COPY_RAW(ds.d_nh, nh, NTYP_MAX, int);
    CUDA_ALLOC_COPY_RAW(ds.d_ijtoh, ijtoh, IJTOH_N1 * IJTOH_N2 * IJTOH_N3, int);
    CUDA_ALLOC_COPY_RAW(ds.d_mill, mill, MILL_DIM1 * MILL_DIM2, int);
    CUDA_ALLOC_COPY_RAW(ds.d_dfftt__nl, dfftt__nl, NL_SIZE, int);
    CUDA_ALLOC_COPY_RAW(ds.d_dfftt__nl_sorted, dfftt__nl_sorted, NL_SIZE, int);
    CUDA_ALLOC_COPY_RAW(ds.d_dfftt__nl_ix, dfftt__nl_ix, NL_SIZE, int);
    return ds;
}

static void free_device_arrays_soa(DeviceArraysSoA& ds) {
    cudaFree(ds.d_rhoc_re); cudaFree(ds.d_rhoc_im);
    cudaFree(ds.d_becphi_re); cudaFree(ds.d_becphi_im);
    cudaFree(ds.d_becpsi_re); cudaFree(ds.d_becpsi_im);
    cudaFree(ds.d_qgm_re); cudaFree(ds.d_qgm_im);
    cudaFree(ds.d_eigts1_re); cudaFree(ds.d_eigts1_im);
    cudaFree(ds.d_eigts2_re); cudaFree(ds.d_eigts2_im);
    cudaFree(ds.d_eigts3_re); cudaFree(ds.d_eigts3_im);
    cudaFree(ds.d_qgm_T_re); cudaFree(ds.d_qgm_T_im);
    cudaFree(ds.d_eigts1_T_re); cudaFree(ds.d_eigts1_T_im);
    cudaFree(ds.d_eigts2_T_re); cudaFree(ds.d_eigts2_T_im);
    cudaFree(ds.d_eigts3_T_re); cudaFree(ds.d_eigts3_T_im);
    cudaFree(ds.d_xkq); cudaFree(ds.d_xk); cudaFree(ds.d_tau);
    cudaFree(ds.d_upf_tvanp); cudaFree(ds.d_nij_type); cudaFree(ds.d_ityp);
    cudaFree(ds.d_ofsbeta); cudaFree(ds.d_nh); cudaFree(ds.d_ijtoh);
    cudaFree(ds.d_mill); cudaFree(ds.d_dfftt__nl);
    cudaFree(ds.d_dfftt__nl_sorted); cudaFree(ds.d_dfftt__nl_ix);
}

// Reset SoA rhoc to initial rhoc values
static void reset_rhoc_soa(DeviceArraysSoA& ds) {
    double* h_re = new double[RHOC_SIZE];
    double* h_im = new double[RHOC_SIZE];
    aos_to_soa(rhoc, h_re, h_im, RHOC_SIZE);
    cudaMemcpy(ds.d_rhoc_re, h_re, RHOC_SIZE * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(ds.d_rhoc_im, h_im, RHOC_SIZE * sizeof(double), cudaMemcpyHostToDevice);
    delete[] h_re;
    delete[] h_im;
}

// Correctness check: download SoA rhoc from device, reassemble AoS, compare
static bool check_correctness_soa(DeviceArraysSoA& ds, const Complex_DP* ref, int n, const char* label) {
    double* h_re = new double[n];
    double* h_im = new double[n];
    cudaMemcpy(h_re, ds.d_rhoc_re, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(h_im, ds.d_rhoc_im, n * sizeof(double), cudaMemcpyDeviceToHost);
    double maxdiff = 0.0;
    for (int i = 0; i < n; i++) {
        double dr = h_re[i] - cuCreal(ref[i]);
        double di = h_im[i] - cuCimag(ref[i]);
        double d = sqrt(dr * dr + di * di);
        if (d > maxdiff) maxdiff = d;
    }
    delete[] h_re;
    delete[] h_im;
    if (maxdiff >= 1.0e-8) {
        printf("%s -- FAILED correctness check. Max diff: %12.4e\n", label, maxdiff);
        return false;
    }
    printf("%s -- PASSED correctness test.\n", label);
    return true;
}

// ============================================================
// Profiling routines
// ============================================================
static void profile_cpu_original() {
    printf("\n(1) CPU baseline\n");
    memcpy(rhoc_out_sim, rhoc, RHOC_SIZE * sizeof(Complex_DP));
    addusxx_g_cpu(rhoc_out_sim, xkq, xk, tau, becphi_c, becpsi_c,
                  nkb, ngms, nat, ntyp, upf_tvanp, nij_type,
                  ityp, ofsbeta, nh, ijtoh, qgm, eigts1, eigts2, eigts3,
                  mill, dfftt__nl);
    check_correctness(rhoc_out_sim, rhoc_out, RHOC_SIZE, "(1) CPU");

    cudaEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventCreate(&start_ev[i]);
        cudaEventCreate(&stop_ev[i]);
    }

    printf("Profiling for %d iterations, with %d warm-up iterations.\n", NUM_ITERS, NUM_WARMUP);
    memcpy(rhoc_out_sim, rhoc, RHOC_SIZE * sizeof(Complex_DP));
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventRecord(start_ev[i], 0);
        addusxx_g_cpu(rhoc_out_sim, xkq, xk, tau, becphi_c, becpsi_c,
                      nkb, ngms, nat, ntyp, upf_tvanp, nij_type,
                      ityp, ofsbeta, nh, ijtoh, qgm, eigts1, eigts2, eigts3,
                      mill, dfftt__nl);
        cudaEventRecord(stop_ev[i], 0);
        cudaEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        cudaEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }
    printf("Measured kernel time: %13.6f ms\n", total / NUM_ITERS);

    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventDestroy(start_ev[i]);
        cudaEventDestroy(stop_ev[i]);
    }
}

static void profile_gpu_baseline(DeviceArrays& da) {
    printf("\n>>> addusxx: GPU, baseline\n");

    // Correctness check
    reset_rhoc_device(da);
    addusxx_g_gpu(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                  da.d_becphi_c, da.d_becpsi_c,
                  nkb, ngms, nat, ntyp,
                  da.d_upf_tvanp, da.d_nij_type, da.d_ityp,
                  da.d_ofsbeta, da.d_nh, da.d_ijtoh,
                  da.d_qgm, da.d_eigts1, da.d_eigts2, da.d_eigts3,
                  da.d_mill, da.d_dfftt__nl,
                  upf_tvanp, nij_type, nh);
    cudaMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), cudaMemcpyDeviceToHost);
    check_correctness(rhoc_out_sim, rhoc_out, RHOC_SIZE, "GPU baseline");

    // Profile
    cudaEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventCreate(&start_ev[i]);
        cudaEventCreate(&stop_ev[i]);
    }

    printf("Profiling for %d iterations, with %d warm-up iterations.\n", NUM_ITERS, NUM_WARMUP);
    reset_rhoc_device(da);
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventRecord(start_ev[i], 0);
        addusxx_g_gpu(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                      da.d_becphi_c, da.d_becpsi_c,
                      nkb, ngms, nat, ntyp,
                      da.d_upf_tvanp, da.d_nij_type, da.d_ityp,
                      da.d_ofsbeta, da.d_nh, da.d_ijtoh,
                      da.d_qgm, da.d_eigts1, da.d_eigts2, da.d_eigts3,
                      da.d_mill, da.d_dfftt__nl,
                      upf_tvanp, nij_type, nh);
        cudaEventRecord(stop_ev[i], 0);
        cudaEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        cudaEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }
    printf("Measured kernel time: %13.6f ms\n", total / NUM_ITERS);

    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventDestroy(start_ev[i]);
        cudaEventDestroy(stop_ev[i]);
    }
}

static void profile_gpu_optimized(DeviceArrays& da) {
    printf("\n>>> addusxx: GPU, optimized\n");

    reset_rhoc_device(da);
    addusxx_g_gpu_optimized(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                            da.d_becphi_c, da.d_becpsi_c,
                            nkb, ngms, nat, ntyp,
                            da.d_upf_tvanp, da.d_nij_type, da.d_ityp,
                            da.d_ofsbeta, da.d_nh, da.d_ijtoh,
                            da.d_qgm_T, da.d_eigts1_T, da.d_eigts2_T, da.d_eigts3_T,
                            da.d_mill, da.d_dfftt__nl_sorted, da.d_dfftt__nl_ix,
                            upf_tvanp, nij_type, nh);
    cudaMemcpy(rhoc_out_sim, da.d_rhoc, RHOC_SIZE * sizeof(Complex_DP), cudaMemcpyDeviceToHost);
    check_correctness(rhoc_out_sim, rhoc_out, RHOC_SIZE, "GPU optimized");

    cudaEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventCreate(&start_ev[i]);
        cudaEventCreate(&stop_ev[i]);
    }

    printf("Profiling for %d iterations, with %d warm-up iterations.\n", NUM_ITERS, NUM_WARMUP);
    reset_rhoc_device(da);
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventRecord(start_ev[i], 0);
        addusxx_g_gpu_optimized(da.d_rhoc, da.d_xkq, da.d_xk, da.d_tau,
                                da.d_becphi_c, da.d_becpsi_c,
                                nkb, ngms, nat, ntyp,
                                da.d_upf_tvanp, da.d_nij_type, da.d_ityp,
                                da.d_ofsbeta, da.d_nh, da.d_ijtoh,
                                da.d_qgm_T, da.d_eigts1_T, da.d_eigts2_T, da.d_eigts3_T,
                                da.d_mill, da.d_dfftt__nl_sorted, da.d_dfftt__nl_ix,
                                upf_tvanp, nij_type, nh);
        cudaEventRecord(stop_ev[i], 0);
        cudaEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        cudaEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }
    printf("Measured kernel time: %13.6f ms\n", total / NUM_ITERS);

    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventDestroy(start_ev[i]);
        cudaEventDestroy(stop_ev[i]);
    }
}

static void profile_gpu_baseline_soa(DeviceArraysSoA& ds) {
    printf("\n>>> addusxx: GPU, baseline SoA\n");

    reset_rhoc_soa(ds);
    addusxx_g_gpu_soa(
        ds.d_rhoc_re, ds.d_rhoc_im,
        ds.d_xkq, ds.d_xk, ds.d_tau,
        ds.d_becphi_re, ds.d_becphi_im,
        ds.d_becpsi_re, ds.d_becpsi_im,
        nkb, ngms, nat, ntyp,
        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp,
        ds.d_ofsbeta, ds.d_nh, ds.d_ijtoh,
        ds.d_qgm_re, ds.d_qgm_im,
        ds.d_eigts1_re, ds.d_eigts1_im,
        ds.d_eigts2_re, ds.d_eigts2_im,
        ds.d_eigts3_re, ds.d_eigts3_im,
        ds.d_mill, ds.d_dfftt__nl,
        upf_tvanp, nij_type, nh);
    check_correctness_soa(ds, rhoc_out, RHOC_SIZE, "GPU baseline SoA");

    cudaEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventCreate(&start_ev[i]);
        cudaEventCreate(&stop_ev[i]);
    }

    printf("Profiling for %d iterations, with %d warm-up iterations.\n", NUM_ITERS, NUM_WARMUP);
    reset_rhoc_soa(ds);
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventRecord(start_ev[i], 0);
        addusxx_g_gpu_soa(
            ds.d_rhoc_re, ds.d_rhoc_im,
            ds.d_xkq, ds.d_xk, ds.d_tau,
            ds.d_becphi_re, ds.d_becphi_im,
            ds.d_becpsi_re, ds.d_becpsi_im,
            nkb, ngms, nat, ntyp,
            ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp,
            ds.d_ofsbeta, ds.d_nh, ds.d_ijtoh,
            ds.d_qgm_re, ds.d_qgm_im,
            ds.d_eigts1_re, ds.d_eigts1_im,
            ds.d_eigts2_re, ds.d_eigts2_im,
            ds.d_eigts3_re, ds.d_eigts3_im,
            ds.d_mill, ds.d_dfftt__nl,
            upf_tvanp, nij_type, nh);
        cudaEventRecord(stop_ev[i], 0);
        cudaEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        cudaEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }
    printf("Measured kernel time: %13.6f ms\n", total / NUM_ITERS);

    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventDestroy(start_ev[i]);
        cudaEventDestroy(stop_ev[i]);
    }
}

static void profile_gpu_optimized_soa(DeviceArraysSoA& ds) {
    printf("\n>>> addusxx: GPU, optimized SoA\n");

    reset_rhoc_soa(ds);
    addusxx_g_gpu_optimized_soa(
        ds.d_rhoc_re, ds.d_rhoc_im,
        ds.d_xkq, ds.d_xk, ds.d_tau,
        ds.d_becphi_re, ds.d_becphi_im,
        ds.d_becpsi_re, ds.d_becpsi_im,
        nkb, ngms, nat, ntyp,
        ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp,
        ds.d_ofsbeta, ds.d_nh, ds.d_ijtoh,
        ds.d_qgm_T_re, ds.d_qgm_T_im,
        ds.d_eigts1_T_re, ds.d_eigts1_T_im,
        ds.d_eigts2_T_re, ds.d_eigts2_T_im,
        ds.d_eigts3_T_re, ds.d_eigts3_T_im,
        ds.d_mill, ds.d_dfftt__nl_sorted, ds.d_dfftt__nl_ix,
        upf_tvanp, nij_type, nh);
    check_correctness_soa(ds, rhoc_out, RHOC_SIZE, "GPU optimized SoA");

    cudaEvent_t start_ev[TOTAL_ITERS], stop_ev[TOTAL_ITERS];
    float ms[TOTAL_ITERS];
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventCreate(&start_ev[i]);
        cudaEventCreate(&stop_ev[i]);
    }

    printf("Profiling for %d iterations, with %d warm-up iterations.\n", NUM_ITERS, NUM_WARMUP);
    reset_rhoc_soa(ds);
    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventRecord(start_ev[i], 0);
        addusxx_g_gpu_optimized_soa(
            ds.d_rhoc_re, ds.d_rhoc_im,
            ds.d_xkq, ds.d_xk, ds.d_tau,
            ds.d_becphi_re, ds.d_becphi_im,
            ds.d_becpsi_re, ds.d_becpsi_im,
            nkb, ngms, nat, ntyp,
            ds.d_upf_tvanp, ds.d_nij_type, ds.d_ityp,
            ds.d_ofsbeta, ds.d_nh, ds.d_ijtoh,
            ds.d_qgm_T_re, ds.d_qgm_T_im,
            ds.d_eigts1_T_re, ds.d_eigts1_T_im,
            ds.d_eigts2_T_re, ds.d_eigts2_T_im,
            ds.d_eigts3_T_re, ds.d_eigts3_T_im,
            ds.d_mill, ds.d_dfftt__nl_sorted, ds.d_dfftt__nl_ix,
            upf_tvanp, nij_type, nh);
        cudaEventRecord(stop_ev[i], 0);
        cudaEventSynchronize(stop_ev[i]);
    }

    float total = 0.0f;
    for (int i = NUM_WARMUP; i < TOTAL_ITERS; i++) {
        cudaEventElapsedTime(&ms[i], start_ev[i], stop_ev[i]);
        total += ms[i];
    }
    printf("Measured kernel time: %13.6f ms\n", total / NUM_ITERS);

    for (int i = 0; i < TOTAL_ITERS; i++) {
        cudaEventDestroy(start_ev[i]);
        cudaEventDestroy(stop_ev[i]);
    }
}

int main() {
    init_addusxx_data();

    // Allocate and copy all inputs to device (not measured)
    DeviceArrays da = allocate_device_arrays();
    DeviceArraysSoA ds = allocate_device_arrays_soa();

    profile_cpu_original();
    profile_gpu_baseline(da);
    profile_gpu_optimized(da);
    profile_gpu_baseline_soa(ds);
    profile_gpu_optimized_soa(ds);

    free_device_arrays(da);
    free_device_arrays_soa(ds);
    return 0;
}
