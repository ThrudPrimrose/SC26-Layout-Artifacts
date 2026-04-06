#include "usxx_kernels_cpu.h"
#include "usxx_kernels_cpu_impl.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <utility>

#include <omp.h>
#include <sys/mman.h>
#include <sys/stat.h>

// ============================================================
// Data loading
// ============================================================
static bool file_exists(const char* fn) { struct stat st; return stat(fn, &st) == 0; }
static long file_size_bytes(const char* fn) { struct stat st; return (stat(fn,&st)==0) ? st.st_size : -1; }

static bool data_load_int(int& v, const char* fn) {
    if (!file_exists(fn)) { printf("Missing: %s\n",fn); return false; }
    FILE*f=fopen(fn,"rb"); fread(&v,4,1,f); fclose(f); return true;
}
static bool data_load_int_array(int* d, int n, const char* fn) {
    if (!file_exists(fn)) { printf("Missing: %s\n",fn); return false; }
    long fs=file_size_bytes(fn);
    if (n != fs/4) { printf("Size mismatch %s: %d vs %ld\n",fn,n,fs); return false; }
    FILE*f=fopen(fn,"rb"); fread(d,4,n,f); fclose(f); return true;
}
static bool data_load_real_array(double* d, int n, const char* fn) {
    if (!file_exists(fn)) { printf("Missing: %s\n",fn); return false; }
    long fs=file_size_bytes(fn);
    if (n != fs/8) { printf("Size mismatch %s: %d vs %ld\n",fn,n,fs); return false; }
    FILE*f=fopen(fn,"rb"); fread(d,8,n,f); fclose(f); return true;
}
static bool data_load_cmplx_array(Complex_DP* d, int n, const char* fn) {
    if (!file_exists(fn)) { printf("Missing: %s\n",fn); return false; }
    long fs=file_size_bytes(fn);
    if (n != fs/16) { printf("Size mismatch %s: %d vs %ld\n",fn,n,fs); return false; }
    FILE*f=fopen(fn,"rb"); fread(d,16,n,f); fclose(f); return true;
}

// ============================================================
// NUMA helpers
// ============================================================
template<typename T>
static T* numa_alloc(size_t count) {
    size_t bytes = count * sizeof(T);
    void* p = mmap(nullptr, bytes, PROT_READ|PROT_WRITE,
                   MAP_PRIVATE|MAP_ANONYMOUS|MAP_NORESERVE, -1, 0);
    if (p == MAP_FAILED) { perror("mmap"); abort(); }
    madvise(p, bytes, MADV_HUGEPAGE);
    return static_cast<T*>(p);
}

template<typename T>
static void first_touch_zero(T* arr, size_t count) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count; i++) arr[i] = T{};
}

template<typename T>
static void first_touch_copy(T* dst, const T* src, size_t count) {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < count; i++) dst[i] = src[i];
}

// Cache flush
static constexpr size_t FLUSH_ELEMS = 256ULL * 1024 * 1024 / sizeof(double);
static double* g_flush = nullptr;

static void init_flush() {
    if (!g_flush) {
        g_flush = numa_alloc<double>(FLUSH_ELEMS);
        first_touch_zero(g_flush, FLUSH_ELEMS);
    }
}

static void flush_caches() {
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < FLUSH_ELEMS; i++) g_flush[i] = 1.0;
}

// ============================================================
// Array dimensions
// ============================================================
static constexpr int RHOC_SIZE   = 120000;
static constexpr int BECPC_SIZE  = 122;
static constexpr int QGM_ROWS    = 55191;
static constexpr int QGM_COLS    = 397;
static constexpr int EIGTS1_DIM1 = 217;
static constexpr int EIGTS1_DIM2 = 10;
static constexpr int EIGTS2_DIM1 = 109;
static constexpr int EIGTS2_DIM2 = 10;
static constexpr int NAT_MAX     = 10;
static constexpr int NTYP_MAX    = 3;
static constexpr int MILL_DIM1   = 3;
static constexpr int MILL_DIM2   = 156121;
static constexpr int NL_SIZE     = 55191;
static constexpr int IJTOH_N1    = 19;
static constexpr int IJTOH_N2    = 19;
static constexpr int IJTOH_N3    = 3;

static constexpr int NUM_ITERS   = 100;
static constexpr int NUM_WARMUP  = 5;

static int file_id = 12;

// ============================================================
// NUMA-placed arrays
// ============================================================
struct NUMAArrays {
    Complex_DP *rhoc, *rhoc_init, *rhoc_ref;
    Complex_DP *becphi_c, *becpsi_c;
    Complex_DP *qgm, *eigts1, *eigts2, *eigts3;
    Complex_DP *eigts1_T, *eigts2_T, *eigts3_T;
    double *rhoc_re, *rhoc_im, *rhoc_re_init, *rhoc_im_init;
    double *becphi_re, *becphi_im, *becpsi_re, *becpsi_im;
    double *qgm_re, *qgm_im;
    double *eigts1_re, *eigts1_im, *eigts2_re, *eigts2_im, *eigts3_re, *eigts3_im;
    double *eigts1_T_re, *eigts1_T_im, *eigts2_T_re, *eigts2_T_im, *eigts3_T_re, *eigts3_T_im;
    DP xkq[3], xk[3], tau[NAT_MAX * 3];
    int nkb, ngms, nat, ntyp;
    int nij_type[NTYP_MAX], ityp[NAT_MAX], ofsbeta[NAT_MAX], nh[NTYP_MAX];
    int upf_tvanp[NTYP_MAX];
    int *ijtoh, *mill, *dfftt__nl;
    int *dfftt__nl_sorted, *dfftt__nl_ix;
};

static void transpose_2d(const Complex_DP* A, Complex_DP* B, int m, int n) {
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            B[i * n + j] = A[j * m + i];
}

static NUMAArrays load_and_place() {
    NUMAArrays a;
    char fn[256];

    #define L_INT(v, name) snprintf(fn,256,"./bin/addusxx_g__" name "_%d.bin",file_id); data_load_int(v,fn);
    #define L_INTA(v,n,name) snprintf(fn,256,"./bin/addusxx_g__" name "_%d.bin",file_id); data_load_int_array(v,n,fn);
    #define L_REAL(v,n,name) snprintf(fn,256,"./bin/addusxx_g__" name "_%d.bin",file_id); data_load_real_array(v,n,fn);

    printf("Loading data...\n");
    L_INT(a.nkb,"nkb"); L_INT(a.ngms,"ngms"); L_INT(a.nat,"nat"); L_INT(a.ntyp,"ntyp");
    L_REAL(a.xkq,3,"xkq"); L_REAL(a.xk,3,"xk"); L_REAL(a.tau,3*NAT_MAX,"tau");
    L_INTA(a.nij_type,NTYP_MAX,"nij_type"); L_INTA(a.ityp,NAT_MAX,"ityp");
    L_INTA(a.ofsbeta,NAT_MAX,"ofsbeta"); L_INTA(a.nh,NTYP_MAX,"nh");
    L_INTA(a.upf_tvanp,NTYP_MAX,"upf_tvanp");

    int qgm_total = QGM_ROWS * QGM_COLS;
    Complex_DP* tmp_qgm = new Complex_DP[qgm_total];
    snprintf(fn,256,"./bin/addusxx_g__qgm_%d.bin",file_id);
    data_load_cmplx_array(tmp_qgm, qgm_total, fn);
    a.qgm = numa_alloc<Complex_DP>(qgm_total);
    first_touch_copy(a.qgm, tmp_qgm, qgm_total);
    a.qgm_re = numa_alloc<double>(qgm_total);
    a.qgm_im = numa_alloc<double>(qgm_total);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < qgm_total; i++) { a.qgm_re[i] = tmp_qgm[i].x; a.qgm_im[i] = tmp_qgm[i].y; }
    delete[] tmp_qgm;

    int e1_total = EIGTS1_DIM1 * EIGTS1_DIM2;
    int e2_total = EIGTS2_DIM1 * EIGTS2_DIM2;
    Complex_DP *tmp_e1 = new Complex_DP[e1_total];
    Complex_DP *tmp_e2 = new Complex_DP[e2_total];
    Complex_DP *tmp_e3 = new Complex_DP[e2_total];
    snprintf(fn,256,"./bin/addusxx_g__eigts1_%d.bin",file_id); data_load_cmplx_array(tmp_e1,e1_total,fn);
    snprintf(fn,256,"./bin/addusxx_g__eigts2_%d.bin",file_id); data_load_cmplx_array(tmp_e2,e2_total,fn);
    snprintf(fn,256,"./bin/addusxx_g__eigts3_%d.bin",file_id); data_load_cmplx_array(tmp_e3,e2_total,fn);

    a.eigts1 = numa_alloc<Complex_DP>(e1_total); first_touch_copy(a.eigts1, tmp_e1, e1_total);
    a.eigts2 = numa_alloc<Complex_DP>(e2_total); first_touch_copy(a.eigts2, tmp_e2, e2_total);
    a.eigts3 = numa_alloc<Complex_DP>(e2_total); first_touch_copy(a.eigts3, tmp_e3, e2_total);

    a.eigts1_T = numa_alloc<Complex_DP>(e1_total); transpose_2d(tmp_e1, a.eigts1_T, EIGTS1_DIM1, EIGTS1_DIM2);
    a.eigts2_T = numa_alloc<Complex_DP>(e2_total); transpose_2d(tmp_e2, a.eigts2_T, EIGTS2_DIM1, EIGTS2_DIM2);
    a.eigts3_T = numa_alloc<Complex_DP>(e2_total); transpose_2d(tmp_e3, a.eigts3_T, EIGTS2_DIM1, EIGTS2_DIM2);

    a.eigts1_re = numa_alloc<double>(e1_total); a.eigts1_im = numa_alloc<double>(e1_total);
    a.eigts2_re = numa_alloc<double>(e2_total); a.eigts2_im = numa_alloc<double>(e2_total);
    a.eigts3_re = numa_alloc<double>(e2_total); a.eigts3_im = numa_alloc<double>(e2_total);
    for (int i=0;i<e1_total;i++){a.eigts1_re[i]=tmp_e1[i].x; a.eigts1_im[i]=tmp_e1[i].y;}
    for (int i=0;i<e2_total;i++){a.eigts2_re[i]=tmp_e2[i].x; a.eigts2_im[i]=tmp_e2[i].y;}
    for (int i=0;i<e2_total;i++){a.eigts3_re[i]=tmp_e3[i].x; a.eigts3_im[i]=tmp_e3[i].y;}
    a.eigts1_T_re = numa_alloc<double>(e1_total); a.eigts1_T_im = numa_alloc<double>(e1_total);
    a.eigts2_T_re = numa_alloc<double>(e2_total); a.eigts2_T_im = numa_alloc<double>(e2_total);
    a.eigts3_T_re = numa_alloc<double>(e2_total); a.eigts3_T_im = numa_alloc<double>(e2_total);
    for (int i=0;i<e1_total;i++){a.eigts1_T_re[i]=a.eigts1_T[i].x; a.eigts1_T_im[i]=a.eigts1_T[i].y;}
    for (int i=0;i<e2_total;i++){a.eigts2_T_re[i]=a.eigts2_T[i].x; a.eigts2_T_im[i]=a.eigts2_T[i].y;}
    for (int i=0;i<e2_total;i++){a.eigts3_T_re[i]=a.eigts3_T[i].x; a.eigts3_T_im[i]=a.eigts3_T[i].y;}

    delete[] tmp_e1; delete[] tmp_e2; delete[] tmp_e3;

    Complex_DP tmp_bp[BECPC_SIZE], tmp_bs[BECPC_SIZE];
    snprintf(fn,256,"./bin/addusxx_g__becphi_c_%d.bin",file_id); data_load_cmplx_array(tmp_bp,BECPC_SIZE,fn);
    snprintf(fn,256,"./bin/addusxx_g__becpsi_c_%d.bin",file_id); data_load_cmplx_array(tmp_bs,BECPC_SIZE,fn);
    a.becphi_c = numa_alloc<Complex_DP>(BECPC_SIZE); memcpy(a.becphi_c, tmp_bp, BECPC_SIZE*sizeof(Complex_DP));
    a.becpsi_c = numa_alloc<Complex_DP>(BECPC_SIZE); memcpy(a.becpsi_c, tmp_bs, BECPC_SIZE*sizeof(Complex_DP));
    a.becphi_re = numa_alloc<double>(BECPC_SIZE); a.becphi_im = numa_alloc<double>(BECPC_SIZE);
    a.becpsi_re = numa_alloc<double>(BECPC_SIZE); a.becpsi_im = numa_alloc<double>(BECPC_SIZE);
    aos_to_soa(tmp_bp, a.becphi_re, a.becphi_im, BECPC_SIZE);
    aos_to_soa(tmp_bs, a.becpsi_re, a.becpsi_im, BECPC_SIZE);

    Complex_DP tmp_rhoc[RHOC_SIZE], tmp_rhoc_out[RHOC_SIZE];
    snprintf(fn,256,"./bin/addusxx_g__rhoc_%d.bin",file_id); data_load_cmplx_array(tmp_rhoc,RHOC_SIZE,fn);
    snprintf(fn,256,"./bin/addusxx_g__rhoc_out_%d.bin",file_id); data_load_cmplx_array(tmp_rhoc_out,RHOC_SIZE,fn);
    a.rhoc_init = numa_alloc<Complex_DP>(RHOC_SIZE); first_touch_copy(a.rhoc_init, tmp_rhoc, RHOC_SIZE);
    a.rhoc_ref  = numa_alloc<Complex_DP>(RHOC_SIZE); first_touch_copy(a.rhoc_ref, tmp_rhoc_out, RHOC_SIZE);
    a.rhoc      = numa_alloc<Complex_DP>(RHOC_SIZE);
    a.rhoc_re_init = numa_alloc<double>(RHOC_SIZE); a.rhoc_im_init = numa_alloc<double>(RHOC_SIZE);
    aos_to_soa(tmp_rhoc, a.rhoc_re_init, a.rhoc_im_init, RHOC_SIZE);
    a.rhoc_re = numa_alloc<double>(RHOC_SIZE);
    a.rhoc_im = numa_alloc<double>(RHOC_SIZE);

    int ijtoh_total = IJTOH_N1*IJTOH_N2*IJTOH_N3;
    int mill_total = MILL_DIM1*MILL_DIM2;
    int* tmp_ijtoh = new int[ijtoh_total]; int* tmp_mill = new int[mill_total]; int* tmp_nl = new int[NL_SIZE];
    snprintf(fn,256,"./bin/addusxx_g__ijtoh_%d.bin",file_id); data_load_int_array(tmp_ijtoh,ijtoh_total,fn);
    snprintf(fn,256,"./bin/addusxx_g__mill_%d.bin",file_id); data_load_int_array(tmp_mill,mill_total,fn);
    snprintf(fn,256,"./bin/addusxx_g__dfftt__nl_%d.bin",file_id); data_load_int_array(tmp_nl,NL_SIZE,fn);

    a.ijtoh = numa_alloc<int>(ijtoh_total); memcpy(a.ijtoh, tmp_ijtoh, ijtoh_total*sizeof(int));
    a.mill  = numa_alloc<int>(mill_total);  first_touch_copy(a.mill, tmp_mill, mill_total);
    a.dfftt__nl = numa_alloc<int>(NL_SIZE); first_touch_copy(a.dfftt__nl, tmp_nl, NL_SIZE);

    struct NLEntry { int val, orig; };
    NLEntry* entries = new NLEntry[a.ngms];
    for (int i=0; i<a.ngms; i++) { entries[i].val=tmp_nl[i]; entries[i].orig=i+1; }
    std::sort(entries, entries+a.ngms, [](const NLEntry&a, const NLEntry&b){return a.val<b.val;});
    a.dfftt__nl_sorted = numa_alloc<int>(NL_SIZE);
    a.dfftt__nl_ix     = numa_alloc<int>(NL_SIZE);
    for (int i=0; i<a.ngms; i++) { a.dfftt__nl_sorted[i]=entries[i].val; a.dfftt__nl_ix[i]=entries[i].orig; }
    delete[] entries;

    delete[] tmp_ijtoh; delete[] tmp_mill; delete[] tmp_nl;

    return a;
}

// ============================================================
// Correctness
// ============================================================
static bool check_aos(const Complex_DP* sim, const Complex_DP* ref, int n, const char* label) {
    double maxd = 0;
    for (int i=0; i<n; i++) { double d=cabs_val(csub(sim[i],ref[i])); if(d>maxd) maxd=d; }
    if (maxd >= 1e-8) { printf("  %s FAILED (maxdiff=%.4e)\n",label,maxd); return false; }
    printf("  %s PASSED\n",label); return true;
}

static bool check_soa(const double* re, const double* im, const Complex_DP* ref, int n, const char* label) {
    double maxd = 0;
    for (int i=0; i<n; i++) {
        double dr=re[i]-ref[i].x, di=im[i]-ref[i].y;
        double d=sqrt(dr*dr+di*di); if(d>maxd) maxd=d;
    }
    if (maxd >= 1e-8) { printf("  %s FAILED (maxdiff=%.4e)\n",label,maxd); return false; }
    printf("  %s PASSED\n",label); return true;
}

// ============================================================
// Profiling
// ============================================================
static void profile(FILE* csv, const char* variant, int blocksize, NUMAArrays& a,
                    std::function<void()> reset_fn,
                    std::function<void()> kernel_fn,
                    std::function<bool()> verify_fn,
                    int nthreads)
{
    omp_set_num_threads(nthreads);
    reset_fn(); kernel_fn();
    if (!verify_fn()) return;
    for (int i = 0; i < NUM_WARMUP; i++) { reset_fn(); kernel_fn(); }
    for (int r = 0; r < NUM_ITERS; r++) {
        reset_fn();
        flush_caches();
        double t0 = omp_get_wtime();
        kernel_fn();
        double t1 = omp_get_wtime();
        double ms = (t1 - t0) * 1000.0;
        fprintf(csv, "%s,%d,%d,%d,%d,%d,%.6f\n", variant, blocksize, a.ngms, RHOC_SIZE, nthreads, r, ms);
    }
}

// ============================================================
// Template dispatch: run all 6 kernels for a given BLOCKSIZE
// ============================================================
template<int BS>
static void run_all_for_blocksize(FILE* csv, NUMAArrays& a, int nt) {
    char label[64];

    snprintf(label, sizeof(label), "baseline_aos");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc, a.rhoc_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_baseline_aos<BS>(a.rhoc, a.xkq, a.xk, a.tau,
                a.becphi_c, a.becpsi_c, a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm, a.eigts1, a.eigts2, a.eigts3, a.mill, a.dfftt__nl); },
        [&]{ return check_aos(a.rhoc, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    snprintf(label, sizeof(label), "eigts_t_aos");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc, a.rhoc_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_eigts_transposed_aos<BS>(a.rhoc, a.xkq, a.xk, a.tau,
                a.becphi_c, a.becpsi_c, a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm, a.eigts1_T, a.eigts2_T, a.eigts3_T, a.mill, a.dfftt__nl); },
        [&]{ return check_aos(a.rhoc, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    snprintf(label, sizeof(label), "sorted_aos");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc, a.rhoc_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_sorted_aos<BS>(a.rhoc, a.xkq, a.xk, a.tau,
                a.becphi_c, a.becpsi_c, a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm, a.eigts1_T, a.eigts2_T, a.eigts3_T, a.mill,
                a.dfftt__nl_sorted, a.dfftt__nl_ix); },
        [&]{ return check_aos(a.rhoc, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    snprintf(label, sizeof(label), "baseline_soa");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc_re, a.rhoc_re_init, RHOC_SIZE);
             first_touch_copy(a.rhoc_im, a.rhoc_im_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_baseline_soa<BS>(a.rhoc_re, a.rhoc_im, a.xkq, a.xk, a.tau,
                a.becphi_re, a.becphi_im, a.becpsi_re, a.becpsi_im,
                a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm_re, a.qgm_im,
                a.eigts1_re, a.eigts1_im, a.eigts2_re, a.eigts2_im,
                a.eigts3_re, a.eigts3_im, a.mill, a.dfftt__nl); },
        [&]{ return check_soa(a.rhoc_re, a.rhoc_im, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    snprintf(label, sizeof(label), "eigts_t_soa");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc_re, a.rhoc_re_init, RHOC_SIZE);
             first_touch_copy(a.rhoc_im, a.rhoc_im_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_eigts_transposed_soa<BS>(a.rhoc_re, a.rhoc_im, a.xkq, a.xk, a.tau,
                a.becphi_re, a.becphi_im, a.becpsi_re, a.becpsi_im,
                a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm_re, a.qgm_im,
                a.eigts1_T_re, a.eigts1_T_im, a.eigts2_T_re, a.eigts2_T_im,
                a.eigts3_T_re, a.eigts3_T_im, a.mill, a.dfftt__nl); },
        [&]{ return check_soa(a.rhoc_re, a.rhoc_im, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    snprintf(label, sizeof(label), "sorted_soa");
    profile(csv, label, BS, a,
        [&]{ first_touch_copy(a.rhoc_re, a.rhoc_re_init, RHOC_SIZE);
             first_touch_copy(a.rhoc_im, a.rhoc_im_init, RHOC_SIZE); },
        [&]{ cpu_addusxx_sorted_soa<BS>(a.rhoc_re, a.rhoc_im, a.xkq, a.xk, a.tau,
                a.becphi_re, a.becphi_im, a.becpsi_re, a.becpsi_im,
                a.nkb, a.ngms, a.nat, a.ntyp,
                a.upf_tvanp, a.nij_type, a.ityp, a.ofsbeta, a.nh, a.ijtoh,
                a.qgm_re, a.qgm_im,
                a.eigts1_T_re, a.eigts1_T_im, a.eigts2_T_re, a.eigts2_T_im,
                a.eigts3_T_re, a.eigts3_T_im, a.mill,
                a.dfftt__nl_sorted, a.dfftt__nl_ix); },
        [&]{ return check_soa(a.rhoc_re, a.rhoc_im, a.rhoc_ref, RHOC_SIZE, label); },
        nt);

    fflush(csv);
}

// ============================================================
// Compile-time iteration over powers of 2
// ============================================================
template<int BS>
struct BlockSweep {
    static void run(FILE* csv, NUMAArrays& a, int nt) {
        BlockSweep<BS / 2>::run(csv, a, nt);   // recurse smaller first
        printf("\n--- blocksize=%d ---\n", BS);
        run_all_for_blocksize<BS>(csv, a, nt);
    }
};

template<>
struct BlockSweep<1> {
    static void run(FILE* csv, NUMAArrays& a, int nt) {
        printf("\n--- blocksize=1 ---\n");
        run_all_for_blocksize<1>(csv, a, nt);
    }
};

// ============================================================
// Main
// ============================================================
int main() {
    init_flush();
    NUMAArrays a = load_and_place();

    int max_threads = omp_get_max_threads();
    std::vector<int> thread_counts;
    if (max_threads >= 96) thread_counts = {96};
    else if (max_threads >= 72) thread_counts = {288};
    else thread_counts = {std::max(1,max_threads/2), max_threads};

    printf("Threads: "); for (int t : thread_counts) printf("%d ",t); printf("\n");

    FILE* csv = fopen("addusxx_cpu_sweep.csv", "w");
    fprintf(csv, "variant,blocksize,ngms,rhoc_size,nthreads,rep,time_ms\n");

    for (int nt : thread_counts) {
        printf("\n=== nthreads=%d ===\n", nt);
        BlockSweep<256>::run(csv, a, nt);
    }

    fclose(csv);
    printf("\nDone. CSV: addusxx_cpu_sweep.csv\n");
    return 0;
}