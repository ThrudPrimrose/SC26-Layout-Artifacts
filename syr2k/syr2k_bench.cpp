// syr2k_bench.cpp — layout/schedule/blocking/tiling explorer (OpenMP integrated)
// Build: g++ -O3 -fopenmp -std=c++17 -DSZ_N=1024 -DSZ_M=1024 -o syr2k_bench syr2k_bench.cpp -lopenblas
// Run:   ./syr2k_bench <mode:1-5> [output.csv]
//  1 = layout permutations (8 variants)
//  2 = loop orderings (6 variants, now parallel where safe)
//  3 = blocked arrays + blocked loops (parallel on bi loop)
//  4 = tiled loops (parallel on tile loops)
//  5 = OpenBLAS syr2k (alpha = 1.0)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <omp.h>
#include <cblas.h>

#ifndef SZ_N
#define SZ_N 1024
#endif
#ifndef SZ_M
#define SZ_M 1024
#endif
#ifndef TI
#define TI 32
#endif
#ifndef TJ
#define TJ 32
#endif
#ifndef TK
#define TK 32
#endif
#ifndef BI
#define BI 32
#endif
#ifndef BJ
#define BJ 32
#endif
#ifndef BK
#define BK 32
#endif
#ifndef NRUNS
#define NRUNS 20
#endif
#ifndef NWARM
#define NWARM 5
#endif

static constexpr int N = SZ_N, M = SZ_M;
static constexpr double alpha = 1.5, beta_val = 1.2;
using Clock = std::chrono::high_resolution_clock;

// Fast parallel init
static void init_rand(double *a, long n) {
    #pragma omp parallel
    {
        uint32_t s = 42u + (uint32_t)omp_get_thread_num() * 2654435761u;
        #pragma omp for schedule(static)
        for (long i = 0; i < n; i++) {
            s ^= s<<13; s ^= s>>17; s ^= s<<5;
            a[i] = (s & 0xFFFFu) * (1.0/65536.0);
        }
    }
}

// Layout helpers (unchanged)
enum Lay { R, C };
template<Lay L> inline int ix(int r, int c, int nr, int nc) {
    return L==R ? r*nc+c : c*nr+r;
}
static void to_rm(const double *src, double *dst, int nr, int nc, Lay lay) {
    if (lay == R) { memcpy(dst, src, (size_t)nr*nc*sizeof(double)); return; }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            dst[i*nc+j] = src[j*nr+i];
}
static void to_lay(const double *rm, double *dst, int nr, int nc, Lay lay) {
    if (lay == R) { memcpy(dst, rm, (size_t)nr*nc*sizeof(double)); return; }
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            dst[j*nr+i] = rm[i*nc+j];
}

// ── MODE 1: layout-permutation kernels (serial) ─────────────────────────
template<Lay LA, Lay LB, Lay LC>
static void kern_lay(double *Cd, const double *Ad, const double *Bd) {
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < M; k++)
            for (int j = 0; j <= i; j++)
                Cd[ix<LC>(i,j,N,N)] += alpha*Ad[ix<LA>(i,k,N,M)]*Bd[ix<LB>(j,k,N,M)]
                                      + alpha*Bd[ix<LB>(i,k,N,M)]*Ad[ix<LA>(j,k,N,M)];
    }
}

// Verification
static double verify(const double *got, const double *ref, int n) {
    double mx = 0;
    for (int i = 0; i < n; i++)
        mx = std::max(mx, std::abs(got[i]-ref[i])/(1.0+std::abs(ref[i])));
    return mx;
}

template<Lay LA, Lay LB, Lay LC>
static void run_layout(const double *A, const double *B, const double *C0,
                       const double *Cref, FILE *fp, const char *tag) {
    double *Ac = new double[N*M], *Bc = new double[N*M];
    double *C0c = new double[N*N], *Cw = new double[N*N], *Cback = new double[N*N];
    to_lay(A, Ac, N, M, LA);
    to_lay(B, Bc, N, M, LB);
    to_lay(C0, C0c, N, N, LC);

    for (int r = 0; r < NRUNS; r++) {
        memcpy(Cw, C0c, (size_t)N*N*sizeof(double));
        auto t0 = Clock::now();
        kern_lay<LA,LB,LC>(Cw, Ac, Bc);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,%d,1,%s,0,0,0,0,0,0,%d,%.2f\n",
                    N,M,tag,r-NWARM,
                    std::chrono::duration<double,std::micro>(t1-t0).count());
    }
    to_rm(Cw, Cback, N, N, LC);
    double e = verify(Cback, Cref, N*N);
    fprintf(stderr, "%s %s: maxerr=%.2e\n", e>1e-6?"FAIL":"OK", tag, e);
    delete[] Ac; delete[] Bc; delete[] C0c; delete[] Cw; delete[] Cback;
}

// ── MODE 2: loop-order kernels (row-major, OpenMP where safe) ───────────
#define RM(a,r,c,cols) a[(r)*(cols)+(c)]

// Safe: parallel on i (each row updated by one thread)
static void kern_ikj(double *C, const double *A, const double *B) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < M; k++)
            for (int j = 0; j <= i; j++)
                RM(C,i,j,N) += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                             + alpha * RM(B,i,k,M) * RM(A,j,k,M);
    }
}

// Safe: parallel on i (row updates independent)
static void kern_ijk(double *C, const double *A, const double *B) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            double sum = 0.0;
            for (int k = 0; k < M; k++)
                sum += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                     + alpha * RM(B,i,k,M) * RM(A,j,k,M);
            RM(C,i,j,N) += sum;
        }
    }
}

// Safe: parallel on j (each column block updated by one thread)
static void kern_jik(double *C, const double *A, const double *B) {
    #pragma omp parallel for
    for (int j = 0; j < N; j++) {
        for (int i = j; i < N; i++) {
            double sum = 0.0;
            for (int k = 0; k < M; k++)
                sum += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                     + alpha * RM(B,i,k,M) * RM(A,j,k,M);
            RM(C,i,j,N) += sum;
        }
    }
}

// Unsafe (races if parallelized naively) – left serial
static void kern_kij(double *C, const double *A, const double *B) {
    for (int k = 0; k < M; k++)
        #pragma omp parallel for
        for (int i = 0; i < N; i++)
            for (int j = 0; j <= i; j++)
                RM(C,i,j,N) += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                             + alpha * RM(B,i,k,M) * RM(A,j,k,M);
}

static void kern_kji(double *C, const double *A, const double *B) {
    for (int k = 0; k < M; k++)
        #pragma omp parallel for
        for (int j = 0; j < N; j++)
            for (int i = j; i < N; i++)
                RM(C,i,j,N) += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                             + alpha * RM(B,i,k,M) * RM(A,j,k,M);
}

static void kern_jki(double *C, const double *A, const double *B) {
    #pragma omp parallel for 
    for (int j = 0; j < N; j++)
        for (int k = 0; k < M; k++)
            for (int i = j; i < N; i++)
                RM(C,i,j,N) += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                             + alpha * RM(B,i,k,M) * RM(A,j,k,M);
}


using KernFn = void(*)(double*, const double*, const double*);

static void run_loop(KernFn fn, const double *A, const double *B, const double *C0,
                     const double *Cref, FILE *fp, const char *tag) {
    double *Cw = new double[N*N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(Cw, C0, (size_t)N*N*sizeof(double));
        auto t0 = Clock::now();
        fn(Cw, A, B);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,%d,2,%s,0,0,0,0,0,0,%d,%.2f\n",
                    N,M,tag,r-NWARM,
                    std::chrono::duration<double,std::micro>(t1-t0).count());
    }
    double e = verify(Cw, Cref, N*N);
    fprintf(stderr, "%s %s: maxerr=%.2e\n", e>1e-6?"FAIL":"OK", tag, e);
    delete[] Cw;
}

// ── MODE 3: blocked arrays + blocked loops (parallel on bi) ─────────────
inline int bidx(int r, int c, int nr, int nc, int br, int bc) {
    int nbr = (nr+br-1)/br, nbc = (nc+bc-1)/bc;
    (void)nbr;
    return ((r/br)*nbc + c/bc)*br*bc + (r%br)*bc + (c%bc);
}
inline int bsz(int nr, int nc, int br, int bc) {
    return ((nr+br-1)/br)*((nc+bc-1)/bc)*br*bc;
}

static void to_blk(const double *rm, double *blk, int nr, int nc, int br, int bc) {
    memset(blk, 0, bsz(nr,nc,br,bc)*sizeof(double));
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            blk[bidx(i,j,nr,nc,br,bc)] = rm[i*nc+j];
}
static void from_blk(const double *blk, double *rm, int nr, int nc, int br, int bc) {
    for (int i = 0; i < nr; i++)
        for (int j = 0; j < nc; j++)
            rm[i*nc+j] = blk[bidx(i,j,nr,nc,br,bc)];
}

static void kern_blocked(double *Cd, const double *Ad, const double *Bd) {
    #pragma omp parallel for collapse(2)
    for (int bi = 0; bi < N; bi += BI) {
        for (int bk = 0; bk < M; bk += BK) {
            for (int bj = 0; bj < std::min(N, bi + BI); bj += BJ) {
                int ie = std::min(bi + BI, N);
                int ke = std::min(bk + BK, M);
                int je = std::min(bj + BJ, N);
                for (int i = bi; i < ie; i++)
                    for (int k = bk; k < ke; k++)
                        for (int j = bj; j < std::min(je, i + 1); j++)
                            Cd[bidx(i,j,N,N,BI,BJ)] +=
                                alpha * Ad[bidx(i,k,N,M,BI,BK)] * Bd[bidx(j,k,N,M,BI,BK)]
                              + alpha * Bd[bidx(i,k,N,M,BI,BK)] * Ad[bidx(j,k,N,M,BI,BK)];
            }
        }
    }
}

static void run_blocked(const double *A, const double *B, const double *C0,
                        const double *Cref, FILE *fp) {
    double *Ab = new double[bsz(N,M,BI,BK)];
    double *Bb = new double[bsz(N,M,BI,BK)];
    double *C0b = new double[bsz(N,N,BI,BJ)];
    double *Cw  = new double[bsz(N,N,BI,BJ)];
    double *Cback = new double[N*N];
    to_blk(A,Ab,N,M,BI,BK); to_blk(B,Bb,N,M,BI,BK); to_blk(C0,C0b,N,N,BI,BJ);

    for (int r = 0; r < NRUNS; r++) {
        memcpy(Cw, C0b, bsz(N,N,BI,BJ)*sizeof(double));
        auto t0 = Clock::now();
        kern_blocked(Cw, Ab, Bb);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,%d,3,blocked,%d,%d,%d,0,0,0,%d,%.2f\n",
                    N,M,BI,BJ,BK,r-NWARM,
                    std::chrono::duration<double,std::micro>(t1-t0).count());
    }
    from_blk(Cw, Cback, N, N, BI, BJ);
    double e = verify(Cback, Cref, N*N);
    fprintf(stderr, "%s blocked(%d,%d,%d): maxerr=%.2e\n", e>1e-6?"FAIL":"OK",BI,BJ,BK,e);
    delete[] Ab; delete[] Bb; delete[] C0b; delete[] Cw; delete[] Cback;
}

// ── MODE 4: tiled loops (row-major, parallel on tile loops) ─────────────
static void kern_tiled(double *C, const double *A, const double *B) {
    #pragma omp parallel for collapse(2)
    for (int ii = 0; ii < N; ii += TI) {
        for (int jj = 0; jj < N; jj += TJ) {
            int ie = std::min(ii + TI, N);
            int je = std::min(jj + TJ, N);
            for (int kk = 0; kk < M; kk += TK) {
                int ke = std::min(kk + TK, M);
                for (int i = ii; i < ie; i++) {
                    for (int k = kk; k < ke; k++) {
                        for (int j = jj; j < std::min(je, i + 1); j++) {
                            RM(C,i,j,N) += alpha * RM(A,i,k,M) * RM(B,j,k,M)
                                         + alpha * RM(B,i,k,M) * RM(A,j,k,M);
                        }
                    }
                }
            }
        }
    }
}

static void run_tiled(const double *A, const double *B, const double *C0,
                      const double *Cref, FILE *fp) {
    double *Cw = new double[N*N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(Cw, C0, (size_t)N*N*sizeof(double));
        auto t0 = Clock::now();
        kern_tiled(Cw, A, B);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,%d,4,tiled,0,0,0,%d,%d,%d,%d,%.2f\n",
                    N,M,TI,TJ,TK,r-NWARM,
                    std::chrono::duration<double,std::micro>(t1-t0).count());
    }
    double e = verify(Cw, Cref, N*N);
    fprintf(stderr, "%s tiled(%d,%d,%d): maxerr=%.2e\n", e>1e-6?"FAIL":"OK",TI,TJ,TK,e);
    delete[] Cw;
}

// ── MODE 5: OpenBLAS syr2k (alpha = 1.0) ────────────────────────────────
static void kern_openblas(double *C, const double *A, const double *B) {
    cblas_dsyr2k(CblasRowMajor, CblasLower, CblasNoTrans,
                 N, M, 1.0, A, M, B, M, beta_val, C, N);
}

static void run_openblas(const double *A, const double *B, const double *C0,
                          const double *Cref, FILE *fp) {
    double *Cw = new double[N*N];
    for (int r = 0; r < NRUNS; r++) {
        memcpy(Cw, C0, (size_t)N*N*sizeof(double));
        auto t0 = Clock::now();
        kern_openblas(Cw, A, B);
        auto t1 = Clock::now();
        if (r >= NWARM)
            fprintf(fp, "%d,%d,5,openblas,0,0,0,0,0,0,%d,%.2f\n",
                    N,M,r-NWARM,
                    std::chrono::duration<double,std::micro>(t1-t0).count());
    }
    double e = verify(Cw, Cref, N*N);
    fprintf(stderr, "%s openblas: maxerr=%.2e\n", e>1e-6?"FAIL":"OK", e);
    delete[] Cw;
}

// ── main ────────────────────────────────────────────────────────────────
int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <mode:1-5> [output.csv]\n", argv[0]);
        return 1;
    }
    int mode = atoi(argv[1]);
    const char *csv = argc > 2 ? argv[2] : "results.csv";

    double *A = new double[N*M], *B = new double[N*M];
    double *C0 = new double[N*N], *Cref = new double[N*N];
    init_rand(A, (long)N*M);
    init_rand(B, (long)N*M);
    init_rand(C0, (long)N*N);

    // Reference (serial ikj)
    fprintf(stderr, "computing reference (N=%d, M=%d)...\n", N, M);
    memcpy(Cref, C0, (size_t)N*N*sizeof(double));
    // Temporarily disable OpenMP for reference to keep it serial
    int saved_threads = omp_get_max_threads();
    omp_set_num_threads(1);
    kern_ikj(Cref, A, B);
    omp_set_num_threads(saved_threads);

    FILE *fp = fopen(csv, "a");
    fseek(fp, 0, SEEK_END);
    if (ftell(fp) == 0)
        fprintf(fp, "n,m,mode,variant,bi,bj,bk,ti,tj,tk,run,time_us\n");

    switch (mode) {
    case 1:
        fprintf(stderr, "mode 1: layout permutations\n");
        run_layout<R,R,R>(A,B,C0,Cref,fp,"RRR");
        run_layout<R,R,C>(A,B,C0,Cref,fp,"RRC");
        run_layout<R,C,R>(A,B,C0,Cref,fp,"RCR");
        run_layout<R,C,C>(A,B,C0,Cref,fp,"RCC");
        run_layout<C,R,R>(A,B,C0,Cref,fp,"CRR");
        run_layout<C,R,C>(A,B,C0,Cref,fp,"CRC");
        run_layout<C,C,R>(A,B,C0,Cref,fp,"CCR");
        run_layout<C,C,C>(A,B,C0,Cref,fp,"CCC");
        break;
    case 2:
        fprintf(stderr, "mode 2: loop orderings (parallel where safe)\n");
        run_loop(kern_ikj, A,B,C0,Cref,fp,"ikj");
        run_loop(kern_ijk, A,B,C0,Cref,fp,"ijk");
        run_loop(kern_kij, A,B,C0,Cref,fp,"kij");   // serial
        run_loop(kern_kji, A,B,C0,Cref,fp,"kji");   // serial
        run_loop(kern_jik, A,B,C0,Cref,fp,"jik");
        run_loop(kern_jki, A,B,C0,Cref,fp,"jki");   // serial
        break;
    case 3:
        fprintf(stderr, "mode 3: blocked arrays BI=%d BJ=%d BK=%d\n", BI,BJ,BK);
        run_blocked(A,B,C0,Cref,fp);
        break;
    case 4:
        fprintf(stderr, "mode 4: tiled loops TI=%d TJ=%d TK=%d\n", TI,TJ,TK);
        run_tiled(A,B,C0,Cref,fp);
        break;
    case 5:
        fprintf(stderr, "mode 5: OpenBLAS SYR2K (alpha=1.0)\n");
        run_openblas(A,B,C0,Cref,fp);
        break;
    default:
        fprintf(stderr, "unknown mode %d\n", mode);
        fclose(fp); return 1;
    }

    fclose(fp);
    delete[] A; delete[] B; delete[] C0; delete[] Cref;
    fprintf(stderr, "done. results appended to %s\n", csv);
}