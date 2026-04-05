/*
 * bench_cpu.cpp -- NUMA-aware CPU z_v_grad_w stencil benchmark
 *
 * Arrays have separate horizontal dimensions:
 *   N_e = edges (direct arrays), N_c = cells (w target), N_v = verts (z_w_v target)
 */
#include "bench_common.h"
#include "icon_data_loader.h"
#include <ctime>
#include <omp.h>
#include <utility>
#include <cassert>

static SchedKind sched_for_par_for(int V) {
  int kV = kern_v(V);
  return (kV <= 2) ? SCHED_JK_OUTER : SCHED_JE_OUTER;
}

/* ================================================================ */
/*  Kernel signature: N_e, N_c, N_v, nlev, nlev_end                  */
/* ================================================================ */
typedef void (*kern_t)(double *, const double *, const double *, const double *,
                       const int *, const double *, const double *,
                       const double *, const double *, const int *,
                       int, int, int, int, int);

/* ---- Unblocked kernels ---- */
template <int V>
static void cpu_par_for(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  if constexpr (V <= 2) {
#pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev_end; jk++)
      for (int je = 0; je < N_e; je++) { STENCIL_BODY(V) }
  } else {
#pragma omp parallel for schedule(static)
    for (int je = 0; je < N_e; je++)
      for (int jk = 0; jk < nlev_end; jk++) { STENCIL_BODY(V) }
  }
}

template <int V>
static void cpu_collapse2(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  if constexpr (V <= 2) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int jk = 0; jk < nlev_end; jk++)
      for (int je = 0; je < N_e; je++) { STENCIL_BODY(V) }
  } else {
#pragma omp parallel for collapse(2) schedule(static)
    for (int je = 0; je < N_e; je++)
      for (int jk = 0; jk < nlev_end; jk++) { STENCIL_BODY(V) }
  }
}

template <int V>
static void cpu_collapse2_je(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
#pragma omp parallel for collapse(2) schedule(static)
  for (int je = 0; je < N_e; je++)
    for (int jk = 0; jk < nlev_end; jk++) { STENCIL_BODY(V) }
}

template <int V>
static void cpu_numa4(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  constexpr int ND = NUMA_DOMAINS;
  int tpd = std::max(1, omp_get_max_threads() / ND);
  int chunk = N_e / ND;
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int d = std::min(tid / tpd, ND - 1);
    int je0 = d * chunk;
    int je1 = (d == ND-1) ? N_e : je0 + chunk;
    int ltid = tid - d * tpd;
    int ln = (d == ND-1) ? (omp_get_num_threads() - d*tpd) : tpd;
    int my0 = je0 + (int)((long long)(je1-je0) * ltid / ln);
    int my1 = je0 + (int)((long long)(je1-je0) * (ltid+1) / ln);
    for (int je = my0; je < my1; je++)
      for (int jk = 0; jk < nlev_end; jk++) { STENCIL_BODY(V) }
  }
}

static kern_t par_tbl[]    = {cpu_par_for<1>,      cpu_par_for<2>,      cpu_par_for<3>,      cpu_par_for<4>};
static kern_t col_tbl[]    = {cpu_collapse2<1>,    cpu_collapse2<2>,    cpu_collapse2<3>,    cpu_collapse2<4>};
static kern_t col_je_tbl[] = {cpu_collapse2_je<1>, cpu_collapse2_je<2>, cpu_collapse2_je<3>, cpu_collapse2_je<4>};
static kern_t numa4_tbl[]  = {cpu_numa4<1>,        cpu_numa4<2>,        cpu_numa4<3>,        cpu_numa4<4>};

/* ---- Blocked kernels ---- */
template <int B>
static void cpu_blocked_for(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  int nblocks = N_e / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++)
    for (int jk = 0; jk < nlev_end; jk++) {
#pragma omp simd
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
    }
}
template <int B>
static void cpu_blocked_col(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  int nblocks = N_e / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++)
    for (int jk = 0; jk < nlev_end; jk++) {
#pragma omp simd
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
    }
}
static kern_t bfor_tbl[] = {cpu_blocked_for<8>, cpu_blocked_for<16>,
                            cpu_blocked_for<32>, cpu_blocked_for<64>,
                            cpu_blocked_for<128>};
static kern_t bcol_tbl[] = {cpu_blocked_col<8>, cpu_blocked_col<16>,
                            cpu_blocked_col<32>, cpu_blocked_col<64>,
                            cpu_blocked_col<128>};

/* ---- Reference (serial) ---- */
template <int V>
static void cpu_ref(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  for (int jk = 0; jk < nlev_end; jk++)
    for (int je = 0; je < N_e; je++) { STENCIL_BODY(V) }
}
static void cpu_ref_v(int V, double *o, const double *vn, const double *id,
                      const double *w, const int *ci, const double *vt,
                      const double *ip, const double *tg, const double *zw,
                      const int *vi, int Ne, int Nc, int Nv, int nl, int nle) {
  int kV = kern_v(V);
  switch (kV) {
  case 1: cpu_ref<1>(o,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle); break;
  case 2: cpu_ref<2>(o,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle); break;
  case 3: cpu_ref<3>(o,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle); break;
  case 4: cpu_ref<4>(o,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle); break;
  }
}
template <int B>
static void cpu_ref_blocked(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end) {
  for (int jb = 0; jb < N_e / B; jb++)
    for (int jk = 0; jk < nlev_end; jk++)
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
}
static kern_t bref_tbl[] = {cpu_ref_blocked<8>, cpu_ref_blocked<16>,
                            cpu_ref_blocked<32>, cpu_ref_blocked<64>,
                            cpu_ref_blocked<128>};

/* ---- Verify ---- */
static bool verify(const double *got, const double *ref, size_t n,
                   int *nf, double *mr) {
  *nf=0; *mr=0;
  for (size_t i=0;i<n;i++) {
    double d=std::abs(got[i]-ref[i]),
           dn=std::max(std::abs(ref[i]),1e-300), r=d/dn;
    if (r>*mr) *mr=r;
    if (d > 1e-12+1e-8*std::abs(ref[i])) (*nf)++;
  }
  return *nf==0;
}

/* ---- Flush ---- */
static constexpr int FN=8192*2, FS=3;
static double *fb0=nullptr, *fb1=nullptr;
static void flush() {
  static bool init=false;
  if (!init) {
    size_t n=(size_t)FN*FN;
    fb0=numa_alloc_unfaulted<double>(n); fb1=numa_alloc_unfaulted<double>(n);
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<n;i++) {
      uint64_t h=splitmix64(12345ULL+i);
      fb0[i]=(double)(h>>11)/(double)(1ULL<<53); fb1[i]=fb0[i];
    }
    init=true;
  }
  double *A=fb0, *B=fb1;
  for (int s=0;s<FS;s++) {
#pragma omp parallel for schedule(static)
    for (int i=1;i<FN-1;i++)
      for (int j=1;j<FN-1;j++)
        B[i*FN+j]=0.25*(A[(i-1)*FN+j]+A[(i+1)*FN+j]+A[i*FN+(j-1)]+A[i*FN+(j+1)]);
    std::swap(A,B);
  }
}

/* ---- helper: benchmark one kernel variant ---- */
static void bench_kernel(FILE *f, kern_t kern, const char *par_label,
                         int V, BenchData &bd, int nlev_end,
                         const char *dl, double *hr) {
  int Ne=bd.N_e, Nc=bd.N_c, Nv=bd.N_v, nl=bd.nlev;
  for (int r=0; r<WARMUP; r++) {
    flush();
    kern(bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
         bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
         bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
  }
  flush();
  { int nf; double mr;
    verify(bd.h_out, hr, bd.sz_e, &nf, &mr);
    printf("V%d %-12s %-18s %s mr=%.2e\n", V, dl, par_label, nf?"FAIL":"OK", mr); }
  for (int r=0; r<NRUNS; r++) {
    flush();
    auto t0=std::chrono::high_resolution_clock::now();
    kern(bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
         bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
         bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
    auto t1=std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,%d,%d,%d,%d,%d,%s,0,%s,%d,%.9f\n",
            V, nlev_end, Ne, Nc, Nv, dl,
            par_label, r, std::chrono::duration<double,std::milli>(t1-t0).count());
    flush();
  }
}

/* ---- Run unblocked ---- */
static void run_unblocked(FILE *f, int V, BenchData &bd, int nlev_end,
                          const char *dl, double *hr) {
  int Vi = kern_v(V) - 1;
  int Ne=bd.N_e, Nc=bd.N_c, Nv=bd.N_v, nl=bd.nlev;

  memset(hr, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(sched_for_par_for(V));
  cpu_ref_v(V, hr, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx, bd.h_z_vt_ie,
            bd.inv_primal, bd.tangent_o, bd.h_z_w_v, bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
  flush();

  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(sched_for_par_for(V));
  bench_kernel(f, par_tbl[Vi], "omp_for", V, bd, nlev_end, dl, hr);

  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(SCHED_COLLAPSE2);
  bench_kernel(f, col_tbl[Vi], "omp_collapse2", V, bd, nlev_end, dl, hr);

  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(SCHED_COLLAPSE2_JE);
  bench_kernel(f, col_je_tbl[Vi], "omp_collapse2_je", V, bd, nlev_end, dl, hr);

  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(SCHED_NUMA4);
  bench_kernel(f, numa4_tbl[Vi], "omp_numa4", V, bd, nlev_end, dl, hr);

  printf("Done: nlev=%d(%d) dist=%-12s V=%d  (N_e=%d N_c=%d N_v=%d)\n",
         nl, nlev_end, dl, V, Ne, Nc, Nv);
}

/* ---- Run blocked ---- */
static void run_blocked(FILE *f, int bi, BenchData &bd, int nlev_end,
                        const char *dl, double *hr) {
  int B = BLOCK_SIZES[bi];
  int Ne=bd.N_e, Nc=bd.N_c, Nv=bd.N_v, nl=bd.nlev;
  memset(hr, 0, bd.sz_e*sizeof(double));
  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bref_tbl[bi](hr, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
               bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
               bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
  flush();
  bd.change_schedule(SCHED_JE_OUTER);
  for (int r=0;r<WARMUP;r++) {
    flush();
    bfor_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
  }
  flush();
  { int nf; double mr;
    verify(bd.h_out, hr, bd.sz_e, &nf, &mr);
    printf("B%d %-12s omp_for  %s mr=%.2e\n", B, dl, nf?"FAIL":"OK", mr); }
  for (int r=0;r<NRUNS;r++) {
    flush();
    auto t0=std::chrono::high_resolution_clock::now();
    bfor_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
    auto t1=std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,0,%d,%d,%d,%d,%s,%d,blocked_omp_for,%d,%.9f\n",
            nlev_end, Ne, Nc, Nv, dl, B, r,
            std::chrono::duration<double,std::milli>(t1-t0).count());
    flush();
  }
  memset(bd.h_out, 0, bd.sz_e*sizeof(double));
  bd.change_schedule(SCHED_COLLAPSE2);
  for (int r=0;r<WARMUP;r++) {
    flush();
    bcol_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
    flush();
  }
  { int nf; double mr;
    verify(bd.h_out, hr, bd.sz_e, &nf, &mr);
    printf("B%d %-12s collapse %s mr=%.2e\n", B, dl, nf?"FAIL":"OK", mr); }
  flush();
  for (int r=0;r<NRUNS;r++) {
    flush();
    auto t0=std::chrono::high_resolution_clock::now();
    bcol_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, Ne, Nc, Nv, nl, nlev_end);
    auto t1=std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,0,%d,%d,%d,%d,%s,%d,blocked_collapse2,%d,%.9f\n",
            nlev_end, Ne, Nc, Nv, dl, B, r,
            std::chrono::duration<double,std::milli>(t1-t0).count());
    flush();
  }
  printf("Done: nlev=%d(%d) dist=%-12s B=%d\n", nl, nlev_end, dl, B);
}

/* ---- Helper: run unblocked V_start..V_end for one dist ---- */
static void run_unblocked_block(
    FILE* fcsv, int N_e, int N_c, int N_v, int nlev, int nlev_end,
    const char* dl, int V_start, int V_end,
    int* cell_logical, int* vert_logical,
    double* icon_inv_dual, double* icon_inv_primal, double* icon_tangent)
{
    BenchData bd;
    bd.alloc(N_e, N_c, N_v, nlev); bd.fill(nlev);
    if (icon_inv_dual)
      for (int je=0;je<N_e;je++) {
        bd.inv_dual[je]=icon_inv_dual[je];
        bd.inv_primal[je]=icon_inv_primal[je];
        bd.tangent_o[je]=icon_tangent[je];
      }
    double *hr = numa_alloc_unfaulted<double>(bd.sz_e);
#pragma omp parallel for schedule(static)
    for (size_t i=0;i<bd.sz_e;i++) hr[i]=0;
    for (int V=V_start; V<=V_end; V++) {
      int kV = kern_v(V);
      bd.set_variant(kV, cell_logical, vert_logical, sched_for_par_for(V));
      run_unblocked(fcsv, V, bd, nlev_end, dl, hr);
      fflush(fcsv);
    }
    numa_dealloc(hr, bd.sz_e);
    bd.free_all();
}

/* ---- main ---- */
int main(int argc, char *argv[]) {
  int icon_step = (argc >= 2) ? atoi(argv[1]) : 9;
  int L1_bytes  = (argc >= 3) ? atoi(argv[2]) : 32768;

  omp_set_max_active_levels(2);

  FILE *fcsv = fopen("z_v_grad_w_cpu.csv", "w");
  fprintf(fcsv, "backend,variant,nlev,N_e,N_c,N_v,cell_dist,blocking,"
                "parallelization,run_id,time_ms\n");

  /* ---- load ICON exact data ---- */
  std::string gp = icon_global_path(icon_step);
  std::string pp = icon_patch_path(icon_step);
  int icon_nproma = icon_read_nproma(gp.c_str());
  if (icon_nproma <= 0)
    fprintf(stderr, "WARNING: could not read nproma from '%s'\n", gp.c_str());
  printf("Loading ICON: %s  (nproma=%d)\n", pp.c_str(), icon_nproma);

  IconEdgeData ied;
  bool have_exact = (icon_nproma > 0) && icon_load_patch(pp.c_str(), icon_nproma, ied);
  if (have_exact)
    printf("ICON: nproma=%d  n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
           ied.nproma, ied.n_edges, ied.n_edges_valid, ied.n_cells, ied.n_verts);
  assert(have_exact);
  if (have_exact) icon_print_locality_metrics(ied);

  /* Dimensions from ICON topology */
  const int N_e = ied.n_edges_valid;  /* 122880 */
  const int N_c = ied.n_cells;        /* 81920  */
  const int N_v = ied.n_verts;        /* 81920  */
  printf("Benchmark dimensions: N_e=%d (valid edges)  N_c=%d (cells)  N_v=%d (verts)\n",
         N_e, N_c, N_v);

  std::mt19937 rng(42);
  int *cell_logical = new int[N_e*2];
  int *vert_logical = new int[N_e*2];

  int nthreads = omp_get_max_threads();
  printf("OMP threads: %d  NUMA domains: %d  L1: %d bytes\n",
         nthreads, NUMA_DOMAINS, L1_bytes);
  for (int bi=0;bi<N_BLOCK_SIZES;bi++)
    printf("  B=%3d L1%%=%.1f%%\n", BLOCK_SIZES[bi], l1_ratio(BLOCK_SIZES[bi],L1_bytes)*100);
  srand((unsigned)time(NULL));
  flush();
  printf("Ready\n\n");

  for (int ni=0; ni<N_NLEVS; ni++) {
    int nlev_base   = NLEVS[ni];
    int nlev_padded = icon_pad_nlev(nlev_base);

    /* Synthetic distributions */
    for (int di=0; di<2; di++) {
      gen_idx_logical(cell_logical, N_e, N_c, (CellDist)di, rng);
      gen_idx_logical(vert_logical, N_e, N_v, (CellDist)di, rng);
      if (ni == 0) print_synthetic_locality(dist_name[di], cell_logical,
                                            vert_logical, N_e, N_c, N_v);
      run_unblocked_block(fcsv, N_e, N_c, N_v, nlev_base, nlev_base,
          dist_name[di], 1, 4, cell_logical, vert_logical,
          nullptr, nullptr, nullptr);
      if (nlev_padded != nlev_base)
        run_unblocked_block(fcsv, N_e, N_c, N_v, nlev_padded, nlev_base,
            dist_name[di], 5, 5, cell_logical, vert_logical,
            nullptr, nullptr, nullptr);
    }

    /* EXACT distribution */
    if (have_exact) {
      int *ecl=new int[N_e*2], *evl=new int[N_e*2];
      for (int i=0;i<N_e*2;i++) { ecl[i]=ied.cell_idx[i]; evl[i]=ied.vert_idx[i]; }
      run_unblocked_block(fcsv, N_e, N_c, N_v, nlev_base, nlev_base,
          "exact", 1, 4, ecl, evl,
          ied.inv_dual.data(), ied.inv_primal.data(), ied.tangent_o.data());
      if (nlev_padded != nlev_base)
        run_unblocked_block(fcsv, N_e, N_c, N_v, nlev_padded, nlev_base,
            "exact", 5, 5, ecl, evl,
            ied.inv_dual.data(), ied.inv_primal.data(), ied.tangent_o.data());
      delete[] ecl; delete[] evl;
    }

    /* Blocked synthetic */
    for (int di=0; di<2; di++) {
      rng.seed(42);
      gen_idx_logical(cell_logical, N_e, N_c, (CellDist)di, rng);
      gen_idx_logical(vert_logical, N_e, N_v, (CellDist)di, rng);
      for (int bi=0;bi<N_BLOCK_SIZES;bi++) {
        int B=BLOCK_SIZES[bi];
        if (N_e % B != 0) { printf("SKIP B=%d !| N_e=%d\n", B, N_e); continue; }
        if (N_c % B != 0) { printf("SKIP B=%d !| N_c=%d\n", B, N_c); continue; }
        BenchData bd; bd.alloc(N_e, N_c, N_v, nlev_base); bd.fill(nlev_base);
        double *hr=numa_alloc_unfaulted<double>(bd.sz_e);
#pragma omp parallel for schedule(static)
        for (size_t i=0;i<bd.sz_e;i++) hr[i]=0;
        bd.set_variant_blocked(B, cell_logical, vert_logical, SCHED_JE_OUTER);
        run_blocked(fcsv, bi, bd, nlev_base, dist_name[di], hr);
        fflush(fcsv);
        numa_dealloc(hr, bd.sz_e); bd.free_all();
      }
    }

    /* Blocked exact */
    if (have_exact) {
      int *ecl=new int[N_e*2], *evl=new int[N_e*2];
      for (int i=0;i<N_e*2;i++) { ecl[i]=ied.cell_idx[i]; evl[i]=ied.vert_idx[i]; }
      for (int bi=0;bi<N_BLOCK_SIZES;bi++) {
        int B=BLOCK_SIZES[bi];
        if (N_e % B != 0) { printf("SKIP B=%d !| N_e=%d\n", B, N_e); continue; }
        if (N_c % B != 0) { printf("SKIP B=%d !| N_c=%d\n", B, N_c); continue; }
        BenchData bd; bd.alloc(N_e, N_c, N_v, nlev_base); bd.fill(nlev_base);
        for (int je=0;je<N_e;je++) {
          bd.inv_dual[je]=ied.inv_dual[je];
          bd.inv_primal[je]=ied.inv_primal[je];
          bd.tangent_o[je]=ied.tangent_o[je];
        }
        double *hr=numa_alloc_unfaulted<double>(bd.sz_e);
#pragma omp parallel for schedule(static)
        for (size_t i=0;i<bd.sz_e;i++) hr[i]=0;
        bd.set_variant_blocked(B, ecl, evl, SCHED_JE_OUTER);
        run_blocked(fcsv, bi, bd, nlev_base, "exact", hr);
        fflush(fcsv);
        numa_dealloc(hr, bd.sz_e); bd.free_all();
      }
      delete[] ecl; delete[] evl;
    }
  }

  numa_dealloc(fb0, (size_t)FN*FN);
  numa_dealloc(fb1, (size_t)FN*FN);
  delete[] cell_logical;
  delete[] vert_logical;
  if (have_exact) ied.free_all();
  fclose(fcsv);
  printf("\nWritten: z_v_grad_w_cpu.csv\n");
  return 0;
}