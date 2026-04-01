/*
 * bench_cpu.cpp -- NUMA-aware CPU z_v_grad_w stencil benchmark
 *
 * Unblocked (V1-V4) + blocked (B=8..128) × 2 OMP schedules × 5 distributions
 *
 * Compile: g++ -O3 -fopenmp -march=native -std=c++17 bench_cpu.cpp -o bench_cpu
 * Run:     ./bench_cpu [timestep] [L1_bytes]
 */
#include "bench_common.h"
#include "icon_data_loader.h"
#include <ctime>
#include <omp.h>
#include <utility>

static SchedKind sched_for_par_for(int V) {
  return (V <= 2) ? SCHED_JK_OUTER : SCHED_JE_OUTER;
}

/* ---- Unblocked kernels ---- */
template <int V>
static void cpu_par_for(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx, int N,
    int nlev) {
  if constexpr (V <= 2) {
#pragma omp parallel for schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int je = 0; je < N; je++) {
        STENCIL_BODY(V)
      }
  } else {
#pragma omp parallel for schedule(static)
    for (int je = 0; je < N; je++)
      for (int jk = 0; jk < nlev; jk++) {
        STENCIL_BODY(V)
      }
  }
}
template <int V>
static void cpu_collapse2(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx, int N,
    int nlev) {
  if constexpr (V <= 2) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int jk = 0; jk < nlev; jk++)
      for (int je = 0; je < N; je++) {
        STENCIL_BODY(V)
      }
  } else {
#pragma omp parallel for collapse(2) schedule(static)
    for (int je = 0; je < N; je++)
      for (int jk = 0; jk < nlev; jk++) {
        STENCIL_BODY(V)
      }
  }
}

typedef void (*kern_t)(double *, const double *, const double *, const double *,
                       const int *, const double *, const double *,
                       const double *, const double *, const int *, int, int);
static kern_t par_tbl[] = {cpu_par_for<1>, cpu_par_for<2>, cpu_par_for<3>,
                           cpu_par_for<4>};
static kern_t col_tbl[] = {cpu_collapse2<1>, cpu_collapse2<2>, cpu_collapse2<3>,
                           cpu_collapse2<4>};

/* ---- Blocked kernels ---- */
template <int B>
static void cpu_blocked_for(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx, int N,
    int nlev) {
  int nblocks = N / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++) {
#pragma omp unroll
    for (int jk = 0; jk < nlev; jk++) {
#pragma omp simd
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
    }
  }
}
template <int B>
static void cpu_blocked_col(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx, int N,
    int nlev) {
  int nblocks = N / B;
#pragma omp parallel for schedule(static)
  for (int jb = 0; jb < nblocks; jb++) {
#pragma omp unroll
    for (int jk = 0; jk < nlev; jk++) {
#pragma omp simd
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
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
static void
cpu_ref(double *__restrict__ out, const double *__restrict__ vn_ie,
        const double *__restrict__ inv_dual, const double *__restrict__ w,
        const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
        const double *__restrict__ inv_primal,
        const double *__restrict__ tangent, const double *__restrict__ z_w_v,
        const int *__restrict__ vert_idx, int N, int nlev) {
  for (int jk = 0; jk < nlev; jk++)
    for (int je = 0; je < N; je++) {
      STENCIL_BODY(V)
    }
}
static void cpu_ref_v(int V, double *o, const double *vn, const double *id,
                      const double *w, const int *ci, const double *vt,
                      const double *ip, const double *tg, const double *zw,
                      const int *vi, int N, int nl) {
  switch (V) {
  case 1:
    cpu_ref<1>(o, vn, id, w, ci, vt, ip, tg, zw, vi, N, nl);
    break;
  case 2:
    cpu_ref<2>(o, vn, id, w, ci, vt, ip, tg, zw, vi, N, nl);
    break;
  case 3:
    cpu_ref<3>(o, vn, id, w, ci, vt, ip, tg, zw, vi, N, nl);
    break;
  case 4:
    cpu_ref<4>(o, vn, id, w, ci, vt, ip, tg, zw, vi, N, nl);
    break;
  }
}
template <int B>
static void cpu_ref_blocked(
    double *__restrict__ out, const double *__restrict__ vn_ie,
    const double *__restrict__ inv_dual, const double *__restrict__ w,
    const int *__restrict__ cell_idx, const double *__restrict__ z_vt_ie,
    const double *__restrict__ inv_primal, const double *__restrict__ tangent,
    const double *__restrict__ z_w_v, const int *__restrict__ vert_idx, int N,
    int nlev) {
  for (int jb = 0; jb < N / B; jb++)
    for (int jk = 0; jk < nlev; jk++)
      for (int jl = 0; jl < B; jl++) {
        int je = jb * B + jl;
        STENCIL_BODY_BLOCKED(B)
      }
}
static kern_t bref_tbl[] = {cpu_ref_blocked<8>, cpu_ref_blocked<16>,
                            cpu_ref_blocked<32>, cpu_ref_blocked<64>,
                            cpu_ref_blocked<128>};

/* ---- Verify ---- */
static bool verify(const double *got, const double *ref, size_t n, int *nf,
                   double *mr) {
  *nf = 0;
  *mr = 0;
  for (size_t i = 0; i < n; i++) {
    double d = std::abs(got[i] - ref[i]),
           dn = std::max(std::abs(ref[i]), 1e-300), r = d / dn;
    if (r > *mr)
      *mr = r;
    if (d > 1e-15 + 1e-12 * std::abs(ref[i]))
      (*nf)++;
  }
  return *nf == 0;
}

/* ---- Flush ---- */
static constexpr int FN = 8192 * 4, FS = 3;
static double *fb0 = nullptr, *fb1 = nullptr;
static void flush() {
  static bool init = false;
  if (!init) {
    size_t n = (size_t)FN * FN;
    fb0 = numa_alloc_unfaulted<double>(n);
    fb1 = numa_alloc_unfaulted<double>(n);
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n; i++) {
      uint64_t h = splitmix64(12345ULL + i);
      fb0[i] = (double)(h >> 11) / (double)(1ULL << 53);
      fb1[i] = fb0[i];
    }
    init = true;
  }
  double *A = fb0, *B = fb1;
  for (int s = 0; s < FS; s++) {
#pragma omp parallel for schedule(static)
    for (int i = 1; i < FN - 1; i++)
      for (int j = 1; j < FN - 1; j++)
        B[i * FN + j] = 0.25 * (A[(i - 1) * FN + j] + A[(i + 1) * FN + j] +
                                A[i * FN + (j - 1)] + A[i * FN + (j + 1)]);
    std::swap(A, B);
  }
  printf("  [flush] A[0]=%.6e\n", A[0]);
}

/* ---- Run unblocked ---- */
static void run_unblocked(FILE *f, int V, int N, int nlev, const char *dl,
                          BenchData &bd, double *hr) {
  bd.change_schedule(sched_for_par_for(V));
  cpu_ref_v(V, hr, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx, bd.h_z_vt_ie,
            bd.inv_primal, bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
  flush();
  for (int r = 0; r < WARMUP; r++) {
    flush();
    par_tbl[V - 1](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                   bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                   bd.h_vidx, N, nlev);
  }
  flush();
  {
    int nf;
    double mr;
    verify(bd.h_out, hr, bd.sz2d, &nf, &mr);
    printf("V%d %-12s omp_for  %s mr=%.2e\n", V, dl, nf ? "FAIL" : "OK", mr);
  }
  for (int r = 0; r < NRUNS; r++) {
    flush();
    auto t0 = std::chrono::high_resolution_clock::now();
    par_tbl[V - 1](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                   bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                   bd.h_vidx, N, nlev);
    auto t1 = std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,%d,%d,%d,%s,0,omp_for,%d,%.9f\n", V, nlev, N, dl, r,
            std::chrono::duration<double, std::milli>(t1 - t0).count());
    flush();
  }
  bd.change_schedule(SCHED_COLLAPSE2);
  for (int r = 0; r < WARMUP; r++) {
    flush();
    col_tbl[V - 1](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                   bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                   bd.h_vidx, N, nlev);
    flush();
  }
  {
    int nf;
    double mr;
    verify(bd.h_out, hr, bd.sz2d, &nf, &mr);
    printf("V%d %-12s collapse %s mr=%.2e\n", V, dl, nf ? "FAIL" : "OK", mr);
  }
  flush();
  for (int r = 0; r < NRUNS; r++) {
    flush();
    auto t0 = std::chrono::high_resolution_clock::now();
    col_tbl[V - 1](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                   bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                   bd.h_vidx, N, nlev);
    auto t1 = std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,%d,%d,%d,%s,0,omp_collapse2,%d,%.9f\n", V, nlev, N, dl, r,
            std::chrono::duration<double, std::milli>(t1 - t0).count());
    flush();
  }
  printf("Done: nlev=%d dist=%-12s V=%d\n", nlev, dl, V);
}

/* ---- Run blocked ---- */
static void run_blocked(FILE *f, int bi, int N, int nlev, const char *dl,
                        BenchData &bd, double *hr) {
  int B = BLOCK_SIZES[bi];
  bref_tbl[bi](hr, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx, bd.h_z_vt_ie,
               bd.inv_primal, bd.tangent_o, bd.h_z_w_v, bd.h_vidx, N, nlev);
  flush();
  bd.change_schedule(SCHED_JE_OUTER);
  for (int r = 0; r < WARMUP; r++) {
    flush();
    bfor_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, N, nlev);
  }
  flush();
  {
    int nf;
    double mr;
    verify(bd.h_out, hr, bd.sz2d, &nf, &mr);
    printf("B%d %-12s omp_for  %s mr=%.2e\n", B, dl, nf ? "FAIL" : "OK", mr);
  }
  for (int r = 0; r < NRUNS; r++) {
    flush();
    auto t0 = std::chrono::high_resolution_clock::now();
    bfor_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, N, nlev);
    auto t1 = std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,0,%d,%d,%s,%d,blocked_omp_for,%d,%.9f\n", nlev, N, dl, B, r,
            std::chrono::duration<double, std::milli>(t1 - t0).count());
    flush();
  }
  bd.change_schedule(SCHED_COLLAPSE2);
  for (int r = 0; r < WARMUP; r++) {
    flush();
    bcol_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, N, nlev);
    flush();
  }
  {
    int nf;
    double mr;
    verify(bd.h_out, hr, bd.sz2d, &nf, &mr);
    printf("B%d %-12s collapse %s mr=%.2e\n", B, dl, nf ? "FAIL" : "OK", mr);
  }
  flush();
  for (int r = 0; r < NRUNS; r++) {
    flush();
    auto t0 = std::chrono::high_resolution_clock::now();
    bcol_tbl[bi](bd.h_out, bd.h_vn_ie, bd.inv_dual, bd.h_w, bd.h_cidx,
                 bd.h_z_vt_ie, bd.inv_primal, bd.tangent_o, bd.h_z_w_v,
                 bd.h_vidx, N, nlev);
    auto t1 = std::chrono::high_resolution_clock::now();
    fprintf(f, "cpu,0,%d,%d,%s,%d,blocked_collapse2,%d,%.9f\n", nlev, N, dl, B,
            r, std::chrono::duration<double, std::milli>(t1 - t0).count());
    flush();
  }
  printf("Done: nlev=%d dist=%-12s B=%d\n", nlev, dl, B);
}

/* ---- main ---- */
int main(int argc, char *argv[]) {
  int icon_step = (argc >= 2) ? atoi(argv[1]) : 9;
  int L1_bytes = (argc >= 3) ? atoi(argv[2]) : 32768;

  FILE *fcsv = fopen("z_v_grad_w_cpu.csv", "w");
  fprintf(fcsv, "backend,variant,nlev,nproma,cell_dist,blocking,"
                "parallelization,run_id,time_ms\n");

  const int N = NPROMA;
  std::mt19937 rng(42);
  VertData vd;
  vd.init(N, rng);
  int *cell_logical = new int[N * 2];

  std::string pp = icon_patch_path(icon_step);
  printf("Loading ICON: %s\n", pp.c_str());
  IconEdgeData ied;
  bool have_exact = icon_load_patch(pp.c_str(), ied);
  if (have_exact)
    printf("ICON: Ne=%d Nc=%d Nv=%d\n", ied.n_edges_valid, ied.n_cells,
           ied.n_verts);

  printf("OMP threads: %d  L1: %d bytes\n", omp_get_max_threads(), L1_bytes);
  for (int bi = 0; bi < N_BLOCK_SIZES; bi++)
    printf("  B=%3d L1%%=%.1f%%\n", BLOCK_SIZES[bi],
           l1_ratio(BLOCK_SIZES[bi], L1_bytes) * 100);
  srand((unsigned)time(NULL));
  flush();
  printf("Ready\n\n");

  for (int ni = 0; ni < N_NLEVS; ni++) {
    int nlev = NLEVS[ni];

    /* -- Unblocked synthetic -- */
    {
      BenchData bd;
      bd.alloc(N, nlev);
      bd.fill(nlev);
      double *hr = numa_alloc_unfaulted<double>(bd.sz2d);
#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < bd.sz2d; i++)
        hr[i] = 0;
      for (int di = 0; di < 4; di++) {
        gen_cell_idx_logical(cell_logical, N, (CellDist)di, rng);
        for (int V = 1; V <= 4; V++) {
          bd.set_variant(V, cell_logical, vd.logical, sched_for_par_for(V));
          run_unblocked(fcsv, V, N, nlev, dist_name[di], bd, hr);
          fflush(fcsv);
        }
      }
      numa_dealloc(hr, bd.sz2d);
      bd.free_all();
    }

    /* -- Unblocked exact -- */
    if (have_exact) {
      int Ne = ied.n_edges_valid;
      BenchData bd;
      bd.alloc(Ne, nlev);
      bd.fill(nlev);
      double *hr = numa_alloc_unfaulted<double>(bd.sz2d);
#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < bd.sz2d; i++)
        hr[i] = 0;
      int *ecl = new int[Ne * 2], *evl = new int[Ne * 2];
      for (int i = 0; i < Ne * 2; i++) {
        ecl[i] = ied.cell_idx[i];
        evl[i] = ied.vert_idx[i];
      }
      for (int je = 0; je < Ne; je++) {
        bd.inv_dual[je] = ied.inv_dual[je];
        bd.inv_primal[je] = ied.inv_primal[je];
        bd.tangent_o[je] = ied.tangent_o[je];
      }
      for (int V = 1; V <= 4; V++) {
        bd.set_variant(V, ecl, evl, sched_for_par_for(V));
        run_unblocked(fcsv, V, Ne, nlev, "exact", bd, hr);
        fflush(fcsv);
      }
      delete[] ecl;
      delete[] evl;
      numa_dealloc(hr, bd.sz2d);
      bd.free_all();
    }

    /* -- Blocked synthetic -- */
    for (int di = 0; di < 4; di++) {
      rng.seed(42);
      gen_cell_idx_logical(cell_logical, N, (CellDist)di, rng);
      VertData vd2;
      vd2.init(N, rng);
      for (int bi = 0; bi < N_BLOCK_SIZES; bi++) {
        int B = BLOCK_SIZES[bi];
        if (N % B != 0)
          continue;
        BenchData bd;
        bd.alloc(N, nlev);
        bd.fill(nlev);
        double *hr = numa_alloc_unfaulted<double>(bd.sz2d);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < bd.sz2d; i++)
          hr[i] = 0;
        bd.set_variant_blocked(B, cell_logical, vd2.logical, SCHED_JE_OUTER);
        run_blocked(fcsv, bi, N, nlev, dist_name[di], bd, hr);
        fflush(fcsv);
        numa_dealloc(hr, bd.sz2d);
        bd.free_all();
      }
      vd2.free_all();
    }

    /* -- Blocked exact -- */
    if (have_exact) {
      int Ne = ied.n_edges_valid;
      int *ecl = new int[Ne * 2], *evl = new int[Ne * 2];
      for (int i = 0; i < Ne * 2; i++) {
        ecl[i] = ied.cell_idx[i];
        evl[i] = ied.vert_idx[i];
      }
      for (int bi = 0; bi < N_BLOCK_SIZES; bi++) {
        int B = BLOCK_SIZES[bi];
        if (Ne % B != 0) {
          printf("SKIP B=%d !| Ne=%d\n", B, Ne);
          continue;
        }
        BenchData bd;
        bd.alloc(Ne, nlev);
        bd.fill(nlev);
        for (int je = 0; je < Ne; je++) {
          bd.inv_dual[je] = ied.inv_dual[je];
          bd.inv_primal[je] = ied.inv_primal[je];
          bd.tangent_o[je] = ied.tangent_o[je];
        }
        double *hr = numa_alloc_unfaulted<double>(bd.sz2d);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < bd.sz2d; i++)
          hr[i] = 0;
        bd.set_variant_blocked(B, ecl, evl, SCHED_JE_OUTER);
        run_blocked(fcsv, bi, Ne, nlev, "exact", bd, hr);
        fflush(fcsv);
        numa_dealloc(hr, bd.sz2d);
        bd.free_all();
      }
      delete[] ecl;
      delete[] evl;
    }
  }

  numa_dealloc(fb0, (size_t)FN * FN);
  numa_dealloc(fb1, (size_t)FN * FN);
  vd.free_all();
  delete[] cell_logical;
  if (have_exact)
    ied.free_all();
  fclose(fcsv);
  printf("\nWritten: z_v_grad_w_cpu.csv\n");
  return 0;
}