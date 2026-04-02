#ifndef ICON_DATA_LOADER_H
#define ICON_DATA_LOADER_H
/*
 * icon_data_loader.h -- Standalone reader for ICON's serialised p_patch files.
 *
 * Reads the text-based serde format produced by the f2dace infrastructure and
 * extracts the edge-level arrays needed by the z_v_grad_w micro-benchmark.
 *
 * Key design decisions:
 *   - nproma is read from the global data file (authoritative source).
 *   - All arrays are flattened to their FULL blocked extent:
 *       edges:    nproma * nblks_e
 *       cells:    nproma * nblks_c
 *       vertices: nproma * nblks_v
 *   - Connectivity values (cell_idx, cell_blk, vertex_idx, vertex_blk) are
 *     Fortran 1-based.  We linearise with:
 *       flat_target = (blk_val - 1) * nproma + (idx_val - 1)
 *   - The serde array stride (size[0]) is used for reading the raw data;
 *     nproma is used for the output flat layout.  These are normally equal
 *     but the distinction is maintained for correctness.
 *   - Out-of-range or invalid references (e.g. boundary padding with
 *     idx<=0 or blk<=0) are clamped to 0.
 */

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

/* ------------------------------------------------------------------ */
/*  Minimal serde helpers (mirrors the project's serde scheme)        */
/* ------------------------------------------------------------------ */
namespace icon_serde {

inline void scroll_space(std::istream &s) {
  while (!s.eof() && s.peek() != EOF && std::isspace(s.peek()))
    s.get();
}

inline std::string read_line(std::istream &s) {
  scroll_space(s);
  std::string line;
  std::getline(s, line);
  return line;
}

inline bool expect_line(std::istream &s, const char *tag) {
  std::string l = read_line(s);
  if (l.find(tag) == std::string::npos) {
    fprintf(stderr, "[icon_serde] expected '%s', got '%s'\n", tag, l.c_str());
    return false;
  }
  return true;
}

inline int read_int(std::istream &s) {
  scroll_space(s);
  int v = 0;
  s >> v;
  return v;
}

inline double read_double(std::istream &s) {
  scroll_space(s);
  long double v = 0;
  s >> v;
  return static_cast<double>(v);
}

inline bool read_bool(std::istream &s) {
  scroll_space(s);
  char c = '0';
  s >> c;
  return c == '1';
}

struct ArrayMeta {
  int rank = 0;
  int size[4] = {};
  int lbound[4] = {};
  int volume() const {
    int v = 1;
    for (int i = 0; i < rank; i++)
      v *= size[i];
    return v;
  }
};

inline ArrayMeta read_meta(std::istream &s) {
  ArrayMeta m;
  expect_line(s, "# rank");
  m.rank = read_int(s);
  assert(m.rank >= 1 && m.rank <= 4);
  expect_line(s, "# size");
  for (int i = 0; i < m.rank; i++)
    m.size[i] = read_int(s);
  expect_line(s, "# lbound");
  for (int i = 0; i < m.rank; i++)
    m.lbound[i] = read_int(s);
  return m;
}

inline void skip_entries(std::istream &s, int count) {
  std::string tok;
  for (int i = 0; i < count; i++) {
    scroll_space(s);
    s >> tok;
  }
}

/* Read an "# alloc"-guarded array block.
   If `out` is non-null, fills it; otherwise skips the data.
   Returns meta (rank==0 if the array was not allocated). */
template <typename T>
ArrayMeta read_alloc_array(std::istream &s, std::vector<T> *out = nullptr) {
  expect_line(s, "# alloc");
  bool alloc = read_bool(s);
  if (!alloc)
    return {};

  ArrayMeta m = read_meta(s);
  expect_line(s, "# entries");
  int vol = m.volume();
  if (out) {
    out->resize(vol);
    for (int i = 0; i < vol; i++) {
      scroll_space(s);
      if constexpr (std::is_same_v<T, int>)
        s >> (*out)[i];
      else {
        long double tmp;
        s >> tmp;
        (*out)[i] = static_cast<T>(tmp);
      }
    }
  } else {
    skip_entries(s, vol);
  }
  return m;
}

/* Read an "# assoc / # missing"-guarded pointer array block. */
template <typename T>
ArrayMeta read_assoc_array(std::istream &s, std::vector<T> *out = nullptr) {
  expect_line(s, "# assoc");
  bool assoc = read_bool(s);
  if (!assoc)
    return {};

  expect_line(s, "# missing");
  read_int(s); /* always 1 */

  ArrayMeta m = read_meta(s);
  expect_line(s, "# entries");
  int vol = m.volume();
  if (out) {
    out->resize(vol);
    for (int i = 0; i < vol; i++) {
      scroll_space(s);
      if constexpr (std::is_same_v<T, int>)
        s >> (*out)[i];
      else {
        long double tmp;
        s >> tmp;
        (*out)[i] = static_cast<T>(tmp);
      }
    }
  } else {
    skip_entries(s, vol);
  }
  return m;
}

/* ---- skip helpers for compound types we don't need ---- */

inline void skip_decomp_info(std::istream &s) {
  expect_line(s, "# owner_mask");
  read_alloc_array<int>(s);
}

inline void skip_grid_cells(std::istream &s) {
  expect_line(s, "# neighbor_idx");
  read_alloc_array<int>(s);
  expect_line(s, "# neighbor_blk");
  read_alloc_array<int>(s);
  expect_line(s, "# edge_idx");
  read_alloc_array<int>(s);
  expect_line(s, "# edge_blk");
  read_alloc_array<int>(s);
  expect_line(s, "# area");
  read_assoc_array<double>(s);
  expect_line(s, "# start_index");
  read_alloc_array<int>(s);
  expect_line(s, "# end_index");
  read_alloc_array<int>(s);
  expect_line(s, "# start_block");
  read_alloc_array<int>(s);
  expect_line(s, "# end_block");
  read_alloc_array<int>(s);
  expect_line(s, "# decomp_info");
  skip_decomp_info(s);
}

inline void skip_grid_vertices(std::istream &s) {
  expect_line(s, "# cell_idx");
  read_alloc_array<int>(s);
  expect_line(s, "# cell_blk");
  read_alloc_array<int>(s);
  expect_line(s, "# edge_idx");
  read_alloc_array<int>(s);
  expect_line(s, "# edge_blk");
  read_alloc_array<int>(s);
  expect_line(s, "# start_index");
  read_alloc_array<int>(s);
  expect_line(s, "# end_index");
  read_alloc_array<int>(s);
  expect_line(s, "# start_block");
  read_alloc_array<int>(s);
  expect_line(s, "# end_block");
  read_alloc_array<int>(s);
}

} // namespace icon_serde

/* ------------------------------------------------------------------ */
/*  Public data structure returned by the loader                      */
/* ------------------------------------------------------------------ */
struct IconEdgeData {
  int nproma = 0;
  int nblks_c = 0;
  int nblks_e = 0;
  int nblks_v = 0;
  int n_edges = 0;       /* nproma * nblks_e  (full blocked extent) */
  int n_cells = 0;       /* nproma * nblks_c */
  int n_verts = 0;       /* nproma * nblks_v */
  int n_edges_valid = 0; /* actual valid edges (from end_block/end_index) */

  /* Flat connectivity: cell_idx[je*2 + n], n=0,1 -> flat 0-based cell id
   *                    (0 <= cell_id < n_cells)
   *                    vert_idx[je*2 + n]        -> flat 0-based vert id
   *                    (0 <= vert_id < n_verts)
   * Length: n_edges * 2 each. */
  std::vector<int> cell_idx;
  std::vector<int> vert_idx;

  /* Flat 1-D edge geometry (length = n_edges) */
  std::vector<double> inv_dual;
  std::vector<double> inv_primal;
  std::vector<double> tangent_o;

  void free_all() {
    cell_idx.clear();    cell_idx.shrink_to_fit();
    vert_idx.clear();    vert_idx.shrink_to_fit();
    inv_dual.clear();    inv_dual.shrink_to_fit();
    inv_primal.clear();  inv_primal.shrink_to_fit();
    tangent_o.clear();   tangent_o.shrink_to_fit();
  }
};

/* ------------------------------------------------------------------ */
/*  Linearisation helpers                                             */
/* ------------------------------------------------------------------ */

/*
 * Linearise a blocked (s0, s1, 2) connectivity array into flat [n_edges * 2].
 *
 * s0, s1 are the serde array dimensions (s0 = allocated first dim, s1 = nblks).
 * nproma is the real ICON nproma (used for both output stride and target
 * linearisation).
 *
 * Stored values (idx_val, blk_val) are Fortran 1-based references into
 * the *target* entity grid (cells or vertices).  The flat 0-based target
 * index is:
 *   target_flat = (blk_val - 1) * nproma + (idx_val - 1)
 *
 * Invalid references (boundary padding: idx<=0 or blk<=0) -> clamped to 0.
 */
static void linearise_connectivity(
    const std::vector<int> &raw_idx, const icon_serde::ArrayMeta &meta_idx,
    const std::vector<int> &raw_blk, int nproma,
    int target_max, /* n_cells or n_verts */
    int n_edges,    /* nproma * nblks_e -- full extent */
    std::vector<int> &out) {
  int s0 = meta_idx.size[0]; /* serde first dim (== nproma normally)  */
  int s1 = meta_idx.size[1]; /* nblks_e                               */

  out.assign(n_edges * 2, 0);

  for (int n = 0; n < 2; n++) {
    for (int jb = 0; jb < s1; jb++) {
      for (int jc = 0; jc < nproma; jc++) {
        int edge_linear = jb * nproma + jc;
        if (edge_linear >= n_edges)
          continue;

        int serde_flat = jc + jb * s0 + n * s0 * s1;
        int idx_val = raw_idx[serde_flat]; /* Fortran 1-based */
        int blk_val = raw_blk[serde_flat]; /* Fortran 1-based */

        int target_flat = (blk_val - 1) * nproma + (idx_val - 1);

        if (idx_val <= 0 || blk_val <= 0 || target_flat < 0 ||
            target_flat >= target_max)
          target_flat = 0;

        out[edge_linear * 2 + n] = target_flat;
      }
    }
  }
}

/*
 * Linearise a 2-D blocked (s0, nblks) double array into flat [nproma * nblks].
 *
 * s0 is the serde first dimension; nproma is used for the output stride.
 */
static void linearise_2d_double(const std::vector<double> &raw,
                                int s0_serde, int nproma, int nblks,
                                std::vector<double> &out) {
  int n_flat = nproma * nblks;
  out.assign(n_flat, 0.0);
  for (int jb = 0; jb < nblks; jb++) {
    for (int jc = 0; jc < nproma; jc++) {
      int serde_idx = jc + jb * s0_serde;
      int flat_idx  = jb * nproma + jc;
      out[flat_idx] = raw[serde_idx];
    }
  }
}

/* ------------------------------------------------------------------ */
/*  Read nproma from the global data file                             */
/*                                                                    */
/*  The global data serde format (from deserialize_global_data):      */
/*    # nflatlev      array (rank/size/lbound/entries, NO alloc guard) */
/*    # i_am_accel_node   scalar bool                                 */
/*    # lextra_diffu      scalar bool                                 */
/*    # nproma            scalar int   <-- what we want               */
/*    ...                                                             */
/* ------------------------------------------------------------------ */
inline int icon_read_nproma(const char *global_path) {
  using namespace icon_serde;

  std::ifstream f(global_path);
  if (!f.is_open()) {
    fprintf(stderr, "[icon_data_loader] cannot open global data '%s'\n",
            global_path);
    return -1;
  }

  /* 1. Skip nflatlev array (no # alloc guard -- bare array) */
  expect_line(f, "# nflatlev");
  {
    ArrayMeta m = read_meta(f);
    expect_line(f, "# entries");
    skip_entries(f, m.volume());
  }

  /* 2. Skip i_am_accel_node (scalar bool) */
  expect_line(f, "# i_am_accel_node");
  read_int(f);

  /* 3. Skip lextra_diffu (scalar bool) */
  expect_line(f, "# lextra_diffu");
  read_int(f);

  /* 4. Read nproma */
  expect_line(f, "# nproma");
  int nproma = read_int(f);

  fprintf(stderr, "[icon_data_loader] global nproma = %d\n", nproma);
  return nproma;
}

/* ------------------------------------------------------------------ */
/*  Main loader                                                       */
/*                                                                    */
/*  nproma must be provided (read it from global data first).         */
/* ------------------------------------------------------------------ */
inline bool icon_load_patch(const char *path, int nproma, IconEdgeData &ed) {
  using namespace icon_serde;

  if (nproma <= 0) {
    fprintf(stderr, "[icon_data_loader] invalid nproma=%d\n", nproma);
    return false;
  }

  std::ifstream f(path);
  if (!f.is_open()) {
    fprintf(stderr, "[icon_data_loader] cannot open '%s'\n", path);
    return false;
  }
  fprintf(stderr, "[icon_data_loader] reading %s  (nproma=%d)\n",
          path, nproma);

  ed.nproma = nproma;

  /* ---- top-level scalars ---- */
  expect_line(f, "# nblks_c");
  ed.nblks_c = read_int(f);
  expect_line(f, "# nblks_e");
  ed.nblks_e = read_int(f);
  expect_line(f, "# nblks_v");
  ed.nblks_v = read_int(f);

  /* ---- derived flat sizes ---- */
  ed.n_edges = nproma * ed.nblks_e;
  ed.n_cells = nproma * ed.nblks_c;
  ed.n_verts = nproma * ed.nblks_v;

  /* ---- skip cells ---- */
  expect_line(f, "# cells");
  skip_grid_cells(f);

  /* ---- read edges ---- */
  expect_line(f, "# edges");

  std::vector<int> raw_cell_idx, raw_cell_blk;
  std::vector<int> raw_vert_idx, raw_vert_blk;
  std::vector<double> raw_inv_primal, raw_inv_dual, raw_tangent;
  std::vector<int> raw_end_index, raw_end_block;

  ArrayMeta m_cidx, m_cblk, m_vidx, m_vblk;
  ArrayMeta m_tang, m_invp, m_invd;

  /* 1 */ expect_line(f, "# cell_idx");
  m_cidx = read_alloc_array<int>(f, &raw_cell_idx);
  /* 2 */ expect_line(f, "# cell_blk");
  m_cblk = read_alloc_array<int>(f, &raw_cell_blk);
  /* 3 */ expect_line(f, "# vertex_idx");
  m_vidx = read_alloc_array<int>(f, &raw_vert_idx);
  /* 4 */ expect_line(f, "# vertex_blk");
  m_vblk = read_alloc_array<int>(f, &raw_vert_blk);
  /* 5 */ expect_line(f, "# tangent_orientation");
  m_tang = read_alloc_array<double>(f, &raw_tangent);
  /* 6 */ expect_line(f, "# quad_idx");
  read_alloc_array<int>(f); /* skip */
  /* 7 */ expect_line(f, "# quad_blk");
  read_alloc_array<int>(f); /* skip */
  /* 8 */ expect_line(f, "# inv_primal_edge_length");
  m_invp = read_alloc_array<double>(f, &raw_inv_primal);
  /* 9 */ expect_line(f, "# inv_dual_edge_length");
  m_invd = read_alloc_array<double>(f, &raw_inv_dual);
  /*10 */ expect_line(f, "# area_edge");
  read_alloc_array<double>(f); /* skip */
  /*11 */ expect_line(f, "# f_e");
  read_alloc_array<double>(f); /* skip */
  /*12 */ expect_line(f, "# fn_e");
  read_alloc_array<double>(f); /* skip */
  /*13 */ expect_line(f, "# ft_e");
  read_alloc_array<double>(f); /* skip */
  /*14 */ expect_line(f, "# start_index");
  read_alloc_array<int>(f); /* skip */
  /*15 */ expect_line(f, "# end_index");
  read_alloc_array<int>(f, &raw_end_index);
  /*16 */ expect_line(f, "# start_block");
  read_alloc_array<int>(f); /* skip */
  /*17 */ expect_line(f, "# end_block");
  read_alloc_array<int>(f, &raw_end_block);

  /* ---- sanity-check serde dims vs nproma ---- */
  {
    int s0 = m_cidx.size[0];
    if (s0 != nproma) {
      fprintf(stderr,
              "[icon_data_loader] WARNING: serde size[0]=%d != nproma=%d\n",
              s0, nproma);
    }
    int s1 = m_cidx.size[1];
    if (s1 != ed.nblks_e) {
      fprintf(stderr,
              "[icon_data_loader] WARNING: serde size[1]=%d != nblks_e=%d\n",
              s1, ed.nblks_e);
    }
  }

  /* ---- compute n_edges_valid (informational) ---- */
  {
    int max_flat = 0;
    int n_entries = (int)raw_end_block.size();
    for (int i = 0; i < n_entries; i++) {
      if (raw_end_block[i] <= 0 || raw_end_index[i] <= 0)
        continue;
      int flat = (raw_end_block[i] - 1) * nproma + raw_end_index[i];
      if (flat > max_flat)
        max_flat = flat;
    }
    ed.n_edges_valid = max_flat;
    if (ed.n_edges_valid <= 0)
      ed.n_edges_valid = ed.n_edges;
    if (ed.n_edges_valid > ed.n_edges)
      ed.n_edges_valid = ed.n_edges;
  }

  fprintf(stderr,
          "[icon_data_loader] nproma=%d  nblks_e=%d  nblks_c=%d  nblks_v=%d\n"
          "                   n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
          ed.nproma, ed.nblks_e, ed.nblks_c, ed.nblks_v,
          ed.n_edges, ed.n_edges_valid, ed.n_cells, ed.n_verts);

  /* ---- linearise connectivity (full nproma * nblks_e extent) ----
   *
   * cell_idx targets cells -> target_max = n_cells = nproma * nblks_c
   * vert_idx targets verts -> target_max = n_verts = nproma * nblks_v
   */
  linearise_connectivity(raw_cell_idx, m_cidx, raw_cell_blk, nproma,
                         ed.n_cells, ed.n_edges, ed.cell_idx);

  linearise_connectivity(raw_vert_idx, m_vidx, raw_vert_blk, nproma,
                         ed.n_verts, ed.n_edges, ed.vert_idx);

  /* ---- linearise 1-D geometry (full nproma * nblks_e) ---- */
  linearise_2d_double(raw_inv_dual,   m_invd.size[0], nproma, ed.nblks_e, ed.inv_dual);
  linearise_2d_double(raw_inv_primal, m_invp.size[0], nproma, ed.nblks_e, ed.inv_primal);
  linearise_2d_double(raw_tangent,    m_tang.size[0], nproma, ed.nblks_e, ed.tangent_o);

  fprintf(stderr,
          "[icon_data_loader] done.  Sample cell_idx[0]={%d,%d}  "
          "vert_idx[0]={%d,%d}\n",
          ed.cell_idx[0], ed.cell_idx[1], ed.vert_idx[0], ed.vert_idx[1]);

  return true;
}

/* ------------------------------------------------------------------ */
/*  Convenience: build paths from ICON_DATA_PATH env + timestep       */
/* ------------------------------------------------------------------ */
inline std::string icon_data_dir() {
  const char *env = std::getenv("ICON_DATA_PATH");
  std::string dir =
      env ? env : "/home/primrose/Work/icon-artifacts/velocity/data_r02b05";
  while (!dir.empty() && dir.back() == '/')
    dir.pop_back();
  return dir;
}

inline std::string icon_patch_path(int timestep = 9) {
  return icon_data_dir() + "/p_patch." + std::to_string(timestep) + ".data";
}

inline std::string icon_global_path(int timestep = 9) {
  return icon_data_dir() + "/global_data." + std::to_string(timestep) + ".data";
}

#endif /* ICON_DATA_LOADER_H */