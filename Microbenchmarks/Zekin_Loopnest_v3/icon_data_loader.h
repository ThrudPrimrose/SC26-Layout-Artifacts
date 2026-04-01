#ifndef ICON_DATA_LOADER_H
#define ICON_DATA_LOADER_H
/*
 * icon_data_loader.h -- Standalone reader for ICON's serialised p_patch files.
 *
 * Reads the text-based serde format produced by the f2dace infrastructure and
 * extracts the edge-level arrays needed by the z_v_grad_w micro-benchmark.
 *
 * Key design decisions:
 *   - We read end_block / end_index from the edges to compute the actual
 *     number of valid edges (the last block is typically partially filled).
 *   - Connectivity values (cell_idx, cell_blk, vertex_idx, vertex_blk) are
 *     Fortran 1-based.  We linearise with:
 *       flat_target = (blk_val - 1) * nproma + (idx_val - 1)
 *     This is correct regardless of the lbounds of the *source* array,
 *     because the stored values are always 1-based block/index references.
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

inline void scroll_space(std::istream& s) {
    while (!s.eof() && s.peek() != EOF && std::isspace(s.peek()))
        s.get();
}

inline std::string read_line(std::istream& s) {
    scroll_space(s);
    std::string line;
    std::getline(s, line);
    return line;
}

inline bool expect_line(std::istream& s, const char* tag) {
    std::string l = read_line(s);
    if (l.find(tag) == std::string::npos) {
        fprintf(stderr, "[icon_serde] expected '%s', got '%s'\n",
                tag, l.c_str());
        return false;
    }
    return true;
}

inline int read_int(std::istream& s) {
    scroll_space(s);
    int v = 0;
    s >> v;
    return v;
}

inline double read_double(std::istream& s) {
    scroll_space(s);
    long double v = 0;
    s >> v;
    return static_cast<double>(v);
}

inline bool read_bool(std::istream& s) {
    scroll_space(s);
    char c = '0';
    s >> c;
    return c == '1';
}

struct ArrayMeta {
    int rank = 0;
    int size[4]   = {};
    int lbound[4] = {};
    int volume() const {
        int v = 1;
        for (int i = 0; i < rank; i++) v *= size[i];
        return v;
    }
};

inline ArrayMeta read_meta(std::istream& s) {
    ArrayMeta m;
    expect_line(s, "# rank");
    m.rank = read_int(s);
    assert(m.rank >= 1 && m.rank <= 4);
    expect_line(s, "# size");
    for (int i = 0; i < m.rank; i++) m.size[i]   = read_int(s);
    expect_line(s, "# lbound");
    for (int i = 0; i < m.rank; i++) m.lbound[i] = read_int(s);
    return m;
}

inline void skip_entries(std::istream& s, int count) {
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
ArrayMeta read_alloc_array(std::istream& s,
                           std::vector<T>* out = nullptr)
{
    expect_line(s, "# alloc");
    bool alloc = read_bool(s);
    if (!alloc) return {};

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
ArrayMeta read_assoc_array(std::istream& s,
                           std::vector<T>* out = nullptr)
{
    expect_line(s, "# assoc");
    bool assoc = read_bool(s);
    if (!assoc) return {};

    expect_line(s, "# missing");
    read_int(s);  /* always 1 */

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

inline void skip_decomp_info(std::istream& s) {
    expect_line(s, "# owner_mask");
    read_alloc_array<int>(s);
}

inline void skip_grid_cells(std::istream& s) {
    expect_line(s, "# neighbor_idx");  read_alloc_array<int>(s);
    expect_line(s, "# neighbor_blk");  read_alloc_array<int>(s);
    expect_line(s, "# edge_idx");      read_alloc_array<int>(s);
    expect_line(s, "# edge_blk");      read_alloc_array<int>(s);
    expect_line(s, "# area");          read_assoc_array<double>(s);
    expect_line(s, "# start_index");   read_alloc_array<int>(s);
    expect_line(s, "# end_index");     read_alloc_array<int>(s);
    expect_line(s, "# start_block");   read_alloc_array<int>(s);
    expect_line(s, "# end_block");     read_alloc_array<int>(s);
    expect_line(s, "# decomp_info");   skip_decomp_info(s);
}

inline void skip_grid_vertices(std::istream& s) {
    expect_line(s, "# cell_idx");      read_alloc_array<int>(s);
    expect_line(s, "# cell_blk");      read_alloc_array<int>(s);
    expect_line(s, "# edge_idx");      read_alloc_array<int>(s);
    expect_line(s, "# edge_blk");      read_alloc_array<int>(s);
    expect_line(s, "# start_index");   read_alloc_array<int>(s);
    expect_line(s, "# end_index");     read_alloc_array<int>(s);
    expect_line(s, "# start_block");   read_alloc_array<int>(s);
    expect_line(s, "# end_block");     read_alloc_array<int>(s);
}

} // namespace icon_serde

/* ------------------------------------------------------------------ */
/*  Public data structure returned by the loader                      */
/* ------------------------------------------------------------------ */
struct IconEdgeData {
    int nproma        = 0;
    int nblks_c       = 0;
    int nblks_e       = 0;
    int nblks_v       = 0;
    int n_edges_valid = 0;   /* actual valid edges (from end_block/end_index) */
    int n_edges_alloc = 0;   /* nproma * nblks_e (includes last-block padding) */
    int n_cells       = 0;   /* nproma * nblks_c */
    int n_verts       = 0;   /* nproma * nblks_v */

    /* Flat connectivity: cell_idx[je*2 + n], n=0,1 -> flat 0-based cell id */
    std::vector<int> cell_idx;
    std::vector<int> vert_idx;

    /* Flat 1-D edge geometry (length = n_edges_valid) */
    std::vector<double> inv_dual;
    std::vector<double> inv_primal;
    std::vector<double> tangent_o;

    void free_all() {
        cell_idx.clear();   cell_idx.shrink_to_fit();
        vert_idx.clear();   vert_idx.shrink_to_fit();
        inv_dual.clear();   inv_dual.shrink_to_fit();
        inv_primal.clear(); inv_primal.shrink_to_fit();
        tangent_o.clear();  tangent_o.shrink_to_fit();
    }
};

/* ------------------------------------------------------------------ */
/*  Linearisation helpers                                             */
/* ------------------------------------------------------------------ */

/*
 * Linearise a blocked (nproma, nblks_e, 2) connectivity pair into
 * flat [n_valid * 2].
 *
 * serde column-major order:
 *   serde_flat(jc, jb, n) = jc + jb * s0 + n * s0 * s1
 *   where s0 = size[0] = nproma, s1 = size[1] = nblks_e
 *
 * Stored values (idx_val, blk_val) are Fortran 1-based references into
 * the *target* entity grid (cells or vertices).  The flat 0-based target
 * index is ALWAYS:
 *   target_flat = (blk_val - 1) * nproma + (idx_val - 1)
 *
 * This does NOT depend on the lbounds of the source array -- the values
 * in cell_idx/cell_blk are 1-based Fortran indices regardless.
 *
 * Invalid references (boundary padding: idx<=0 or blk<=0) -> clamped to 0.
 */
static void linearise_connectivity(
    const std::vector<int>& raw_idx,
    const icon_serde::ArrayMeta& meta_idx,
    const std::vector<int>& raw_blk,
    int nproma,
    int target_max,        /* n_cells or n_verts (for bounds check) */
    int n_valid,           /* only linearise this many edges        */
    std::vector<int>& out)
{
    int s0 = meta_idx.size[0];   /* nproma        */
    int s1 = meta_idx.size[1];   /* nblks_e       */

    out.resize(n_valid * 2);

    for (int n = 0; n < 2; n++) {
        for (int jb = 0; jb < s1; jb++) {
            for (int jc = 0; jc < s0; jc++) {
                int edge_linear = jb * s0 + jc;
                if (edge_linear >= n_valid) continue;

                int serde_flat = jc + jb * s0 + n * s0 * s1;
                int idx_val = raw_idx[serde_flat];   /* Fortran 1-based */
                int blk_val = raw_blk[serde_flat];   /* Fortran 1-based */

                /* Convert to flat 0-based.  1-based -> 0-based:
                   flat = (blk - 1) * nproma + (idx - 1) */
                int target_flat = (blk_val - 1) * nproma + (idx_val - 1);

                /* Clamp invalid / out-of-range to 0 */
                if (idx_val <= 0 || blk_val <= 0 ||
                    target_flat < 0 || target_flat >= target_max)
                    target_flat = 0;

                out[edge_linear * 2 + n] = target_flat;
            }
        }
    }
}

/*
 * Linearise a 2-D blocked (nproma, nblks_e) double array,
 * taking only the first n_valid elements.
 */
static void linearise_2d_double(
    const std::vector<double>& raw,
    int n_valid,
    std::vector<double>& out)
{
    out.resize(n_valid);
    /* serde column-major with dim-0 = nproma fastest is identical to
       our flat linearisation order, so just copy the first n_valid. */
    for (int i = 0; i < n_valid; i++)
        out[i] = raw[i];
}

/* ------------------------------------------------------------------ */
/*  Main loader                                                       */
/* ------------------------------------------------------------------ */
inline bool icon_load_patch(const char* path, IconEdgeData& ed) {
    using namespace icon_serde;

    std::ifstream f(path);
    if (!f.is_open()) {
        fprintf(stderr, "[icon_data_loader] cannot open '%s'\n", path);
        return false;
    }
    fprintf(stderr, "[icon_data_loader] reading %s ...\n", path);

    /* ---- top-level scalars ---- */
    expect_line(f, "# nblks_c");  ed.nblks_c = read_int(f);
    expect_line(f, "# nblks_e");  ed.nblks_e = read_int(f);
    expect_line(f, "# nblks_v");  ed.nblks_v = read_int(f);

    /* ---- skip cells ---- */
    expect_line(f, "# cells");
    skip_grid_cells(f);

    /* ---- read edges ---- */
    expect_line(f, "# edges");

    /*
     * Field order in serde for t_grid_edges (from the deserialize function):
     *  1. cell_idx            alloc  (nproma, nblks_e, 2)  int
     *  2. cell_blk            alloc  (nproma, nblks_e, 2)  int
     *  3. vertex_idx          alloc  (nproma, nblks_e, 2)  int
     *  4. vertex_blk          alloc  (nproma, nblks_e, 2)  int
     *  5. tangent_orientation alloc  (nproma, nblks_e)     double
     *  6. quad_idx            alloc  (nproma, nblks_e, 4)  int    [skip]
     *  7. quad_blk            alloc  (nproma, nblks_e, 4)  int    [skip]
     *  8. inv_primal_edge_len alloc  (nproma, nblks_e)     double
     *  9. inv_dual_edge_len   alloc  (nproma, nblks_e)     double
     * 10. area_edge           alloc  (nproma, nblks_e)     double [skip]
     * 11. f_e                 alloc  (nproma, nblks_e)     double [skip]
     * 12. fn_e                alloc  (nproma, nblks_e)     double [skip]
     * 13. ft_e                alloc  (nproma, nblks_e)     double [skip]
     * 14. start_index         alloc  (nlevs)               int    [skip]
     * 15. end_index           alloc  (nlevs)               int    [READ]
     * 16. start_block         alloc  (nlevs)               int    [skip]
     * 17. end_block           alloc  (nlevs)               int    [READ]
     */

    std::vector<int>    raw_cell_idx, raw_cell_blk;
    std::vector<int>    raw_vert_idx, raw_vert_blk;
    std::vector<double> raw_inv_primal, raw_inv_dual, raw_tangent;
    std::vector<int>    raw_end_index, raw_end_block;

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
            read_alloc_array<int>(f);  /* skip */
    /* 7 */ expect_line(f, "# quad_blk");
            read_alloc_array<int>(f);  /* skip */
    /* 8 */ expect_line(f, "# inv_primal_edge_length");
            m_invp = read_alloc_array<double>(f, &raw_inv_primal);
    /* 9 */ expect_line(f, "# inv_dual_edge_length");
            m_invd = read_alloc_array<double>(f, &raw_inv_dual);
    /*10 */ expect_line(f, "# area_edge");
            read_alloc_array<double>(f);  /* skip */
    /*11 */ expect_line(f, "# f_e");
            read_alloc_array<double>(f);  /* skip */
    /*12 */ expect_line(f, "# fn_e");
            read_alloc_array<double>(f);  /* skip */
    /*13 */ expect_line(f, "# ft_e");
            read_alloc_array<double>(f);  /* skip */
    /*14 */ expect_line(f, "# start_index");
            read_alloc_array<int>(f);  /* skip */
    /*15 */ expect_line(f, "# end_index");
            read_alloc_array<int>(f, &raw_end_index);
    /*16 */ expect_line(f, "# start_block");
            read_alloc_array<int>(f);  /* skip */
    /*17 */ expect_line(f, "# end_block");
            read_alloc_array<int>(f, &raw_end_block);

    /* ---- derive sizes ---- */
    ed.nproma        = m_cidx.size[0];
    ed.n_edges_alloc = ed.nproma * ed.nblks_e;
    ed.n_cells       = ed.nproma * ed.nblks_c;
    ed.n_verts       = ed.nproma * ed.nblks_v;

    /* ---- compute actual valid edge count from end_block / end_index ----
     *
     * end_block and end_index are 1-D arrays indexed by refin_ctrl level.
     * The LAST entry corresponds to the outermost halo (or the final
     * valid category for a single-rank run).  The actual number of valid
     * edges is:
     *   n_valid = (max_end_block - 1) * nproma + max_end_index
     *
     * We take the MAX across all entries to handle empty categories
     * (end_index=0 / end_block=0).
     */
    {
        int max_flat = 0;
        int n_entries = (int)raw_end_block.size();
        fprintf(stderr, "[icon_data_loader] end_block/end_index have %d "
                "refin_ctrl entries\n", n_entries);
        for (int i = 0; i < n_entries; i++) {
            if (raw_end_block[i] <= 0 || raw_end_index[i] <= 0) continue;
            /* n_valid up to and including entry (end_block[i], end_index[i]):
               = (end_block[i] - 1) * nproma + end_index[i]
               because end_block is 1-based and end_index is 1-based
               (end_index is within-block, so end_index=nproma means
                the entire block is valid). */
            int flat = (raw_end_block[i] - 1) * ed.nproma + raw_end_index[i];
            if (flat > max_flat) max_flat = flat;
        }
        ed.n_edges_valid = max_flat;

        if (ed.n_edges_valid <= 0) {
            fprintf(stderr, "[icon_data_loader] WARNING: n_edges_valid=%d, "
                    "falling back to n_edges_alloc=%d\n",
                    ed.n_edges_valid, ed.n_edges_alloc);
            ed.n_edges_valid = ed.n_edges_alloc;
        }
        if (ed.n_edges_valid > ed.n_edges_alloc) {
            fprintf(stderr, "[icon_data_loader] WARNING: n_edges_valid=%d > "
                    "n_edges_alloc=%d, clamping\n",
                    ed.n_edges_valid, ed.n_edges_alloc);
            ed.n_edges_valid = ed.n_edges_alloc;
        }
    }

    fprintf(stderr, "[icon_data_loader] nproma=%d  nblks_e=%d  "
            "n_edges_valid=%d  n_edges_alloc=%d  n_cells=%d  n_verts=%d\n",
            ed.nproma, ed.nblks_e,
            ed.n_edges_valid, ed.n_edges_alloc,
            ed.n_cells, ed.n_verts);

    /* ---- linearise connectivity (only valid edges) ---- */
    linearise_connectivity(raw_cell_idx, m_cidx, raw_cell_blk,
                           ed.nproma, ed.n_cells, ed.n_edges_valid,
                           ed.cell_idx);

    linearise_connectivity(raw_vert_idx, m_vidx, raw_vert_blk,
                           ed.nproma, ed.n_verts, ed.n_edges_valid,
                           ed.vert_idx);

    /* ---- linearise 1-D geometry (only valid edges) ---- */
    linearise_2d_double(raw_inv_dual,   ed.n_edges_valid, ed.inv_dual);
    linearise_2d_double(raw_inv_primal, ed.n_edges_valid, ed.inv_primal);
    linearise_2d_double(raw_tangent,    ed.n_edges_valid, ed.tangent_o);

    fprintf(stderr, "[icon_data_loader] done.  Sample cell_idx[0]={%d,%d}  "
            "vert_idx[0]={%d,%d}\n",
            ed.cell_idx[0], ed.cell_idx[1],
            ed.vert_idx[0], ed.vert_idx[1]);

    return true;
}

/* ------------------------------------------------------------------ */
/*  Convenience: build path from ICON_DATA_PATH env + timestep        */
/*                                                                    */
/*  ICON_DATA_PATH  defaults to:                                      */
/*    /home/primrose/Work/icon-artifacts/velocity/data_r02b05         */
/*  Timestep defaults to 9.                                           */
/* ------------------------------------------------------------------ */
inline std::string icon_patch_path(int timestep = 9) {
    const char* env = std::getenv("ICON_DATA_PATH");
    std::string dir = env ? env
                          : "/home/primrose/Work/icon-artifacts/velocity/data_r02b05";
    /* strip trailing slash */
    while (!dir.empty() && dir.back() == '/') dir.pop_back();
    return dir + "/p_patch." + std::to_string(timestep) + ".data";
}

#endif /* ICON_DATA_LOADER_H */