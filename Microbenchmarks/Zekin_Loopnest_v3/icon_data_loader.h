#ifndef ICON_DATA_LOADER_H
#define ICON_DATA_LOADER_H
/*
 * icon_data_loader.h -- Standalone reader for ICON's serialised p_patch files.
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
/*  Minimal serde helpers                                             */
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
inline int read_int(std::istream &s) { scroll_space(s); int v=0; s>>v; return v; }
inline double read_double(std::istream &s) { scroll_space(s); long double v=0; s>>v; return (double)v; }
inline bool read_bool(std::istream &s) { scroll_space(s); char c='0'; s>>c; return c=='1'; }

struct ArrayMeta {
  int rank=0, size[4]={}, lbound[4]={};
  int volume() const { int v=1; for(int i=0;i<rank;i++) v*=size[i]; return v; }
};
inline ArrayMeta read_meta(std::istream &s) {
  ArrayMeta m;
  expect_line(s,"# rank"); m.rank=read_int(s);
  assert(m.rank>=1 && m.rank<=4);
  expect_line(s,"# size");   for(int i=0;i<m.rank;i++) m.size[i]=read_int(s);
  expect_line(s,"# lbound"); for(int i=0;i<m.rank;i++) m.lbound[i]=read_int(s);
  return m;
}
inline void skip_entries(std::istream &s, int count) {
  std::string tok; for(int i=0;i<count;i++){scroll_space(s);s>>tok;}
}

/* Print array dimensions */
inline void print_meta(const char *section, const char *name, const ArrayMeta &m, const char *dtype) {
  fprintf(stderr, "[icon_sizes] %s/%-28s  %s  rank=%d  (", section, name, dtype, m.rank);
  for (int i=0;i<m.rank;i++) fprintf(stderr, "%s%d", i?", ":"", m.size[i]);
  long long bytes = (long long)m.volume() * (strcmp(dtype,"int")==0 ? 4 : 8);
  fprintf(stderr, ")  elements=%d  %.1f MB\n", m.volume(), bytes/1e6);
}

template<typename T>
ArrayMeta read_alloc_array(std::istream &s, std::vector<T> *out=nullptr) {
  expect_line(s,"# alloc"); bool alloc=read_bool(s);
  if(!alloc) return {};
  ArrayMeta m=read_meta(s); expect_line(s,"# entries"); int vol=m.volume();
  if(out){ out->resize(vol);
    for(int i=0;i<vol;i++){scroll_space(s);
      if constexpr(std::is_same_v<T,int>) s>>(*out)[i];
      else{long double tmp;s>>tmp;(*out)[i]=(T)tmp;}
    }
  } else skip_entries(s,vol);
  return m;
}
template<typename T>
ArrayMeta read_assoc_array(std::istream &s, std::vector<T> *out=nullptr) {
  expect_line(s,"# assoc"); bool assoc=read_bool(s);
  if(!assoc) return {};
  expect_line(s,"# missing"); read_int(s);
  ArrayMeta m=read_meta(s); expect_line(s,"# entries"); int vol=m.volume();
  if(out){ out->resize(vol);
    for(int i=0;i<vol;i++){scroll_space(s);
      if constexpr(std::is_same_v<T,int>) s>>(*out)[i];
      else{long double tmp;s>>tmp;(*out)[i]=(T)tmp;}
    }
  } else skip_entries(s,vol);
  return m;
}

inline void skip_decomp_info(std::istream &s) { expect_line(s,"# owner_mask"); read_alloc_array<int>(s); }
inline void skip_grid_cells(std::istream &s) {
  expect_line(s,"# neighbor_idx"); read_alloc_array<int>(s);
  expect_line(s,"# neighbor_blk"); read_alloc_array<int>(s);
  expect_line(s,"# edge_idx");     read_alloc_array<int>(s);
  expect_line(s,"# edge_blk");     read_alloc_array<int>(s);
  expect_line(s,"# area");         read_assoc_array<double>(s);
  expect_line(s,"# start_index");  read_alloc_array<int>(s);
  expect_line(s,"# end_index");    read_alloc_array<int>(s);
  expect_line(s,"# start_block");  read_alloc_array<int>(s);
  expect_line(s,"# end_block");    read_alloc_array<int>(s);
  expect_line(s,"# decomp_info");  skip_decomp_info(s);
}
inline void skip_grid_vertices(std::istream &s) {
  expect_line(s,"# cell_idx");    read_alloc_array<int>(s);
  expect_line(s,"# cell_blk");    read_alloc_array<int>(s);
  expect_line(s,"# edge_idx");    read_alloc_array<int>(s);
  expect_line(s,"# edge_blk");    read_alloc_array<int>(s);
  expect_line(s,"# start_index"); read_alloc_array<int>(s);
  expect_line(s,"# end_index");   read_alloc_array<int>(s);
  expect_line(s,"# start_block"); read_alloc_array<int>(s);
  expect_line(s,"# end_block");   read_alloc_array<int>(s);
}

} // namespace icon_serde

/* ------------------------------------------------------------------ */
struct IconEdgeData {
  int nproma=0, nblks_c=0, nblks_e=0, nblks_v=0;
  int n_edges=0, n_cells=0, n_verts=0, n_edges_valid=0;
  std::vector<int> cell_idx, vert_idx;
  std::vector<int> start_idx, end_idx, start_blk, end_blk;
  std::vector<double> inv_dual, inv_primal, tangent_o;
  void free_all() {
    cell_idx.clear(); cell_idx.shrink_to_fit();
    vert_idx.clear(); vert_idx.shrink_to_fit();
    start_idx.clear(); start_idx.shrink_to_fit();
    end_idx.clear(); end_idx.shrink_to_fit();
    start_blk.clear(); start_blk.shrink_to_fit();
    end_blk.clear(); end_blk.shrink_to_fit();
    inv_dual.clear(); inv_dual.shrink_to_fit();
    inv_primal.clear(); inv_primal.shrink_to_fit();
    tangent_o.clear(); tangent_o.shrink_to_fit();
  }
};

/* ------------------------------------------------------------------ */
static void linearise_connectivity(
    const std::vector<int>&raw_idx, const icon_serde::ArrayMeta&meta_idx,
    const std::vector<int>&raw_blk, int nproma,
    int target_max, int n_edges, std::vector<int>&out) {
  int s0=meta_idx.size[0], s1=meta_idx.size[1];
  out.assign(n_edges*2, 0);
  for(int n=0;n<2;n++) for(int jb=0;jb<s1;jb++) for(int jc=0;jc<nproma;jc++){
    int el=jb*nproma+jc; if(el>=n_edges) continue;
    int sf=jc+jb*s0+n*s0*s1;
    int iv=raw_idx[sf], bv=raw_blk[sf];
    int tf=(bv-1)*nproma+(iv-1);
    if(iv<=0||bv<=0||tf<0||tf>=target_max) tf=0;
    out[el*2+n]=tf;
  }
}
static void linearise_2d_double(const std::vector<double>&raw,
    int s0_serde, int nproma, int nblks, std::vector<double>&out) {
  out.assign(nproma*nblks, 0.0);
  for(int jb=0;jb<nblks;jb++) for(int jc=0;jc<nproma;jc++)
    out[jb*nproma+jc]=raw[jc+jb*s0_serde];
}

/* ------------------------------------------------------------------ */
inline int icon_read_nproma(const char *global_path) {
  using namespace icon_serde;
  std::ifstream f(global_path);
  if(!f.is_open()){fprintf(stderr,"[icon_data_loader] cannot open '%s'\n",global_path);return -1;}
  expect_line(f,"# nflatlev");
  {ArrayMeta m=read_meta(f);expect_line(f,"# entries");skip_entries(f,m.volume());}
  expect_line(f,"# i_am_accel_node"); read_int(f);
  expect_line(f,"# lextra_diffu");    read_int(f);
  expect_line(f,"# nproma");
  int nproma=read_int(f);
  fprintf(stderr,"[icon_data_loader] global nproma = %d\n",nproma);
  return nproma;
}

inline bool icon_load_patch(const char *path, int nproma, IconEdgeData &ed) {
  using namespace icon_serde;
  if(nproma<=0){fprintf(stderr,"[icon_data_loader] invalid nproma=%d\n",nproma);return false;}
  std::ifstream f(path);
  if(!f.is_open()){fprintf(stderr,"[icon_data_loader] cannot open '%s'\n",path);return false;}
  fprintf(stderr,"[icon_data_loader] reading %s  (nproma=%d)\n",path,nproma);
  ed.nproma=nproma;
  expect_line(f,"# nblks_c"); ed.nblks_c=read_int(f);
  expect_line(f,"# nblks_e"); ed.nblks_e=read_int(f);
  expect_line(f,"# nblks_v"); ed.nblks_v=read_int(f);
  ed.n_edges=nproma*ed.nblks_e; ed.n_cells=nproma*ed.nblks_c; ed.n_verts=nproma*ed.nblks_v;

  fprintf(stderr,"[icon_sizes] Grid topology:\n");
  fprintf(stderr,"[icon_sizes]   nproma=%d  nblks_c=%d  nblks_e=%d  nblks_v=%d\n",
          nproma, ed.nblks_c, ed.nblks_e, ed.nblks_v);
  fprintf(stderr,"[icon_sizes]   n_edges = nproma*nblks_e = %d\n", ed.n_edges);
  fprintf(stderr,"[icon_sizes]   n_cells = nproma*nblks_c = %d\n", ed.n_cells);
  fprintf(stderr,"[icon_sizes]   n_verts = nproma*nblks_v = %d\n", ed.n_verts);

  expect_line(f,"# cells"); skip_grid_cells(f);
  expect_line(f,"# edges");

  std::vector<int> rc_idx,rc_blk,rv_idx,rv_blk;
  std::vector<double> r_invp,r_invd,r_tang;
  ArrayMeta m_ci,m_cb,m_vi,m_vb,m_tg,m_qi,m_qb,m_ip,m_id,m_ae,m_fe,m_fn,m_ft;
  ArrayMeta m_si,m_ei,m_sb,m_eb;

  fprintf(stderr,"\n[icon_sizes] Edge arrays (serde dimensions):\n");

  expect_line(f,"# cell_idx");               m_ci=read_alloc_array<int>(f,&rc_idx);
  print_meta("edges","cell_idx",m_ci,"int");
  expect_line(f,"# cell_blk");               m_cb=read_alloc_array<int>(f,&rc_blk);
  print_meta("edges","cell_blk",m_cb,"int");
  expect_line(f,"# vertex_idx");             m_vi=read_alloc_array<int>(f,&rv_idx);
  print_meta("edges","vertex_idx",m_vi,"int");
  expect_line(f,"# vertex_blk");             m_vb=read_alloc_array<int>(f,&rv_blk);
  print_meta("edges","vertex_blk",m_vb,"int");
  expect_line(f,"# tangent_orientation");    m_tg=read_alloc_array<double>(f,&r_tang);
  print_meta("edges","tangent_orientation",m_tg,"dbl");
  expect_line(f,"# quad_idx");               m_qi=read_alloc_array<int>(f);
  print_meta("edges","quad_idx",m_qi,"int");
  expect_line(f,"# quad_blk");               m_qb=read_alloc_array<int>(f);
  print_meta("edges","quad_blk",m_qb,"int");
  expect_line(f,"# inv_primal_edge_length"); m_ip=read_alloc_array<double>(f,&r_invp);
  print_meta("edges","inv_primal_edge_length",m_ip,"dbl");
  expect_line(f,"# inv_dual_edge_length");   m_id=read_alloc_array<double>(f,&r_invd);
  print_meta("edges","inv_dual_edge_length",m_id,"dbl");
  expect_line(f,"# area_edge");              m_ae=read_alloc_array<double>(f);
  print_meta("edges","area_edge",m_ae,"dbl");
  expect_line(f,"# f_e");                    m_fe=read_alloc_array<double>(f);
  print_meta("edges","f_e",m_fe,"dbl");
  expect_line(f,"# fn_e");                   m_fn=read_alloc_array<double>(f);
  print_meta("edges","fn_e",m_fn,"dbl");
  expect_line(f,"# ft_e");                   m_ft=read_alloc_array<double>(f);
  print_meta("edges","ft_e",m_ft,"dbl");
  expect_line(f,"# start_index");            m_si=read_alloc_array<int>(f,&ed.start_idx);
  print_meta("edges","start_index",m_si,"int");
  expect_line(f,"# end_index");              m_ei=read_alloc_array<int>(f,&ed.end_idx);
  print_meta("edges","end_index",m_ei,"int");
  expect_line(f,"# start_block");            m_sb=read_alloc_array<int>(f,&ed.start_blk);
  print_meta("edges","start_block",m_sb,"int");
  expect_line(f,"# end_block");              m_eb=read_alloc_array<int>(f,&ed.end_blk);
  print_meta("edges","end_block",m_eb,"int");

  /* linearised / collapsed dimensions */
  fprintf(stderr,"\n[icon_sizes] Linearised arrays used in stencil:\n");
  fprintf(stderr,"[icon_sizes]   cell_idx     [n_edges*2]     = [%d]       (%d -> [0,%d))\n",
          ed.n_edges*2, ed.n_edges, ed.n_cells);
  fprintf(stderr,"[icon_sizes]   vert_idx     [n_edges*2]     = [%d]       (%d -> [0,%d))\n",
          ed.n_edges*2, ed.n_edges, ed.n_verts);
  fprintf(stderr,"[icon_sizes]   inv_dual     [n_edges]       = [%d]\n", ed.n_edges);
  fprintf(stderr,"[icon_sizes]   inv_primal   [n_edges]       = [%d]\n", ed.n_edges);
  fprintf(stderr,"[icon_sizes]   tangent_o    [n_edges]       = [%d]\n", ed.n_edges);
  fprintf(stderr,"[icon_sizes]   vn_ie        [n_edges*nlev]  (runtime, per-level)\n");
  fprintf(stderr,"[icon_sizes]   w            [n_cells*nlev]  (indirect target, n_cells=%d)\n", ed.n_cells);
  fprintf(stderr,"[icon_sizes]   z_w_v        [n_verts*nlev]  (indirect target, n_verts=%d)\n", ed.n_verts);
  fprintf(stderr,"[icon_sizes]   out          [n_edges*nlev]  (output)\n");

  {int s0=m_ci.size[0]; if(s0!=nproma) fprintf(stderr,"[icon_data_loader] WARNING: serde size[0]=%d != nproma=%d\n",s0,nproma);}
  {int mx=0; for(int i=0;i<(int)ed.end_blk.size();i++){
    if(ed.end_blk[i]<=0||ed.end_idx[i]<=0) continue;
    int fl=(ed.end_blk[i]-1)*nproma+ed.end_idx[i]; if(fl>mx) mx=fl;}
    ed.n_edges_valid=mx; if(ed.n_edges_valid<=0) ed.n_edges_valid=ed.n_edges;
    if(ed.n_edges_valid>ed.n_edges) ed.n_edges_valid=ed.n_edges;}

  fprintf(stderr,"\n[icon_sizes] n_edges_valid=%d (%.1f%% of n_edges=%d)\n",
          ed.n_edges_valid, 100.0*ed.n_edges_valid/ed.n_edges, ed.n_edges);

  fprintf(stderr,"\n[icon_sizes] Edge range table (%d levels):\n", (int)ed.end_idx.size());
  fprintf(stderr,"[icon_sizes]   %4s  %14s  %14s  %10s  %8s\n",
          "lvl", "start(blk,idx)", "end(blk,idx)", "linear", "n_edges");
  for (int i=0; i<(int)ed.end_idx.size(); i++) {
    int s_lin = (ed.start_blk[i]>0 && ed.start_idx[i]>0)
              ? (ed.start_blk[i]-1)*nproma + (ed.start_idx[i]-1) : -1;
    int e_lin = (ed.end_blk[i]>0 && ed.end_idx[i]>0)
              ? (ed.end_blk[i]-1)*nproma + (ed.end_idx[i]-1) : -1;
    int count = (s_lin>=0 && e_lin>=0) ? e_lin-s_lin+1 : 0;
    fprintf(stderr,"[icon_sizes]   %4d  (%4d,%6d)  (%4d,%6d)  [%6d,%6d]  %8d\n",
            i, ed.start_blk[i], ed.start_idx[i], ed.end_blk[i], ed.end_idx[i],
            s_lin, e_lin, count);
  }

  linearise_connectivity(rc_idx,m_ci,rc_blk,nproma,ed.n_cells,ed.n_edges,ed.cell_idx);
  linearise_connectivity(rv_idx,m_vi,rv_blk,nproma,ed.n_verts,ed.n_edges,ed.vert_idx);
  linearise_2d_double(r_invd,m_id.size[0],nproma,ed.nblks_e,ed.inv_dual);
  linearise_2d_double(r_invp,m_ip.size[0],nproma,ed.nblks_e,ed.inv_primal);
  linearise_2d_double(r_tang,m_tg.size[0],nproma,ed.nblks_e,ed.tangent_o);

  fprintf(stderr,"[icon_data_loader] done.  Sample cell_idx[0]={%d,%d}  vert_idx[0]={%d,%d}\n",
    ed.cell_idx[0],ed.cell_idx[1],ed.vert_idx[0],ed.vert_idx[1]);
  return true;
}

/* ------------------------------------------------------------------ */
/*  Indirection locality metrics                                      */
/*                                                                     */
/*  Key metrics for cache performance:                                 */
/*    consecutive: |idx[je+1] - idx[je]|  → cache reuse between edges  */
/*    pair:        |ci0 - ci1|            → reuse within one edge      */
/*    unique:      # distinct idx values  → indirect working set       */
/* ------------------------------------------------------------------ */
inline void print_idx_locality(const char *label, const char *idx_name,
                               const int *logical, int N, int idx_max) {
  static constexpr int BUCKETS[] = {1,2,4,8,16,32,64,128,256,1024,4096,16384,65536};
  static constexpr int NB = sizeof(BUCKETS)/sizeof(BUCKETS[0]);

  /* --- consecutive distance: |idx[(je+1)*2+n] - idx[je*2+n]| --- */
  double consec_sum = 0;
  long long consec_same_cl = 0, consec_within_page = 0;
  int consec_max = 0;
  long long consec_hist[NB+1] = {};
  for (int je = 0; je < N-1; je++) {
    for (int n = 0; n < 2; n++) {
      int d = std::abs(logical[(je+1)*2+n] - logical[je*2+n]);
      consec_sum += d;
      if (d < 8) consec_same_cl++;
      if (d < 512) consec_within_page++;
      if (d > consec_max) consec_max = d;
      int b = NB;
      for (int i = 0; i < NB; i++) { if (d < BUCKETS[i]) { b = i; break; } }
      consec_hist[b]++;
    }
  }
  long long ctotal = (long long)(N-1) * 2;

  /* --- pair distance: |ci0 - ci1| for same edge --- */
  double pair_sum = 0;
  long long pair_same_cl = 0;
  int pair_max = 0;
  for (int je = 0; je < N; je++) {
    int d = std::abs(logical[je*2+0] - logical[je*2+1]);
    pair_sum += d;
    if (d < 8) pair_same_cl++;
    if (d > pair_max) pair_max = d;
  }

  /* --- unique count (working set) via bitset --- */
  int unique = 0;
  {
    std::vector<bool> seen(idx_max, false);
    for (int je = 0; je < N; je++)
      for (int n = 0; n < 2; n++) {
        int idx = logical[je*2+n];
        if (idx >= 0 && idx < idx_max && !seen[idx]) { seen[idx] = true; unique++; }
      }
  }

  fprintf(stderr, "[locality] %s %s (N=%d, idx_max=%d, unique=%d = %.1f%%):\n",
          label, idx_name, N, idx_max, unique, 100.0*unique/idx_max);
  fprintf(stderr, "[locality]   consecutive |idx[je+1]-idx[je]|: avg=%5.1f  max=%d  same_cl(<8)=%5.1f%%  within_page(<512)=%5.1f%%\n",
          consec_sum/ctotal, consec_max, 100.0*consec_same_cl/ctotal, 100.0*consec_within_page/ctotal);
  fprintf(stderr, "[locality]   pair |ci0-ci1|:                   avg=%5.1f  max=%d  same_cl(<8)=%5.1f%%\n",
          pair_sum/N, pair_max, 100.0*pair_same_cl/N);
  fprintf(stderr, "[locality]   consecutive histogram:\n");
  fprintf(stderr, "[locality]     ");
  long long cum = 0;
  for (int i = 0; i < NB; i++) {
    cum += consec_hist[i];
    fprintf(stderr, "<%d:%.0f%% ", BUCKETS[i], 100.0*cum/ctotal);
  }
  fprintf(stderr, "\n");
}

inline void icon_print_locality_metrics(const IconEdgeData &ed) {
  fprintf(stderr, "\n");
  print_idx_locality("exact", "cell_idx", ed.cell_idx.data(), ed.n_edges, ed.n_cells);
  print_idx_locality("exact", "vert_idx", ed.vert_idx.data(), ed.n_edges, ed.n_verts);
  fprintf(stderr, "\n");
}

inline void print_synthetic_locality(const char *dist_label,
                                     const int *cell_logical,
                                     const int *vert_logical,
                                     int N, int cell_max, int vert_max) {
  print_idx_locality(dist_label, "cell_idx", cell_logical, N, cell_max);
  if (vert_logical)
    print_idx_locality(dist_label, "vert_idx", vert_logical, N, vert_max);
}

/* ------------------------------------------------------------------ */
inline std::string icon_data_dir() {
  const char*env=std::getenv("ICON_DATA_PATH");
  std::string dir=env?env:"/home/primrose/Work/icon-artifacts/velocity/data_r02b05";
  while(!dir.empty()&&dir.back()=='/') dir.pop_back();
  return dir;
}
inline std::string icon_patch_path(int ts=9)  {return icon_data_dir()+"/p_patch."+std::to_string(ts)+".data";}
inline std::string icon_global_path(int ts=9) {return icon_data_dir()+"/global_data.t0."+std::to_string(ts)+".data";}

inline int icon_pad_nlev(int nlev, int align=32) { return nlev; }

/*
 * Variant mapping:
 *   V1 = je-first, basic IN    (kernel V=1)
 *   V2 = je-first, optimised IN (kernel V=2)
 *   V3 = jk-first, basic IN    (kernel V=3)
 *   V4 = jk-first, optimised IN (kernel V=4)
 *   V5 = jk-first, optimised IN, padded nlev (kernel V=4, nlev=pad(nlev_end))
 *   V6 = flat 1D, jk-first, padded nlev      (kernel V=4 layout)
 *   V7 = flat 1D with TX, jk-first, unpadded  (kernel V=4 layout)
 */
inline int kern_v(int V) { return (V >= 5) ? 4 : V; }

#endif /* ICON_DATA_LOADER_H */