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
  std::vector<double> inv_dual, inv_primal, tangent_o;
  void free_all() {
    cell_idx.clear(); cell_idx.shrink_to_fit();
    vert_idx.clear(); vert_idx.shrink_to_fit();
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
  expect_line(f,"# cells"); skip_grid_cells(f);
  expect_line(f,"# edges");
  std::vector<int> rc_idx,rc_blk,rv_idx,rv_blk; std::vector<int> re_idx,re_blk;
  std::vector<double> r_invp,r_invd,r_tang;
  ArrayMeta m_ci,m_cb,m_vi,m_tg,m_ip,m_id;
  expect_line(f,"# cell_idx");               m_ci=read_alloc_array<int>(f,&rc_idx);
  expect_line(f,"# cell_blk");               m_cb=read_alloc_array<int>(f,&rc_blk);
  expect_line(f,"# vertex_idx");             m_vi=read_alloc_array<int>(f,&rv_idx);
  expect_line(f,"# vertex_blk");             read_alloc_array<int>(f,&rv_blk);
  expect_line(f,"# tangent_orientation");    m_tg=read_alloc_array<double>(f,&r_tang);
  expect_line(f,"# quad_idx");               read_alloc_array<int>(f);
  expect_line(f,"# quad_blk");               read_alloc_array<int>(f);
  expect_line(f,"# inv_primal_edge_length"); m_ip=read_alloc_array<double>(f,&r_invp);
  expect_line(f,"# inv_dual_edge_length");   m_id=read_alloc_array<double>(f,&r_invd);
  expect_line(f,"# area_edge");              read_alloc_array<double>(f);
  expect_line(f,"# f_e");                    read_alloc_array<double>(f);
  expect_line(f,"# fn_e");                   read_alloc_array<double>(f);
  expect_line(f,"# ft_e");                   read_alloc_array<double>(f);
  expect_line(f,"# start_index");            read_alloc_array<int>(f);
  expect_line(f,"# end_index");              read_alloc_array<int>(f,&re_idx);
  expect_line(f,"# start_block");            read_alloc_array<int>(f);
  expect_line(f,"# end_block");              read_alloc_array<int>(f,&re_blk);
  {int s0=m_ci.size[0]; if(s0!=nproma) fprintf(stderr,"[icon_data_loader] WARNING: serde size[0]=%d != nproma=%d\n",s0,nproma);}
  {int mx=0; for(int i=0;i<(int)re_blk.size();i++){
    if(re_blk[i]<=0||re_idx[i]<=0) continue;
    int fl=(re_blk[i]-1)*nproma+re_idx[i]; if(fl>mx) mx=fl;}
    ed.n_edges_valid=mx; if(ed.n_edges_valid<=0) ed.n_edges_valid=ed.n_edges;
    if(ed.n_edges_valid>ed.n_edges) ed.n_edges_valid=ed.n_edges;}
  fprintf(stderr,"[icon_data_loader] nproma=%d  nblks_e=%d  nblks_c=%d  nblks_v=%d\n"
    "                   n_edges=%d (valid=%d)  n_cells=%d  n_verts=%d\n",
    ed.nproma,ed.nblks_e,ed.nblks_c,ed.nblks_v,ed.n_edges,ed.n_edges_valid,ed.n_cells,ed.n_verts);
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
inline std::string icon_data_dir() {
  const char*env=std::getenv("ICON_DATA_PATH");
  std::string dir=env?env:"/home/primrose/Work/icon-artifacts/velocity/data_r02b05";
  while(!dir.empty()&&dir.back()=='/') dir.pop_back();
  return dir;
}
inline std::string icon_patch_path(int ts=9)  {return icon_data_dir()+"/p_patch."+std::to_string(ts)+".data";}
inline std::string icon_global_path(int ts=9) {return icon_data_dir()+"/global_data.t0."+std::to_string(ts)+".data";}

/* ---- Padding & variant helpers ---- */
inline int icon_pad_nlev(int nlev, int align=32) { return ((nlev+align-1)/align)*align; }

/*
 * Variant mapping:
 *   V1 = je-first, basic IN    (kernel V=1)
 *   V2 = je-first, optimised IN (kernel V=2)
 *   V3 = jk-first, basic IN    (kernel V=3)
 *   V4 = jk-first, optimised IN (kernel V=4)
 *   V5 = jk-first, optimised IN, padded nlev (kernel V=4, nlev=pad(nlev_end))
 */
inline int kern_v(int V) { return (V==5) ? 4 : V; }

#endif /* ICON_DATA_LOADER_H */