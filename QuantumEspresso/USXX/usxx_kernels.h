#pragma once
#include "types.h"

// ============================================================
// Fortran column-major indexing helpers
// ============================================================
#define IDX2(row, col, nrows) ((col) * (nrows) + (row))
#define IDX3(i, j, k, n1, n2) ((k) * (n1) * (n2) + (j) * (n1) + (i))

#define EIGTS1_IDX(ig, na) ((na) * 217 + ((ig) + 108))
#define EIGTS2_IDX(ig, na) ((na) * 109 + ((ig) + 54))
#define EIGTS3_IDX(ig, na) ((na) * 109 + ((ig) + 54))

#define EIGTS1T_IDX(na, ig) (((ig) + 108) * 10 + (na))
#define EIGTS2T_IDX(na, ig) (((ig) + 54) * 10 + (na))
#define EIGTS3T_IDX(na, ig) (((ig) + 54) * 10 + (na))

// ============================================================
// CPU baseline
// ============================================================
void addusxx_g_cpu(
    Complex_DP* rhoc, const DP* xkq, const DP* xk, const DP* tau,
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const Complex_DP* qgm, const Complex_DP* eigts1,
    const Complex_DP* eigts2, const Complex_DP* eigts3,
    const int* mill, const int* dfftt__nl);

// ============================================================
// GPU baseline AoS — faithful to Fortran: writes per atom
// Scratch: d_eigqts (nat)
// ============================================================
void addusxx_g_gpu(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1,
    const Complex_DP* d_eigts2, const Complex_DP* d_eigts3,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    Complex_DP* d_eigqts);

// ============================================================
// GPU eigts-transposed AoS — natural g-vector order,
// original qgm (coalesced), transposed eigts, register accum.
// Scratch: d_eigqts (nat)
// ============================================================
void addusxx_g_gpu_eigts_transposed(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    Complex_DP* d_eigqts);

// ============================================================
// GPU shared-bec AoS — eigts-transposed + becphi/becpsi in
// shared memory. Same signature as eigts_transposed.
// Scratch: d_eigqts (nat)
// ============================================================
void addusxx_g_gpu_shared_bec(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    Complex_DP* d_eigqts);

// ============================================================
// GPU optimized AoS — fused compute+scatter, sorted writes,
// transposed layouts. No intermediate buffer.
// Scratch: d_eigqts (nat)
// ============================================================
void addusxx_g_gpu_optimized(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const Complex_DP* d_qgm_T, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl_sorted, const int* d_dfftt__nl_ix,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    Complex_DP* d_eigqts);

// ============================================================
// GPU baseline SoA — faithful to Fortran: writes per atom
// Scratch: d_eigqts_re/im (nat each)
// ============================================================
void addusxx_g_gpu_soa(
    double* d_rhoc_re, double* d_rhoc_im,
    const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_becphi_re, const double* d_becphi_im,
    const double* d_becpsi_re, const double* d_becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const double* d_qgm_re, const double* d_qgm_im,
    const double* d_eigts1_re, const double* d_eigts1_im,
    const double* d_eigts2_re, const double* d_eigts2_im,
    const double* d_eigts3_re, const double* d_eigts3_im,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    double* d_eigqts_re, double* d_eigqts_im);

// ============================================================
// GPU eigts-transposed SoA — natural order, original qgm,
// transposed eigts, register accum.
// Scratch: d_eigqts_re/im (nat each)
// ============================================================
void addusxx_g_gpu_eigts_transposed_soa(
    double* d_rhoc_re, double* d_rhoc_im,
    const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_becphi_re, const double* d_becphi_im,
    const double* d_becpsi_re, const double* d_becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const double* d_qgm_re, const double* d_qgm_im,
    const double* d_eigts1_T_re, const double* d_eigts1_T_im,
    const double* d_eigts2_T_re, const double* d_eigts2_T_im,
    const double* d_eigts3_T_re, const double* d_eigts3_T_im,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    double* d_eigqts_re, double* d_eigqts_im);

// ============================================================
// GPU shared-bec SoA — eigts-transposed + shared memory bec
// Scratch: d_eigqts_re/im (nat each)
// ============================================================
void addusxx_g_gpu_shared_bec_soa(
    double* d_rhoc_re, double* d_rhoc_im,
    const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_becphi_re, const double* d_becphi_im,
    const double* d_becpsi_re, const double* d_becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const double* d_qgm_re, const double* d_qgm_im,
    const double* d_eigts1_T_re, const double* d_eigts1_T_im,
    const double* d_eigts2_T_re, const double* d_eigts2_T_im,
    const double* d_eigts3_T_re, const double* d_eigts3_T_im,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    double* d_eigqts_re, double* d_eigqts_im);

// ============================================================
// GPU optimized SoA — fused compute+scatter, sorted writes
// Scratch: d_eigqts_re/im (nat each)
// ============================================================
void addusxx_g_gpu_optimized_soa(
    double* d_rhoc_re, double* d_rhoc_im,
    const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_becphi_re, const double* d_becphi_im,
    const double* d_becpsi_re, const double* d_becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_upf_tvanp, const int* d_nij_type, const int* d_ityp,
    const int* d_ofsbeta, const int* d_nh, const int* d_ijtoh,
    const double* d_qgm_T_re, const double* d_qgm_T_im,
    const double* d_eigts1_T_re, const double* d_eigts1_T_im,
    const double* d_eigts2_T_re, const double* d_eigts2_T_im,
    const double* d_eigts3_T_re, const double* d_eigts3_T_im,
    const int* d_mill, const int* d_dfftt__nl_sorted, const int* d_dfftt__nl_ix,
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh,
    int tblock_size, int coarsen,
    double* d_eigqts_re, double* d_eigqts_im);

// ============================================================
// Utility
// ============================================================
inline void aos_to_soa(const Complex_DP* aos, double* re, double* im, int n) {
    for (int i = 0; i < n; i++) { re[i] = creal_val(aos[i]); im[i] = cimag_val(aos[i]); }
}
inline void soa_to_aos(const double* re, const double* im, Complex_DP* aos, int n) {
    for (int i = 0; i < n; i++) { aos[i] = make_cmplx(re[i], im[i]); }
}