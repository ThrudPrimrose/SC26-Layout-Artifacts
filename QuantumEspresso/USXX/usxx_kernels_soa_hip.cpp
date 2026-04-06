#include "hip/hip_runtime.h"
#include "usxx_kernels.h"
#include <cstdio>
#include <cmath>

#define SOA_MUL_RE(ar, ai, br, bi) ((ar)*(br) - (ai)*(bi))
#define SOA_MUL_IM(ar, ai, br, bi) ((ar)*(bi) + (ai)*(br))
#define SOA_MULCONJ_RE(ar, ai, br, bi) ((ar)*(br) + (ai)*(bi))
#define SOA_MULCONJ_IM(ar, ai, br, bi) ((ai)*(br) - (ar)*(bi))

static constexpr int QGM_NROWS  = 55191;
static constexpr int QGMT_NROWS = 397;
static constexpr int IJTOH_N1   = 19;
static constexpr int IJTOH_N2   = 19;

__global__ void kernel_eigqts_soa(
    double* eigqts_re, double* eigqts_im,
    const double* xkq, const double* xk, const double* tau, int nat)
{
    int na = blockIdx.x * blockDim.x + threadIdx.x;
    if (na >= nat) return;
    double arg = TPI * ((xk[0] - xkq[0]) * tau[na * 3 + 0]
                      + (xk[1] - xkq[1]) * tau[na * 3 + 1]
                      + (xk[2] - xkq[2]) * tau[na * 3 + 2]);
    eigqts_re[na] = cos(arg);
    eigqts_im[na] = -sin(arg);
}

__device__ inline void compute_sf_soa(
    double& sr, double& si, int na,
    const double* eigqts_re, const double* eigqts_im,
    const double* e1_re, const double* e1_im,
    const double* e2_re, const double* e2_im,
    const double* e3_re, const double* e3_im,
    int e1_idx, int e2_idx, int e3_idx)
{
    sr = eigqts_re[na]; si = eigqts_im[na];
    double tr, ti;
    tr = sr; ti = si;
    sr = SOA_MUL_RE(tr, ti, e1_re[e1_idx], e1_im[e1_idx]);
    si = SOA_MUL_IM(tr, ti, e1_re[e1_idx], e1_im[e1_idx]);
    tr = sr; ti = si;
    sr = SOA_MUL_RE(tr, ti, e2_re[e2_idx], e2_im[e2_idx]);
    si = SOA_MUL_IM(tr, ti, e2_re[e2_idx], e2_im[e2_idx]);
    tr = sr; ti = si;
    sr = SOA_MUL_RE(tr, ti, e3_re[e3_idx], e3_im[e3_idx]);
    si = SOA_MUL_IM(tr, ti, e3_re[e3_idx], e3_im[e3_idx]);
}

// ============================================================
// GPU baseline SoA — faithful to Fortran: writes per atom
// ============================================================
__global__ void kernel_addusxx_baseline_soa(
    double* __restrict__ rhoc_re, double* __restrict__ rhoc_im,
    const double* __restrict__ becphi_re, const double* __restrict__ becphi_im,
    const double* __restrict__ becpsi_re, const double* __restrict__ becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const double* __restrict__ qgm_re, const double* __restrict__ qgm_im,
    const double* __restrict__ eigts1_re, const double* __restrict__ eigts1_im,
    const double* __restrict__ eigts2_re, const double* __restrict__ eigts2_im,
    const double* __restrict__ eigts3_re, const double* __restrict__ eigts3_im,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const double* __restrict__ eigqts_re, const double* __restrict__ eigqts_im,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];
        int nl_idx = dfftt__nl[gi] - 1;

        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            int ijkb0 = ofsbeta[na];
            double a2r = 0.0, a2i = 0.0;

            for (int ih = 0; ih < nh_nt; ih++) {
                double a1r = 0.0, a1i = 0.0;
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), IJTOH_N1, IJTOH_N2)];
                    int idx = (nij + ijtoh_val - 1) * QGM_NROWS + gi;
                    double qr = qgm_re[idx], qi = qgm_im[idx];
                    double br = becpsi_re[ijkb0 + jh], bi = becpsi_im[ijkb0 + jh];
                    a1r += SOA_MUL_RE(qr, qi, br, bi);
                    a1i += SOA_MUL_IM(qr, qi, br, bi);
                }
                double pr = becphi_re[ijkb0 + ih], pi = becphi_im[ijkb0 + ih];
                a2r += SOA_MULCONJ_RE(a1r, a1i, pr, pi);
                a2i += SOA_MULCONJ_IM(a1r, a1i, pr, pi);
            }

            double sr, si;
            compute_sf_soa(sr, si, na,
                           eigqts_re, eigqts_im,
                           eigts1_re, eigts1_im,
                           eigts2_re, eigts2_im,
                           eigts3_re, eigts3_im,
                           EIGTS1_IDX(m1, na), EIGTS2_IDX(m2, na), EIGTS3_IDX(m3, na));

            double rr = SOA_MUL_RE(a2r, a2i, sr, si);
            double ri = SOA_MUL_IM(a2r, a2i, sr, si);

            // Write per atom — matches Fortran
            rhoc_re[nl_idx] += rr;
            rhoc_im[nl_idx] += ri;
        }
    }
}

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
    double* d_eigqts_re, double* d_eigqts_im)
{
    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_baseline_soa<<<blocks, tblock_size>>>(
            d_rhoc_re, d_rhoc_im,
            d_becphi_re, d_becphi_im,
            d_becpsi_re, d_becpsi_im,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_re, d_qgm_im,
            d_eigts1_re, d_eigts1_im,
            d_eigts2_re, d_eigts2_im,
            d_eigts3_re, d_eigts3_im,
            d_mill, d_dfftt__nl,
            d_eigqts_re, d_eigqts_im,
            coarsen);
    }
}

// ============================================================
// GPU eigts-transposed SoA — natural order, original qgm,
// transposed eigts, register accumulation
// ============================================================
__global__ void kernel_addusxx_eigts_transposed_soa(
    double* __restrict__ rhoc_re, double* __restrict__ rhoc_im,
    const double* __restrict__ becphi_re, const double* __restrict__ becphi_im,
    const double* __restrict__ becpsi_re, const double* __restrict__ becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const double* __restrict__ qgm_re, const double* __restrict__ qgm_im,
    const double* __restrict__ eigts1_T_re, const double* __restrict__ eigts1_T_im,
    const double* __restrict__ eigts2_T_re, const double* __restrict__ eigts2_T_im,
    const double* __restrict__ eigts3_T_re, const double* __restrict__ eigts3_T_im,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const double* __restrict__ eigqts_re, const double* __restrict__ eigqts_im,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];

        double acc_r = 0.0, acc_i = 0.0;

        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            int ijkb0 = ofsbeta[na];
            double a2r = 0.0, a2i = 0.0;

            for (int ih = 0; ih < nh_nt; ih++) {
                double a1r = 0.0, a1i = 0.0;
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), IJTOH_N1, IJTOH_N2)];
                    // Original qgm layout — coalesced
                    int idx = (nij + ijtoh_val - 1) * QGM_NROWS + gi;
                    double qr = qgm_re[idx], qi = qgm_im[idx];
                    double br = becpsi_re[ijkb0 + jh], bi = becpsi_im[ijkb0 + jh];
                    a1r += SOA_MUL_RE(qr, qi, br, bi);
                    a1i += SOA_MUL_IM(qr, qi, br, bi);
                }
                double pr = becphi_re[ijkb0 + ih], pi = becphi_im[ijkb0 + ih];
                a2r += SOA_MULCONJ_RE(a1r, a1i, pr, pi);
                a2i += SOA_MULCONJ_IM(a1r, a1i, pr, pi);
            }

            // Transposed eigts
            double sr = eigqts_re[na], si = eigqts_im[na];
            double tr, ti;
            int idx;
            idx = EIGTS1T_IDX(na, m1);
            tr = sr; ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts1_T_re[idx], eigts1_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts1_T_re[idx], eigts1_T_im[idx]);
            idx = EIGTS2T_IDX(na, m2);
            tr = sr; ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts2_T_re[idx], eigts2_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts2_T_re[idx], eigts2_T_im[idx]);
            idx = EIGTS3T_IDX(na, m3);
            tr = sr; ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts3_T_re[idx], eigts3_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts3_T_re[idx], eigts3_T_im[idx]);

            acc_r += SOA_MUL_RE(a2r, a2i, sr, si);
            acc_i += SOA_MUL_IM(a2r, a2i, sr, si);
        }

        int nl_idx = dfftt__nl[gi] - 1;
        rhoc_re[nl_idx] += acc_r;
        rhoc_im[nl_idx] += acc_i;
    }
}

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
    double* d_eigqts_re, double* d_eigqts_im)
{
    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_eigts_transposed_soa<<<blocks, tblock_size>>>(
            d_rhoc_re, d_rhoc_im,
            d_becphi_re, d_becphi_im,
            d_becpsi_re, d_becpsi_im,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_re, d_qgm_im,
            d_eigts1_T_re, d_eigts1_T_im,
            d_eigts2_T_re, d_eigts2_T_im,
            d_eigts3_T_re, d_eigts3_T_im,
            d_mill, d_dfftt__nl,
            d_eigqts_re, d_eigqts_im,
            coarsen);
    }
}

// ============================================================
// GPU shared-bec SoA — eigts-transposed + shared memory bec
// ============================================================
__global__ void kernel_addusxx_shared_bec_soa(
    double* __restrict__ rhoc_re, double* __restrict__ rhoc_im,
    const double* __restrict__ becphi_re, const double* __restrict__ becphi_im,
    const double* __restrict__ becpsi_re, const double* __restrict__ becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const double* __restrict__ qgm_re, const double* __restrict__ qgm_im,
    const double* __restrict__ eigts1_T_re, const double* __restrict__ eigts1_T_im,
    const double* __restrict__ eigts2_T_re, const double* __restrict__ eigts2_T_im,
    const double* __restrict__ eigts3_T_re, const double* __restrict__ eigts3_T_im,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const double* __restrict__ eigqts_re, const double* __restrict__ eigqts_im,
    int coarsen)
{
    __shared__ double s_bpsi_re[20], s_bpsi_im[20];
    __shared__ double s_bphi_re[20], s_bphi_im[20];

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    double acc_r[8], acc_i[8];
    for (int c = 0; c < coarsen; c++) { acc_r[c] = 0; acc_i[c] = 0; }

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;
        int ijkb0 = ofsbeta[na];

        if (threadIdx.x < nh_nt) {
            s_bpsi_re[threadIdx.x] = becpsi_re[ijkb0 + threadIdx.x];
            s_bpsi_im[threadIdx.x] = becpsi_im[ijkb0 + threadIdx.x];
            s_bphi_re[threadIdx.x] = becphi_re[ijkb0 + threadIdx.x];
            s_bphi_im[threadIdx.x] = becphi_im[ijkb0 + threadIdx.x];
        }
        __syncthreads();

        for (int c = 0; c < coarsen; c++) {
            int gi = g_start + c;
            if (gi >= ngms) break;

            double a2r = 0, a2i = 0;
            for (int ih = 0; ih < nh_nt; ih++) {
                double a1r = 0, a1i = 0;
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based-1), IJTOH_N1, IJTOH_N2)];
                    int idx = (nij + ijtoh_val - 1) * QGM_NROWS + gi;
                    double qr = qgm_re[idx], qi = qgm_im[idx];
                    a1r += SOA_MUL_RE(qr, qi, s_bpsi_re[jh], s_bpsi_im[jh]);
                    a1i += SOA_MUL_IM(qr, qi, s_bpsi_re[jh], s_bpsi_im[jh]);
                }
                a2r += SOA_MULCONJ_RE(a1r, a1i, s_bphi_re[ih], s_bphi_im[ih]);
                a2i += SOA_MULCONJ_IM(a1r, a1i, s_bphi_re[ih], s_bphi_im[ih]);
            }

            int m1 = mill[gi*3+0], m2 = mill[gi*3+1], m3 = mill[gi*3+2];
            int e1 = EIGTS1T_IDX(na,m1), e2 = EIGTS2T_IDX(na,m2), e3 = EIGTS3T_IDX(na,m3);
            double sr = eigqts_re[na], si = eigqts_im[na], tr, ti;
            tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]); si=SOA_MUL_IM(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]);
            tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]); si=SOA_MUL_IM(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]);
            tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]); si=SOA_MUL_IM(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]);

            acc_r[c] += SOA_MUL_RE(a2r, a2i, sr, si);
            acc_i[c] += SOA_MUL_IM(a2r, a2i, sr, si);
        }

        __syncthreads();
    }

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) break;
        int nl_idx = dfftt__nl[gi] - 1;
        rhoc_re[nl_idx] += acc_r[c];
        rhoc_im[nl_idx] += acc_i[c];
    }
}

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
    double* d_eigqts_re, double* d_eigqts_im)
{
    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_shared_bec_soa<<<blocks, tblock_size>>>(
            d_rhoc_re, d_rhoc_im,
            d_becphi_re, d_becphi_im,
            d_becpsi_re, d_becpsi_im,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_re, d_qgm_im,
            d_eigts1_T_re, d_eigts1_T_im,
            d_eigts2_T_re, d_eigts2_T_im,
            d_eigts3_T_re, d_eigts3_T_im,
            d_mill, d_dfftt__nl,
            d_eigqts_re, d_eigqts_im,
            coarsen);
    }
}

// ============================================================
// GPU optimized SoA — FUSED compute+scatter
// Sorted iteration: sequential writes, random reads via orig_g
// ============================================================
__global__ void kernel_addusxx_optimized_fused_soa(
    double* __restrict__ rhoc_re, double* __restrict__ rhoc_im,
    const double* __restrict__ becphi_re, const double* __restrict__ becphi_im,
    const double* __restrict__ becpsi_re, const double* __restrict__ becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const double* __restrict__ qgm_T_re, const double* __restrict__ qgm_T_im,
    const double* __restrict__ eigts1_T_re, const double* __restrict__ eigts1_T_im,
    const double* __restrict__ eigts2_T_re, const double* __restrict__ eigts2_T_im,
    const double* __restrict__ eigts3_T_re, const double* __restrict__ eigts3_T_im,
    const int* __restrict__ mill,
    const int* __restrict__ dfftt__nl_sorted,
    const int* __restrict__ dfftt__nl_ix,
    const double* __restrict__ eigqts_re, const double* __restrict__ eigqts_im,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        int nl_idx = dfftt__nl_sorted[gi] - 1;   // sequential write target
        int orig_g = dfftt__nl_ix[gi] - 1;        // original g-vector for reads

        int m1 = mill[orig_g * 3 + 0];
        int m2 = mill[orig_g * 3 + 1];
        int m3 = mill[orig_g * 3 + 2];

        double acc_r = 0.0, acc_i = 0.0;

        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            int ijkb0 = ofsbeta[na];
            double a2r = 0.0, a2i = 0.0;

            for (int ih = 0; ih < nh_nt; ih++) {
                double a1r = 0.0, a1i = 0.0;
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), IJTOH_N1, IJTOH_N2)];
                    int idx = orig_g * QGMT_NROWS + (nij + ijtoh_val - 1);
                    double qr = qgm_T_re[idx], qi = qgm_T_im[idx];
                    double br = becpsi_re[ijkb0 + jh], bi = becpsi_im[ijkb0 + jh];
                    a1r += SOA_MUL_RE(qr, qi, br, bi);
                    a1i += SOA_MUL_IM(qr, qi, br, bi);
                }
                double pr = becphi_re[ijkb0 + ih], pi = becphi_im[ijkb0 + ih];
                a2r += SOA_MULCONJ_RE(a1r, a1i, pr, pi);
                a2i += SOA_MULCONJ_IM(a1r, a1i, pr, pi);
            }

            double sr, si;
            compute_sf_soa(sr, si, na,
                           eigqts_re, eigqts_im,
                           eigts1_T_re, eigts1_T_im,
                           eigts2_T_re, eigts2_T_im,
                           eigts3_T_re, eigts3_T_im,
                           EIGTS1T_IDX(na, m1), EIGTS2T_IDX(na, m2), EIGTS3T_IDX(na, m3));

            acc_r += SOA_MUL_RE(a2r, a2i, sr, si);
            acc_i += SOA_MUL_IM(a2r, a2i, sr, si);
        }

        // Single sequential write
        rhoc_re[nl_idx] += acc_r;
        rhoc_im[nl_idx] += acc_i;
    }
}

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
    double* d_eigqts_re, double* d_eigqts_im)
{
    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_optimized_fused_soa<<<blocks, tblock_size>>>(
            d_rhoc_re, d_rhoc_im,
            d_becphi_re, d_becphi_im,
            d_becpsi_re, d_becpsi_im,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_T_re, d_qgm_T_im,
            d_eigts1_T_re, d_eigts1_T_im,
            d_eigts2_T_re, d_eigts2_T_im,
            d_eigts3_T_re, d_eigts3_T_im,
            d_mill, d_dfftt__nl_sorted, d_dfftt__nl_ix,
            d_eigqts_re, d_eigqts_im,
            coarsen);
    }
}