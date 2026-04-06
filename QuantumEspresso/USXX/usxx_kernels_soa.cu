#include "usxx_kernels.h"
#include <cstdio>
#include <cmath>

// ============================================================
// SoA complex arithmetic (inline device helpers)
// ============================================================
// mul: (a_re + i*a_im)*(b_re + i*b_im)
//    = (a_re*b_re - a_im*b_im) + i*(a_re*b_im + a_im*b_re)
#define SOA_MUL_RE(ar, ai, br, bi) ((ar)*(br) - (ai)*(bi))
#define SOA_MUL_IM(ar, ai, br, bi) ((ar)*(bi) + (ai)*(br))

// conj(a) = (a_re, -a_im)
// mul by conj: (a_re + i*a_im)*(b_re - i*b_im)
//            = (a_re*b_re + a_im*b_im) + i*(a_im*b_re - a_re*b_im)
#define SOA_MULCONJ_RE(ar, ai, br, bi) ((ar)*(br) + (ai)*(bi))
#define SOA_MULCONJ_IM(ar, ai, br, bi) ((ai)*(br) - (ar)*(bi))

static constexpr int QGM_NROWS  = 55191;
static constexpr int QGMT_NROWS = 397;
static constexpr int IJTOH_N1   = 19;
static constexpr int IJTOH_N2   = 19;

// ============================================================
// eigqts kernel (SoA output)
// ============================================================
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

// ============================================================
// GPU baseline — SoA
// ============================================================
__global__ void kernel_addusxx_baseline_soa(
    double* rhoc_re, double* rhoc_im,
    const double* becphi_re, const double* becphi_im,
    const double* becpsi_re, const double* becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* ityp, const int* ofsbeta, const int* ijtoh,
    const double* qgm_re, const double* qgm_im,
    const double* eigts1_re, const double* eigts1_im,
    const double* eigts2_re, const double* eigts2_im,
    const double* eigts3_re, const double* eigts3_im,
    const int* mill, const int* dfftt__nl,
    const double* eigqts_re, const double* eigqts_im)
{
    int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= ngms) return;

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
                double br = becpsi_re[ijkb0 + jh - 1], bi = becpsi_im[ijkb0 + jh - 1];
                a1r += SOA_MUL_RE(qr, qi, br, bi);
                a1i += SOA_MUL_IM(qr, qi, br, bi);
            }

            // aux2 += aux1 * conj(becphi(ikb))
            double pr = becphi_re[ijkb0 + ih - 1], pi = becphi_im[ijkb0 + ih - 1];
            a2r += SOA_MULCONJ_RE(a1r, a1i, pr, pi);
            a2i += SOA_MULCONJ_IM(a1r, a1i, pr, pi);
        }

        // Multiply by eigqts * eigts1 * eigts2 * eigts3
        int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];

        double sr = eigqts_re[na], si = eigqts_im[na];

        {   // *= eigts1
            int idx = EIGTS1_IDX(m1, na);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts1_re[idx], eigts1_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts1_re[idx], eigts1_im[idx]);
        }
        {   // *= eigts2
            int idx = EIGTS2_IDX(m2, na);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts2_re[idx], eigts2_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts2_re[idx], eigts2_im[idx]);
        }
        {   // *= eigts3
            int idx = EIGTS3_IDX(m3, na);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts3_re[idx], eigts3_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts3_re[idx], eigts3_im[idx]);
        }

        // aux2 *= sf
        double rr = SOA_MUL_RE(a2r, a2i, sr, si);
        double ri = SOA_MUL_IM(a2r, a2i, sr, si);

        int nl_idx = dfftt__nl[gi] - 1;
        atomicAdd(&rhoc_re[nl_idx], rr);
        atomicAdd(&rhoc_im[nl_idx], ri);
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
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh)
{
    double *d_eigqts_re, *d_eigqts_im;
    cudaMalloc(&d_eigqts_re, nat * sizeof(double));
    cudaMalloc(&d_eigqts_im, nat * sizeof(double));

    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int threads = 256;
    int blocks = (ngms + threads - 1) / threads;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_baseline_soa<<<blocks, threads>>>(
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
            d_eigqts_re, d_eigqts_im);
    }

    cudaFree(d_eigqts_re);
    cudaFree(d_eigqts_im);
}

// ============================================================
// GPU optimized — SoA (transposed layouts + sorted scatter)
// ============================================================
__global__ void kernel_addusxx_opt_soa_compute(
    double* aux2_re, double* aux2_im, // (nat, ngms) col-major
    const double* becphi_re, const double* becphi_im,
    const double* becpsi_re, const double* becpsi_im,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* ityp, const int* ofsbeta, const int* ijtoh,
    const double* qgm_T_re, const double* qgm_T_im,
    const double* eigts1_T_re, const double* eigts1_T_im,
    const double* eigts2_T_re, const double* eigts2_T_im,
    const double* eigts3_T_re, const double* eigts3_T_im,
    const int* mill,
    const double* eigqts_re, const double* eigqts_im)
{
    int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= ngms) return;

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;

        int ijkb0 = ofsbeta[na];
        double a2r = 0.0, a2i = 0.0;

        for (int ih = 0; ih < nh_nt; ih++) {
            double a1r = 0.0, a1i = 0.0;

            #pragma unroll 4
            for (int jh = 0; jh < nh_nt; jh++) {
                int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), IJTOH_N1, IJTOH_N2)];
                int idx = gi * QGMT_NROWS + (nij + ijtoh_val - 1);
                double qr = qgm_T_re[idx], qi = qgm_T_im[idx];
                double br = becpsi_re[ijkb0 + jh - 1], bi = becpsi_im[ijkb0 + jh - 1];
                a1r += SOA_MUL_RE(qr, qi, br, bi);
                a1i += SOA_MUL_IM(qr, qi, br, bi);
            }

            double pr = becphi_re[ijkb0 + ih - 1], pi = becphi_im[ijkb0 + ih - 1];
            a2r += SOA_MULCONJ_RE(a1r, a1i, pr, pi);
            a2i += SOA_MULCONJ_IM(a1r, a1i, pr, pi);
        }

        int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];

        double sr = eigqts_re[na], si = eigqts_im[na];

        {
            int idx = EIGTS1T_IDX(na, m1);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts1_T_re[idx], eigts1_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts1_T_re[idx], eigts1_T_im[idx]);
        }
        {
            int idx = EIGTS2T_IDX(na, m2);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts2_T_re[idx], eigts2_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts2_T_re[idx], eigts2_T_im[idx]);
        }
        {
            int idx = EIGTS3T_IDX(na, m3);
            double tr = sr, ti = si;
            sr = SOA_MUL_RE(tr, ti, eigts3_T_re[idx], eigts3_T_im[idx]);
            si = SOA_MUL_IM(tr, ti, eigts3_T_re[idx], eigts3_T_im[idx]);
        }

        double rr = SOA_MUL_RE(a2r, a2i, sr, si);
        double ri = SOA_MUL_IM(a2r, a2i, sr, si);

        int out_idx = gi * nat + na;
        aux2_re[out_idx] = rr;
        aux2_im[out_idx] = ri;
    }
}

__global__ void kernel_addusxx_opt_soa_scatter(
    double* rhoc_re, double* rhoc_im,
    const double* aux2_re, const double* aux2_im,
    int ngms, int nat, int nt_1based,
    const int* ityp,
    const int* dfftt__nl_sorted, const int* dfftt__nl_ix)
{
    int gi = blockIdx.x * blockDim.x + threadIdx.x;
    if (gi >= ngms) return;

    int nl_idx = dfftt__nl_sorted[gi] - 1;
    int ix = dfftt__nl_ix[gi] - 1;

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;
        int src = ix * nat + na;
        atomicAdd(&rhoc_re[nl_idx], aux2_re[src]);
        atomicAdd(&rhoc_im[nl_idx], aux2_im[src]);
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
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh)
{
    double *d_eigqts_re, *d_eigqts_im;
    cudaMalloc(&d_eigqts_re, nat * sizeof(double));
    cudaMalloc(&d_eigqts_im, nat * sizeof(double));

    double *d_aux2_re, *d_aux2_im;
    size_t aux2_bytes = (size_t)nat * ngms * sizeof(double);
    cudaMalloc(&d_aux2_re, aux2_bytes);
    cudaMalloc(&d_aux2_im, aux2_bytes);

    kernel_eigqts_soa<<<(nat + 255) / 256, 256>>>(
        d_eigqts_re, d_eigqts_im, d_xkq, d_xk, d_tau, nat);

    int threads = 256;
    int blocks = (ngms + threads - 1) / threads;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;

        cudaMemset(d_aux2_re, 0, aux2_bytes);
        cudaMemset(d_aux2_im, 0, aux2_bytes);

        kernel_addusxx_opt_soa_compute<<<blocks, threads>>>(
            d_aux2_re, d_aux2_im,
            d_becphi_re, d_becphi_im,
            d_becpsi_re, d_becpsi_im,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_T_re, d_qgm_T_im,
            d_eigts1_T_re, d_eigts1_T_im,
            d_eigts2_T_re, d_eigts2_T_im,
            d_eigts3_T_re, d_eigts3_T_im,
            d_mill,
            d_eigqts_re, d_eigqts_im);

        kernel_addusxx_opt_soa_scatter<<<blocks, threads>>>(
            d_rhoc_re, d_rhoc_im,
            d_aux2_re, d_aux2_im,
            ngms, nat, nt + 1,
            d_ityp,
            d_dfftt__nl_sorted, d_dfftt__nl_ix);
    }

    cudaFree(d_aux2_re);
    cudaFree(d_aux2_im);
    cudaFree(d_eigqts_re);
    cudaFree(d_eigqts_im);
}
