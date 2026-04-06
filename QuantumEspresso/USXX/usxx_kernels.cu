#include "usxx_kernels.h"
#include <cstdio>
#include <cmath>
#include <cstring>

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
    const int* mill, const int* dfftt__nl)
{
    Complex_DP* eigqts = new Complex_DP[nat];
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0] - xkq[0]) * tau[na * 3 + 0]
                      + (xk[1] - xkq[1]) * tau[na * 3 + 1]
                      + (xk[2] - xkq[2]) * tau[na * 3 + 2]);
        eigqts[na] = make_cmplx(cos(arg), -sin(arg));
    }

    constexpr int blocksize = 256;
    int numblock = (ngms + blocksize - 1) / blocksize;
    constexpr int qgm_nrows = 55191;
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            Complex_DP* aux1 = new Complex_DP[blocksize];
            Complex_DP* aux2 = new Complex_DP[blocksize];

            #pragma omp for
            for (int iblock = 0; iblock < numblock; iblock++) {
                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int offset = iblock * blocksize;
                    int rbs = (ngms - offset < blocksize) ? (ngms - offset) : blocksize;
                    int ijkb0 = ofsbeta[na];

                    for (int s = 0; s < rbs; s++) aux2[s] = make_cmplx(0.0, 0.0);
                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) aux1[s] = make_cmplx(0.0, 0.0);
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, ijtoh_n1, ijtoh_n2)];
                            int qgm_col = nij + ijtoh_val - 1;
                            Complex_DP bpsi = becpsi_c[ijkb0 + jh];
                            for (int s = 0; s < rbs; s++)
                                aux1[s] = cadd(aux1[s], cmul(qgm[qgm_col * qgm_nrows + offset + s], bpsi));
                        }
                        Complex_DP cbphi = cconj(becphi_c[ijkb0 + ih]);
                        for (int s = 0; s < rbs; s++)
                            aux2[s] = cadd(aux2[s], cmul(aux1[s], cbphi));
                    }
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int m1 = mill[g * 3 + 0], m2 = mill[g * 3 + 1], m3 = mill[g * 3 + 2];
                        Complex_DP sf = cmul(eigqts[na],
                                       cmul(eigts1[EIGTS1_IDX(m1, na)],
                                       cmul(eigts2[EIGTS2_IDX(m2, na)],
                                            eigts3[EIGTS3_IDX(m3, na)])));
                        aux2[s] = cmul(aux2[s], sf);
                    }
                    for (int s = 0; s < rbs; s++) {
                        int nl_idx = dfftt__nl[offset + s] - 1;
                        rhoc[nl_idx] = cadd(rhoc[nl_idx], aux2[s]);
                    }
                }
            }
            delete[] aux1;
            delete[] aux2;
        }
    }
    delete[] eigqts;
}

// ============================================================
// Shared eigqts kernel
// ============================================================
__global__ void kernel_eigqts(Complex_DP* eigqts, const DP* xkq, const DP* xk,
                              const DP* tau, int nat) {
    int na = blockIdx.x * blockDim.x + threadIdx.x;
    if (na >= nat) return;
    DP arg = TPI * ((xk[0] - xkq[0]) * tau[na * 3 + 0]
                  + (xk[1] - xkq[1]) * tau[na * 3 + 1]
                  + (xk[2] - xkq[2]) * tau[na * 3 + 2]);
    eigqts[na] = make_cmplx(cos(arg), -sin(arg));
}

// ============================================================
// GPU baseline AoS — coarsened, no atomics
// ============================================================
__global__ void kernel_addusxx_baseline(
    Complex_DP* __restrict__ rhoc,
    const Complex_DP* __restrict__ becphi_c,
    const Complex_DP* __restrict__ becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const Complex_DP* __restrict__ qgm,
    const Complex_DP* __restrict__ eigts1,
    const Complex_DP* __restrict__ eigts2,
    const Complex_DP* __restrict__ eigts3,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const Complex_DP* __restrict__ eigqts,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    constexpr int qgm_nrows = 55191;
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        Complex_DP accum = make_cmplx(0.0, 0.0);

        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            int ijkb0 = ofsbeta[na];

            Complex_DP aux2_val = make_cmplx(0.0, 0.0);
            for (int ih = 0; ih < nh_nt; ih++) {
                Complex_DP aux1_val = make_cmplx(0.0, 0.0);
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), ijtoh_n1, ijtoh_n2)];
                    aux1_val = cadd(aux1_val,
                        cmul(qgm[(nij + ijtoh_val - 1) * qgm_nrows + gi],
                             becpsi_c[ijkb0 + jh]));
                }
                aux2_val = cadd(aux2_val,
                    cmul(aux1_val, cconj(becphi_c[ijkb0 + ih])));
            }

            int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];
            Complex_DP sf = cmul(eigqts[na],
                           cmul(eigts1[EIGTS1_IDX(m1, na)],
                           cmul(eigts2[EIGTS2_IDX(m2, na)],
                                eigts3[EIGTS3_IDX(m3, na)])));
            accum = cadd(accum, cmul(aux2_val, sf));
        }

        int nl_idx = dfftt__nl[gi] - 1;
        rhoc[nl_idx] = cadd(rhoc[nl_idx], accum);
    }
}

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
    int tblock_size, int coarsen)
{
    Complex_DP* d_eigqts;
    cudaMalloc(&d_eigqts, nat * sizeof(Complex_DP));
    kernel_eigqts<<<(nat + 255) / 256, 256>>>(d_eigqts, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        kernel_addusxx_baseline<<<blocks, tblock_size>>>(
            d_rhoc, d_becphi_c, d_becpsi_c,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm, d_eigts1, d_eigts2, d_eigts3,
            d_mill, d_dfftt__nl, d_eigqts,
            coarsen);
    }
    cudaFree(d_eigqts);
}

// ============================================================
// GPU optimized AoS — compute with coarsening
// ============================================================
__global__ void kernel_addusxx_optimized_compute(
    Complex_DP* __restrict__ aux2_array,
    const Complex_DP* __restrict__ becphi_c,
    const Complex_DP* __restrict__ becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const Complex_DP* __restrict__ qgm_T,
    const Complex_DP* __restrict__ eigts1_T,
    const Complex_DP* __restrict__ eigts2_T,
    const Complex_DP* __restrict__ eigts3_T,
    const int* __restrict__ mill,
    const Complex_DP* __restrict__ eigqts,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;
    constexpr int qgmT_nrows = 397;
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            int ijkb0 = ofsbeta[na];
            Complex_DP aux2_val = make_cmplx(0.0, 0.0);

            for (int ih = 0; ih < nh_nt; ih++) {
                Complex_DP aux1_val = make_cmplx(0.0, 0.0);
                #pragma unroll 4
                for (int jh = 0; jh < nh_nt; jh++) {
                    int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), ijtoh_n1, ijtoh_n2)];
                    aux1_val = cadd(aux1_val,
                        cmul(qgm_T[gi * qgmT_nrows + (nij + ijtoh_val - 1)],
                             becpsi_c[ijkb0 + jh]));
                }
                aux2_val = cadd(aux2_val,
                    cmul(aux1_val, cconj(becphi_c[ijkb0 + ih])));
            }

            int m1 = mill[gi * 3 + 0], m2 = mill[gi * 3 + 1], m3 = mill[gi * 3 + 2];
            Complex_DP sf = cmul(eigqts[na],
                           cmul(eigts1_T[EIGTS1T_IDX(na, m1)],
                           cmul(eigts2_T[EIGTS2T_IDX(na, m2)],
                                eigts3_T[EIGTS3T_IDX(na, m3)])));
            aux2_array[gi * nat + na] = cmul(aux2_val, sf);
        }
    }
}

// Scatter — no atomics
__global__ void kernel_addusxx_optimized_scatter(
    Complex_DP* __restrict__ rhoc,
    const Complex_DP* __restrict__ aux2_array,
    int ngms, int nat, int nt_1based,
    const int* __restrict__ ityp,
    const int* __restrict__ dfftt__nl_sorted,
    const int* __restrict__ dfftt__nl_ix,
    int coarsen)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int g_start = tid * coarsen;

    for (int c = 0; c < coarsen; c++) {
        int gi = g_start + c;
        if (gi >= ngms) return;

        int nl_idx = dfftt__nl_sorted[gi] - 1;
        int ix = dfftt__nl_ix[gi] - 1;

        Complex_DP accum = make_cmplx(0.0, 0.0);
        for (int na = 0; na < nat; na++) {
            if (ityp[na] != nt_1based) continue;
            accum = cadd(accum, aux2_array[ix * nat + na]);
        }
        rhoc[nl_idx] = cadd(rhoc[nl_idx], accum);
    }
}

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
    int tblock_size, int coarsen)
{
    Complex_DP* d_eigqts;
    cudaMalloc(&d_eigqts, nat * sizeof(Complex_DP));
    Complex_DP* d_aux2_array;
    cudaMalloc(&d_aux2_array, (size_t)nat * ngms * sizeof(Complex_DP));

    kernel_eigqts<<<(nat + 255) / 256, 256>>>(d_eigqts, d_xkq, d_xk, d_tau, nat);

    int n_logical = (ngms + coarsen - 1) / coarsen;
    int blocks = (n_logical + tblock_size - 1) / tblock_size;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        cudaMemset(d_aux2_array, 0, (size_t)nat * ngms * sizeof(Complex_DP));

        kernel_addusxx_optimized_compute<<<blocks, tblock_size>>>(
            d_aux2_array, d_becphi_c, d_becpsi_c,
            nkb, ngms, nat,
            h_nij_type[nt], h_nh[nt], nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_T, d_eigts1_T, d_eigts2_T, d_eigts3_T,
            d_mill, d_eigqts, coarsen);

        kernel_addusxx_optimized_scatter<<<blocks, tblock_size>>>(
            d_rhoc, d_aux2_array,
            ngms, nat, nt + 1, d_ityp,
            d_dfftt__nl_sorted, d_dfftt__nl_ix, coarsen);
    }

    cudaFree(d_aux2_array);
    cudaFree(d_eigqts);
}
