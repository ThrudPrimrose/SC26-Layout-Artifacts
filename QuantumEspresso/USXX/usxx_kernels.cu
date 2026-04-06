#include "usxx_kernels.h"
#include <cstdio>
#include <cmath>
#include <cstring>

// ============================================================
// CPU baseline (mirrors Fortran addusxx_g)
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
    // All indices are 0-based in C; Fortran arrays loaded column-major.
    // Fortran 1-based offsets into qgm second dim etc. are adjusted.

    Complex_DP* eigqts = new Complex_DP[nat];
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0] - xkq[0]) * tau[na * 3 + 0]
                      + (xk[1] - xkq[1]) * tau[na * 3 + 1]
                      + (xk[2] - xkq[2]) * tau[na * 3 + 2]);
        eigqts[na] = make_cmplx(cos(arg), -sin(arg));
    }

    constexpr int blocksize = 256;
    int numblock = (ngms + blocksize - 1) / blocksize;

    // qgm is Fortran (ngms_dim, nqgm_dim) col-major => qgm[col * ngms_dim + row]
    // ngms_dim = 55191 (first dim of qgm allocation)
    constexpr int qgm_nrows = 55191;
    // ijtoh is Fortran (19,19,3) col-major
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt]; // Fortran 1-based offset
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            Complex_DP* aux1 = new Complex_DP[blocksize];
            Complex_DP* aux2 = new Complex_DP[blocksize];

            #pragma omp for
            for (int iblock = 0; iblock < numblock; iblock++) {
                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue; // ityp is 1-based type index

                    int offset = iblock * blocksize;
                    int rbs = (ngms - offset < blocksize) ? (ngms - offset) : blocksize;

                    int ijkb0 = ofsbeta[na]; // Fortran 1-based

                    for (int s = 0; s < rbs; s++) aux2[s] = make_cmplx(0.0, 0.0);

                    for (int ih = 0; ih < nh_nt; ih++) {
                        int ikb = ijkb0 + ih; // 0-based into becphi_c: Fortran ikb = ijkb0+ih, C index = ikb-1
                        for (int s = 0; s < rbs; s++) aux1[s] = make_cmplx(0.0, 0.0);

                        for (int jh = 0; jh < nh_nt; jh++) {
                            int jkb = ijkb0 + jh;
                            // ijtoh(ih+1, jh+1, nt+1) in Fortran => IDX3(ih, jh, nt, 19, 19)
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, ijtoh_n1, ijtoh_n2)]; // 1-based
                            // qgm column index (0-based): nij + ijtoh_val - 1
                            int qgm_col = nij + ijtoh_val - 1; // nij is 1-based offset in Fortran, ijtoh_val is 1-based
                            // Actually in Fortran: qgm(g, nij+ijtoh(ih,jh,nt))
                            // nij and ijtoh are both Fortran integers used as 1-based column indices
                            // In C 0-based: column = (nij + ijtoh_val) - 1
                            Complex_DP bpsi = becpsi_c[jkb - 1]; // jkb is 1-based
                            for (int s = 0; s < rbs; s++) {
                                int g = offset + s; // 0-based g-vector index
                                // qgm(g+1, nij+ijtoh_val) => qgm[(nij+ijtoh_val-1) * qgm_nrows + g]
                                aux1[s] = cadd(aux1[s], cmul(qgm[qgm_col * qgm_nrows + g], bpsi));
                            }
                        }

                        Complex_DP cbphi = cconj(becphi_c[ikb - 1]); // ikb is 1-based
                        for (int s = 0; s < rbs; s++) {
                            aux2[s] = cadd(aux2[s], cmul(aux1[s], cbphi));
                        }
                    }

                    // Multiply by eigqts and structure factors
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        // mill is Fortran (3, 156121) col-major: mill(dim+1, g+1) => mill[g*3 + dim]
                        int m1 = mill[g * 3 + 0]; // these are actual integer values (signed)
                        int m2 = mill[g * 3 + 1];
                        int m3 = mill[g * 3 + 2];
                        Complex_DP sf = cmul(eigqts[na],
                                       cmul(eigts1[EIGTS1_IDX(m1, na)],
                                       cmul(eigts2[EIGTS2_IDX(m2, na)],
                                            eigts3[EIGTS3_IDX(m3, na)])));
                        aux2[s] = cmul(aux2[s], sf);
                    }

                    // Scatter into rhoc
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int nl_idx = dfftt__nl[g] - 1; // Fortran 1-based index
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
// GPU baseline CUDA kernel
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

__global__ void kernel_addusxx_baseline(
    Complex_DP* rhoc,
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* ityp, const int* ofsbeta, const int* ijtoh,
    const Complex_DP* qgm, const Complex_DP* eigts1,
    const Complex_DP* eigts2, const Complex_DP* eigts3,
    const int* mill, const int* dfftt__nl,
    const Complex_DP* eigqts)
{
    int ngms_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (ngms_i >= ngms) return;

    constexpr int qgm_nrows = 55191;
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;

        int ijkb0 = ofsbeta[na]; // 1-based

        Complex_DP aux2_val = make_cmplx(0.0, 0.0);

        for (int ih = 0; ih < nh_nt; ih++) {
            Complex_DP aux1_val = make_cmplx(0.0, 0.0);

            #pragma unroll 4
            for (int jh = 0; jh < nh_nt; jh++) {
                int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), ijtoh_n1, ijtoh_n2)];
                int qgm_col = nij + ijtoh_val - 1;
                aux1_val = cadd(aux1_val,
                    cmul(qgm[qgm_col * qgm_nrows + ngms_i],
                         becpsi_c[ijkb0 + jh - 1]));
            }

            aux2_val = cadd(aux2_val,
                cmul(aux1_val, cconj(becphi_c[ijkb0 + ih - 1])));
        }

        int m1 = mill[ngms_i * 3 + 0];
        int m2 = mill[ngms_i * 3 + 1];
        int m3 = mill[ngms_i * 3 + 2];
        Complex_DP sf = cmul(eigqts[na],
                       cmul(eigts1[EIGTS1_IDX(m1, na)],
                       cmul(eigts2[EIGTS2_IDX(m2, na)],
                            eigts3[EIGTS3_IDX(m3, na)])));
        aux2_val = cmul(aux2_val, sf);

        int nl_idx = dfftt__nl[ngms_i] - 1;
        // Atomic add (real and imag parts)
        atomicAdd(&(reinterpret_cast<double*>(rhoc))[nl_idx * 2], cuCreal(aux2_val));
        atomicAdd(&(reinterpret_cast<double*>(rhoc))[nl_idx * 2 + 1], cuCimag(aux2_val));
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
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh)
{
    Complex_DP* d_eigqts;
    cudaMalloc(&d_eigqts, nat * sizeof(Complex_DP));

    kernel_eigqts<<<(nat + 255) / 256, 256>>>(d_eigqts, d_xkq, d_xk, d_tau, nat);

    int threads = 256;
    int blocks = (ngms + threads - 1) / threads;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        int nij = h_nij_type[nt];
        int nh_nt = h_nh[nt];

        kernel_addusxx_baseline<<<blocks, threads>>>(
            d_rhoc, d_becphi_c, d_becpsi_c,
            nkb, ngms, nat,
            nij, nh_nt, nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm, d_eigts1, d_eigts2, d_eigts3,
            d_mill, d_dfftt__nl, d_eigqts);
    }

    cudaFree(d_eigqts);
}

// ============================================================
// GPU optimized CUDA kernel (transposed layouts)
// ============================================================
__global__ void kernel_addusxx_optimized_compute(
    Complex_DP* aux2_array, // (nat, ngms) col-major
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* ityp, const int* ofsbeta, const int* ijtoh,
    const Complex_DP* qgm_T, const Complex_DP* eigts1_T,
    const Complex_DP* eigts2_T, const Complex_DP* eigts3_T,
    const int* mill,
    const Complex_DP* eigqts)
{
    int ngms_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (ngms_i >= ngms) return;

    // qgm_T is Fortran (397, 55191) col-major: qgm_T(col, row) => qgm_T[row * 397 + col]
    constexpr int qgmT_nrows = 397;
    constexpr int ijtoh_n1 = 19, ijtoh_n2 = 19;

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;

        int ijkb0 = ofsbeta[na];
        Complex_DP aux2_val = make_cmplx(0.0, 0.0);

        for (int ih = 0; ih < nh_nt; ih++) {
            Complex_DP aux1_val = make_cmplx(0.0, 0.0);

            #pragma unroll 4
            for (int jh = 0; jh < nh_nt; jh++) {
                int ijtoh_val = ijtoh[IDX3(ih, jh, (nt_1based - 1), ijtoh_n1, ijtoh_n2)];
                // qgm_T(nij+ijtoh_val, ngms_i+1) => qgm_T[ngms_i * qgmT_nrows + (nij+ijtoh_val-1)]
                int qgmT_idx = ngms_i * qgmT_nrows + (nij + ijtoh_val - 1);
                aux1_val = cadd(aux1_val,
                    cmul(qgm_T[qgmT_idx], becpsi_c[ijkb0 + jh - 1]));
            }

            aux2_val = cadd(aux2_val,
                cmul(aux1_val, cconj(becphi_c[ijkb0 + ih - 1])));
        }

        int m1 = mill[ngms_i * 3 + 0];
        int m2 = mill[ngms_i * 3 + 1];
        int m3 = mill[ngms_i * 3 + 2];
        Complex_DP sf = cmul(eigqts[na],
                       cmul(eigts1_T[EIGTS1T_IDX(na, m1)],
                       cmul(eigts2_T[EIGTS2T_IDX(na, m2)],
                            eigts3_T[EIGTS3T_IDX(na, m3)])));
        aux2_val = cmul(aux2_val, sf);

        // aux2_array is (nat, ngms) col-major: aux2_array[ngms_i * nat + na]
        aux2_array[ngms_i * nat + na] = aux2_val;
    }
}

__global__ void kernel_addusxx_optimized_scatter(
    Complex_DP* rhoc,
    const Complex_DP* aux2_array,
    int ngms, int nat, int nt_1based,
    const int* ityp,
    const int* dfftt__nl_sorted, const int* dfftt__nl_ix)
{
    int ngms_i = blockIdx.x * blockDim.x + threadIdx.x;
    if (ngms_i >= ngms) return;

    int nl_idx = dfftt__nl_sorted[ngms_i] - 1; // 1-based to 0-based
    int ix = dfftt__nl_ix[ngms_i] - 1;          // 1-based to 0-based

    for (int na = 0; na < nat; na++) {
        if (ityp[na] != nt_1based) continue;
        Complex_DP val = aux2_array[ix * nat + na];
        atomicAdd(&(reinterpret_cast<double*>(rhoc))[nl_idx * 2], cuCreal(val));
        atomicAdd(&(reinterpret_cast<double*>(rhoc))[nl_idx * 2 + 1], cuCimag(val));
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
    const int* h_upf_tvanp, const int* h_nij_type, const int* h_nh)
{
    Complex_DP* d_eigqts;
    cudaMalloc(&d_eigqts, nat * sizeof(Complex_DP));

    Complex_DP* d_aux2_array;
    cudaMalloc(&d_aux2_array, (size_t)nat * ngms * sizeof(Complex_DP));

    kernel_eigqts<<<(nat + 255) / 256, 256>>>(d_eigqts, d_xkq, d_xk, d_tau, nat);

    int threads = 256;
    int blocks = (ngms + threads - 1) / threads;

    for (int nt = 0; nt < ntyp; nt++) {
        if (h_upf_tvanp[nt] != 1) continue;
        int nij = h_nij_type[nt];
        int nh_nt = h_nh[nt];

        // Zero aux2_array
        cudaMemset(d_aux2_array, 0, (size_t)nat * ngms * sizeof(Complex_DP));

        kernel_addusxx_optimized_compute<<<blocks, threads>>>(
            d_aux2_array,
            d_becphi_c, d_becpsi_c,
            nkb, ngms, nat,
            nij, nh_nt, nt + 1,
            d_ityp, d_ofsbeta, d_ijtoh,
            d_qgm_T, d_eigts1_T, d_eigts2_T, d_eigts3_T,
            d_mill, d_eigqts);

        kernel_addusxx_optimized_scatter<<<blocks, threads>>>(
            d_rhoc, d_aux2_array,
            ngms, nat, nt + 1,
            d_ityp,
            d_dfftt__nl_sorted, d_dfftt__nl_ix);
    }

    cudaFree(d_aux2_array);
    cudaFree(d_eigqts);
}
