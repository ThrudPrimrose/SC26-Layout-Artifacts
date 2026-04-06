#include "usxx_kernels_cpu.h"
#include <omp.h>

static constexpr int QGM_NROWS  = 55191;
static constexpr int IJTOH_N1   = 19;
static constexpr int IJTOH_N2   = 19;
static constexpr int BLOCKSIZE  = 256;

// ============================================================
// Shared: compute eigqts
// ============================================================
static void compute_eigqts(Complex_DP* eigqts, const DP* xkq, const DP* xk,
                           const DP* tau, int nat) {
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0]-xkq[0])*tau[na*3+0]
                      + (xk[1]-xkq[1])*tau[na*3+1]
                      + (xk[2]-xkq[2])*tau[na*3+2]);
        eigqts[na] = make_cmplx(cos(arg), -sin(arg));
    }
}

// ============================================================
// Baseline AoS: original layout, per-atom write, cache-blocked
// Matches Fortran CPU kernel structure exactly
// ============================================================
void cpu_addusxx_baseline_aos(
    Complex_DP* rhoc, const DP* xkq, const DP* xk, const DP* tau,
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const Complex_DP* qgm, const Complex_DP* eigts1,
    const Complex_DP* eigts2, const Complex_DP* eigts3,
    const int* mill, const int* dfftt__nl)
{
    Complex_DP eigqts[10];
    compute_eigqts(eigqts, xkq, xk, tau, nat);

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            Complex_DP aux1[BLOCKSIZE], aux2[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    for (int s = 0; s < rbs; s++) aux2[s] = make_cmplx(0.0, 0.0);
                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) aux1[s] = make_cmplx(0.0, 0.0);
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int qgm_col = nij + ijtoh_val - 1;
                            Complex_DP bpsi = becpsi_c[ijkb0 + jh];
                            for (int s = 0; s < rbs; s++)
                                aux1[s] = cadd(aux1[s], cmul(qgm[qgm_col*QGM_NROWS + offset+s], bpsi));
                        }
                        Complex_DP cbphi = cconj(becphi_c[ijkb0 + ih]);
                        for (int s = 0; s < rbs; s++)
                            aux2[s] = cadd(aux2[s], cmul(aux1[s], cbphi));
                    }
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int m1 = mill[g*3+0], m2 = mill[g*3+1], m3 = mill[g*3+2];
                        Complex_DP sf = cmul(eigqts[na],
                                       cmul(eigts1[EIGTS1_IDX(m1,na)],
                                       cmul(eigts2[EIGTS2_IDX(m2,na)],
                                            eigts3[EIGTS3_IDX(m3,na)])));
                        aux2[s] = cmul(aux2[s], sf);
                    }
                    for (int s = 0; s < rbs; s++) {
                        int nl_idx = dfftt__nl[offset+s] - 1;
                        rhoc[nl_idx] = cadd(rhoc[nl_idx], aux2[s]);
                    }
                }
            }
        }
    }
}

// ============================================================
// Eigts-transposed AoS: original qgm, transposed eigts,
// register accumulation across atoms, cache-blocked
// ============================================================
void cpu_addusxx_eigts_transposed_aos(
    Complex_DP* rhoc, const DP* xkq, const DP* xk, const DP* tau,
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const Complex_DP* qgm, const Complex_DP* eigts1_T,
    const Complex_DP* eigts2_T, const Complex_DP* eigts3_T,
    const int* mill, const int* dfftt__nl)
{
    Complex_DP eigqts[10];
    compute_eigqts(eigqts, xkq, xk, tau, nat);

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            Complex_DP aux1[BLOCKSIZE], accum[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                // Zero accumulators for this block
                for (int s = 0; s < rbs; s++) accum[s] = make_cmplx(0.0, 0.0);

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    Complex_DP aux2_local[BLOCKSIZE];
                    for (int s = 0; s < rbs; s++) aux2_local[s] = make_cmplx(0.0, 0.0);

                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) aux1[s] = make_cmplx(0.0, 0.0);
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int qgm_col = nij + ijtoh_val - 1;
                            Complex_DP bpsi = becpsi_c[ijkb0 + jh];
                            for (int s = 0; s < rbs; s++)
                                aux1[s] = cadd(aux1[s], cmul(qgm[qgm_col*QGM_NROWS + offset+s], bpsi));
                        }
                        Complex_DP cbphi = cconj(becphi_c[ijkb0 + ih]);
                        for (int s = 0; s < rbs; s++)
                            aux2_local[s] = cadd(aux2_local[s], cmul(aux1[s], cbphi));
                    }
                    // Transposed eigts + accumulate across atoms
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int m1 = mill[g*3+0], m2 = mill[g*3+1], m3 = mill[g*3+2];
                        Complex_DP sf = cmul(eigqts[na],
                                       cmul(eigts1_T[EIGTS1T_IDX(na,m1)],
                                       cmul(eigts2_T[EIGTS2T_IDX(na,m2)],
                                            eigts3_T[EIGTS3T_IDX(na,m3)])));
                        accum[s] = cadd(accum[s], cmul(aux2_local[s], sf));
                    }
                }
                // Single write per block element
                for (int s = 0; s < rbs; s++) {
                    int nl_idx = dfftt__nl[offset+s] - 1;
                    rhoc[nl_idx] = cadd(rhoc[nl_idx], accum[s]);
                }
            }
        }
    }
}

// ============================================================
// Sorted AoS: iterate in sorted dfftt__nl order,
// transposed eigts, register accumulation
// ============================================================
void cpu_addusxx_sorted_aos(
    Complex_DP* rhoc, const DP* xkq, const DP* xk, const DP* tau,
    const Complex_DP* becphi_c, const Complex_DP* becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const Complex_DP* qgm, const Complex_DP* eigts1_T,
    const Complex_DP* eigts2_T, const Complex_DP* eigts3_T,
    const int* mill, const int* dfftt__nl_sorted, const int* dfftt__nl_ix)
{
    Complex_DP eigqts[10];
    compute_eigqts(eigqts, xkq, xk, tau, nat);

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            Complex_DP aux1[BLOCKSIZE], accum[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                for (int s = 0; s < rbs; s++) accum[s] = make_cmplx(0.0, 0.0);

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    Complex_DP aux2_local[BLOCKSIZE];
                    for (int s = 0; s < rbs; s++) aux2_local[s] = make_cmplx(0.0, 0.0);

                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) aux1[s] = make_cmplx(0.0, 0.0);
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int qgm_col = nij + ijtoh_val - 1;
                            Complex_DP bpsi = becpsi_c[ijkb0 + jh];
                            for (int s = 0; s < rbs; s++) {
                                int orig_g = dfftt__nl_ix[offset+s] - 1;
                                aux1[s] = cadd(aux1[s], cmul(qgm[qgm_col*QGM_NROWS + orig_g], bpsi));
                            }
                        }
                        Complex_DP cbphi = cconj(becphi_c[ijkb0 + ih]);
                        for (int s = 0; s < rbs; s++)
                            aux2_local[s] = cadd(aux2_local[s], cmul(aux1[s], cbphi));
                    }
                    for (int s = 0; s < rbs; s++) {
                        int orig_g = dfftt__nl_ix[offset+s] - 1;
                        int m1 = mill[orig_g*3+0], m2 = mill[orig_g*3+1], m3 = mill[orig_g*3+2];
                        Complex_DP sf = cmul(eigqts[na],
                                       cmul(eigts1_T[EIGTS1T_IDX(na,m1)],
                                       cmul(eigts2_T[EIGTS2T_IDX(na,m2)],
                                            eigts3_T[EIGTS3T_IDX(na,m3)])));
                        accum[s] = cadd(accum[s], cmul(aux2_local[s], sf));
                    }
                }
                // Sequential writes to rhoc
                for (int s = 0; s < rbs; s++) {
                    int nl_idx = dfftt__nl_sorted[offset+s] - 1;
                    rhoc[nl_idx] = cadd(rhoc[nl_idx], accum[s]);
                }
            }
        }
    }
}

// ============================================================
// Baseline SoA
// ============================================================
void cpu_addusxx_baseline_soa(
    double* rhoc_re, double* rhoc_im,
    const DP* xkq, const DP* xk, const DP* tau,
    const double* becphi_re, const double* becphi_im,
    const double* becpsi_re, const double* becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const double* qgm_re, const double* qgm_im,
    const double* eigts1_re, const double* eigts1_im,
    const double* eigts2_re, const double* eigts2_im,
    const double* eigts3_re, const double* eigts3_im,
    const int* mill, const int* dfftt__nl)
{
    double eqr[10], eqi[10];
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0]-xkq[0])*tau[na*3+0]
                      + (xk[1]-xkq[1])*tau[na*3+1]
                      + (xk[2]-xkq[2])*tau[na*3+2]);
        eqr[na] = cos(arg); eqi[na] = -sin(arg);
    }

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            double a1r[BLOCKSIZE], a1i[BLOCKSIZE], a2r[BLOCKSIZE], a2i[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    for (int s = 0; s < rbs; s++) { a2r[s] = 0; a2i[s] = 0; }
                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) { a1r[s] = 0; a1i[s] = 0; }
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int idx_base = (nij + ijtoh_val - 1) * QGM_NROWS + offset;
                            double br = becpsi_re[ijkb0+jh], bi = becpsi_im[ijkb0+jh];
                            for (int s = 0; s < rbs; s++) {
                                double qr = qgm_re[idx_base+s], qi = qgm_im[idx_base+s];
                                a1r[s] += SOA_MUL_RE(qr,qi,br,bi);
                                a1i[s] += SOA_MUL_IM(qr,qi,br,bi);
                            }
                        }
                        double pr = becphi_re[ijkb0+ih], pi = becphi_im[ijkb0+ih];
                        for (int s = 0; s < rbs; s++) {
                            a2r[s] += SOA_MULCONJ_RE(a1r[s],a1i[s],pr,pi);
                            a2i[s] += SOA_MULCONJ_IM(a1r[s],a1i[s],pr,pi);
                        }
                    }
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int m1 = mill[g*3+0], m2 = mill[g*3+1], m3 = mill[g*3+2];
                        int e1 = EIGTS1_IDX(m1,na), e2 = EIGTS2_IDX(m2,na), e3 = EIGTS3_IDX(m3,na);
                        double sr = eqr[na], si = eqi[na], tr, ti;
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts1_re[e1],eigts1_im[e1]); si=SOA_MUL_IM(tr,ti,eigts1_re[e1],eigts1_im[e1]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts2_re[e2],eigts2_im[e2]); si=SOA_MUL_IM(tr,ti,eigts2_re[e2],eigts2_im[e2]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts3_re[e3],eigts3_im[e3]); si=SOA_MUL_IM(tr,ti,eigts3_re[e3],eigts3_im[e3]);
                        double rr = SOA_MUL_RE(a2r[s],a2i[s],sr,si);
                        double ri = SOA_MUL_IM(a2r[s],a2i[s],sr,si);
                        int nl_idx = dfftt__nl[g] - 1;
                        rhoc_re[nl_idx] += rr;
                        rhoc_im[nl_idx] += ri;
                    }
                }
            }
        }
    }
}

// ============================================================
// Eigts-transposed SoA
// ============================================================
void cpu_addusxx_eigts_transposed_soa(
    double* rhoc_re, double* rhoc_im,
    const DP* xkq, const DP* xk, const DP* tau,
    const double* becphi_re, const double* becphi_im,
    const double* becpsi_re, const double* becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const double* qgm_re, const double* qgm_im,
    const double* eigts1_T_re, const double* eigts1_T_im,
    const double* eigts2_T_re, const double* eigts2_T_im,
    const double* eigts3_T_re, const double* eigts3_T_im,
    const int* mill, const int* dfftt__nl)
{
    double eqr[10], eqi[10];
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0]-xkq[0])*tau[na*3+0]
                      + (xk[1]-xkq[1])*tau[na*3+1]
                      + (xk[2]-xkq[2])*tau[na*3+2]);
        eqr[na] = cos(arg); eqi[na] = -sin(arg);
    }

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            double a1r[BLOCKSIZE], a1i[BLOCKSIZE];
            double acc_r[BLOCKSIZE], acc_i[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                for (int s = 0; s < rbs; s++) { acc_r[s] = 0; acc_i[s] = 0; }

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    double a2r[BLOCKSIZE], a2i[BLOCKSIZE];
                    for (int s = 0; s < rbs; s++) { a2r[s] = 0; a2i[s] = 0; }

                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) { a1r[s] = 0; a1i[s] = 0; }
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int idx_base = (nij + ijtoh_val - 1) * QGM_NROWS + offset;
                            double br = becpsi_re[ijkb0+jh], bi = becpsi_im[ijkb0+jh];
                            for (int s = 0; s < rbs; s++) {
                                double qr = qgm_re[idx_base+s], qi = qgm_im[idx_base+s];
                                a1r[s] += SOA_MUL_RE(qr,qi,br,bi);
                                a1i[s] += SOA_MUL_IM(qr,qi,br,bi);
                            }
                        }
                        double pr = becphi_re[ijkb0+ih], pi = becphi_im[ijkb0+ih];
                        for (int s = 0; s < rbs; s++) {
                            a2r[s] += SOA_MULCONJ_RE(a1r[s],a1i[s],pr,pi);
                            a2i[s] += SOA_MULCONJ_IM(a1r[s],a1i[s],pr,pi);
                        }
                    }
                    // Transposed eigts + accumulate
                    for (int s = 0; s < rbs; s++) {
                        int g = offset + s;
                        int m1 = mill[g*3+0], m2 = mill[g*3+1], m3 = mill[g*3+2];
                        int e1 = EIGTS1T_IDX(na,m1), e2 = EIGTS2T_IDX(na,m2), e3 = EIGTS3T_IDX(na,m3);
                        double sr = eqr[na], si = eqi[na], tr, ti;
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]); si=SOA_MUL_IM(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]); si=SOA_MUL_IM(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]); si=SOA_MUL_IM(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]);
                        acc_r[s] += SOA_MUL_RE(a2r[s],a2i[s],sr,si);
                        acc_i[s] += SOA_MUL_IM(a2r[s],a2i[s],sr,si);
                    }
                }
                for (int s = 0; s < rbs; s++) {
                    int nl_idx = dfftt__nl[offset+s] - 1;
                    rhoc_re[nl_idx] += acc_r[s];
                    rhoc_im[nl_idx] += acc_i[s];
                }
            }
        }
    }
}

// ============================================================
// Sorted SoA
// ============================================================
void cpu_addusxx_sorted_soa(
    double* rhoc_re, double* rhoc_im,
    const DP* xkq, const DP* xk, const DP* tau,
    const double* becphi_re, const double* becphi_im,
    const double* becpsi_re, const double* becpsi_im,
    int nkb, int ngms, int nat, int ntyp,
    const int* upf_tvanp, const int* nij_type, const int* ityp,
    const int* ofsbeta, const int* nh, const int* ijtoh,
    const double* qgm_re, const double* qgm_im,
    const double* eigts1_T_re, const double* eigts1_T_im,
    const double* eigts2_T_re, const double* eigts2_T_im,
    const double* eigts3_T_re, const double* eigts3_T_im,
    const int* mill, const int* dfftt__nl_sorted, const int* dfftt__nl_ix)
{
    double eqr[10], eqi[10];
    for (int na = 0; na < nat; na++) {
        DP arg = TPI * ((xk[0]-xkq[0])*tau[na*3+0]
                      + (xk[1]-xkq[1])*tau[na*3+1]
                      + (xk[2]-xkq[2])*tau[na*3+2]);
        eqr[na] = cos(arg); eqi[na] = -sin(arg);
    }

    int numblock = (ngms + BLOCKSIZE - 1) / BLOCKSIZE;

    for (int nt = 0; nt < ntyp; nt++) {
        if (upf_tvanp[nt] != 1) continue;
        int nij = nij_type[nt];
        int nh_nt = nh[nt];

        #pragma omp parallel
        {
            double a1r[BLOCKSIZE], a1i[BLOCKSIZE];
            double acc_r[BLOCKSIZE], acc_i[BLOCKSIZE];

            #pragma omp for schedule(static)
            for (int iblock = 0; iblock < numblock; iblock++) {
                int offset = iblock * BLOCKSIZE;
                int rbs = (ngms - offset < BLOCKSIZE) ? (ngms - offset) : BLOCKSIZE;

                for (int s = 0; s < rbs; s++) { acc_r[s] = 0; acc_i[s] = 0; }

                for (int na = 0; na < nat; na++) {
                    if (ityp[na] != (nt + 1)) continue;
                    int ijkb0 = ofsbeta[na];

                    double a2r[BLOCKSIZE], a2i[BLOCKSIZE];
                    for (int s = 0; s < rbs; s++) { a2r[s] = 0; a2i[s] = 0; }

                    for (int ih = 0; ih < nh_nt; ih++) {
                        for (int s = 0; s < rbs; s++) { a1r[s] = 0; a1i[s] = 0; }
                        for (int jh = 0; jh < nh_nt; jh++) {
                            int ijtoh_val = ijtoh[IDX3(ih, jh, nt, IJTOH_N1, IJTOH_N2)];
                            int qgm_col = nij + ijtoh_val - 1;
                            double br = becpsi_re[ijkb0+jh], bi = becpsi_im[ijkb0+jh];
                            for (int s = 0; s < rbs; s++) {
                                int orig_g = dfftt__nl_ix[offset+s] - 1;
                                double qr = qgm_re[qgm_col*QGM_NROWS + orig_g];
                                double qi = qgm_im[qgm_col*QGM_NROWS + orig_g];
                                a1r[s] += SOA_MUL_RE(qr,qi,br,bi);
                                a1i[s] += SOA_MUL_IM(qr,qi,br,bi);
                            }
                        }
                        double pr = becphi_re[ijkb0+ih], pi = becphi_im[ijkb0+ih];
                        for (int s = 0; s < rbs; s++) {
                            a2r[s] += SOA_MULCONJ_RE(a1r[s],a1i[s],pr,pi);
                            a2i[s] += SOA_MULCONJ_IM(a1r[s],a1i[s],pr,pi);
                        }
                    }
                    for (int s = 0; s < rbs; s++) {
                        int orig_g = dfftt__nl_ix[offset+s] - 1;
                        int m1 = mill[orig_g*3+0], m2 = mill[orig_g*3+1], m3 = mill[orig_g*3+2];
                        int e1 = EIGTS1T_IDX(na,m1), e2 = EIGTS2T_IDX(na,m2), e3 = EIGTS3T_IDX(na,m3);
                        double sr = eqr[na], si = eqi[na], tr, ti;
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]); si=SOA_MUL_IM(tr,ti,eigts1_T_re[e1],eigts1_T_im[e1]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]); si=SOA_MUL_IM(tr,ti,eigts2_T_re[e2],eigts2_T_im[e2]);
                        tr=sr; ti=si; sr=SOA_MUL_RE(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]); si=SOA_MUL_IM(tr,ti,eigts3_T_re[e3],eigts3_T_im[e3]);
                        acc_r[s] += SOA_MUL_RE(a2r[s],a2i[s],sr,si);
                        acc_i[s] += SOA_MUL_IM(a2r[s],a2i[s],sr,si);
                    }
                }
                // Sequential writes
                for (int s = 0; s < rbs; s++) {
                    int nl_idx = dfftt__nl_sorted[offset+s] - 1;
                    rhoc_re[nl_idx] += acc_r[s];
                    rhoc_im[nl_idx] += acc_i[s];
                }
            }
        }
    }
}