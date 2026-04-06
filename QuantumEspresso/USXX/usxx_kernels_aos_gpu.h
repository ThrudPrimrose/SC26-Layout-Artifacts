#pragma once
#include "usxx_kernels_gpu.h"
#include <cmath>

__global__ void kernel_eigqts(Complex_DP* eigqts, const DP* xkq, const DP* xk,
                              const DP* tau, int nat) {
    int na = blockIdx.x * blockDim.x + threadIdx.x;
    if (na >= nat) return;
    DP arg = TPI * ((xk[0]-xkq[0])*tau[na*3+0]
                  + (xk[1]-xkq[1])*tau[na*3+1]
                  + (xk[2]-xkq[2])*tau[na*3+2]);
    eigqts[na] = make_cmplx(cos(arg), -sin(arg));
}

// Body appears twice in source (unrolled / not), only one compiles per instantiation.
#define COARSEN_LOOP_BEGIN \
    if constexpr (UNROLL_C) { _Pragma("unroll") \
        for (int c = 0; c < COARSEN; c++) {
#define COARSEN_LOOP_MID \
        } \
    } else { \
        for (int c = 0; c < COARSEN; c++) {
#define COARSEN_LOOP_END \
        } \
    }

// ============================================================
// (1) Baseline AoS
// ============================================================
template<int COARSEN, bool UNROLL_C>
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
    const Complex_DP* __restrict__ eigqts)
{
    constexpr int R=55191, N1=19, N2=19;
    int g_start = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=g_start+c; if(gi<ngms){ \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; int nl=dfftt__nl[gi]-1; \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt_1based)continue; int ij=ofsbeta[na]; \
    Complex_DP a2=make_cmplx(0,0); \
    for(int ih=0;ih<nh_nt;ih++){Complex_DP a1=make_cmplx(0,0); \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijtoh[IDX3(ih,jh,nt_1based-1,N1,N2)]-1; \
    a1=cadd(a1,cmul(qgm[col*R+gi],becpsi_c[ij+jh]));} \
    a2=cadd(a2,cmul(a1,cconj(becphi_c[ij+ih])));} \
    Complex_DP sf=cmul(eigqts[na],cmul(eigts1[EIGTS1_IDX(m1,na)], \
    cmul(eigts2[EIGTS2_IDX(m2,na)],eigts3[EIGTS3_IDX(m3,na)]))); \
    rhoc[nl]=cadd(rhoc[nl],cmul(a2,sf));}}
    COARSEN_LOOP_BEGIN B COARSEN_LOOP_MID B COARSEN_LOOP_END
#undef B
}

// ============================================================
// (2) Eigts-transposed AoS
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_eigts_t(
    Complex_DP* __restrict__ rhoc,
    const Complex_DP* __restrict__ becphi_c,
    const Complex_DP* __restrict__ becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const Complex_DP* __restrict__ qgm,
    const Complex_DP* __restrict__ eigts1_T,
    const Complex_DP* __restrict__ eigts2_T,
    const Complex_DP* __restrict__ eigts3_T,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const Complex_DP* __restrict__ eigqts)
{
    constexpr int R=55191, N1=19, N2=19;
    int g_start = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=g_start+c; if(gi<ngms){ \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; Complex_DP acc=make_cmplx(0,0); \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt_1based)continue; int ij=ofsbeta[na]; \
    Complex_DP a2=make_cmplx(0,0); \
    for(int ih=0;ih<nh_nt;ih++){Complex_DP a1=make_cmplx(0,0); \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijtoh[IDX3(ih,jh,nt_1based-1,N1,N2)]-1; \
    a1=cadd(a1,cmul(qgm[col*R+gi],becpsi_c[ij+jh]));} \
    a2=cadd(a2,cmul(a1,cconj(becphi_c[ij+ih])));} \
    Complex_DP sf=cmul(eigqts[na],cmul(eigts1_T[EIGTS1T_IDX(na,m1)], \
    cmul(eigts2_T[EIGTS2T_IDX(na,m2)],eigts3_T[EIGTS3T_IDX(na,m3)]))); \
    acc=cadd(acc,cmul(a2,sf));} \
    int nl=dfftt__nl[gi]-1; rhoc[nl]=cadd(rhoc[nl],acc);}
    COARSEN_LOOP_BEGIN B COARSEN_LOOP_MID B COARSEN_LOOP_END
#undef B
}

// ============================================================
// (3) Shared-bec AoS
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_shared_bec(
    Complex_DP* __restrict__ rhoc,
    const Complex_DP* __restrict__ becphi_c,
    const Complex_DP* __restrict__ becpsi_c,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt_1based,
    const int* __restrict__ ityp, const int* __restrict__ ofsbeta,
    const int* __restrict__ ijtoh,
    const Complex_DP* __restrict__ qgm,
    const Complex_DP* __restrict__ eigts1_T,
    const Complex_DP* __restrict__ eigts2_T,
    const Complex_DP* __restrict__ eigts3_T,
    const int* __restrict__ mill, const int* __restrict__ dfftt__nl,
    const Complex_DP* __restrict__ eigqts)
{
    __shared__ Complex_DP s_bpsi[20], s_bphi[20];
    constexpr int R=55191, N1=19, N2=19;
    int g_start = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
    Complex_DP acc[COARSEN];

#define BINIT acc[c]=make_cmplx(0,0);
    COARSEN_LOOP_BEGIN BINIT COARSEN_LOOP_MID BINIT COARSEN_LOOP_END
#undef BINIT

    for(int na=0;na<nat;na++){
        if(ityp[na]!=nt_1based) continue;
        int ij=ofsbeta[na];
        if(threadIdx.x<nh_nt){
            s_bpsi[threadIdx.x]=becpsi_c[ij+threadIdx.x];
            s_bphi[threadIdx.x]=becphi_c[ij+threadIdx.x];
        }
        __syncthreads();

#define BCOMP \
    int gi=g_start+c; if(gi<ngms){ \
    Complex_DP a2=make_cmplx(0,0); \
    for(int ih=0;ih<nh_nt;ih++){Complex_DP a1=make_cmplx(0,0); \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijtoh[IDX3(ih,jh,nt_1based-1,N1,N2)]-1; \
    a1=cadd(a1,cmul(qgm[col*R+gi],s_bpsi[jh]));} \
    a2=cadd(a2,cmul(a1,cconj(s_bphi[ih])));} \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; \
    Complex_DP sf=cmul(eigqts[na],cmul(eigts1_T[EIGTS1T_IDX(na,m1)], \
    cmul(eigts2_T[EIGTS2T_IDX(na,m2)],eigts3_T[EIGTS3T_IDX(na,m3)]))); \
    acc[c]=cadd(acc[c],cmul(a2,sf));}
        COARSEN_LOOP_BEGIN BCOMP COARSEN_LOOP_MID BCOMP COARSEN_LOOP_END
#undef BCOMP

        __syncthreads();
    }

#define BWRITE int gi=g_start+c; if(gi<ngms){int nl=dfftt__nl[gi]-1; rhoc[nl]=cadd(rhoc[nl],acc[c]);}
    COARSEN_LOOP_BEGIN BWRITE COARSEN_LOOP_MID BWRITE COARSEN_LOOP_END
#undef BWRITE
}

// ============================================================
// (4) Optimized AoS — sorted fused
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_optimized(
    Complex_DP* __restrict__ rhoc,
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
    const int* __restrict__ dfftt__nl_sorted,
    const int* __restrict__ dfftt__nl_ix,
    const Complex_DP* __restrict__ eigqts)
{
    constexpr int QT=397, N1=19, N2=19;
    int g_start = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=g_start+c; if(gi<ngms){ \
    int nl=dfftt__nl_sorted[gi]-1, og=dfftt__nl_ix[gi]-1; \
    int m1=mill[og*3],m2=mill[og*3+1],m3=mill[og*3+2]; Complex_DP acc=make_cmplx(0,0); \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt_1based)continue; int ij=ofsbeta[na]; \
    Complex_DP a2=make_cmplx(0,0); \
    for(int ih=0;ih<nh_nt;ih++){Complex_DP a1=make_cmplx(0,0); \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijtoh[IDX3(ih,jh,nt_1based-1,N1,N2)]-1; \
    a1=cadd(a1,cmul(qgm_T[og*QT+col],becpsi_c[ij+jh]));} \
    a2=cadd(a2,cmul(a1,cconj(becphi_c[ij+ih])));} \
    Complex_DP sf=cmul(eigqts[na],cmul(eigts1_T[EIGTS1T_IDX(na,m1)], \
    cmul(eigts2_T[EIGTS2T_IDX(na,m2)],eigts3_T[EIGTS3T_IDX(na,m3)]))); \
    acc=cadd(acc,cmul(a2,sf));} \
    rhoc[nl]=cadd(rhoc[nl],acc);}
    COARSEN_LOOP_BEGIN B COARSEN_LOOP_MID B COARSEN_LOOP_END
#undef B
}

#undef COARSEN_LOOP_BEGIN
#undef COARSEN_LOOP_MID
#undef COARSEN_LOOP_END

// ============================================================
// AoS host wrappers
// ============================================================
inline void addusxx_g_gpu(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofsbeta, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1,
    const Complex_DP* d_eigts2, const Complex_DP* d_eigts3,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, Complex_DP* d_eigqts)
{
    kernel_eigqts<<<(nat+tbs-1)/tbs,tbs>>>(d_eigqts,d_xkq,d_xk,d_tau,nat);
    int nl=(ngms+coarsen-1)/coarsen, bl=(nl+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_baseline,bl,tbs,coarsen,unroll,
            d_rhoc,d_becphi_c,d_becpsi_c,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofsbeta,d_ijtoh,
            d_qgm,d_eigts1,d_eigts2,d_eigts3,d_mill,d_dfftt__nl,d_eigqts);}
}

inline void addusxx_g_gpu_eigts_t(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofsbeta, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, Complex_DP* d_eigqts)
{
    kernel_eigqts<<<(nat+tbs-1)/tbs,tbs>>>(d_eigqts,d_xkq,d_xk,d_tau,nat);
    int nl=(ngms+coarsen-1)/coarsen, bl=(nl+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_eigts_t,bl,tbs,coarsen,unroll,
            d_rhoc,d_becphi_c,d_becpsi_c,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofsbeta,d_ijtoh,
            d_qgm,d_eigts1_T,d_eigts2_T,d_eigts3_T,d_mill,d_dfftt__nl,d_eigqts);}
}

inline void addusxx_g_gpu_shared_bec(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofsbeta, const int* d_ijtoh,
    const Complex_DP* d_qgm, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, Complex_DP* d_eigqts)
{
    kernel_eigqts<<<(nat+tbs-1)/tbs,tbs>>>(d_eigqts,d_xkq,d_xk,d_tau,nat);
    int nl=(ngms+coarsen-1)/coarsen, bl=(nl+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_shared_bec,bl,tbs,coarsen,unroll,
            d_rhoc,d_becphi_c,d_becpsi_c,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofsbeta,d_ijtoh,
            d_qgm,d_eigts1_T,d_eigts2_T,d_eigts3_T,d_mill,d_dfftt__nl,d_eigqts);}
}

inline void addusxx_g_gpu_optimized(
    Complex_DP* d_rhoc, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const Complex_DP* d_becphi_c, const Complex_DP* d_becpsi_c,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofsbeta, const int* d_ijtoh,
    const Complex_DP* d_qgm_T, const Complex_DP* d_eigts1_T,
    const Complex_DP* d_eigts2_T, const Complex_DP* d_eigts3_T,
    const int* d_mill, const int* d_dfftt__nl_sorted, const int* d_dfftt__nl_ix,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, Complex_DP* d_eigqts)
{
    kernel_eigqts<<<(nat+tbs-1)/tbs,tbs>>>(d_eigqts,d_xkq,d_xk,d_tau,nat);
    int nl=(ngms+coarsen-1)/coarsen, bl=(nl+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_optimized,bl,tbs,coarsen,unroll,
            d_rhoc,d_becphi_c,d_becpsi_c,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofsbeta,d_ijtoh,
            d_qgm_T,d_eigts1_T,d_eigts2_T,d_eigts3_T,
            d_mill,d_dfftt__nl_sorted,d_dfftt__nl_ix,d_eigqts);}
}