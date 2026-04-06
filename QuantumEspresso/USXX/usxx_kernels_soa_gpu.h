#pragma once
#include "usxx_kernels_gpu.h"

static constexpr int QGM_R  = 55191;
static constexpr int QGMT_R = 397;
static constexpr int IJ_N1  = 19, IJ_N2 = 19;

__global__ void kernel_eigqts_soa(
    double* er, double* ei, const double* xkq, const double* xk,
    const double* tau, int nat) {
    int na = blockIdx.x * blockDim.x + threadIdx.x;
    if (na >= nat) return;
    double arg = TPI * ((xk[0]-xkq[0])*tau[na*3]
                      + (xk[1]-xkq[1])*tau[na*3+1]
                      + (xk[2]-xkq[2])*tau[na*3+2]);
    er[na] = cos(arg); ei[na] = -sin(arg);
}

#define CFOR_BEGIN \
    if constexpr (UNROLL_C) { _Pragma("unroll") \
        for (int c = 0; c < COARSEN; c++) {
#define CFOR_MID \
        } \
    } else { \
        for (int c = 0; c < COARSEN; c++) {
#define CFOR_END \
        } \
    }

// ============================================================
// (1) Baseline SoA
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_baseline_soa(
    double* __restrict__ yr, double* __restrict__ yi,
    const double* __restrict__ pr, const double* __restrict__ pi,
    const double* __restrict__ sr, const double* __restrict__ si,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt1,
    const int* __restrict__ ityp, const int* __restrict__ ofs,
    const int* __restrict__ ijt,
    const double* __restrict__ qr, const double* __restrict__ qi,
    const double* __restrict__ e1r, const double* __restrict__ e1i,
    const double* __restrict__ e2r, const double* __restrict__ e2i,
    const double* __restrict__ e3r, const double* __restrict__ e3i,
    const int* __restrict__ mill, const int* __restrict__ nl,
    const double* __restrict__ eqr, const double* __restrict__ eqi)
{
    int gs = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=gs+c; if(gi<ngms){ \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; int nli=nl[gi]-1; \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt1)continue; int ij=ofs[na]; \
    double a2r=0,a2i=0; \
    for(int ih=0;ih<nh_nt;ih++){double a1r=0,a1i=0; \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijt[IDX3(ih,jh,nt1-1,IJ_N1,IJ_N2)]-1; \
    int ix=col*QGM_R+gi; \
    a1r+=SOA_MUL_RE(qr[ix],qi[ix],sr[ij+jh],si[ij+jh]); \
    a1i+=SOA_MUL_IM(qr[ix],qi[ix],sr[ij+jh],si[ij+jh]);} \
    a2r+=SOA_MULCONJ_RE(a1r,a1i,pr[ij+ih],pi[ij+ih]); \
    a2i+=SOA_MULCONJ_IM(a1r,a1i,pr[ij+ih],pi[ij+ih]);} \
    int i1=EIGTS1_IDX(m1,na),i2=EIGTS2_IDX(m2,na),i3=EIGTS3_IDX(m3,na); \
    double xr=eqr[na],xi_=eqi[na],tr,ti; \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e1r[i1],e1i[i1]);xi_=SOA_MUL_IM(tr,ti,e1r[i1],e1i[i1]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e2r[i2],e2i[i2]);xi_=SOA_MUL_IM(tr,ti,e2r[i2],e2i[i2]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e3r[i3],e3i[i3]);xi_=SOA_MUL_IM(tr,ti,e3r[i3],e3i[i3]); \
    yr[nli]+=SOA_MUL_RE(a2r,a2i,xr,xi_); yi[nli]+=SOA_MUL_IM(a2r,a2i,xr,xi_);}}
    CFOR_BEGIN B CFOR_MID B CFOR_END
#undef B
}

// ============================================================
// (2) Eigts-transposed SoA
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_eigts_t_soa(
    double* __restrict__ yr, double* __restrict__ yi,
    const double* __restrict__ pr, const double* __restrict__ pi,
    const double* __restrict__ sr, const double* __restrict__ si,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt1,
    const int* __restrict__ ityp, const int* __restrict__ ofs,
    const int* __restrict__ ijt,
    const double* __restrict__ qr, const double* __restrict__ qi,
    const double* __restrict__ e1r, const double* __restrict__ e1i,
    const double* __restrict__ e2r, const double* __restrict__ e2i,
    const double* __restrict__ e3r, const double* __restrict__ e3i,
    const int* __restrict__ mill, const int* __restrict__ nl,
    const double* __restrict__ eqr, const double* __restrict__ eqi)
{
    int gs = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=gs+c; if(gi<ngms){ \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; double acr=0,aci=0; \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt1)continue; int ij=ofs[na]; \
    double a2r=0,a2i=0; \
    for(int ih=0;ih<nh_nt;ih++){double a1r=0,a1i=0; \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijt[IDX3(ih,jh,nt1-1,IJ_N1,IJ_N2)]-1; \
    int ix=col*QGM_R+gi; \
    a1r+=SOA_MUL_RE(qr[ix],qi[ix],sr[ij+jh],si[ij+jh]); \
    a1i+=SOA_MUL_IM(qr[ix],qi[ix],sr[ij+jh],si[ij+jh]);} \
    a2r+=SOA_MULCONJ_RE(a1r,a1i,pr[ij+ih],pi[ij+ih]); \
    a2i+=SOA_MULCONJ_IM(a1r,a1i,pr[ij+ih],pi[ij+ih]);} \
    int i1=EIGTS1T_IDX(na,m1),i2=EIGTS2T_IDX(na,m2),i3=EIGTS3T_IDX(na,m3); \
    double xr=eqr[na],xi_=eqi[na],tr,ti; \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e1r[i1],e1i[i1]);xi_=SOA_MUL_IM(tr,ti,e1r[i1],e1i[i1]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e2r[i2],e2i[i2]);xi_=SOA_MUL_IM(tr,ti,e2r[i2],e2i[i2]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e3r[i3],e3i[i3]);xi_=SOA_MUL_IM(tr,ti,e3r[i3],e3i[i3]); \
    acr+=SOA_MUL_RE(a2r,a2i,xr,xi_); aci+=SOA_MUL_IM(a2r,a2i,xr,xi_);} \
    int nli=nl[gi]-1; yr[nli]+=acr; yi[nli]+=aci;}
    CFOR_BEGIN B CFOR_MID B CFOR_END
#undef B
}

// ============================================================
// (3) Shared-bec SoA
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_shared_bec_soa(
    double* __restrict__ yr, double* __restrict__ yi,
    const double* __restrict__ pr, const double* __restrict__ pi,
    const double* __restrict__ sr, const double* __restrict__ si,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt1,
    const int* __restrict__ ityp, const int* __restrict__ ofs,
    const int* __restrict__ ijt,
    const double* __restrict__ qr, const double* __restrict__ qi,
    const double* __restrict__ e1r, const double* __restrict__ e1i,
    const double* __restrict__ e2r, const double* __restrict__ e2i,
    const double* __restrict__ e3r, const double* __restrict__ e3i,
    const int* __restrict__ mill, const int* __restrict__ nl,
    const double* __restrict__ eqr, const double* __restrict__ eqi)
{
    __shared__ double spr[20],spi[20],ssr[20],ssi[20];
    int gs = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
    double acr[COARSEN], aci[COARSEN];

#define BINIT acr[c]=0; aci[c]=0;
    CFOR_BEGIN BINIT CFOR_MID BINIT CFOR_END
#undef BINIT

    for(int na=0;na<nat;na++){
        if(ityp[na]!=nt1) continue;
        int ij=ofs[na];
        if(threadIdx.x<nh_nt){
            ssr[threadIdx.x]=sr[ij+threadIdx.x]; ssi[threadIdx.x]=si[ij+threadIdx.x];
            spr[threadIdx.x]=pr[ij+threadIdx.x]; spi[threadIdx.x]=pi[ij+threadIdx.x];
        }
        __syncthreads();

#define BCOMP \
    int gi=gs+c; if(gi<ngms){ \
    double a2r=0,a2i=0; \
    for(int ih=0;ih<nh_nt;ih++){double a1r=0,a1i=0; \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijt[IDX3(ih,jh,nt1-1,IJ_N1,IJ_N2)]-1; \
    int ix=col*QGM_R+gi; \
    a1r+=SOA_MUL_RE(qr[ix],qi[ix],ssr[jh],ssi[jh]); \
    a1i+=SOA_MUL_IM(qr[ix],qi[ix],ssr[jh],ssi[jh]);} \
    a2r+=SOA_MULCONJ_RE(a1r,a1i,spr[ih],spi[ih]); \
    a2i+=SOA_MULCONJ_IM(a1r,a1i,spr[ih],spi[ih]);} \
    int m1=mill[gi*3],m2=mill[gi*3+1],m3=mill[gi*3+2]; \
    int i1=EIGTS1T_IDX(na,m1),i2=EIGTS2T_IDX(na,m2),i3=EIGTS3T_IDX(na,m3); \
    double xr=eqr[na],xi_=eqi[na],tr,ti; \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e1r[i1],e1i[i1]);xi_=SOA_MUL_IM(tr,ti,e1r[i1],e1i[i1]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e2r[i2],e2i[i2]);xi_=SOA_MUL_IM(tr,ti,e2r[i2],e2i[i2]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e3r[i3],e3i[i3]);xi_=SOA_MUL_IM(tr,ti,e3r[i3],e3i[i3]); \
    acr[c]+=SOA_MUL_RE(a2r,a2i,xr,xi_); aci[c]+=SOA_MUL_IM(a2r,a2i,xr,xi_);}
        CFOR_BEGIN BCOMP CFOR_MID BCOMP CFOR_END
#undef BCOMP

        __syncthreads();
    }

#define BWRITE int gi=gs+c; if(gi<ngms){int i=nl[gi]-1; yr[i]+=acr[c]; yi[i]+=aci[c];}
    CFOR_BEGIN BWRITE CFOR_MID BWRITE CFOR_END
#undef BWRITE
}

// ============================================================
// (4) Optimized SoA — sorted fused
// ============================================================
template<int COARSEN, bool UNROLL_C>
__global__ void kernel_addusxx_optimized_soa(
    double* __restrict__ yr, double* __restrict__ yi,
    const double* __restrict__ pr, const double* __restrict__ pi,
    const double* __restrict__ sr, const double* __restrict__ si,
    int nkb, int ngms, int nat,
    int nij, int nh_nt, int nt1,
    const int* __restrict__ ityp, const int* __restrict__ ofs,
    const int* __restrict__ ijt,
    const double* __restrict__ qTr, const double* __restrict__ qTi,
    const double* __restrict__ e1r, const double* __restrict__ e1i,
    const double* __restrict__ e2r, const double* __restrict__ e2i,
    const double* __restrict__ e3r, const double* __restrict__ e3i,
    const int* __restrict__ mill,
    const int* __restrict__ nls, const int* __restrict__ nlix,
    const double* __restrict__ eqr, const double* __restrict__ eqi)
{
    int gs = (blockIdx.x*blockDim.x+threadIdx.x)*COARSEN;
#define B \
    int gi=gs+c; if(gi<ngms){ \
    int nli=nls[gi]-1, og=nlix[gi]-1; \
    int m1=mill[og*3],m2=mill[og*3+1],m3=mill[og*3+2]; double acr=0,aci=0; \
    for(int na=0;na<nat;na++){if(ityp[na]!=nt1)continue; int ij=ofs[na]; \
    double a2r=0,a2i=0; \
    for(int ih=0;ih<nh_nt;ih++){double a1r=0,a1i=0; \
    for(int jh=0;jh<nh_nt;jh++){int col=nij+ijt[IDX3(ih,jh,nt1-1,IJ_N1,IJ_N2)]-1; \
    int ix=og*QGMT_R+col; \
    a1r+=SOA_MUL_RE(qTr[ix],qTi[ix],sr[ij+jh],si[ij+jh]); \
    a1i+=SOA_MUL_IM(qTr[ix],qTi[ix],sr[ij+jh],si[ij+jh]);} \
    a2r+=SOA_MULCONJ_RE(a1r,a1i,pr[ij+ih],pi[ij+ih]); \
    a2i+=SOA_MULCONJ_IM(a1r,a1i,pr[ij+ih],pi[ij+ih]);} \
    int i1=EIGTS1T_IDX(na,m1),i2=EIGTS2T_IDX(na,m2),i3=EIGTS3T_IDX(na,m3); \
    double xr=eqr[na],xi_=eqi[na],tr,ti; \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e1r[i1],e1i[i1]);xi_=SOA_MUL_IM(tr,ti,e1r[i1],e1i[i1]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e2r[i2],e2i[i2]);xi_=SOA_MUL_IM(tr,ti,e2r[i2],e2i[i2]); \
    tr=xr;ti=xi_;xr=SOA_MUL_RE(tr,ti,e3r[i3],e3i[i3]);xi_=SOA_MUL_IM(tr,ti,e3r[i3],e3i[i3]); \
    acr+=SOA_MUL_RE(a2r,a2i,xr,xi_); aci+=SOA_MUL_IM(a2r,a2i,xr,xi_);} \
    yr[nli]+=acr; yi[nli]+=aci;}
    CFOR_BEGIN B CFOR_MID B CFOR_END
#undef B
}

#undef CFOR_BEGIN
#undef CFOR_MID
#undef CFOR_END

// ============================================================
// SoA host wrappers
// ============================================================
inline void addusxx_g_gpu_soa(
    double* d_yr, double* d_yi, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_pr, const double* d_pi, const double* d_sr, const double* d_si,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofs, const int* d_ijt,
    const double* d_qr, const double* d_qi,
    const double* d_e1r, const double* d_e1i,
    const double* d_e2r, const double* d_e2i,
    const double* d_e3r, const double* d_e3i,
    const int* d_mill, const int* d_nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, double* d_eqr, double* d_eqi)
{
    kernel_eigqts_soa<<<(nat+tbs-1)/tbs,tbs>>>(d_eqr,d_eqi,d_xkq,d_xk,d_tau,nat);
    int nl_=(ngms+coarsen-1)/coarsen, bl=(nl_+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_baseline_soa,bl,tbs,coarsen,unroll,
            d_yr,d_yi,d_pr,d_pi,d_sr,d_si,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofs,d_ijt,
            d_qr,d_qi,d_e1r,d_e1i,d_e2r,d_e2i,d_e3r,d_e3i,
            d_mill,d_nl,d_eqr,d_eqi);}
}

inline void addusxx_g_gpu_eigts_t_soa(
    double* d_yr, double* d_yi, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_pr, const double* d_pi, const double* d_sr, const double* d_si,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofs, const int* d_ijt,
    const double* d_qr, const double* d_qi,
    const double* d_e1Tr, const double* d_e1Ti,
    const double* d_e2Tr, const double* d_e2Ti,
    const double* d_e3Tr, const double* d_e3Ti,
    const int* d_mill, const int* d_nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, double* d_eqr, double* d_eqi)
{
    kernel_eigqts_soa<<<(nat+tbs-1)/tbs,tbs>>>(d_eqr,d_eqi,d_xkq,d_xk,d_tau,nat);
    int nl_=(ngms+coarsen-1)/coarsen, bl=(nl_+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_eigts_t_soa,bl,tbs,coarsen,unroll,
            d_yr,d_yi,d_pr,d_pi,d_sr,d_si,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofs,d_ijt,
            d_qr,d_qi,d_e1Tr,d_e1Ti,d_e2Tr,d_e2Ti,d_e3Tr,d_e3Ti,
            d_mill,d_nl,d_eqr,d_eqi);}
}

inline void addusxx_g_gpu_shared_bec_soa(
    double* d_yr, double* d_yi, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_pr, const double* d_pi, const double* d_sr, const double* d_si,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofs, const int* d_ijt,
    const double* d_qr, const double* d_qi,
    const double* d_e1Tr, const double* d_e1Ti,
    const double* d_e2Tr, const double* d_e2Ti,
    const double* d_e3Tr, const double* d_e3Ti,
    const int* d_mill, const int* d_nl,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, double* d_eqr, double* d_eqi)
{
    kernel_eigqts_soa<<<(nat+tbs-1)/tbs,tbs>>>(d_eqr,d_eqi,d_xkq,d_xk,d_tau,nat);
    int nl_=(ngms+coarsen-1)/coarsen, bl=(nl_+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_shared_bec_soa,bl,tbs,coarsen,unroll,
            d_yr,d_yi,d_pr,d_pi,d_sr,d_si,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofs,d_ijt,
            d_qr,d_qi,d_e1Tr,d_e1Ti,d_e2Tr,d_e2Ti,d_e3Tr,d_e3Ti,
            d_mill,d_nl,d_eqr,d_eqi);}
}

inline void addusxx_g_gpu_optimized_soa(
    double* d_yr, double* d_yi, const DP* d_xkq, const DP* d_xk, const DP* d_tau,
    const double* d_pr, const double* d_pi, const double* d_sr, const double* d_si,
    int nkb, int ngms, int nat, int ntyp,
    const int* d_ityp, const int* d_ofs, const int* d_ijt,
    const double* d_qTr, const double* d_qTi,
    const double* d_e1Tr, const double* d_e1Ti,
    const double* d_e2Tr, const double* d_e2Ti,
    const double* d_e3Tr, const double* d_e3Ti,
    const int* d_mill, const int* d_nls, const int* d_nlix,
    const int* h_tvanp, const int* h_nij, const int* h_nh,
    int tbs, int coarsen, bool unroll, double* d_eqr, double* d_eqi)
{
    kernel_eigqts_soa<<<(nat+tbs-1)/tbs,tbs>>>(d_eqr,d_eqi,d_xkq,d_xk,d_tau,nat);
    int nl_=(ngms+coarsen-1)/coarsen, bl=(nl_+tbs-1)/tbs;
    for(int nt=0;nt<ntyp;nt++){if(h_tvanp[nt]!=1)continue;
        DISPATCH_KERNEL(kernel_addusxx_optimized_soa,bl,tbs,coarsen,unroll,
            d_yr,d_yi,d_pr,d_pi,d_sr,d_si,nkb,ngms,nat,
            h_nij[nt],h_nh[nt],nt+1,d_ityp,d_ofs,d_ijt,
            d_qTr,d_qTi,d_e1Tr,d_e1Ti,d_e2Tr,d_e2Ti,d_e3Tr,d_e3Ti,
            d_mill,d_nls,d_nlix,d_eqr,d_eqi);}
}