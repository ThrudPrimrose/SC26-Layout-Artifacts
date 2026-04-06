#pragma once
#include "types.h"

// ============================================================
// Indexing helpers (unchanged)
// ============================================================
#define IDX3(i, j, k, n1, n2) ((k)*(n1)*(n2) + (j)*(n1) + (i))
#define EIGTS1_IDX(ig, na)  ((na) * 217 + ((ig) + 108))
#define EIGTS2_IDX(ig, na)  ((na) * 109 + ((ig) + 54))
#define EIGTS3_IDX(ig, na)  ((na) * 109 + ((ig) + 54))
#define EIGTS1T_IDX(na, ig) (((ig) + 108) * 10 + (na))
#define EIGTS2T_IDX(na, ig) (((ig) + 54) * 10 + (na))
#define EIGTS3T_IDX(na, ig) (((ig) + 54) * 10 + (na))

// ============================================================
// SoA arithmetic
// ============================================================
#define SOA_MUL_RE(ar,ai,br,bi)     ((ar)*(br)-(ai)*(bi))
#define SOA_MUL_IM(ar,ai,br,bi)     ((ar)*(bi)+(ai)*(br))
#define SOA_MULCONJ_RE(ar,ai,br,bi) ((ar)*(br)+(ai)*(bi))
#define SOA_MULCONJ_IM(ar,ai,br,bi) ((ai)*(br)-(ar)*(bi))

// ============================================================
// No device lambdas — kernels use if constexpr with duplicated
// loop header to apply #pragma unroll conditionally.
// ============================================================

// ============================================================
// Dispatch macro: runtime (coarsen, unroll) → template instantiation
// kern must be template<int COARSEN, bool UNROLL_C>
// ============================================================
#define DISPATCH_KERNEL(kern, grid, tbs, coarsen, unroll, ...) do { \
    if (unroll) { \
        switch(coarsen) { \
        case 1: kern<1,true><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 2: kern<2,true><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 4: kern<4,true><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 8: kern<8,true><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 16: kern<16,true><<<grid,tbs>>>(__VA_ARGS__); break; \
        } \
    } else { \
        switch(coarsen) { \
        case 1: kern<1,false><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 2: kern<2,false><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 4: kern<4,false><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 8: kern<8,false><<<grid,tbs>>>(__VA_ARGS__); break; \
        case 16: kern<16,false><<<grid,tbs>>>(__VA_ARGS__); break; \
        } \
    } \
} while(0)

// ============================================================
// CPU baseline (non-template, no GPU)
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
// Utility
// ============================================================
inline void aos_to_soa(const Complex_DP* aos, double* re, double* im, int n) {
    for (int i = 0; i < n; i++) { re[i] = creal_val(aos[i]); im[i] = cimag_val(aos[i]); }
}