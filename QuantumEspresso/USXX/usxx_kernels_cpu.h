#pragma once
#include <cmath>
#include <cstring>

// ============================================================
// Types (no GPU headers)
// ============================================================
using DP = double;
struct Complex_DP { double x, y; };

inline Complex_DP make_cmplx(double r, double i) { Complex_DP c; c.x = r; c.y = i; return c; }
inline double creal_val(Complex_DP a) { return a.x; }
inline double cimag_val(Complex_DP a) { return a.y; }
inline Complex_DP cmul(Complex_DP a, Complex_DP b) {
    return {a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x};
}
inline Complex_DP cadd(Complex_DP a, Complex_DP b) { return {a.x+b.x, a.y+b.y}; }
inline Complex_DP csub(Complex_DP a, Complex_DP b) { return {a.x-b.x, a.y-b.y}; }
inline Complex_DP cconj(Complex_DP a) { return {a.x, -a.y}; }
inline double cabs_val(Complex_DP a) { return std::sqrt(a.x*a.x + a.y*a.y); }

static constexpr double PI  = 3.14159265358979323846;
static constexpr double TPI = 2.0 * PI;

// ============================================================
// Indexing helpers
// ============================================================
#define IDX3(i, j, k, n1, n2) ((k) * (n1) * (n2) + (j) * (n1) + (i))

#define EIGTS1_IDX(ig, na) ((na) * 217 + ((ig) + 108))
#define EIGTS2_IDX(ig, na) ((na) * 109 + ((ig) + 54))
#define EIGTS3_IDX(ig, na) ((na) * 109 + ((ig) + 54))

#define EIGTS1T_IDX(na, ig) (((ig) + 108) * 10 + (na))
#define EIGTS2T_IDX(na, ig) (((ig) + 54) * 10 + (na))
#define EIGTS3T_IDX(na, ig) (((ig) + 54) * 10 + (na))

// ============================================================
// SoA arithmetic
// ============================================================
#define SOA_MUL_RE(ar, ai, br, bi) ((ar)*(br) - (ai)*(bi))
#define SOA_MUL_IM(ar, ai, br, bi) ((ar)*(bi) + (ai)*(br))
#define SOA_MULCONJ_RE(ar, ai, br, bi) ((ar)*(br) + (ai)*(bi))
#define SOA_MULCONJ_IM(ar, ai, br, bi) ((ai)*(br) - (ar)*(bi))

// ============================================================
// CPU kernel declarations
// ============================================================
// Templated see impl.h
// Utility
inline void aos_to_soa(const Complex_DP* aos, double* re, double* im, int n) {
    for (int i = 0; i < n; i++) { re[i] = aos[i].x; im[i] = aos[i].y; }
}