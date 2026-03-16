// copy_cast.h — mixed-type CopyND specialization for precision lowering.
//
// DaCe's CopyND<T>::ConstDst::Copy requires matching src/dst types.
// When arrays are excluded from lowering (stay double) but downstream
// scalars are lowered (float), the codegen emits CopyND<float>::Copy(double*,...).
// This specialization handles the implicit cast.
#pragma once

#include <dace/copy.h>

namespace dace {

template <>
struct CopyND<float, 1, false, 1> {
    template <int DST_STRIDE, int... OTHER_DSTDIMS>
    struct ConstDst {
        // Same-type copy
        template <typename... Args>
        static DACE_HDFI void Copy(const float *src, float *dst,
                                   const int &src_stride,
                                   const Args &...src_otherdims) {
            dst[0] = src[0];
        }

        // Mixed-type copy (double -> float)
        template <typename... Args>
        static DACE_HDFI void Copy(const double *src, float *dst,
                                   const int &src_stride,
                                   const Args &...src_otherdims) {
            dst[0] = static_cast<float>(src[0]);
        }
    };

    template <int SRC_STRIDE, int... OTHER_SRCDIMS>
    struct ConstSrc {
        // Same-type copy
        template <typename... Args>
        static DACE_HDFI void Copy(const float *src, float *dst,
                                   const int &dst_stride,
                                   const Args &...dst_otherdims) {
            dst[0] = src[0];
        }

        // Mixed-type copy (float -> double)
        template <typename... Args>
        static DACE_HDFI void Copy(const float *src, double *dst,
                                   const int &dst_stride,
                                   const Args &...dst_otherdims) {
            dst[0] = static_cast<double>(src[0]);
        }
    };
};

template <>
struct CopyND<double, 1, false, 1> {
    template <int DST_STRIDE, int... OTHER_DSTDIMS>
    struct ConstDst {
        // Same-type copy
        template <typename... Args>
        static DACE_HDFI void Copy(const double *src, double *dst,
                                   const int &src_stride,
                                   const Args &...src_otherdims) {
            dst[0] = src[0];
        }

        // Mixed-type copy (float -> double)
        template <typename... Args>
        static DACE_HDFI void Copy(const float *src, double *dst,
                                   const int &src_stride,
                                   const Args &...src_otherdims) {
            dst[0] = static_cast<double>(src[0]);
        }
    };

    template <int SRC_STRIDE, int... OTHER_SRCDIMS>
    struct ConstSrc {
        // Same-type copy
        template <typename... Args>
        static DACE_HDFI void Copy(const double *src, double *dst,
                                   const int &dst_stride,
                                   const Args &...dst_otherdims) {
            dst[0] = src[0];
        }

        // Mixed-type copy (double -> float)
        template <typename... Args>
        static DACE_HDFI void Copy(const double *src, float *dst,
                                   const int &dst_stride,
                                   const Args &...dst_otherdims) {
            dst[0] = static_cast<float>(src[0]);
        }
    };
};

} // namespace dace
