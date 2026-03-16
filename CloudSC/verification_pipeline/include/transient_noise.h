// transient_noise.h — Multiplicative noise injection for transient sensitivity analysis.
//
// Usage: after a write to a target transient array element, call:
//   array[idx] = transient_noise::perturb(array[idx], eps);
//
// This adds relative noise: value * (1 + eps * uniform(-1, 1))
// Thread-safe via thread_local RNG.
#pragma once

#include <random>
#include <cmath>

namespace transient_noise {

inline thread_local std::mt19937 rng{std::random_device{}()};
inline thread_local std::uniform_real_distribution<double> dist{-1.0, 1.0};

template <typename T>
static inline T perturb(T value, double eps) {
    return static_cast<T>(static_cast<double>(value) * (1.0 + eps * dist(rng)));
}

} // namespace transient_noise
