// sensitivity.h — forward-difference sensitivity sweep for CPU kernels.
//
// For each input field, perturbs it by random noise scaled by eps, re-runs
// the kernel, measures the RMS change in output fields vs unperturbed baseline.
// Reports dimensionless amplification: delta_out_rms / (eps * in_rms)
//
// Ported from velocity/include/sensitivity.cuh (CPU-only, no CUDA).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace sensitivity {

struct Field {
  double *ptr;
  int n;
  std::string name;
};

inline void run(const std::function<void()> &restore,
                const std::function<void()> &execute,
                const std::vector<Field> &inputs,
                const std::vector<Field> &outputs, double eps = 1e-5,
                int n_samples = 10) {
  // Baseline run
  restore();
  execute();

  long long out_n = 0;
  std::vector<std::vector<double>> baseline(outputs.size());
  for (size_t i = 0; i < outputs.size(); ++i) {
    baseline[i].resize(outputs[i].n);
    std::memcpy(baseline[i].data(), outputs[i].ptr,
                outputs[i].n * sizeof(double));
    out_n += outputs[i].n;
  }

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  struct Result {
    double amplification;
    double var_out;
    std::string name;
  };
  std::vector<Result> results;
  results.reserve(inputs.size());

  for (const auto &f : inputs) {
    if (!f.ptr || f.n <= 0)
      continue;

    double total_in_var = 0.0;
    double total_out_var = 0.0;

    for (int s = 0; s < n_samples; ++s) {
      // Restore all arrays to clean state
      restore();

      // Perturb this one input field
      double in_ss = 0.0;
      for (int i = 0; i < f.n; i++) {
        double noise = eps * dist(rng);
        double delta = f.ptr[i] * noise;
        f.ptr[i] += delta;
        in_ss += delta * delta;
      }

      execute();

      double out_ss = 0.0;
      for (size_t i = 0; i < outputs.size(); ++i) {
        for (int j = 0; j < outputs[i].n; ++j) {
          double d = outputs[i].ptr[j] - baseline[i][j];
          out_ss += d * d;
        }
      }

      total_in_var += in_ss / f.n;
      total_out_var += out_ss / (out_n > 0 ? out_n : 1);
    }

    double mean_in_var = total_in_var / n_samples;
    double mean_out_var = total_out_var / n_samples;
    double amplification =
        (mean_in_var > 0.0) ? std::sqrt(mean_out_var / mean_in_var) : 0.0;
    results.push_back({amplification, mean_out_var, f.name});
  }

  std::sort(results.begin(), results.end(),
            [](const auto &a, const auto &b) {
              return a.var_out > b.var_out;
            });

  std::cerr << "\n=== SENSITIVITY ANALYSIS (eps=" << std::scientific << eps
            << ", " << std::defaultfloat << n_samples << " samples) ===\n";
  std::cerr << "Outputs:";
  for (const auto &o : outputs)
    std::cerr << " " << o.name;
  std::cerr << "\n";
  std::cerr << std::left << std::setw(24) << "Input" << std::right
            << std::setw(16) << "Var(out)" << std::setw(16) << "RMS(out)"
            << std::setw(16) << "Amplification\n";
  std::cerr << std::string(72, '-') << "\n";
  for (const auto &r : results)
    std::cerr << std::left << std::setw(24) << r.name << std::right
              << std::setw(16) << std::scientific << std::setprecision(4)
              << r.var_out
              << std::setw(16) << std::sqrt(r.var_out)
              << std::setw(16) << r.amplification << "\n";
  std::cerr << "=== END SENSITIVITY ===\n\n";
}

} // namespace sensitivity
