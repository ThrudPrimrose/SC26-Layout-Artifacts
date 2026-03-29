
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <algorithm>
#include <cstring>

int main(int argc, char** argv) {
    int threads = (argc > 1) ? atoi(argv[1]) : 0;
    if (threads > 0) omp_set_num_threads(threads);

    const long long SZ = 64LL * 1024 * 1024;  // 64M doubles = 512 MB
    const size_t bytes = SZ * sizeof(double);
    double* a = (double*)aligned_alloc(64, bytes);
    double* b = (double*)aligned_alloc(64, bytes);

    // Init (parallel, first-touch)
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < SZ; i++) { a[i] = 1.0; b[i] = 0.0; }

    // Warmup
    #pragma omp parallel for schedule(static)
    for (long long i = 0; i < SZ; i++) b[i] = a[i];

    const int REPS = 100;
    double best = 1e30;
    for (int r = 0; r < REPS; r++) {
        double t0 = omp_get_wtime();
        #pragma omp parallel for schedule(static)
        for (long long i = 0; i < SZ; i++) b[i] = a[i];
        double t1 = omp_get_wtime();
        best = std::min(best, t1 - t0);
    }
    // bytes moved: read a + write b = 2 * bytes
    printf("%.1f\n", 2.0 * bytes / best / 1e9);
    free(a); free(b);
    return 0;
}
