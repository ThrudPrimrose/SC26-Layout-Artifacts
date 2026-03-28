/*
 * collapse_demo.c -- 20x20 array, collapse(2), each thread writes its ID
 *
 * Compile: gcc -O2 -fopenmp -o collapse_demo collapse_demo.c
 * Run:     OMP_NUM_THREADS=4 ./collapse_demo
 */
#include <stdio.h>
#include <omp.h>

#define N 20

int main() {
    int A[N][N];

    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < N; j++)
        for (int i = 0; i < N; i++)
            A[i][j] = omp_get_thread_num();

    printf("    ");
    for (int j = 0; j < N; j++) printf("%3d", j);
    printf("\n    ");
    for (int j = 0; j < N; j++) printf("---");
    printf("\n");

    for (int i = 0; i < N; i++) {
        printf("%2d |", i);
        for (int j = 0; j < N; j++)
            printf("%3d", A[i][j]);
        printf("\n");
    }
    return 0;
}