// transpose_hiptensor.cpp — hipTensor 2.x permutation benchmark
// Compile: hipcc -O3 -std=c++17 -o transpose_hiptensor transpose_hiptensor.cpp -lhiptensor
#include <cstdio>
#include <cstdlib>
#include <hip/hip_runtime.h>
#include <hiptensor/hiptensor.h>

#define HIP_CHECK(x) do{hipError_t e=(x);if(e){fprintf(stderr,"HIP %d: %s\n",__LINE__,hipGetErrorString(e));exit(1);}}while(0)
#define HT_CHECK(x) do{hiptensorStatus_t s=(x);if(s!=HIPTENSOR_STATUS_SUCCESS){fprintf(stderr,"hipTensor %d: status %d\n",__LINE__,(int)s);exit(1);}}while(0)

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr,"Usage: %s <N> <csv> [WARMUP=5] [REPS=100]\n",argv[0]); return 1; }
    int N = atoi(argv[1]);
    const char* csv = argv[2];
    int WARMUP = (argc>3) ? atoi(argv[3]) : 5;
    int REPS   = (argc>4) ? atoi(argv[4]) : 100;

    size_t bytes = (size_t)N * N * sizeof(double);
    double *hA = (double*)malloc(bytes);
    for (size_t i = 0; i < (size_t)N*N; i++) hA[i] = (double)i / N;

    double *dA, *dB;
    HIP_CHECK(hipMalloc(&dA, bytes));
    HIP_CHECK(hipMalloc(&dB, bytes));
    HIP_CHECK(hipMemcpy(dA, hA, bytes, hipMemcpyHostToDevice));

    // hipTensor setup
    hiptensorHandle_t* handle;
    HT_CHECK(hiptensorCreate(&handle));

    int64_t extA[] = {N, N};
    int64_t strA[] = {N, 1};   // row-major
    int64_t extB[] = {N, N};
    int64_t strB[] = {1, N};   // transposed

    int32_t modeA[] = {'r', 'c'};
    int32_t modeB[] = {'c', 'r'};

    hiptensorTensorDescriptor_t* descA;
    hiptensorTensorDescriptor_t* descB;
    HT_CHECK(hiptensorCreateTensorDescriptor(handle, &descA, 2, extA, strA,
             HIP_R_64F, HIPTENSOR_OP_IDENTITY));
    HT_CHECK(hiptensorCreateTensorDescriptor(handle, &descB, 2, extB, strB,
             HIP_R_64F, HIPTENSOR_OP_IDENTITY));

    hiptensorOperationDescriptor_t* opDesc;
    HT_CHECK(hiptensorCreatePermutation(handle, &opDesc,
             descA, modeA, HIPTENSOR_OP_IDENTITY,
             descB, modeB,
             HIPTENSOR_COMPUTE_64F));

    hiptensorPlanPreference_t* planPref;
    HT_CHECK(hiptensorCreatePlanPreference(handle, &planPref,
             HIPTENSOR_ALGO_DEFAULT, HIPTENSOR_JIT_MODE_NONE));

    hiptensorPlan_t* plan;
    HT_CHECK(hiptensorCreatePlan(handle, &plan, opDesc, planPref, 0));

    double alpha = 1.0;

    // Warmup
    for (int i = 0; i < WARMUP; i++)
        HT_CHECK(hiptensorPermute(handle, plan, &alpha, dA, dB, 0));
    HIP_CHECK(hipDeviceSynchronize());

    // Timed runs
    hipEvent_t* ev = (hipEvent_t*)malloc((REPS+1)*sizeof(hipEvent_t));
    for (int i = 0; i <= REPS; i++) HIP_CHECK(hipEventCreate(&ev[i]));
    for (int i = 0; i < REPS; i++) {
        HIP_CHECK(hipEventRecord(ev[i]));
        HT_CHECK(hiptensorPermute(handle, plan, &alpha, dA, dB, 0));
    }
    HIP_CHECK(hipEventRecord(ev[REPS]));
    HIP_CHECK(hipEventSynchronize(ev[REPS]));

    float* ims = (float*)malloc(REPS * sizeof(float));
    float tot = 0;
    for (int i = 0; i < REPS; i++) {
        HIP_CHECK(hipEventElapsedTime(&ims[i], ev[i], ev[i+1]));
        tot += ims[i];
    }

    // Checksum
    HIP_CHECK(hipMemcpy(hA, dB, bytes, hipMemcpyDeviceToHost));
    double cksum = 0;
    for (size_t i = 0; i < (size_t)N*N; i++) cksum += hA[i];

    double bpi = 2.0 * N * N * 8;
    double avg_ms = tot / REPS;
    printf("hiptensor N=%d | %.4f ms  %.1f GB/s  cksum=%.6e\n",
           N, avg_ms, bpi/(avg_ms/1000.0)/1e9, cksum);

    FILE* f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++) {
            double gbs = bpi / (ims[i]/1000.0) / 1e9;
            fprintf(f, "hiptensor,%d,0,0,0,0,0,0,%d,%.6f,%.3f,%.6e\n",
                    N, i, ims[i]/1000.0, gbs, cksum);
        }
        fclose(f);
    }

    hiptensorDestroyPlan(plan);
    hiptensorDestroyPlanPreference(planPref);
    hiptensorDestroyOperationDescriptor(opDesc);
    hiptensorDestroyTensorDescriptor(descA);
    hiptensorDestroyTensorDescriptor(descB);
    hiptensorDestroy(handle);
    for (int i = 0; i <= REPS; i++) HIP_CHECK(hipEventDestroy(ev[i]));
    free(ev); free(ims); free(hA);
    HIP_CHECK(hipFree(dA)); HIP_CHECK(hipFree(dB));
    return 0;
}
