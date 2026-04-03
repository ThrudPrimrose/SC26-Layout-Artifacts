// transpose_cutensor.cpp — cuTENSOR 2.x permutation: row-major + blocked (FP32)
// Compile: nvcc -O3 -std=c++17 -o transpose_cutensor transpose_cutensor.cpp -lcutensor
// variant=0: A[N,N] {r,c} -> B[N,N] {c,r}            (row-major)
// variant=1: A[NB,NB,SB,SB] {a,b,r,c} -> {b,a,c,r}   (blocked layout)
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cuda_runtime.h>
#include <cutensor.h>

#define CC(x) do{cudaError_t e=(x);if(e){fprintf(stderr,"CUDA %d: %s\n",__LINE__,cudaGetErrorString(e));exit(1);}}while(0)
#define CT(x) do{cutensorStatus_t s=(x);if(s){fprintf(stderr,"CT %d: %s\n",__LINE__,cutensorGetErrorString(s));exit(1);}}while(0)

// ---------- verification helpers ----------

static inline float orig_val(int row, int col, int N) {
    return (float)((size_t)row * N + col) / N;
}

static int verify_var0(const float* hB, int N, int num_checks) {
    int errs = 0;
    for (int k = 0; k < num_checks; k++) {
        int i = (int)((long)k * 997 % N);
        int j = (int)((long)k * 1013 % N);
        float expected = orig_val(i, j, N);
        float got      = hB[(size_t)j * N + i];
        if (fabsf(got - expected) > 1e-4f) {
            if (errs < 5)
                fprintf(stderr, "  MISMATCH var0: B[%d*%d+%d]=%.6f  expected=%.6f\n",
                        j, N, i, got, expected);
            errs++;
        }
    }
    return errs;
}

static int verify_var1(const float* hB, int N, int SB, int num_checks) {
    int NB = N / SB;
    int errs = 0;
    for (int k = 0; k < num_checks; k++) {
        int a = (int)((long)k * 997 % NB);
        int b = (int)((long)k * 1013 % NB);
        int r = (int)((long)k * 1021 % SB);
        int c = (int)((long)k * 1031 % SB);
        float expected = orig_val(a * SB + r, b * SB + c, N);
        size_t idx = (size_t)b * NB * SB * SB + a * SB * SB + c * SB + r;
        float got = hB[idx];
        if (fabsf(got - expected) > 1e-4f) {
            if (errs < 5)
                fprintf(stderr, "  MISMATCH var1: B[b=%d,a=%d,c=%d,r=%d]=%.6f  expected=%.6f\n",
                        b, a, c, r, got, expected);
            errs++;
        }
    }
    return errs;
}

// ---------- main ----------

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,"Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=5] [REPS=100]\n",argv[0]);
        return 1;
    }
    int N=atoi(argv[1]), VAR=atoi(argv[2]);
    const char* csv=argv[3];
    int SB=(argc>4)?atoi(argv[4]):32, WU=(argc>5)?atoi(argv[5]):5, REPS=(argc>6)?atoi(argv[6]):100;

    size_t elems=(size_t)N*N, bytes=elems*sizeof(float);

    /* Print cuTENSOR version for diagnostics */
    size_t ct_ver = cutensorGetVersion();
    fprintf(stderr, "  cuTENSOR version: %zu  (major=%zu minor=%zu patch=%zu)\n",
            ct_ver, ct_ver/10000, (ct_ver/100)%100, ct_ver%100);

    float*hA=(float*)malloc(bytes);
    for(size_t i=0;i<elems;i++) hA[i]=(float)i/N;

    float*hB=nullptr;
    if(VAR==1){
        if(N%SB){fprintf(stderr,"N%%SB!=0\n");return 1;}
        hB=(float*)malloc(bytes);
        int NB=N/SB;
        for(int r=0;r<N;r++) for(int c=0;c<N;c++)
            hB[(r/SB*NB+c/SB)*SB*SB+(r%SB)*SB+c%SB]=hA[r*N+c];
    }

    float*dA,*dB;
    CC(cudaMalloc(&dA,bytes)); CC(cudaMalloc(&dB,bytes));
    CC(cudaMemcpy(dA,VAR==1?hB:hA,bytes,cudaMemcpyHostToDevice));

    cutensorHandle_t handle;          CT(cutensorCreate(&handle));
    cutensorTensorDescriptor_t dscA, dscB;
    cutensorOperationDescriptor_t opDesc;
    cutensorPlanPreference_t pref;
    cutensorPlan_t plan;

    /* ═══════════════════════════════════════════════════════════════════
     * cuTENSOR 2.x convention:
     *   - NULL strides → column-major (Fortran) order, i.e. stride[0]=1
     *   - To describe C-order (row-major) data, reverse the dimension
     *     order and mode labels.  This matches NVIDIA's own samples.
     *
     * Alignment: cudaMalloc guarantees ≥256-byte; 128 is safe.
     * ═══════════════════════════════════════════════════════════════════ */
    const uint32_t ALIGN = 128;

    const char* vname;
    if(VAR==0){
        vname="cutensor";

        /* C-order: A[row][col] at offset row*N + col
         * Column-major {N,N} NULL strides → {1, N}
         * dim0='c' dim1='r' → element(c,r) at c + r*N = row*N+col  ✓
         *
         * Transposed B[col][row] at offset col*N + row
         * dim0='r' dim1='c' → element(r,c) at r + c*N = col*N+row  ✓ */
        int64_t ext[] = {(int64_t)N, (int64_t)N};
        int32_t mA[]  = {'c', 'r'};
        int32_t mB[]  = {'r', 'c'};

        CT(cutensorCreateTensorDescriptor(handle, &dscA, 2, ext, NULL,
                                          CUTENSOR_R_32F, ALIGN));
        CT(cutensorCreateTensorDescriptor(handle, &dscB, 2, ext, NULL,
                                          CUTENSOR_R_32F, ALIGN));
        CT(cutensorCreatePermutation(handle, &opDesc,
                                     dscA, mA, CUTENSOR_OP_IDENTITY,
                                     dscB, mB,
                                     CUTENSOR_COMPUTE_DESC_32F));
    } else {
        vname="cutensor_blk";
        int NB=N/SB;

        /* C-order: A[a][b][r][c] strides {NB·SB², SB², SB, 1}
         * Column-major reversed: extents {SB,SB,NB,NB}  strides {1,SB,SB²,NB·SB²}
         * modes {'c','r','b','a'}  → offset c+r·SB+b·SB²+a·NB·SB²  ✓
         *
         * Transposed B[b][a][c][r] strides {NB·SB², SB², SB, 1}
         * Column-major reversed: same extents, modes {'r','c','a','b'}
         * → offset r+c·SB+a·SB²+b·NB·SB²  ✓ */
        int64_t ext4[] = {(int64_t)SB, (int64_t)SB, (int64_t)NB, (int64_t)NB};
        int32_t mA4[]  = {'c', 'r', 'b', 'a'};
        int32_t mB4[]  = {'r', 'c', 'a', 'b'};

        CT(cutensorCreateTensorDescriptor(handle, &dscA, 4, ext4, NULL,
                                          CUTENSOR_R_32F, ALIGN));
        CT(cutensorCreateTensorDescriptor(handle, &dscB, 4, ext4, NULL,
                                          CUTENSOR_R_32F, ALIGN));
        CT(cutensorCreatePermutation(handle, &opDesc,
                                     dscA, mA4, CUTENSOR_OP_IDENTITY,
                                     dscB, mB4,
                                     CUTENSOR_COMPUTE_DESC_32F));
    }

    CT(cutensorCreatePlanPreference(handle,&pref,CUTENSOR_ALGO_DEFAULT,CUTENSOR_JIT_MODE_NONE));
    CT(cutensorCreatePlan(handle,&plan,opDesc,pref,0));

    float alpha=1.0f;
    for(int i=0;i<WU;i++) CT(cutensorPermute(handle,plan,&alpha,dA,dB,0));
    CC(cudaDeviceSynchronize());

    cudaEvent_t*ev=(cudaEvent_t*)malloc((REPS+1)*sizeof(cudaEvent_t));
    for(int i=0;i<=REPS;i++) CC(cudaEventCreate(&ev[i]));
    for(int i=0;i<REPS;i++){
        CC(cudaEventRecord(ev[i]));
        CT(cutensorPermute(handle,plan,&alpha,dA,dB,0));
    }
    CC(cudaEventRecord(ev[REPS])); CC(cudaEventSynchronize(ev[REPS]));

    float*ims=(float*)malloc(REPS*sizeof(float)); float tot=0;
    for(int i=0;i<REPS;i++){CC(cudaEventElapsedTime(&ims[i],ev[i],ev[i+1]));tot+=ims[i];}

    float* hOut=(float*)malloc(bytes);
    CC(cudaMemcpy(hOut,dB,bytes,cudaMemcpyDeviceToHost));
    double cksum=0; for(size_t i=0;i<elems;i++) cksum+=hOut[i];

    int ncheck = (N < 1000) ? N*N : 10000;
    int errs;
    if (VAR == 0)
        errs = verify_var0(hOut, N, ncheck);
    else
        errs = verify_var1(hOut, N, SB, ncheck);

    if (errs)
        fprintf(stderr, "VERIFY FAIL: %d / %d mismatches\n", errs, ncheck);
    else
        fprintf(stderr, "VERIFY OK: %d element checks passed\n", ncheck);

    double bpi = 2.0 * elems * sizeof(float);
    float avg_ms = tot / REPS;
    float gbps = (float)(bpi / (avg_ms / 1000.0) / 1e9);

    printf("%s N=%d SB=%d | %.4f ms  %.1f GB/s  cksum=%.6e  verify=%s\n",
           vname, N, SB, avg_ms, gbps, cksum, errs ? "FAIL" : "OK");

    FILE*f=fopen(csv,"a");
    if(f){for(int i=0;i<REPS;i++){
        float gbs=(float)(bpi/(ims[i]/1000.0)/1e9);
        fprintf(f,"%s,%d,0,0,0,0,%d,0,%d,%.6f,%.3f,%.6e\n",vname,N,SB,i,ims[i]/1000.0,gbs,cksum);
    }fclose(f);}

    cutensorDestroyPlan(plan);
    cutensorDestroyPlanPreference(pref);
    cutensorDestroyOperationDescriptor(opDesc);
    cutensorDestroyTensorDescriptor(dscA);
    cutensorDestroyTensorDescriptor(dscB);
    cutensorDestroy(handle);
    for(int i=0;i<=REPS;i++) CC(cudaEventDestroy(ev[i]));
    free(ev);free(ims);free(hA);free(hOut);if(hB)free(hB);
    CC(cudaFree(dA));CC(cudaFree(dB));
    return 0;
}