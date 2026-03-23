// transpose_hiptensor.cpp — hipTensor 1.x permutation: row-major + blocked
// Compile: hipcc -O3 -std=c++17 -o transpose_hiptensor transpose_hiptensor.cpp -lhiptensor
// variant=0: A[N,N] {r,c} -> B[N,N] {c,r}            (row-major)
// variant=1: A[NB,NB,SB,SB] {a,b,r,c} -> {b,a,c,r}   (blocked layout)
#include <cstdio>
#include <cstdlib>
#include <hip/hip_runtime.h>
#include <hiptensor/hiptensor.hpp>
#include <hiptensor/hiptensor_types.hpp>

#define HC(x) do{hipError_t e=(x);if(e){fprintf(stderr,"HIP %d: %s\n",__LINE__,hipGetErrorString(e));exit(1);}}while(0)
#define HT(x) do{hiptensorStatus_t s=(x);if(s!=HIPTENSOR_STATUS_SUCCESS){fprintf(stderr,"HT %d: %d\n",__LINE__,(int)s);exit(1);}}while(0)

void verify_v1(float* t, int N) {
    // Check row-major transpose: output[t][i][j] should equal input[j][i]
    // Input was initialized as hA[i*N+j] = i*N + j
    bool ok = true;
    // Check a few elements
    int test_cases[][2] = {{0,1},{1,0},{2,1},{1,2},{N/2, N/2+1}};
    for (auto& rc : test_cases) {
        int r = rc[0], c = rc[1];
        float expected = (float)(c * N + r);
        float actual = t[r * N + c];
        if (actual != expected) {
            printf("verify_v1 error: at (%d,%d) expected %.0f got %.0f\n",
                   r, c, expected, actual);
            ok = false;
        }
    }
    if (ok)
        printf("verify_v1 passed\n");
    else
        exit(1);
}

void verify_v2(float* t, int N, int SB) {
    int NB = N / SB;
    bool ok = true;

    auto idx = [NB, SB](int b, int a, int c, int r) -> size_t {
        return b * NB * SB * SB + a * SB * SB + c * SB + r;
    };

    auto expected = [NB, SB](int b, int a, int c, int r) -> float {
        return a * NB * SB * SB + b * SB * SB + r * SB + c;
    };

    struct Test { int b, a, c, r; };
    Test tests[] = {
        {0, 0, 0, 0},                         // top‑left corner
        {0, 1, 0, 1},                         // (b=0,a=1,c=0,r=1)
        {1, 1, 1, 1},                         // centre of block (1,1)
        {2, 1, 2, 1},                         // (b=2,a=1,c=2,r=1)
        {NB-1, NB-1, SB-1, SB-1}              // bottom‑right corner
    };

    for (auto& tst : tests) {
        float exp_val = expected(tst.b, tst.a, tst.c, tst.r);
        float act_val = t[idx(tst.b, tst.a, tst.c, tst.r)];
        if (act_val != exp_val) {
            printf("verify_v2 error: at (b=%d,a=%d,c=%d,r=%d) expected %.0f got %.0f\n",
                   tst.b, tst.a, tst.c, tst.r, exp_val, act_val);
            ok = false;
        }
    }

    if (ok)
        printf("verify_v2 passed\n");
    else
        exit(1);
}

// Print the output as a normal N×N row‑major matrix
void print_output_as_matrix(float* t, int N, int SB) {
    int NB = N / SB;
    printf("Output matrix (row‑major view):\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int a = i / SB;          // block row index
            int r = i % SB;          // row inside block
            int b = j / SB;          // block column index
            int c = j % SB;          // column inside block
            size_t idx = b * NB * SB * SB + a * SB * SB + c * SB + r;
            printf("%6.0f ", t[idx]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,"Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=5] [REPS=100]\n",argv[0]);
        return 1;
    }
    int N=atoi(argv[1]), VAR=atoi(argv[2]);
    const char* csv=argv[3];
    int SB=(argc>4)?atoi(argv[4]):32, WU=(argc>5)?atoi(argv[5]):5, REPS=(argc>6)?atoi(argv[6]):100;

    size_t elems=(size_t)N*N, bytes=elems*sizeof(float);
    float*hA=(float*)malloc(bytes);
    for(int i=0;i<N;i++) {
        for (int j=0; j<N; j++){
            hA[i*N + j] = (float)(i*N + j);
        }
    }
    //for (int i = 0; i < N*N; i++){
    //    printf("%8.0f ", hA[i]);
    //}
    float*dA,*dB;
    HC(hipMalloc(&dA,bytes)); HC(hipMalloc(&dB,bytes));
    HC(hipMemcpy(dA,hA,bytes,hipMemcpyHostToDevice));

    // hipTensor 1.x: handle is pointer-to-pointer
    hiptensorHandle_t* handle;
    HT(hiptensorCreate(&handle));

    // hipTensor 1.x: descriptors are stack-allocated structs
    hiptensorTensorDescriptor_t dscA, dscB;

    const char* vname;
    int32_t mA2[]={'r','c'}, mB2[]={'c','r'};
    int32_t mA4[]={'a','b','r','c'}, mB4[]={'b','a','c','r'};
    int32_t *mA, *mB;

    if(VAR==0){
        vname="hiptensor";
        int64_t ext[]={N,N}, sA[]={N,1}, sB[]={N,1};
        mA=mA2; mB=mB2;
        // hipTensor 1.x: hiptensorInitTensorDescriptor takes (handle, &desc, ...)
        // desc is passed by pointer to stack struct, not pointer-to-pointer
        HT(hiptensorInitTensorDescriptor(handle,&dscA,2,ext,sA,HIP_R_32F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorInitTensorDescriptor(handle,&dscB,2,ext,sB,HIP_R_32F,HIPTENSOR_OP_IDENTITY));
    } else {
        vname="hiptensor_blk";
        int NB=N/SB;
        int64_t ext4[]={NB,NB,SB,SB};
        int64_t sA4[]={(int64_t)NB*SB*SB,(int64_t)SB*SB,(int64_t)SB,1};
        int64_t sB4[]={(int64_t)NB*SB*SB,(int64_t)SB*SB,(int64_t)SB,1};
        mA=mA4; mB=mB4;
        HT(hiptensorInitTensorDescriptor(handle,&dscA,4,ext4,sA4,HIP_R_32F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorInitTensorDescriptor(handle,&dscB,4,ext4,sB4,HIP_R_32F,HIPTENSOR_OP_IDENTITY));
    }

    // hipTensor 1.x: hiptensorPermutation is a direct call, no plan needed
    // hiptensorPermutation(handle, alpha, A, &descA, modeA, B, &descB, modeB, typeCompute, stream)
    float alpha=1.0;
    for(int i=0;i<WU;i++)
        HT(hiptensorPermutation(handle,&alpha,dA,&dscA,mA,dB,&dscB,mB,HIP_R_32F,0));
    HC(hipDeviceSynchronize());

    hipEvent_t*ev=(hipEvent_t*)malloc((REPS+1)*sizeof(hipEvent_t));
    for(int i=0;i<=REPS;i++) HC(hipEventCreate(&ev[i]));
    for(int i=0;i<REPS;i++){
        HC(hipEventRecord(ev[i]));
        HT(hiptensorPermutation(handle,&alpha,dA,&dscA,mA,dB,&dscB,mB,HIP_R_32F,0));
    }
    HC(hipEventRecord(ev[REPS])); HC(hipEventSynchronize(ev[REPS]));

    float*ims=(float*)malloc(REPS*sizeof(float)); float tot=0;
    for(int i=0;i<REPS;i++){HC(hipEventElapsedTime(&ims[i],ev[i],ev[i+1]));tot+=ims[i];}

    HC(hipMemcpy(hA,dB,bytes,hipMemcpyDeviceToHost));
    float cksum=0; for(size_t i=0;i<elems;i++) cksum+=hA[i];
    float bpi=2.0*elems*4, avg_ms=tot/REPS;
    printf("%s N=%d SB=%d | %.4f ms  %.1f GB/s  cksum=%.6e\n",
           vname,N,SB,avg_ms,bpi/(avg_ms/1000.0)/1e9,cksum);

    //for (int i = 0; i < N*N; i++){
    //    printf("%8.0f ", hA[i]);
    //}
    //printf("\n");

    //if (VAR==0){
    //    verify_v1(hA, N);
    //} else {
    //    verify_v2(hA, N, SB);
    //}

    FILE*f=fopen(csv,"a");
    if(f){for(int i=0;i<REPS;i++){
        float gbs=bpi/(ims[i]/1000.0)/1e9;
        fprintf(f,"%s,%d,0,0,0,0,%d,0,%d,%.6f,%.3f,%.6e\n",vname,N,SB,i,ims[i]/1000.0,gbs,cksum);
    }fclose(f);}

    // hipTensor 1.x: only need to destroy the handle (descriptors are stack-allocated)
    hiptensorDestroy(handle);
    for(int i=0;i<=REPS;i++) HC(hipEventDestroy(ev[i]));
    free(ev);free(ims);free(hA);
    HC(hipFree(dA));HC(hipFree(dB));
    return 0;
}