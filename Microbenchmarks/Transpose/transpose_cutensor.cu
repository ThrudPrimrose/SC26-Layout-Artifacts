// transpose_cutensor.cu — cuTENSOR 2.x permutation: row-major + blocked
// Compile: nvcc -O3 -std=c++17 -o transpose_cutensor transpose_cutensor.cu -lcutensor
// variant=0: A[N,N] {r,c} -> B[N,N] {c,r}            (row-major)
// variant=1: A[NB,NB,SB,SB] {a,b,r,c} -> {b,a,c,r}   (blocked layout)
#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <cutensor.h>

#define CC(x) do{cudaError_t e=(x);if(e){fprintf(stderr,"CUDA %d: %s\n",__LINE__,cudaGetErrorString(e));exit(1);}}while(0)
#define CT(x) do{cutensorStatus_t s=(x);if(s){fprintf(stderr,"CT %d: %s\n",__LINE__,cutensorGetErrorString(s));exit(1);}}while(0)

int main(int argc, char** argv) {
    if (argc < 4) {
        fprintf(stderr,"Usage: %s <N> <variant> <csv> [SB=32] [WARMUP=5] [REPS=100]\n",argv[0]);
        return 1;
    }
    int N=atoi(argv[1]), VAR=atoi(argv[2]);
    const char* csv=argv[3];
    int SB=(argc>4)?atoi(argv[4]):32, WU=(argc>5)?atoi(argv[5]):5, REPS=(argc>6)?atoi(argv[6]):100;

    size_t elems=(size_t)N*N, bytes=elems*sizeof(double);
    double*hA=(double*)malloc(bytes);
    for(size_t i=0;i<elems;i++) hA[i]=(double)i/N;

    double*hB=nullptr;
    if(VAR==1){
        if(N%SB){fprintf(stderr,"N%%SB!=0\n");return 1;}
        hB=(double*)malloc(bytes);
        int NB=N/SB;
        for(int r=0;r<N;r++) for(int c=0;c<N;c++)
            hB[(r/SB*NB+c/SB)*SB*SB+(r%SB)*SB+c%SB]=hA[r*N+c];
    }

    double*dA,*dB;
    CC(cudaMalloc(&dA,bytes)); CC(cudaMalloc(&dB,bytes));
    CC(cudaMemcpy(dA,VAR==1?hB:hA,bytes,cudaMemcpyHostToDevice));

    // cuTENSOR 2.x: all opaque types are value types (pointer typedefs),
    // declare as values and pass by address to Create functions.
    cutensorHandle_t handle;             CT(cutensorCreate(&handle));
    cutensorTensorDescriptor_t dscA, dscB;
    cutensorOperationDescriptor_t opDesc;
    cutensorPlanPreference_t pref;
    cutensorPlan_t plan;

    const char* vname;
    if(VAR==0){
        vname="cutensor";
        int64_t ext[]={N,N}, sA[]={N,1}, sB[]={1,N};
        int32_t mA[]={'r','c'}, mB[]={'c','r'};
        CT(cutensorCreateTensorDescriptor(handle,&dscA,2,ext,sA,CUTENSOR_R_64F,CUTENSOR_OP_IDENTITY));
        CT(cutensorCreateTensorDescriptor(handle,&dscB,2,ext,sB,CUTENSOR_R_64F,CUTENSOR_OP_IDENTITY));
        CT(cutensorCreatePermutation(handle,&opDesc,dscA,mA,CUTENSOR_OP_IDENTITY,dscB,mB,CUTENSOR_COMPUTE_DESC_64F));
    } else {
        vname="cutensor_blk";
        int NB=N/SB;
        int64_t ext4[]={NB,NB,SB,SB};
        int64_t sA4[]={(int64_t)NB*SB*SB,(int64_t)SB*SB,(int64_t)SB,1};
        int64_t sB4[]={(int64_t)SB*SB,(int64_t)NB*SB*SB,1,(int64_t)SB};
        int32_t mA4[]={'a','b','r','c'}, mB4[]={'b','a','c','r'};
        CT(cutensorCreateTensorDescriptor(handle,&dscA,4,ext4,sA4,CUTENSOR_R_64F,CUTENSOR_OP_IDENTITY));
        CT(cutensorCreateTensorDescriptor(handle,&dscB,4,ext4,sB4,CUTENSOR_R_64F,CUTENSOR_OP_IDENTITY));
        CT(cutensorCreatePermutation(handle,&opDesc,dscA,mA4,CUTENSOR_OP_IDENTITY,dscB,mB4,CUTENSOR_COMPUTE_DESC_64F));
    }

    CT(cutensorCreatePlanPreference(handle,&pref,CUTENSOR_ALGO_DEFAULT,CUTENSOR_JIT_MODE_NONE));
    CT(cutensorCreatePlan(handle,&plan,opDesc,pref,0));

    double alpha=1.0;
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

    CC(cudaMemcpy(hA,dB,bytes,cudaMemcpyDeviceToHost));
    double cksum=0; for(size_t i=0;i<elems;i++) cksum+=hA[i];
    double bpi=2.0*elems*8, avg_ms=tot/REPS;
    printf("%s N=%d SB=%d | %.4f ms  %.1f GB/s  cksum=%.6e\n",
           vname,N,SB,avg_ms,bpi/(avg_ms/1000.0)/1e9,cksum);

    FILE*f=fopen(csv,"a");
    if(f){for(int i=0;i<REPS;i++){
        double gbs=bpi/(ims[i]/1000.0)/1e9;
        fprintf(f,"%s,%d,0,0,0,0,%d,0,%d,%.6f,%.3f,%.6e\n",vname,N,SB,i,ims[i]/1000.0,gbs,cksum);
    }fclose(f);}

    cutensorDestroyPlan(plan);
    cutensorDestroyPlanPreference(pref);
    cutensorDestroyOperationDescriptor(opDesc);
    cutensorDestroyTensorDescriptor(dscA);
    cutensorDestroyTensorDescriptor(dscB);
    cutensorDestroy(handle);
    for(int i=0;i<=REPS;i++) CC(cudaEventDestroy(ev[i]));
    free(ev);free(ims);free(hA);if(hB)free(hB);
    CC(cudaFree(dA));CC(cudaFree(dB));
    return 0;
}
