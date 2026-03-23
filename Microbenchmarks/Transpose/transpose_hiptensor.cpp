// transpose_hiptensor.cpp — hipTensor 2.x permutation: row-major + blocked
// Compile: hipcc -O3 -std=c++17 -o transpose_hiptensor transpose_hiptensor.cpp -lhiptensor
// variant=0: A[N,N] {r,c} -> B[N,N] {c,r}            (row-major)
// variant=1: A[NB,NB,SB,SB] {a,b,r,c} -> {b,a,c,r}   (blocked layout)
#include <cstdio>
#include <cstdlib>
#include <hip/hip_runtime.h>
#include <hiptensor/hiptensor.h>

#define HC(x) do{hipError_t e=(x);if(e){fprintf(stderr,"HIP %d: %s\n",__LINE__,hipGetErrorString(e));exit(1);}}while(0)
#define HT(x) do{hiptensorStatus_t s=(x);if(s!=HIPTENSOR_STATUS_SUCCESS){fprintf(stderr,"HT %d: %d\n",__LINE__,(int)s);exit(1);}}while(0)

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

    // Blocked: reorder host data
    double*hB=nullptr;
    if(VAR==1){
        if(N%SB){fprintf(stderr,"N%%SB!=0\n");return 1;}
        hB=(double*)malloc(bytes);
        int NB=N/SB;
        for(int r=0;r<N;r++) for(int c=0;c<N;c++)
            hB[(r/SB*NB+c/SB)*SB*SB+(r%SB)*SB+c%SB]=hA[r*N+c];
    }

    double*dA,*dB;
    HC(hipMalloc(&dA,bytes)); HC(hipMalloc(&dB,bytes));
    HC(hipMemcpy(dA,VAR==1?hB:hA,bytes,hipMemcpyHostToDevice));

    hiptensorHandle_t* handle; HT(hiptensorCreate(&handle));
    hiptensorTensorDescriptor_t *dscA, *dscB;
    hiptensorOperationDescriptor_t* opDesc;

    const char* vname;
    if(VAR==0){
        vname="hiptensor";
        int64_t ext[]={N,N}, sA[]={N,1}, sB[]={1,N};
        int32_t mA[]={'r','c'}, mB[]={'c','r'};
        HT(hiptensorCreateTensorDescriptor(handle,&dscA,2,ext,sA,HIP_R_64F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorCreateTensorDescriptor(handle,&dscB,2,ext,sB,HIP_R_64F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorCreatePermutation(handle,&opDesc,dscA,mA,HIPTENSOR_OP_IDENTITY,dscB,mB,HIPTENSOR_COMPUTE_64F));
    } else {
        vname="hiptensor_blk";
        int NB=N/SB;
        int64_t ext4[]={NB,NB,SB,SB};
        int64_t sA4[]={(int64_t)NB*SB*SB,(int64_t)SB*SB,(int64_t)SB,1};
        int64_t sB4[]={(int64_t)SB*SB,(int64_t)NB*SB*SB,1,(int64_t)SB};
        int32_t mA4[]={'a','b','r','c'}, mB4[]={'b','a','c','r'};
        HT(hiptensorCreateTensorDescriptor(handle,&dscA,4,ext4,sA4,HIP_R_64F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorCreateTensorDescriptor(handle,&dscB,4,ext4,sB4,HIP_R_64F,HIPTENSOR_OP_IDENTITY));
        HT(hiptensorCreatePermutation(handle,&opDesc,dscA,mA4,HIPTENSOR_OP_IDENTITY,dscB,mB4,HIPTENSOR_COMPUTE_64F));
    }

    hiptensorPlanPreference_t* pref;
    HT(hiptensorCreatePlanPreference(handle,&pref,HIPTENSOR_ALGO_DEFAULT,HIPTENSOR_JIT_MODE_NONE));
    hiptensorPlan_t* plan;
    HT(hiptensorCreatePlan(handle,&plan,opDesc,pref,0));

    double alpha=1.0;
    for(int i=0;i<WU;i++) HT(hiptensorPermute(handle,plan,&alpha,dA,dB,0));
    HC(hipDeviceSynchronize());

    hipEvent_t*ev=(hipEvent_t*)malloc((REPS+1)*sizeof(hipEvent_t));
    for(int i=0;i<=REPS;i++) HC(hipEventCreate(&ev[i]));
    for(int i=0;i<REPS;i++){
        HC(hipEventRecord(ev[i]));
        HT(hiptensorPermute(handle,plan,&alpha,dA,dB,0));
    }
    HC(hipEventRecord(ev[REPS])); HC(hipEventSynchronize(ev[REPS]));

    float*ims=(float*)malloc(REPS*sizeof(float)); float tot=0;
    for(int i=0;i<REPS;i++){HC(hipEventElapsedTime(&ims[i],ev[i],ev[i+1]));tot+=ims[i];}

    HC(hipMemcpy(hA,dB,bytes,hipMemcpyDeviceToHost));
    double cksum=0; for(size_t i=0;i<elems;i++) cksum+=hA[i];
    double bpi=2.0*elems*8, avg_ms=tot/REPS;
    printf("%s N=%d SB=%d | %.4f ms  %.1f GB/s  cksum=%.6e\n",
           vname,N,SB,avg_ms,bpi/(avg_ms/1000.0)/1e9,cksum);

    FILE*f=fopen(csv,"a");
    if(f){for(int i=0;i<REPS;i++){
        double gbs=bpi/(ims[i]/1000.0)/1e9;
        fprintf(f,"%s,%d,0,0,0,0,%d,0,%d,%.6f,%.3f,%.6e\n",vname,N,SB,i,ims[i]/1000.0,gbs,cksum);
    }fclose(f);}

    hiptensorDestroyPlan(plan); hiptensorDestroyPlanPreference(pref);
    hiptensorDestroyOperationDescriptor(opDesc);
    hiptensorDestroyTensorDescriptor(dscA); hiptensorDestroyTensorDescriptor(dscB);
    hiptensorDestroy(handle);
    for(int i=0;i<=REPS;i++) HC(hipEventDestroy(ev[i]));
    free(ev);free(ims);free(hA);if(hB)free(hB);
    HC(hipFree(dA));HC(hipFree(dB));
    return 0;
}