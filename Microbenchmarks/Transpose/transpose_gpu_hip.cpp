#include "hip/hip_runtime.h"
#include <cstdio>
#include <cstdlib>
#include <hip/hip_runtime.h>

#define CHECK(x) do { hipError_t e=(x); if(e){fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,hipGetErrorString(e));exit(1);} } while(0)

// ── V0: Naive global memory transpose, each thread does TX×TY elements ──
template<int TX, int TY>
__global__ void tr_naive(const double* __restrict__ in, double* __restrict__ out, int N) {
    int c0 = blockIdx.x * blockDim.x * TX + threadIdx.x * TX;
    int r0 = blockIdx.y * blockDim.y * TY + threadIdx.y * TY;
    #pragma unroll
    for (int dr = 0; dr < TY; dr++) {
        int r = r0 + dr; if (r >= N) continue;
        #pragma unroll
        for (int dc = 0; dc < TX; dc++) {
            int c = c0 + dc; if (c >= N) continue;
            out[c * N + r] = in[r * N + c];
        }
    }
}

// ── V1: Blocked storage transpose (SB×SB blocks), naive element access ──
// Blocked layout: element (r,c) stored at block (r/SB, c/SB), offset (r%SB, c%SB)
// Linear index: (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB  where NB = N/SB
template<int TX, int TY, int SB>
__global__ void tr_blocked(const double* __restrict__ in, double* __restrict__ out, int N) {
    int c0 = blockIdx.x * blockDim.x * TX + threadIdx.x * TX;
    int r0 = blockIdx.y * blockDim.y * TY + threadIdx.y * TY;
    int NB = N / SB;
    #pragma unroll
    for (int dr = 0; dr < TY; dr++) {
        int r = r0 + dr; if (r >= N) continue;
        #pragma unroll
        for (int dc = 0; dc < TX; dc++) {
            int c = c0 + dc; if (c >= N) continue;
            // src in blocked layout
            int si = (r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB;
            // dst in blocked layout (transposed: row=c, col=r)
            int di = (c/SB * NB + r/SB) * SB*SB + (c%SB)*SB + r%SB;
            out[di] = in[si];
        }
    }
}

// ── V2: Shared memory transpose from row-major layout ──
template<int TX, int TY>
__global__ void tr_smem(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    extern __shared__ double sm[];
    int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, tot = blockDim.x * blockDim.y;
    // Load tile row-major into smem
    for (int k = tid; k < BH * BW; k += tot) {
        int lr = k / BW, lc = k % BW;
        int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.0;
    }
    __syncthreads();
    // Write transposed: smem[lc][lr] → out[gc_out * N + gr_out]
    for (int k = tid; k < BW * BH; k += tot) {
        int lc = k / BH, lr = k % BH;  // iterate column-major for coalesced writes
        int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + lc];
    }
}

// ── V3: Shared memory transpose from blocked layout ──
template<int TX, int TY, int SB>
__global__ void tr_smem_blk(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    extern __shared__ double sm[];
    int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, tot = blockDim.x * blockDim.y;
    int NB = N / SB;
    for (int k = tid; k < BH * BW; k += tot) {
        int lr = k / BW, lc = k % BW;
        int gr = r0 + lr, gc = c0 + lc;
        int si = (gr/SB * NB + gc/SB) * SB*SB + (gr%SB)*SB + gc%SB;
        sm[lr * BW + lc] = (gr < N && gc < N) ? in[si] : 0.0;
    }
    __syncthreads();
    for (int k = tid; k < BW * BH; k += tot) {
        int lc = k / BH, lr = k % BH;
        int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) {
            int di = (gc/SB * NB + gr/SB) * SB*SB + (gc%SB)*SB + gr%SB;
            out[di] = sm[lr * BW + lc];
        }
    }
}

// ── V4: Shared memory + padding (row-major) ──
template<int TX, int TY, int PAD>
__global__ void tr_smem_pad(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    const int SW = BW + PAD;
    extern __shared__ double sm[];
    int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, tot = blockDim.x * blockDim.y;
    for (int k = tid; k < BH * BW; k += tot) {
        int lr = k / BW, lc = k % BW;
        int gr = r0 + lr, gc = c0 + lc;
        sm[lr * SW + lc] = (gr < N && gc < N) ? in[gr * N + gc] : 0.0;
    }
    __syncthreads();
    for (int k = tid; k < BW * BH; k += tot) {
        int lc = k / BH, lr = k % BH;
        int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * SW + lc];
    }
}

// ── V5: Shared memory + XOR swizzle (row-major) ──
template<int TX, int TY>
__global__ void tr_smem_swiz(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    extern __shared__ double sm[];
    int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, tot = blockDim.x * blockDim.y;
    for (int k = tid; k < BH * BW; k += tot) {
        int lr = k / BW, lc = k % BW;
        int gr = r0 + lr, gc = c0 + lc;
        sm[lr * BW + (lc ^ lr)] = (gr < N && gc < N) ? in[gr * N + gc] : 0.0;
    }
    __syncthreads();
    for (int k = tid; k < BW * BH; k += tot) {
        int lc = k / BH, lr = k % BH;
        int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) out[gc * N + gr] = sm[lr * BW + (lc ^ lr)];
    }
}

// ── V6: Blocked layout + shared memory + XOR swizzle ──
template<int TX, int TY, int SB>
__global__ void tr_blk_swiz(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    extern __shared__ double sm[];
    int c0 = blockIdx.x * BW, r0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, tot = blockDim.x * blockDim.y;
    int NB = N / SB;
    for (int k = tid; k < BH * BW; k += tot) {
        int lr = k / BW, lc = k % BW;
        int gr = r0 + lr, gc = c0 + lc;
        int si = (gr/SB * NB + gc/SB) * SB*SB + (gr%SB)*SB + gc%SB;
        sm[lr * BW + (lc ^ lr)] = (gr < N && gc < N) ? in[si] : 0.0;
    }
    __syncthreads();
    for (int k = tid; k < BW * BH; k += tot) {
        int lc = k / BH, lr = k % BH;
        int gr = r0 + lr, gc = c0 + lc;
        if (gr < N && gc < N) {
            int di = (gc/SB * NB + gr/SB) * SB*SB + (gc%SB)*SB + gr%SB;
            out[di] = sm[lr * BW + (lc ^ lr)];
        }
    }
}

// ── Dispatch ──
static const char* V_NAMES[] = {"naive","blocked","smem","smem_blk","smem_pad","smem_swiz","blk_swiz"};

// 2-param kernels (TX,TY)
#define D2_BODY(KERN) \
    if(tx==1&&ty==1) KERN<1,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1) KERN<2,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2) KERN<1,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2) KERN<2,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1) KERN<4,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4) KERN<1,4><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4) KERN<4,4><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2) KERN<4,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4) KERN<2,4><<<g,b,sm>>>(in,out,N); \
    else {fprintf(stderr,"Bad TX=%d TY=%d\n",tx,ty);exit(1);}

void d_naive(int tx,int ty,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){D2_BODY(tr_naive)}
void d_smem(int tx,int ty,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){D2_BODY(tr_smem)}
void d_smem_swiz(int tx,int ty,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){D2_BODY(tr_smem_swiz)}

// 3-param kernels (TX,TY,SB)
#define D3_BODY(KERN,SB_VAL) \
    if(tx==1&&ty==1) KERN<1,1,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1) KERN<2,1,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2) KERN<1,2,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2) KERN<2,2,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1) KERN<4,1,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4) KERN<1,4,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4) KERN<4,4,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2) KERN<4,2,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4) KERN<2,4,SB_VAL><<<g,b,sm>>>(in,out,N); \
    else {fprintf(stderr,"Bad TX=%d TY=%d\n",tx,ty);exit(1);}

void d_blocked(int tx,int ty,int sb,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){
    if(sb==16){D3_BODY(tr_blocked,16)} else if(sb==32){D3_BODY(tr_blocked,32)}
    else if(sb==64){D3_BODY(tr_blocked,64)} else{fprintf(stderr,"Bad SB=%d\n",sb);exit(1);}
}
void d_smem_blk(int tx,int ty,int sb,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){
    if(sb==16){D3_BODY(tr_smem_blk,16)} else if(sb==32){D3_BODY(tr_smem_blk,32)}
    else if(sb==64){D3_BODY(tr_smem_blk,64)} else{fprintf(stderr,"Bad SB=%d\n",sb);exit(1);}
}
void d_blk_swiz(int tx,int ty,int sb,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){
    if(sb==16){D3_BODY(tr_blk_swiz,16)} else if(sb==32){D3_BODY(tr_blk_swiz,32)}
    else if(sb==64){D3_BODY(tr_blk_swiz,64)} else{fprintf(stderr,"Bad SB=%d\n",sb);exit(1);}
}

// smem_pad: (TX,TY,PAD)
#define DPAD_BODY(PAD_VAL) \
    if(tx==1&&ty==1) tr_smem_pad<1,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1) tr_smem_pad<2,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2) tr_smem_pad<1,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2) tr_smem_pad<2,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1) tr_smem_pad<4,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4) tr_smem_pad<1,4,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4) tr_smem_pad<4,4,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2) tr_smem_pad<4,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4) tr_smem_pad<2,4,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else {fprintf(stderr,"Bad TX=%d TY=%d\n",tx,ty);exit(1);}

void d_smem_pad(int tx,int ty,int pad,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){
    if(pad==1){DPAD_BODY(1)} else if(pad==2){DPAD_BODY(2)}
    else{fprintf(stderr,"Bad PAD=%d\n",pad);exit(1);}
}

int main(int argc, char** argv) {
    if (argc < 8) {
        fprintf(stderr,"Usage: %s <N> <variant> <BX> <BY> <TX> <TY> <csv> [SB=32] [PAD=1] [WARMUP=5] [REPS=100]\n",argv[0]);
        fprintf(stderr,"  variant: 0=naive 1=blocked 2=smem 3=smem_blk 4=smem_pad 5=smem_swiz 6=blk_swiz\n");
        return 1;
    }
    int N=atoi(argv[1]), VAR=atoi(argv[2]);
    int BX=atoi(argv[3]), BY=atoi(argv[4]), TX=atoi(argv[5]), TY=atoi(argv[6]);
    const char* csv=argv[7];
    int SB=(argc>8)?atoi(argv[8]):32;
    int PAD=(argc>9)?atoi(argv[9]):1;
    int WARMUP=(argc>10)?atoi(argv[10]):5;
    int REPS=(argc>11)?atoi(argv[11]):100;

    size_t bytes = (size_t)N * N * sizeof(double);
    double* h = (double*)malloc(bytes);
    for (size_t i = 0; i < (size_t)N*N; i++) h[i] = (double)i / N;

    double *dI, *dO;
    CHECK(hipMalloc(&dI, bytes)); CHECK(hipMalloc(&dO, bytes));

    // For blocked variants, convert to blocked layout on host
    double* hb = nullptr;
    if (VAR == 1 || VAR == 3 || VAR == 6) {
        hb = (double*)malloc(bytes);
        int NB = N / SB;
        for (int r = 0; r < N; r++)
            for (int c = 0; c < N; c++)
                hb[(r/SB * NB + c/SB) * SB*SB + (r%SB)*SB + c%SB] = h[r * N + c];
        CHECK(hipMemcpy(dI, hb, bytes, hipMemcpyHostToDevice));
    } else {
        CHECK(hipMemcpy(dI, h, bytes, hipMemcpyHostToDevice));
    }

    dim3 block(BX, BY);
    dim3 grid((N + BX*TX - 1)/(BX*TX), (N + BY*TY - 1)/(BY*TY));
    int BW = BX*TX, BH = BY*TY;
    size_t sm_np = BH * BW * sizeof(double);
    size_t sm_pd = BH * (BW + PAD) * sizeof(double);

    auto launch = [&](const double* in, double* out) {
        switch(VAR) {
            case 0: d_naive(TX,TY,grid,block,0,in,out,N); break;
            case 1: d_blocked(TX,TY,SB,grid,block,0,in,out,N); break;
            case 2: d_smem(TX,TY,grid,block,sm_np,in,out,N); break;
            case 3: d_smem_blk(TX,TY,SB,grid,block,sm_np,in,out,N); break;
            case 4: d_smem_pad(TX,TY,PAD,grid,block,sm_pd,in,out,N); break;
            case 5: d_smem_swiz(TX,TY,grid,block,sm_np,in,out,N); break;
            case 6: d_blk_swiz(TX,TY,SB,grid,block,sm_np,in,out,N); break;
        }
    };

    // Warmup
    for (int i = 0; i < WARMUP; i++) launch(dI, dO);
    CHECK(hipDeviceSynchronize());

    // Timed runs
    hipEvent_t* ev = (hipEvent_t*)malloc((REPS+1)*sizeof(hipEvent_t));
    for (int i = 0; i <= REPS; i++) CHECK(hipEventCreate(&ev[i]));
    for (int i = 0; i < REPS; i++) {
        CHECK(hipEventRecord(ev[i]));
        launch(dI, dO);
    }
    CHECK(hipEventRecord(ev[REPS]));
    CHECK(hipEventSynchronize(ev[REPS]));

    float* ims = (float*)malloc(REPS * sizeof(float));
    float tot = 0;
    for (int i = 0; i < REPS; i++) {
        CHECK(hipEventElapsedTime(&ims[i], ev[i], ev[i+1]));
        tot += ims[i];
    }

    // Checksum
    CHECK(hipMemcpy(h, dO, bytes, hipMemcpyDeviceToHost));
    double cksum = 0;
    for (size_t i = 0; i < (size_t)N*N; i++) cksum += h[i];

    // 2 * N^2 * 8 bytes (read + write)
    double bpi = 2.0 * N * N * 8;
    double med_ms = tot / REPS;
    double gbps = bpi / (med_ms / 1000.0) / 1e9;

    printf("%s N=%d BX=%d BY=%d TX=%d TY=%d SB=%d PAD=%d | %.4f ms  %.1f GB/s  cksum=%.6e\n",
           V_NAMES[VAR], N, BX, BY, TX, TY, SB, PAD, med_ms, gbps, cksum);

    FILE* f = fopen(csv, "a");
    if (f) {
        for (int i = 0; i < REPS; i++) {
            double gbs = bpi / (ims[i]/1000.0) / 1e9;
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%.3f,%.6e\n",
                    V_NAMES[VAR], N, BX, BY, TX, TY, SB, PAD, i, ims[i]/1000.0, gbs, cksum);
        }
        fclose(f);
    }

    for (int i = 0; i <= REPS; i++) CHECK(hipEventDestroy(ev[i]));
    free(ev); free(ims); free(h); if(hb) free(hb);
    CHECK(hipFree(dI)); CHECK(hipFree(dO));
    return 0;
}
