#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

#define CHECK(x) do { cudaError_t e = (x); if(e) { fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1); } } while(0)

// ── Stencil traits ──
// ST=0: 5-point cross,  halo=1, 5  FLOPs (4 add + 1 mul)
// ST=1: 3x3 box,        halo=1, 9  FLOPs (8 add + 1 div)
// ST=2: 4x4 box,        halo=2, 16 FLOPs (15 add + 1 div), offsets [-1..+2]
template<int ST> struct Traits;
template<> struct Traits<0> { static constexpr int H=1, F=5; };
template<> struct Traits<1> { static constexpr int H=1, F=9; };
template<> struct Traits<2> { static constexpr int H=2, F=16; };

template<int ST>
__device__ __forceinline__ double apply(const double* p, int s) {
    if constexpr (ST == 0)
        return 0.2 * (p[0] + p[-1] + p[1] + p[-s] + p[s]);
    else if constexpr (ST == 1) {
        double r = 0;
        #pragma unroll
        for (int di = -1; di <= 1; di++)
            #pragma unroll
            for (int dj = -1; dj <= 1; dj++)
                r += p[di * s + dj];
        return r / 9.0;
    } else {
        double r = 0;
        #pragma unroll
        for (int di = -1; di <= 2; di++)
            #pragma unroll
            for (int dj = -1; dj <= 2; dj++)
                r += p[di * s + dj];
        return r / 16.0;
    }
}

// ── V0: Global memory ──
template<int TX, int TY, int ST>
__global__ void kern_global(const double* __restrict__ in, double* __restrict__ out, int N) {
    constexpr int H = Traits<ST>::H;
    int j0 = blockIdx.x * (blockDim.x * TX) + threadIdx.x * TX + H;
    int i0 = blockIdx.y * (blockDim.y * TY) + threadIdx.y * TY + H;
    int stride = N + 2 * H;
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        int i = i0 + di;
        if (i > N + H - 1) continue;
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int j = j0 + dj;
            if (j > N + H - 1) continue;
            out[i * stride + j] = apply<ST>(&in[i * stride + j], stride);
        }
    }
}

// ── V1: Shared memory, naive ──
template<int TX, int TY, int ST>
__global__ void kern_smem(const double* __restrict__ in, double* __restrict__ out, int N) {
    constexpr int H = Traits<ST>::H;
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    const int stride = N + 2 * H;
    extern __shared__ double sm[];
    int SW = BW + 2 * H;
    int gj0 = blockIdx.x * BW, gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, total = blockDim.x * blockDim.y;
    int elems = (BH + 2 * H) * SW;
    for (int k = tid; k < elems; k += total) {
        int si = k / SW, sj = k % SW;
        int gi = gi0 + si, gj = gj0 + sj;
        sm[k] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int li = threadIdx.y * TY + H + di, lj = threadIdx.x * TX + H + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi >= H && gi < N + H && gj >= H && gj < N + H)
                out[gi * stride + gj] = apply<ST>(&sm[li * SW + lj], SW);
        }
    }
}

// ── V2: Shared memory with padding ──
template<int TX, int TY, int ST, int PAD>
__global__ void kern_smem_pad(const double* __restrict__ in, double* __restrict__ out, int N) {
    constexpr int H = Traits<ST>::H;
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    const int stride = N + 2 * H;
    extern __shared__ double sm[];
    int SW = BW + 2 * H + PAD;
    int load_w = BW + 2 * H;
    int gj0 = blockIdx.x * BW, gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, total = blockDim.x * blockDim.y;
    int elems = (BH + 2 * H) * load_w;
    for (int k = tid; k < elems; k += total) {
        int si = k / load_w, sj = k % load_w;
        int gi = gi0 + si, gj = gj0 + sj;
        sm[si * SW + sj] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int li = threadIdx.y * TY + H + di, lj = threadIdx.x * TX + H + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi >= H && gi < N + H && gj >= H && gj < N + H)
                out[gi * stride + gj] = apply<ST>(&sm[li * SW + lj], SW);
        }
    }
}

// ── V3: Shared memory, swapped loop order ──
template<int TX, int TY, int ST>
__global__ void kern_smem_swap(const double* __restrict__ in, double* __restrict__ out, int N) {
    constexpr int H = Traits<ST>::H;
    const int BW = blockDim.x * TX, BH = blockDim.y * TY;
    const int stride = N + 2 * H;
    extern __shared__ double sm[];
    int SW = BW + 2 * H;
    int gj0 = blockIdx.x * BW, gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x, total = blockDim.x * blockDim.y;
    int elems = (BH + 2 * H) * SW;
    for (int k = tid; k < elems; k += total) {
        int si = k / SW, sj = k % SW;
        int gi = gi0 + si, gj = gj0 + sj;
        sm[k] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    #pragma unroll
    for (int dj = 0; dj < TX; dj++) {
        #pragma unroll
        for (int di = 0; di < TY; di++) {
            int li = threadIdx.y * TY + H + di, lj = threadIdx.x * TX + H + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi >= H && gi < N + H && gj >= H && gj < N + H)
                out[gi * stride + gj] = apply<ST>(&sm[li * SW + lj], SW);
        }
    }
}

// ── Dispatch: separate function per kernel to avoid nvcc template-template + <<<>>> parse bug ──
#define DISPATCH_BODY(KERN) \
    if(tx==1&&ty==1&&st==0) KERN<1,1,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==1&&st==1) KERN<1,1,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==1&&st==2) KERN<1,1,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==0) KERN<1,2,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==1) KERN<1,2,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==2) KERN<1,2,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==0) KERN<1,4,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==1) KERN<1,4,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==2) KERN<1,4,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==0) KERN<2,1,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==1) KERN<2,1,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==2) KERN<2,1,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==0) KERN<2,2,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==1) KERN<2,2,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==2) KERN<2,2,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==0) KERN<2,4,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==1) KERN<2,4,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==2) KERN<2,4,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==0) KERN<4,1,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==1) KERN<4,1,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==2) KERN<4,1,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==0) KERN<4,2,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==1) KERN<4,2,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==2) KERN<4,2,2><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==0) KERN<4,4,0><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==1) KERN<4,4,1><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==2) KERN<4,4,2><<<g,b,sm>>>(in,out,N); \
    else { fprintf(stderr,"Bad TX=%d TY=%d ST=%d\n",tx,ty,st); exit(1); }

void d_global(int tx,int ty,int st,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){ DISPATCH_BODY(kern_global) }
void d_smem  (int tx,int ty,int st,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){ DISPATCH_BODY(kern_smem) }
void d_swap  (int tx,int ty,int st,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){ DISPATCH_BODY(kern_smem_swap) }

// smem_pad has 4 template args — expand PAD dimension too
#define DISPATCH_PAD_BODY(PAD_VAL) \
    if(tx==1&&ty==1&&st==0) kern_smem_pad<1,1,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==1&&st==1) kern_smem_pad<1,1,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==1&&st==2) kern_smem_pad<1,1,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==0) kern_smem_pad<1,2,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==1) kern_smem_pad<1,2,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==2&&st==2) kern_smem_pad<1,2,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==0) kern_smem_pad<1,4,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==1) kern_smem_pad<1,4,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==1&&ty==4&&st==2) kern_smem_pad<1,4,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==0) kern_smem_pad<2,1,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==1) kern_smem_pad<2,1,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==1&&st==2) kern_smem_pad<2,1,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==0) kern_smem_pad<2,2,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==1) kern_smem_pad<2,2,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==2&&st==2) kern_smem_pad<2,2,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==0) kern_smem_pad<2,4,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==1) kern_smem_pad<2,4,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==2&&ty==4&&st==2) kern_smem_pad<2,4,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==0) kern_smem_pad<4,1,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==1) kern_smem_pad<4,1,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==1&&st==2) kern_smem_pad<4,1,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==0) kern_smem_pad<4,2,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==1) kern_smem_pad<4,2,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==2&&st==2) kern_smem_pad<4,2,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==0) kern_smem_pad<4,4,0,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==1) kern_smem_pad<4,4,1,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else if(tx==4&&ty==4&&st==2) kern_smem_pad<4,4,2,PAD_VAL><<<g,b,sm>>>(in,out,N); \
    else { fprintf(stderr,"Bad TX=%d TY=%d ST=%d PAD=%d\n",tx,ty,st,pad); exit(1); }

void d_pad(int tx,int ty,int st,int pad,dim3 g,dim3 b,size_t sm,const double*in,double*out,int N){
    if(pad==1){ DISPATCH_PAD_BODY(1) }
    else if(pad==2){ DISPATCH_PAD_BODY(2) }
    else { fprintf(stderr,"Bad PAD=%d\n",pad); exit(1); }
}

static const int HALO_TAB[] = {1, 1, 2};
static const int FLOP_TAB[] = {5, 9, 16};
static const char* ST_NAMES[] = {"5pt","box3","box4"};
static const char* V_NAMES[]  = {"global","smem","smem_pad","smem_swap"};

int main(int argc, char** argv) {
    if (argc < 9) {
        fprintf(stderr, "Usage: %s <N> <tsteps> <stencil> <variant> <BX> <BY> <TX> <TY> [PAD] [csv]\n", argv[0]);
        fprintf(stderr, "  stencil: 0=5pt, 1=box3x3, 2=box4x4\n");
        fprintf(stderr, "  variant: 0=global, 1=smem, 2=smem_pad, 3=smem_swap\n");
        return 1;
    }
    int N   = atoi(argv[1]), tsteps = atoi(argv[2]);
    int ST  = atoi(argv[3]), VAR    = atoi(argv[4]);
    int BX  = atoi(argv[5]), BY     = atoi(argv[6]);
    int TX  = atoi(argv[7]), TY     = atoi(argv[8]);
    int PAD = (argc > 9)  ? atoi(argv[9])  : 1;
    const char* csv = (argc > 10) ? argv[10] : nullptr;

    int H = HALO_TAB[ST];
    size_t S = N + 2 * H;
    size_t bytes = S * S * sizeof(double);

    double* h = (double*)malloc(bytes);
    for (size_t i = 0; i < S; i++)
        for (size_t j = 0; j < S; j++)
            h[i * S + j] = (double)(i * (j + 2)) / N;

    double *dA, *dB;
    CHECK(cudaMalloc(&dA, bytes));
    CHECK(cudaMalloc(&dB, bytes));
    CHECK(cudaMemcpy(dA, h, bytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dB, h, bytes, cudaMemcpyHostToDevice));

    dim3 block(BX, BY);
    dim3 grid((N + BX*TX - 1) / (BX*TX), (N + BY*TY - 1) / (BY*TY));

    int tw = BX*TX + 2*H, th = BY*TY + 2*H;
    size_t sm_np = th * tw * sizeof(double);
    size_t sm_pd = th * (tw + PAD) * sizeof(double);

    auto launch = [&](const double* in, double* out) {
        switch(VAR) {
            case 0: d_global(TX,TY,ST,grid,block,0,in,out,N); break;
            case 1: d_smem(TX,TY,ST,grid,block,sm_np,in,out,N); break;
            case 2: d_pad(TX,TY,ST,PAD,grid,block,sm_pd,in,out,N); break;
            case 3: d_swap(TX,TY,ST,grid,block,sm_np,in,out,N); break;
        }
    };

    for (int t = 0; t < 3; t++) launch((t&1)?dB:dA, (t&1)?dA:dB);
    CHECK(cudaDeviceSynchronize());
    CHECK(cudaMemcpy(dA, h, bytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(dB, h, bytes, cudaMemcpyHostToDevice));

    cudaEvent_t* ev = (cudaEvent_t*)malloc((tsteps+1)*sizeof(cudaEvent_t));
    for (int t = 0; t <= tsteps; t++) CHECK(cudaEventCreate(&ev[t]));
    for (int t = 0; t < tsteps; t++) {
        CHECK(cudaEventRecord(ev[t]));
        launch((t&1)?dB:dA, (t&1)?dA:dB);
    }
    CHECK(cudaEventRecord(ev[tsteps]));
    CHECK(cudaEventSynchronize(ev[tsteps]));

    float* ims = (float*)malloc(tsteps * sizeof(float));
    float tot = 0;
    for (int t = 0; t < tsteps; t++) {
        CHECK(cudaEventElapsedTime(&ims[t], ev[t], ev[t+1]));
        tot += ims[t];
    }
    double secs = tot / 1000.0;

    double* res = (tsteps&1) ? dB : dA;
    CHECK(cudaMemcpy(h, res, bytes, cudaMemcpyDeviceToHost));
    double cksum = 0;
    for (size_t i = H; i < (size_t)(N+H); i++)
        for (size_t j = H; j < (size_t)(N+H); j++)
            cksum += h[i * S + j];

    double fpi = (double)N * N * FLOP_TAB[ST];
    double gflops = fpi * tsteps / secs / 1e9;

    printf("%s %s BX=%d BY=%d TX=%d TY=%d PAD=%d | %.4f s  %.2f GFLOP/s  cksum=%.6e\n",
           ST_NAMES[ST], V_NAMES[VAR], BX, BY, TX, TY, PAD, secs, gflops, cksum);

    if (csv) {
        FILE* f = fopen(csv, "a");
        if (f) {
            for (int t = 0; t < tsteps; t++) {
                double is = ims[t]/1000.0;
                fprintf(f, "%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%.6f,%.3f,%.6e\n",
                        ST_NAMES[ST], V_NAMES[VAR], N, tsteps, BX, BY, TX, TY, PAD,
                        t, is, fpi/is/1e9, cksum);
            }
            fclose(f);
        }
    }

    for (int t = 0; t <= tsteps; t++) CHECK(cudaEventDestroy(ev[t]));
    free(ev); free(ims);
    CHECK(cudaFree(dA)); CHECK(cudaFree(dB)); free(h);
    return 0;
}
