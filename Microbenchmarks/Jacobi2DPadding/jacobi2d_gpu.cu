#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cuda_runtime.h>

#define CHECK(x) do { cudaError_t e = (x); if(e) { fprintf(stderr,"CUDA %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e)); exit(1); } } while(0)

// V0: Global memory only, each thread computes TX x TY points
template<int TX, int TY>
__global__ void jacobi_global(const double* __restrict__ in, double* __restrict__ out, int N) {
    int j0 = blockIdx.x * (blockDim.x * TX) + threadIdx.x * TX + 1;
    int i0 = blockIdx.y * (blockDim.y * TY) + threadIdx.y * TY + 1;
    int stride = N + 2;
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        int i = i0 + di;
        if (i > N) continue;
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int j = j0 + dj;
            if (j > N) continue;
            int idx = i * stride + j;
            out[idx] = 0.2 * (in[idx] + in[idx-1] + in[idx+1] + in[idx-stride] + in[idx+stride]);
        }
    }
}

// V1: Shared memory, naive
template<int TX, int TY>
__global__ void jacobi_smem(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX;
    const int BH = blockDim.y * TY;
    extern __shared__ double smem[];  // (BH+2) x (BW+2)
    int SW = BW + 2;
    int gj0 = blockIdx.x * BW;  // global col start (0-indexed in padded grid)
    int gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    int total = blockDim.x * blockDim.y;
    int stride = N + 2;
    // Cooperative load: tile of (BH+2) x (BW+2)
    int elems = (BH + 2) * SW;
    for (int k = tid; k < elems; k += total) {
        int si = k / SW, sj = k % SW;
        int gi = gi0 + si, gj = gj0 + sj;
        smem[si * SW + sj] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    int lj0 = threadIdx.x * TX + 1;
    int li0 = threadIdx.y * TY + 1;
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int li = li0 + di, lj = lj0 + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi <= N && gj <= N) {
                int sidx = li * SW + lj;
                out[gi * stride + gj] = 0.2 * (smem[sidx] + smem[sidx-1] + smem[sidx+1] + smem[sidx-SW] + smem[sidx+SW]);
            }
        }
    }
}

// V2: Shared memory with padding to reduce bank conflicts
template<int TX, int TY, int PAD>
__global__ void jacobi_smem_pad(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX;
    const int BH = blockDim.y * TY;
    extern __shared__ double smem[];  // (BH+2) x (BW+2+PAD)
    int SW = BW + 2 + PAD;
    int gj0 = blockIdx.x * BW;
    int gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    int total = blockDim.x * blockDim.y;
    int stride = N + 2;
    int load_w = BW + 2;  // actual data width
    int elems = (BH + 2) * load_w;
    for (int k = tid; k < elems; k += total) {
        int si = k / load_w, sj = k % load_w;
        int gi = gi0 + si, gj = gj0 + sj;
        smem[si * SW + sj] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    int lj0 = threadIdx.x * TX + 1;
    int li0 = threadIdx.y * TY + 1;
    #pragma unroll
    for (int di = 0; di < TY; di++) {
        #pragma unroll
        for (int dj = 0; dj < TX; dj++) {
            int li = li0 + di, lj = lj0 + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi <= N && gj <= N) {
                int sidx = li * SW + lj;
                out[gi * stride + gj] = 0.2 * (smem[sidx] + smem[sidx-1] + smem[sidx+1] + smem[sidx-SW] + smem[sidx+SW]);
            }
        }
    }
}

// V3: Shared memory, swapped compute loop order (iterate j in outer, i in inner)
template<int TX, int TY>
__global__ void jacobi_smem_swap(const double* __restrict__ in, double* __restrict__ out, int N) {
    const int BW = blockDim.x * TX;
    const int BH = blockDim.y * TY;
    extern __shared__ double smem[];
    int SW = BW + 2;
    int gj0 = blockIdx.x * BW;
    int gi0 = blockIdx.y * BH;
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    int total = blockDim.x * blockDim.y;
    int stride = N + 2;
    int elems = (BH + 2) * SW;
    for (int k = tid; k < elems; k += total) {
        int si = k / SW, sj = k % SW;
        int gi = gi0 + si, gj = gj0 + sj;
        smem[si * SW + sj] = (gi < stride && gj < stride) ? in[gi * stride + gj] : 0.0;
    }
    __syncthreads();
    int lj0 = threadIdx.x * TX + 1;
    int li0 = threadIdx.y * TY + 1;
    // Swapped: j-outer, i-inner -> different bank access pattern
    #pragma unroll
    for (int dj = 0; dj < TX; dj++) {
        #pragma unroll
        for (int di = 0; di < TY; di++) {
            int li = li0 + di, lj = lj0 + dj;
            int gi = gi0 + li, gj = gj0 + lj;
            if (gi <= N && gj <= N) {
                int sidx = li * SW + lj;
                out[gi * stride + gj] = 0.2 * (smem[sidx] + smem[sidx-1] + smem[sidx+1] + smem[sidx-SW] + smem[sidx+SW]);
            }
        }
    }
}

// Dispatch helpers — instantiate for all TX,TY combos we care about
typedef void (*kernel_fn)(const double*, double*, int);
typedef void (*kernel_pad_fn)(const double*, double*, int);

// We use a macro to stamp out the switch. TX,TY in {1,2,4}.
#define DISPATCH_TXTY(KERN, tx, ty, grid, block, smem_bytes, in, out, N) \
    if(tx==1&&ty==1) KERN<1,1><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==2&&ty==1) KERN<2,1><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==1&&ty==2) KERN<1,2><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==2&&ty==2) KERN<2,2><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==4&&ty==1) KERN<4,1><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==1&&ty==4) KERN<1,4><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==4&&ty==4) KERN<4,4><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==2&&ty==4) KERN<2,4><<<grid,block,smem_bytes>>>(in,out,N); \
    else if(tx==4&&ty==2) KERN<4,2><<<grid,block,smem_bytes>>>(in,out,N); \
    else { fprintf(stderr,"Unsupported TX=%d TY=%d\n",tx,ty); exit(1); }

#define DISPATCH_TXTY_PAD(tx, ty, pad, grid, block, smem_bytes, in, out, N) \
    if(pad==1) { \
        if(tx==1&&ty==1) jacobi_smem_pad<1,1,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==1) jacobi_smem_pad<2,1,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==1&&ty==2) jacobi_smem_pad<1,2,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==2) jacobi_smem_pad<2,2,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==1) jacobi_smem_pad<4,1,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==1&&ty==4) jacobi_smem_pad<1,4,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==4) jacobi_smem_pad<4,4,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==4) jacobi_smem_pad<2,4,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==2) jacobi_smem_pad<4,2,1><<<grid,block,smem_bytes>>>(in,out,N); \
        else { fprintf(stderr,"Unsupported TX=%d TY=%d\n",tx,ty); exit(1); } \
    } else if(pad==2) { \
        if(tx==1&&ty==1) jacobi_smem_pad<1,1,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==1) jacobi_smem_pad<2,1,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==1&&ty==2) jacobi_smem_pad<1,2,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==2) jacobi_smem_pad<2,2,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==1) jacobi_smem_pad<4,1,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==1&&ty==4) jacobi_smem_pad<1,4,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==4) jacobi_smem_pad<4,4,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==2&&ty==4) jacobi_smem_pad<2,4,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else if(tx==4&&ty==2) jacobi_smem_pad<4,2,2><<<grid,block,smem_bytes>>>(in,out,N); \
        else { fprintf(stderr,"Unsupported TX=%d TY=%d\n",tx,ty); exit(1); } \
    } else { fprintf(stderr,"Unsupported PAD=%d\n",pad); exit(1); }

int main(int argc, char** argv) {
    if (argc < 8) {
        fprintf(stderr, "Usage: %s <N> <tsteps> <variant> <BX> <BY> <TX> <TY> [PAD] [csv_file]\n", argv[0]);
        fprintf(stderr, "  variant: 0=global, 1=smem, 2=smem_pad, 3=smem_swap\n");
        return 1;
    }
    int N       = atoi(argv[1]);
    int tsteps  = atoi(argv[2]);
    int variant = atoi(argv[3]);
    int BX      = atoi(argv[4]);
    int BY      = atoi(argv[5]);
    int TX      = atoi(argv[6]);
    int TY      = atoi(argv[7]);
    int PAD     = (argc > 8) ? atoi(argv[8]) : 1;
    const char* csv = (argc > 9) ? argv[9] : nullptr;

    size_t S = N + 2;
    size_t bytes = S * S * sizeof(double);

    // Init host
    double* h = (double*)malloc(bytes);
    for (size_t i = 0; i < S; i++)
        for (size_t j = 0; j < S; j++)
            h[i * S + j] = (double)(i * (j + 2)) / N;

    double *d_A, *d_B;
    CHECK(cudaMalloc(&d_A, bytes));
    CHECK(cudaMalloc(&d_B, bytes));
    CHECK(cudaMemcpy(d_A, h, bytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_B, h, bytes, cudaMemcpyHostToDevice));  // init B too for boundary

    dim3 block(BX, BY);
    dim3 grid((N + BX * TX - 1) / (BX * TX), (N + BY * TY - 1) / (BY * TY));

    int tile_w = BX * TX + 2;
    int tile_h = BY * TY + 2;
    size_t smem_nopad = tile_h * tile_w * sizeof(double);
    size_t smem_pad   = tile_h * (tile_w + PAD) * sizeof(double);

    // Warmup
    for (int t = 0; t < 3; t++) {
        double *in = (t&1) ? d_B : d_A, *out = (t&1) ? d_A : d_B;
        switch(variant) {
            case 0: DISPATCH_TXTY(jacobi_global, TX, TY, grid, block, 0, in, out, N); break;
            case 1: DISPATCH_TXTY(jacobi_smem, TX, TY, grid, block, smem_nopad, in, out, N); break;
            case 2: DISPATCH_TXTY_PAD(TX, TY, PAD, grid, block, smem_pad, in, out, N); break;
            case 3: DISPATCH_TXTY(jacobi_smem_swap, TX, TY, grid, block, smem_nopad, in, out, N); break;
        }
    }
    CHECK(cudaDeviceSynchronize());

    // Re-init
    CHECK(cudaMemcpy(d_A, h, bytes, cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(d_B, h, bytes, cudaMemcpyHostToDevice));

    cudaEvent_t t0, t1;
    CHECK(cudaEventCreate(&t0));
    CHECK(cudaEventCreate(&t1));
    CHECK(cudaEventRecord(t0));

    for (int t = 0; t < tsteps; t++) {
        double *in = (t&1) ? d_B : d_A, *out = (t&1) ? d_A : d_B;
        switch(variant) {
            case 0: DISPATCH_TXTY(jacobi_global, TX, TY, grid, block, 0, in, out, N); break;
            case 1: DISPATCH_TXTY(jacobi_smem, TX, TY, grid, block, smem_nopad, in, out, N); break;
            case 2: DISPATCH_TXTY_PAD(TX, TY, PAD, grid, block, smem_pad, in, out, N); break;
            case 3: DISPATCH_TXTY(jacobi_smem_swap, TX, TY, grid, block, smem_nopad, in, out, N); break;
        }
    }
    CHECK(cudaEventRecord(t1));
    CHECK(cudaEventSynchronize(t1));
    float ms;
    CHECK(cudaEventElapsedTime(&ms, t0, t1));
    double secs = ms / 1000.0;

    // Checksum
    double* res = (tsteps & 1) ? d_B : d_A;
    CHECK(cudaMemcpy(h, res, bytes, cudaMemcpyDeviceToHost));
    double cksum = 0;
    for (size_t i = 1; i <= (size_t)N; i++)
        for (size_t j = 1; j <= (size_t)N; j++)
            cksum += h[i * S + j];

    double gflops = (double)N * N * tsteps * 5.0 / secs / 1e9;
    const char* vnames[] = {"global","smem","smem_pad","smem_swap"};

    printf("%s BX=%d BY=%d TX=%d TY=%d PAD=%d | %.4f s  %.2f GFLOP/s  cksum=%.6e\n",
           vnames[variant], BX, BY, TX, TY, PAD, secs, gflops, cksum);

    if (csv) {
        FILE* f = fopen(csv, "a");
        if (f) {
            fprintf(f, "%s,%d,%d,%d,%d,%d,%d,%d,%.6f,%.3f,%.6e\n",
                    vnames[variant], N, tsteps, BX, BY, TX, TY, PAD, secs, gflops, cksum);
            fclose(f);
        }
    }

    CHECK(cudaFree(d_A));
    CHECK(cudaFree(d_B));
    free(h);
    return 0;
}
