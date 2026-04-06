/*
 * bench_gpu_profile.cu -- 6-kernel GPU profiling for nsight compute
 *
 * Runs EXACT distribution only, 6 specific kernel configs:
 *   K1: V1 + 1x4_96x1   (best orange, new)
 *   K2: V4 + 1x4_96x1   (same config, transformed layout)
 *   K3: V5 + 1x4_96x1   (same config, transformed+padded)
 *   K4: V5 + 2x1_16x16  (best blue, old results)
 *   K5: V6 + v6_128      (best blue, new results)
 *   K6: V1 + 1x8_256x1   (best orange, old results)
 *
 * Compile: nvcc -O3 -arch=sm_80 -std=c++17 -Xcompiler -fopenmp bench_gpu_profile.cu -o bench_gpu_profile
 * Profile: ncu --set full -o profile_report ./bench_gpu_profile
 */
#include "bench_common.h"
#include "icon_data_loader.h"
#include <cassert>

#if __HIP_PLATFORM_AMD__
#include "hip/hip_runtime.h"
#endif

#define CUDA_CHECK(call) do { \
    cudaError_t e=(call); \
    if(e!=cudaSuccess){fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);} \
} while(0)
#define CUDA_LAUNCH_CHECK() do { \
    cudaError_t e=cudaGetLastError(); \
    if(e!=cudaSuccess){fprintf(stderr,"CUDA launch %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(e));exit(1);} \
} while(0)

/* ================================================================ */
/*  GPU kernels (same as bench_gpu.cu)                               */
/* ================================================================ */
template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_je_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    const int je_stride=(int)gridDim.x*BX*TX;
    const int jk_stride=(int)gridDim.y*BY*TY;
    for(int je_base0=((int)blockIdx.x*BX+(int)threadIdx.x)*TX;je_base0<N_e;je_base0+=je_stride){
        int ci0_a[TX],ci1_a[TX],vi0_a[TX],vi1_a[TX];
        double id_a[TX],ip_a[TX],tg_a[TX];
        #pragma unroll
        for(int tx=0;tx<TX;tx++){int je=je_base0+tx;if(je<N_e){
            ci0_a[tx]=cell_idx[IN<V>(je,0,N_e)];ci1_a[tx]=cell_idx[IN<V>(je,1,N_e)];
            vi0_a[tx]=vert_idx[IN<V>(je,0,N_e)];vi1_a[tx]=vert_idx[IN<V>(je,1,N_e)];
            id_a[tx]=inv_dual[je];ip_a[tx]=inv_primal[je];tg_a[tx]=tangent[je];
        }}
        for(int jk_base0=((int)blockIdx.y*BY+(int)threadIdx.y)*TY;jk_base0<nlev_end;jk_base0+=jk_stride){
            #pragma unroll
            for(int ty=0;ty<TY;ty++){int jk=jk_base0+ty;if(jk>=nlev_end)continue;
            #pragma unroll
            for(int tx=0;tx<TX;tx++){int je=je_base0+tx;if(je>=N_e)continue;
                int c2d=IC<V>(je,jk,N_e,nlev);
                out[c2d]=vn_ie[c2d]*id_a[tx]*
                    (w[IC<V>(ci0_a[tx],jk,N_c,nlev)]-w[IC<V>(ci1_a[tx],jk,N_c,nlev)])
                  +z_vt_ie[c2d]*ip_a[tx]*tg_a[tx]*
                    (z_w_v[IC<V>(vi0_a[tx],jk,N_v,nlev)]-z_w_v[IC<V>(vi1_a[tx],jk,N_v,nlev)]);
            }}
        }
    }
}

template<int TX, int TY, int BX, int BY, int V>
__global__ void gpu_kernel_jk_first(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    const int jk_stride=(int)gridDim.x*BX*TX;
    const int je_stride=(int)gridDim.y*BY*TY;
    for(int je_base0=((int)blockIdx.y*BY+(int)threadIdx.y)*TY;je_base0<N_e;je_base0+=je_stride){
        int ci0_a[TY],ci1_a[TY],vi0_a[TY],vi1_a[TY];
        double id_a[TY],ip_a[TY],tg_a[TY];
        #pragma unroll
        for(int ty=0;ty<TY;ty++){int je=je_base0+ty;if(je<N_e){
            ci0_a[ty]=cell_idx[IN<V>(je,0,N_e)];ci1_a[ty]=cell_idx[IN<V>(je,1,N_e)];
            vi0_a[ty]=vert_idx[IN<V>(je,0,N_e)];vi1_a[ty]=vert_idx[IN<V>(je,1,N_e)];
            id_a[ty]=inv_dual[je];ip_a[ty]=inv_primal[je];tg_a[ty]=tangent[je];
        }}
        for(int jk_base0=((int)blockIdx.x*BX+(int)threadIdx.x)*TX;jk_base0<nlev_end;jk_base0+=jk_stride){
            #pragma unroll
            for(int ty=0;ty<TY;ty++){int je=je_base0+ty;if(je>=N_e)continue;
            #pragma unroll
            for(int tx=0;tx<TX;tx++){int jk=jk_base0+tx;if(jk>=nlev_end)continue;
                int c2d=IC<V>(je,jk,N_e,nlev);
                out[c2d]=vn_ie[c2d]*id_a[ty]*
                    (w[IC<V>(ci0_a[ty],jk,N_c,nlev)]-w[IC<V>(ci1_a[ty],jk,N_c,nlev)])
                  +z_vt_ie[c2d]*ip_a[ty]*tg_a[ty]*
                    (z_w_v[IC<V>(vi0_a[ty],jk,N_v,nlev)]-z_w_v[IC<V>(vi1_a[ty],jk,N_v,nlev)]);
            }}
        }
    }
}

template<unsigned BX, int V>
__global__ void gpu_kernel_v6_flat(
    double* __restrict__ out,
    const double* __restrict__ vn_ie, const double* __restrict__ inv_dual,
    const double* __restrict__ w, const int* __restrict__ cell_idx,
    const double* __restrict__ z_vt_ie, const double* __restrict__ inv_primal,
    const double* __restrict__ tangent, const double* __restrict__ z_w_v,
    const int* __restrict__ vert_idx,
    unsigned N_e, unsigned N_c, unsigned N_v,
    unsigned nlev, unsigned nlev_end, unsigned nlev_padded)
{
    unsigned bid=(unsigned)blockIdx.y*gridDim.x+(unsigned)blockIdx.x;
    unsigned gid=bid*BX+threadIdx.x;
    unsigned jk=gid%nlev_padded, je=gid/nlev_padded;
    if(je>=N_e||jk>=nlev_end)return;
    int c2d=IC<V>((int)je,(int)jk,(int)N_e,(int)nlev);
    int ci0=cell_idx[IN<V>((int)je,0,(int)N_e)],ci1=cell_idx[IN<V>((int)je,1,(int)N_e)];
    int vi0=vert_idx[IN<V>((int)je,0,(int)N_e)],vi1=vert_idx[IN<V>((int)je,1,(int)N_e)];
    out[c2d]=vn_ie[c2d]*inv_dual[je]*
        (w[IC<V>(ci0,(int)jk,(int)N_c,(int)nlev)]-w[IC<V>(ci1,(int)jk,(int)N_c,(int)nlev)])
      +z_vt_ie[c2d]*inv_primal[je]*tangent[je]*
        (z_w_v[IC<V>(vi0,(int)jk,(int)N_v,(int)nlev)]-z_w_v[IC<V>(vi1,(int)jk,(int)N_v,(int)nlev)]);
}

/* ================================================================ */
/*  Verify                                                           */
/* ================================================================ */
template<int V>
static void cpu_reference(double* out,
    const double* vn_ie, const double* inv_dual, const double* w,
    const int* cell_idx, const double* z_vt_ie, const double* inv_primal,
    const double* tangent, const double* z_w_v, const int* vert_idx,
    int N_e, int N_c, int N_v, int nlev, int nlev_end)
{
    for(int jk=0;jk<nlev_end;jk++)
        for(int je=0;je<N_e;je++){STENCIL_BODY(V)}
}

static bool verify(const double* got, const double* ref, size_t n) {
    int nf=0; double mr=0;
    for(size_t i=0;i<n;i++){
        double d=std::abs(got[i]-ref[i]),dn=std::max(std::abs(ref[i]),1e-300);
        if(d/dn>mr) mr=d/dn;
        if(d>1e-12+1e-8*std::abs(ref[i])) nf++;
    }
    printf("  verify: %s (max_rel=%.2e, fails=%d/%zu)\n", nf?"FAIL":"OK", mr, nf, n);
    return nf==0;
}

/* ================================================================ */
/*  Named kernel launchers — one per profiling target                */
/* ================================================================ */
#define CLAMP65535(x) ((x)>65535u ? 65535u : (x))

/* K1: V1 + 1x4_96x1  (je-first, original layout) */
static void launch_K1(double* out, const double* vn, const double* id,
    const double* w, const int* ci, const double* vt, const double* ip,
    const double* tg, const double* zw, const int* vi,
    int Ne, int Nc, int Nv, int nl, int nle) {
    dim3 blk(96,1);
    dim3 grd(CLAMP65535(((unsigned)Ne+96*1-1)/(96*1)),
             CLAMP65535(((unsigned)nle+1*4-1)/(1*4)));
    gpu_kernel_je_first<1,4,96,1,1><<<grd,blk>>>(out,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle);
    CUDA_LAUNCH_CHECK();
}

/* K2: V4 + 1x4_96x1  (jk-first, same thread config) */
static void launch_K2(double* out, const double* vn, const double* id,
    const double* w, const int* ci, const double* vt, const double* ip,
    const double* tg, const double* zw, const int* vi,
    int Ne, int Nc, int Nv, int nl, int nle) {
    dim3 blk(96,1);
    dim3 grd(CLAMP65535(((unsigned)nle+96*1-1)/(96*1)),
             CLAMP65535(((unsigned)Ne+1*4-1)/(1*4)));
    gpu_kernel_jk_first<1,4,96,1,4><<<grd,blk>>>(out,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle);
    CUDA_LAUNCH_CHECK();
}

/* K3: V5 + 1x4_96x1  (jk-first, padded nlev) */
static void launch_K3(double* out, const double* vn, const double* id,
    const double* w, const int* ci, const double* vt, const double* ip,
    const double* tg, const double* zw, const int* vi,
    int Ne, int Nc, int Nv, int nl, int nle) {
    dim3 blk(96,1);
    dim3 grd(CLAMP65535(((unsigned)nle+96*1-1)/(96*1)),
             CLAMP65535(((unsigned)Ne+1*4-1)/(1*4)));
    gpu_kernel_jk_first<1,4,96,1,4><<<grd,blk>>>(out,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle);
    CUDA_LAUNCH_CHECK();
}

/* K4: V5 + 2x1_16x16  (best blue, old results) */
static void launch_K4(double* out, const double* vn, const double* id,
    const double* w, const int* ci, const double* vt, const double* ip,
    const double* tg, const double* zw, const int* vi,
    int Ne, int Nc, int Nv, int nl, int nle) {
    dim3 blk(16,16);
    dim3 grd(CLAMP65535(((unsigned)nle+16*2-1)/(16*2)),
             CLAMP65535(((unsigned)Ne+16*1-1)/(16*1)));
    gpu_kernel_jk_first<2,1,16,16,4><<<grd,blk>>>(out,vn,id,w,ci,vt,ip,tg,zw,vi,Ne,Nc,Nv,nl,nle);
    CUDA_LAUNCH_CHECK();
}


/* ================================================================ */
/*  Kernel descriptor                                                */
/* ================================================================ */
struct KernelDesc {
    const char* name;
    int V;              /* layout variant for set_variant */
    bool padded;        /* use padded nlev? */
    bool is_v6;         /* V6 flat kernel? */
};

static constexpr KernelDesc KERNELS[] = {
    {"K1_V1_1x4_96x1",     1, false, false},
    {"K2_V4_1x4_96x1",     4, false, false},
    {"K3_V5_1x4_96x1",     4, true,  false},
    {"K4_V5_2x1_16x16",    4, true,  false},
};
static constexpr int N_KERNELS = 4;
static constexpr int PROFILE_RUNS = 3;

/* ================================================================ */
/*  main                                                             */
/* ================================================================ */
int main(int argc, char* argv[]) {
    int icon_step = (argc>=2) ? atoi(argv[1]) : 9;

    std::string gp=icon_global_path(icon_step), pp=icon_patch_path(icon_step);
    int icon_nproma=icon_read_nproma(gp.c_str());
    IconEdgeData ied;
    bool ok=(icon_nproma>0)&&icon_load_patch(pp.c_str(),icon_nproma,ied);
    assert(ok);

    const int N_e=ied.n_edges_valid, N_c=ied.n_cells, N_v=ied.n_verts;
    const int nlev_base=90, nlev_padded=icon_pad_nlev(nlev_base);
    printf("N_e=%d  N_c=%d  N_v=%d  nlev=%d  nlev_padded=%d\n",N_e,N_c,N_v,nlev_base,nlev_padded);

    cudaDeviceProp prop; CUDA_CHECK(cudaGetDeviceProperties(&prop,0));
    printf("GPU: %s  SM=%d\n\n",prop.name,prop.multiProcessorCount);

    /* exact distribution */
    int *ecl=new int[N_e*2], *evl=new int[N_e*2];
    for(int i=0;i<N_e*2;i++){ecl[i]=ied.cell_idx[i]; evl[i]=ied.vert_idx[i];}

    /* allocate 2 BenchData: base and padded */
    BenchData bd_base, bd_pad;
    bd_base.alloc(N_e,N_c,N_v,nlev_base); bd_base.fill(nlev_base);
    bd_pad.alloc(N_e,N_c,N_v,nlev_padded); bd_pad.fill(nlev_padded);
    for(int je=0;je<N_e;je++){
        bd_base.inv_dual[je]=bd_pad.inv_dual[je]=ied.inv_dual[je];
        bd_base.inv_primal[je]=bd_pad.inv_primal[je]=ied.inv_primal[je];
        bd_base.tangent_o[je]=bd_pad.tangent_o[je]=ied.tangent_o[je];
    }

    /* device arrays — max of base and padded sizes */
    size_t max_e=std::max(bd_base.sz_e,bd_pad.sz_e);
    double *d_vn,*d_w,*d_vt,*d_zw,*d_out,*d_id,*d_ip,*d_tg; int *d_ci,*d_vi;
    CUDA_CHECK(cudaMalloc(&d_vn,max_e*8)); CUDA_CHECK(cudaMalloc(&d_w,bd_base.sz_c*8));
    CUDA_CHECK(cudaMalloc(&d_vt,max_e*8)); CUDA_CHECK(cudaMalloc(&d_zw,bd_base.sz_v*8));
    CUDA_CHECK(cudaMalloc(&d_out,max_e*8));
    CUDA_CHECK(cudaMalloc(&d_id,N_e*8)); CUDA_CHECK(cudaMalloc(&d_ip,N_e*8)); CUDA_CHECK(cudaMalloc(&d_tg,N_e*8));
    CUDA_CHECK(cudaMalloc(&d_ci,N_e*2*4)); CUDA_CHECK(cudaMalloc(&d_vi,N_e*2*4));

    double *h_ref_base=new double[bd_base.sz_e];
    double *h_ref_pad=new double[bd_pad.sz_e];
    double *h_gpu_out=new double[max_e];

    cudaEvent_t ev0,ev1;
    CUDA_CHECK(cudaEventCreate(&ev0)); CUDA_CHECK(cudaEventCreate(&ev1));

    for(int ki=0;ki<N_KERNELS;ki++){
        const KernelDesc& kd=KERNELS[ki];
        BenchData& bd = kd.padded ? bd_pad : bd_base;
        int nlev = kd.padded ? nlev_padded : nlev_base;
        int nlev_end = nlev_base;
        size_t sz_e = bd.sz_e;

        printf("━━━ %s (V=%d, nlev=%d, nlev_end=%d) ━━━\n", kd.name, kd.V, nlev, nlev_end);

        /* set layout */
        bd.set_variant(kd.V, ecl, evl);

        /* CPU reference */
        double* h_ref = kd.padded ? h_ref_pad : h_ref_base;
        memset(h_ref,0,sz_e*sizeof(double));
        if(kd.V<=2) cpu_reference<1>(h_ref,bd.h_vn_ie,bd.inv_dual,bd.h_w,bd.h_cidx,
            bd.h_z_vt_ie,bd.inv_primal,bd.tangent_o,bd.h_z_w_v,bd.h_vidx,N_e,N_c,N_v,nlev,nlev_end);
        else cpu_reference<4>(h_ref,bd.h_vn_ie,bd.inv_dual,bd.h_w,bd.h_cidx,
            bd.h_z_vt_ie,bd.inv_primal,bd.tangent_o,bd.h_z_w_v,bd.h_vidx,N_e,N_c,N_v,nlev,nlev_end);

        /* upload */
        CUDA_CHECK(cudaMemcpy(d_vn,bd.h_vn_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_w,bd.h_w,bd.sz_c*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_vt,bd.h_z_vt_ie,bd.sz_e*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_zw,bd.h_z_w_v,bd.sz_v*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_ci,bd.h_cidx,N_e*2*4,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_vi,bd.h_vidx,N_e*2*4,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_id,bd.inv_dual,N_e*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_ip,bd.inv_primal,N_e*8,cudaMemcpyHostToDevice));
        CUDA_CHECK(cudaMemcpy(d_tg,bd.tangent_o,N_e*8,cudaMemcpyHostToDevice));

        /* warmup + verify */
        CUDA_CHECK(cudaMemset(d_out,0,sz_e*8));
        switch(ki){
        case 0: launch_K1(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
        case 1: launch_K2(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
        case 2: launch_K3(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
        case 3: launch_K4(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
        }
        CUDA_CHECK(cudaDeviceSynchronize());
        CUDA_CHECK(cudaMemcpy(h_gpu_out,d_out,sz_e*8,cudaMemcpyDeviceToHost));
        verify(h_gpu_out,h_ref,sz_e);

        /* timed runs (also captured by ncu) */
        for(int r=0;r<PROFILE_RUNS;r++){
            CUDA_CHECK(cudaMemset(d_out,0,sz_e*8));
            CUDA_CHECK(cudaEventRecord(ev0));
            switch(ki){
            case 0: launch_K1(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
            case 1: launch_K2(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
            case 2: launch_K3(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
            case 3: launch_K4(d_out,d_vn,d_id,d_w,d_ci,d_vt,d_ip,d_tg,d_zw,d_vi,N_e,N_c,N_v,nlev,nlev_end); break;
            }
            CUDA_CHECK(cudaEventRecord(ev1));
            CUDA_CHECK(cudaEventSynchronize(ev1));
            float ms=0; CUDA_CHECK(cudaEventElapsedTime(&ms,ev0,ev1));
            printf("  run %d: %.3f ms\n",r,ms);
        }
        printf("\n");
    }

    CUDA_CHECK(cudaEventDestroy(ev0)); CUDA_CHECK(cudaEventDestroy(ev1));
    cudaFree(d_vn);cudaFree(d_w);cudaFree(d_vt);cudaFree(d_zw);cudaFree(d_out);
    cudaFree(d_id);cudaFree(d_ip);cudaFree(d_tg);cudaFree(d_ci);cudaFree(d_vi);
    delete[]h_ref_base;delete[]h_ref_pad;delete[]h_gpu_out;
    delete[]ecl;delete[]evl;
    bd_base.free_all(); bd_pad.free_all(); ied.free_all();
    printf("Done. Use: ncu --set full -o profile ./bench_gpu_profile\n");
    return 0;
}
