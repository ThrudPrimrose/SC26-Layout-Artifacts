#!/usr/bin/env python3
"""Matrix transpose GPU sweep: 7 kernel variants + library (cuTENSOR/hipTensor), with roofline."""
import subprocess, sys, os, csv
from collections import defaultdict
import numpy as np

AMD = os.environ.get("BEVERIN", "0") == "1"

AMD_FLAGS = ("-O3 -ffast-math -fPIC -Wall -Wextra "
             "-Wno-unused-parameter -munsafe-fp-atomics "
             "-ffp-contract=fast -Wno-ignored-attributes -Wno-unused-result -std=c++17")

# ── Spack CUDA env ──
def setup_cuda_env(spec="cuda@12.9"):
    if AMD is False:
        try:
            prefix = subprocess.check_output(["spack","location","-i",spec], text=True).strip()
            ld = os.environ.get("LD_LIBRARY_PATH", "")
            os.environ["LD_LIBRARY_PATH"] = f"{prefix}/lib64:{ld}"
            print(f"  CUDA: {prefix}")
        except (FileNotFoundError, subprocess.CalledProcessError):
            print("  [WARN] spack not found, using system CUDA")

setup_cuda_env()

BINARY      = "./transpose_gpu"
BINARY_LIB  = "./transpose_lib"
CSV_RAW     = "transpose_raw.csv"
CSV_AGG     = "transpose_results.csv"
N           = 8192*2
REPS        = 100
WARMUP      = 5

V_NAMES = {0:"naive", 1:"blocked", 2:"smem", 3:"smem_blk", 4:"smem_pad", 5:"smem_swiz", 6:"blk_swiz"}
LIB_NAMES = {"cutensor", "cutensor_blk", "hiptensor", "hiptensor_blk"}

CONFIGS = [
    (64, 8, 1, 1), (32, 8, 1, 1), (64, 2, 1, 1),
    (64, 8, 2, 1), (32, 8, 2, 1),
    (64, 8, 1, 2), (32, 8, 1, 2),
    (32, 4, 2, 2), (64, 2, 2, 2),
    (64, 2, 4, 1), (32, 4, 4, 1),
]
VARIANTS  = [0, 1, 2, 3, 4, 5, 6]
SB_VALS   = [8, 64] if AMD else [16, 32]
PAD_VALS  = [1]
LIB_SB_VALS = [8, 32, 64] if AMD else [16, 32]  # block sizes for library blocked variant

# ── Roofline ──
def get_amd_roofline():
    name, peak_fp64, peak_bw = "MI300A", 61300.0, 5300.0
    import tempfile
    src = r'''
#include <hip/hip_runtime.h>
#include <cstdio>
__global__ void kc(const double*s,double*d,long long n){
    long long i=blockIdx.x*(long long)blockDim.x+threadIdx.x;if(i<n)d[i]=s[i];}
int main(){
    long long SZ=128LL*1024*1024; size_t bytes=SZ*8;
    double *a,*b; hipMalloc(&a,bytes); hipMalloc(&b,bytes);
    hipEvent_t s,e; hipEventCreate(&s); hipEventCreate(&e);
    kc<<<(SZ+255)/256,256>>>(a,b,SZ); hipDeviceSynchronize();
    hipEventRecord(s);
    for(int i=0;i<20;i++) kc<<<(SZ+255)/256,256>>>(a,b,SZ);
    hipEventRecord(e); hipEventSynchronize(e);
    float ms; hipEventElapsedTime(&ms,s,e);
    printf("%.1f\n",20.0*2*SZ*8/(ms/1000.0)/1e9);
    hipFree(a); hipFree(b);
}'''
    with tempfile.NamedTemporaryFile(suffix=".cpp", mode="w", delete=False) as f:
        f.write(src); tmp = f.name
    bw_bin = tmp.replace(".cpp", "")
    subprocess.run(f"hipcc {AMD_FLAGS} -o {bw_bin} {tmp}", shell=True, check=True)
    r = subprocess.run(bw_bin, capture_output=True, text=True, check=True)
    emp_bw = float(r.stdout.strip())
    os.remove(tmp); os.remove(bw_bin)
    return dict(name=name, peak_fp64=peak_fp64, peak_bw=peak_bw, emp_bw=emp_bw)

def get_cuda_roofline():
    # GH200 Grace Hopper (HBM3, 96 GB) — NVIDIA DA-11356-002_10
    # Peak BW = 2 * 3143 MHz * (5120-bit / 8) / 1000 = 4023 GB/s
    name, peak_fp64, peak_bw = "GH200", 34000.0, 4023.0

    import tempfile
    src = r'''
#include <cuda_runtime.h>
#include <cstdio>
__global__ void kc(const double*s,double*d,long long n){
    long long i=blockIdx.x*(long long)blockDim.x+threadIdx.x;if(i<n)d[i]=s[i];}
int main(){
    long long SZ=128LL*1024*1024; size_t bytes=SZ*8;
    double *a,*b; cudaMalloc(&a,bytes); cudaMalloc(&b,bytes);
    cudaEvent_t s,e; cudaEventCreate(&s); cudaEventCreate(&e);
    kc<<<(SZ+255)/256,256>>>(a,b,SZ); cudaDeviceSynchronize();
    cudaEventRecord(s);
    for(int i=0;i<20;i++) kc<<<(SZ+255)/256,256>>>(a,b,SZ);
    cudaEventRecord(e); cudaEventSynchronize(e);
    float ms; cudaEventElapsedTime(&ms,s,e);
    printf("%.1f\n",20.0*2*SZ*8/(ms/1000.0)/1e9);
    cudaFree(a); cudaFree(b);
}'''
    with tempfile.NamedTemporaryFile(suffix=".cu", mode="w", delete=False) as f:
        f.write(src); tmp = f.name
    bw_bin = tmp.replace(".cu", "")
    subprocess.run(f"nvcc -O3 -arch=sm_90 -o {bw_bin} {tmp}", shell=True, check=True)
    r = subprocess.run(bw_bin, capture_output=True, text=True, check=True)
    emp_bw = float(r.stdout.strip())
    os.remove(tmp); os.remove(bw_bin)
    return dict(name=name, peak_fp64=peak_fp64, peak_bw=peak_bw, emp_bw=emp_bw)

def get_roofline():
    return get_amd_roofline() if AMD else get_cuda_roofline()

# ── Compile ──
def compile():
    if not AMD:
        cmd = f"nvcc -O3 -std=c++17 -arch=native -o {BINARY} transpose_gpu.cu"
    else:
        cmd = f"hipcc {AMD_FLAGS} -o {BINARY} transpose_gpu_hip.cpp"
    print(f"Compiling kernels: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def compile_lib():
    """Compile cuTENSOR/hipTensor benchmark. Returns True on success."""
    if not AMD:
        cmd = f"nvcc -O3 -std=c++17 -arch=native -o {BINARY_LIB} transpose_cutensor.cu -lcutensor"
        lib = "cuTENSOR"
    else:
        cmd = f"hipcc {AMD_FLAGS} -o {BINARY_LIB} transpose_gpu_hip.cpp -lhiptensor"
        lib = "hipTensor"
    print(f"Compiling {lib}: {cmd}")
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [WARN] {lib} compile failed: {r.stderr.strip()[:120]}")
        print(f"  Library benchmark will be skipped.")
        return False
    return True

# ── Sweep ──
def sweep():
    open(CSV_RAW, "w").close()
    for bx,by,tx,ty in CONFIGS:
        for v in VARIANTS:
            sb_list = SB_VALS if v in (1,3,6) else [32]
            pad_list = PAD_VALS if v == 4 else [1]
            for sb in sb_list:
                for pad in pad_list:
                    if v in (1,3,6) and N % sb != 0: continue
                    cmd = [BINARY, str(N), str(v), str(bx), str(by), str(tx), str(ty),
                           CSV_RAW, str(sb), str(pad), str(WARMUP), str(REPS)]
                    try:
                        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                        if r.returncode != 0:
                            print(f"  FAIL {V_NAMES[v]} {bx}x{by} TX={tx} TY={ty}: {r.stderr.strip()[:80]}")
                            continue
                        print(f"  {r.stdout.strip()}")
                    except subprocess.TimeoutExpired:
                        print(f"  TIMEOUT {V_NAMES[v]} {bx}x{by}")

def sweep_lib():
    """Run library benchmark: row-major (var=0) + blocked (var=1) for each SB."""
    if not os.path.exists(BINARY_LIB):
        return
    lib = "hipTensor" if AMD else "cuTENSOR"
    print(f"\n  ── Library: {lib} ──")
    # Variant 0: row-major transpose
    cmd = [BINARY_LIB, str(N), "0", CSV_RAW, "32", str(WARMUP), str(REPS)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0:
            print(f"  {r.stdout.strip()}")
        else:
            print(f"  FAIL {lib} row-major: {r.stderr.strip()[:100]}")
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {lib} row-major")

    # Variant 1: blocked transpose for each SB
    for sb in LIB_SB_VALS:
        if N % sb != 0: continue
        cmd = [BINARY_LIB, str(N), "1", CSV_RAW, str(sb), str(WARMUP), str(REPS)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0:
                print(f"  {r.stdout.strip()}")
            else:
                print(f"  FAIL {lib} blocked SB={sb}: {r.stderr.strip()[:100]}")
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT {lib} blocked SB={sb}")

# ── Aggregate ──
def aggregate(roof):
    groups = defaultdict(list)
    with open(CSV_RAW) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 12: continue
            key = tuple(p[:8])
            groups[key].append((float(p[9]), float(p[10])))

    rows = []
    for key, iters in groups.items():
        gbs = np.array([g for _,g in iters])
        ts  = np.array([t for t,_ in iters])
        med_gbs = np.median(gbs)
        rows.append(dict(
            variant=key[0], N=key[1], BX=key[2], BY=key[3], TX=key[4], TY=key[5],
            SB=key[6], PAD=key[7],
            med_time_ms=f"{np.median(ts)*1000:.4f}", std_time_ms=f"{np.std(ts)*1000:.4f}",
            min_gbs=f"{np.min(gbs):.1f}", med_gbs=f"{med_gbs:.1f}", max_gbs=f"{np.max(gbs):.1f}",
            pct_peak_bw=f"{100*med_gbs/roof['emp_bw']:.1f}",
        ))

    fields = ["variant","N","BX","BY","TX","TY","SB","PAD",
              "med_time_ms","std_time_ms","min_gbs","med_gbs","max_gbs","pct_peak_bw"]
    with open(CSV_AGG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    return rows

# ── Report ──
def report(roof, rows):
    print(f"\n{'='*88}")
    print(f" GPU: {roof['name']}   Emp BW: {roof['emp_bw']:.0f} GB/s   Theo BW: {roof['peak_bw']:.0f} GB/s")
    print(f" Transpose: N={N}  bytes/call = {2*N*N*8/1e9:.2f} GB (read+write)")
    print(f"{'='*88}")

    lib_rows  = [r for r in rows if r["variant"] in LIB_NAMES]
    kern_rows = [r for r in rows if r["variant"] not in LIB_NAMES]

    # Find best library row-major as reference
    lib_rm = [r for r in lib_rows if "blk" not in r["variant"]]
    lib_bk = [r for r in lib_rows if "blk" in r["variant"]]
    ref_gbs = float(lib_rm[0]["med_gbs"]) if lib_rm else None

    if lib_rm:
        lr = lib_rm[0]
        print(f"\n >>> Library (row-major): {lr['variant']}  {float(lr['med_gbs']):.1f} GB/s ({lr['pct_peak_bw']}% peak) <<<")
    if lib_bk:
        for lr in sorted(lib_bk, key=lambda x: -float(x["med_gbs"])):
            print(f" >>> Library (blocked SB={lr['SB']}): {lr['variant']}  {float(lr['med_gbs']):.1f} GB/s ({lr['pct_peak_bw']}% peak) <<<")

    print(f"\n {'var':<14} {'BX':>3} {'BY':>3} {'TX':>2} {'TY':>2} {'SB':>3} {'P':>1}"
          f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8} {'%pkBW':>6}  {'vs lib':>6}")
    print(f" {'-'*76}")

    for r in sorted(rows, key=lambda x: -float(x["med_gbs"])):
        med = float(r["med_gbs"])
        vs = f"{100*med/ref_gbs:.0f}%" if ref_gbs else "—"
        star = " ★" if r["variant"] in LIB_NAMES else ""
        print(f" {r['variant']:<14} {r['BX']:>3} {r['BY']:>3} {r['TX']:>2} {r['TY']:>2} {r['SB']:>3} {r['PAD']:>1}"
              f"  {med:8.1f} {float(r['min_gbs']):8.1f}"
              f" {float(r['max_gbs']):8.1f}  {r['pct_peak_bw']:>5}%  {vs:>5}{star}")

    if kern_rows:
        best = max(kern_rows, key=lambda x: float(x["med_gbs"]))
        bk = float(best["med_gbs"])
        extra = f" ({100*bk/ref_gbs:.0f}% of library)" if ref_gbs else ""
        print(f"\n Best kernel: {best['variant']} BX={best['BX']} BY={best['BY']}"
              f" TX={best['TX']} TY={best['TY']} SB={best['SB']}"
              f" -> {best['med_gbs']} GB/s ({best['pct_peak_bw']}% peak){extra}")

# ── Main ──
if __name__ == "__main__":
    print("── Roofline ──")
    roof = get_roofline()
    print(f"  Emp BW: {roof['emp_bw']:.0f} GB/s")

    need_compile = not os.path.exists(BINARY) or "--compile" in sys.argv
    if need_compile:
        compile()

    has_lib = False
    if need_compile or not os.path.exists(BINARY_LIB):
        has_lib = compile_lib()
    else:
        has_lib = os.path.exists(BINARY_LIB)

    print(f"\n── Sweep: N={N} reps={REPS} warmup={WARMUP} ──")
    #sweep()
    #if has_lib:
    #    sweep_lib()

    rows = aggregate(roof)
    if rows:
        report(roof, rows)
    print(f"\nPer-iteration: {CSV_RAW}\nAggregated:    {CSV_AGG}")