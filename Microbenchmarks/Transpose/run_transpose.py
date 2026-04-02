#!/usr/bin/env python3
"""Matrix transpose GPU sweep: 7 kernel variants + library (cuTENSOR/hipTensor), with roofline."""
import subprocess, sys, os, csv
from collections import defaultdict
import numpy as np
from pathlib import Path

AMD = os.environ.get("BEVERIN", "0") == "1"

AMD_FLAGS = ("-O3 -ffast-math -fPIC -Wall -Wextra  -DHIP_PLATFORM_AMD=1 -D__HIP_PLATFORM_AMD__=1  "
             "-Wno-unused-parameter -munsafe-fp-atomics  --offload-arch=gfx942 "
             "-ffp-contract=fast -Wno-ignored-attributes -Wno-unused-result -std=c++17")

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".cache")
os.makedirs(CACHE_DIR, exist_ok=True)


BINARY      = "./transpose_gpu"
BINARY_LIB  = "./transpose_lib"
CSV_RAW     = "transpose_raw.csv"
CSV_AGG     = "transpose_results.csv"
N           = 16384
REPS        = 100
WARMUP      = 5

V_NAMES = {0:"naive", 1:"blocked", 2:"smem", 3:"smem_blk", 4:"smem_pad", 5:"smem_swiz", 6:"smem_blk_swiz", 7:"smem_pad_blk"}
LIB_NAMES = {"cutensor", "cutensor_blk", "hiptensor", "hiptensor_blk"}

CONFIGS = [
    # existing ...

    # Higher TY coarsening (TY=4 was best, try TY=8)
    (32,  8, 1, 8), (32, 16, 1, 8), (32,  8, 2, 4), (32, 16, 2, 4),

    # TX=4 with BY=8,16 (TX=4 mostly tested with small BY)
    (32,  8, 4, 1), (32, 16, 4, 1), (32,  8, 4, 2), (32, 16, 4, 2),

    (32,  8, 2, 2), (32,  8, 4, 4),

    # BY=32 (untested entirely)
    (32, 32, 1, 1), (32, 32, 1, 2), (32, 32, 2, 1), (32, 32, 1, 4), (32, 32, 2, 2),

    # BX=16 (untested — smaller blocks, more occupancy)
    (16,  8, 1, 1), (16, 16, 1, 1), (16,  8, 1, 2), (16, 16, 1, 2),
    (16,  8, 2, 1), (16, 16, 2, 2), (16,  8, 1, 4), (16, 16, 1, 4),

    # BX=64 with higher TY (mostly tested with TY=1,2)
    (64,  8, 1, 4), (64, 16, 1, 2), (64, 16, 1, 4), (64, 16, 2, 1),
    (64, 16, 2, 2), 
]
VARIANTS    = [0, 1, 2, 3, 4, 5, 6, 7]
SB_VALS     = [32, 64, 128, 256] if AMD else [8, 32, 64, 128, 256]
PAD_VALS    = [1]
LIB_SB_VALS = [32, 64, 128, 256] if AMD else [8, 32, 64, 128]

# ── Roofline ──
AMD_ROOFLINE_SRC = r'''
#include <hip/hip_runtime.h>
#include <cstdio>
#include <cstdlib>
#define CK(x) do { hipError_t _e=(x); if(_e){ \
    fprintf(stderr,"HIP error %s:%d: %s\n",__FILE__,__LINE__,hipGetErrorString(_e)); \
    exit(1); } } while(0)
__global__ void kc(const double* s, double* d, long long n) {
    long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i < n) d[i] = s[i];
}
int main() {
    const long long SZ = 128LL * 1024 * 1024;
    const size_t bytes = SZ * sizeof(double);
    double *a, *b;
    CK(hipMalloc(&a, bytes)); CK(hipMalloc(&b, bytes));
    hipEvent_t s, e; CK(hipEventCreate(&s)); CK(hipEventCreate(&e));
    const int threads = 256;
    const int blocks  = (int)((SZ + threads - 1) / threads);
    kc<<<blocks, threads>>>(a, b, SZ);
    CK(hipGetLastError()); CK(hipDeviceSynchronize());
    CK(hipEventRecord(s));
    for (int i = 0; i < 20; i++) { kc<<<blocks, threads>>>(a, b, SZ); CK(hipGetLastError()); }
    CK(hipEventRecord(e)); CK(hipEventSynchronize(e));
    float ms; CK(hipEventElapsedTime(&ms, s, e));
    if (ms < 1e-3f) { fprintf(stderr, "timing error: ms=%.6f\n", ms); return 1; }
    printf("%.1f\n", 20.0 * 2 * SZ * sizeof(double) / (ms / 1000.0) / 1e9);
    CK(hipFree(a)); CK(hipFree(b));
    return 0;
}'''

CUDA_ROOFLINE_SRC = r'''
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#define CK(x) do { cudaError_t _e=(x); if(_e){ \
    fprintf(stderr,"CUDA error %s:%d: %s\n",__FILE__,__LINE__,cudaGetErrorString(_e)); \
    exit(1); } } while(0)
__global__ void kc(const double* s, double* d, long long n) {
    long long i = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (i < n) d[i] = s[i];
}
int main() {
    const long long SZ = 128LL * 1024 * 1024;
    const size_t bytes = SZ * sizeof(double);
    double *a, *b;
    CK(cudaMalloc(&a, bytes)); CK(cudaMalloc(&b, bytes));
    cudaEvent_t s, e; CK(cudaEventCreate(&s)); CK(cudaEventCreate(&e));
    const int threads = 256;
    const int blocks  = (int)((SZ + threads - 1) / threads);
    kc<<<blocks, threads>>>(a, b, SZ);
    CK(cudaGetLastError()); CK(cudaDeviceSynchronize());
    CK(cudaEventRecord(s));
    for (int i = 0; i < 20; i++) { kc<<<blocks, threads>>>(a, b, SZ); CK(cudaGetLastError()); }
    CK(cudaEventRecord(e)); CK(cudaEventSynchronize(e));
    float ms; CK(cudaEventElapsedTime(&ms, s, e));
    if (ms < 1e-3f) { fprintf(stderr, "timing error: ms=%.6f\n", ms); return 1; }
    printf("%.1f\n", 20.0 * 2 * SZ * sizeof(double) / (ms / 1000.0) / 1e9);
    CK(cudaFree(a)); CK(cudaFree(b));
    return 0;
}'''

def get_roofline():
    if AMD:
        name, peak_fp64, peak_bw = "MI300A", 61300.0, 5300.0
        suffix, src = ".cpp", AMD_ROOFLINE_SRC
        compile_cmd = lambda obj, src: f"hipcc {AMD_FLAGS} -o {obj} {src}"
    else:
        name, peak_fp64, peak_bw = "GH200", 34000.0, 4023.0
        suffix, src = ".cu", CUDA_ROOFLINE_SRC
        compile_cmd = lambda obj, src: f"nvcc -O3 -arch=sm_90 -o {obj} {src}"

    src_path = os.path.join(CACHE_DIR, f"roofline{suffix}")
    bw_bin   = os.path.join(CACHE_DIR, "roofline")

    with open(src_path, "w") as f:
        f.write(src)

    comp = subprocess.run(compile_cmd(bw_bin, src_path), shell=True,
                          capture_output=True, text=True)
    if comp.returncode != 0:
        print(f"  [ERROR] Roofline compile failed:\n{comp.stderr.strip()}")
        sys.exit(1)

    r = subprocess.run(bw_bin, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  [ERROR] Roofline binary failed:\n{r.stderr.strip()}")
        sys.exit(1)

    emp_bw = float(r.stdout.strip())
    return dict(name=name, peak_fp64=peak_fp64, peak_bw=peak_bw, emp_bw=emp_bw)

# ── Compile ──
def compile():
    if not AMD:
        cmd = f"nvcc -O3 -std=c++17 -arch=native -o {BINARY} transpose_gpu.cu"
    else:
        cmd = f"hipcc {AMD_FLAGS} -o {BINARY} transpose_gpu_hip.cpp"
    print(f"Compiling kernels: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def compile_lib():
    # Skip if already compiled
    if Path(BINARY_LIB).exists():
        print(f"[INFO] {BINARY_LIB} already exists, skipping compilation.")
        return True

    if not AMD:
        cmd = f"nvcc -O3 -std=c++17 -arch=native -o {BINARY_LIB} transpose_cutensor.cu -lcutensor"
        lib = "cuTENSOR"
    else:
        cmd = (f"hipcc {AMD_FLAGS} "
               f"-o {BINARY_LIB} transpose_hiptensor.cpp -lhiptensor")
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
    for bx, by, tx, ty in CONFIGS:
        for v in VARIANTS:
            sb_list  = SB_VALS  if v in (1, 3, 6) else [32]
            pad_list = PAD_VALS if v == 4          else [1]
            for sb in sb_list:
                for pad in pad_list:
                    if v in (1, 3, 6) and N % sb != 0:
                        continue
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
    if not os.path.exists(BINARY_LIB):
        return
    lib = "hipTensor" if AMD else "cuTENSOR"
    print(f"\n  ── Library: {lib} ──")
    # variant 0: row-major
    cmd = [BINARY_LIB, str(N), "0", CSV_RAW, "32", str(WARMUP), str(REPS)]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0: print(f"  {r.stdout.strip()}")
        else: print(f"  FAIL {lib} row-major: {r.stderr.strip()[:100]}")
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT {lib} row-major")
    # variant 1: blocked for each SB
    for sb in LIB_SB_VALS:
        if N % sb != 0:
            continue
        cmd = [BINARY_LIB, str(N), "1", CSV_RAW, str(sb), str(WARMUP), str(REPS)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            if r.returncode == 0: print(f"  {r.stdout.strip()}")
            else: print(f"  FAIL {lib} blocked SB={sb}: {r.stderr.strip()[:100]}")
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT {lib} blocked SB={sb}")

# ── Aggregate ──
def aggregate(roof):
    if not os.path.exists(CSV_RAW):
        print(f"  [ERROR] {CSV_RAW} not found — sweep produced no output.")
        sys.exit(1)

    groups = defaultdict(list)
    with open(CSV_RAW) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 12:
                continue
            key = tuple(p[:8])
            groups[key].append((float(p[9]), float(p[10])))

    if not groups:
        print(f"  [ERROR] {CSV_RAW} is empty — no results to aggregate.")
        sys.exit(1)

    rows = []
    for key, iters in groups.items():
        gbs = np.array([g for _, g in iters])
        ts  = np.array([t for t, _ in iters])
        med_gbs = np.median(gbs)
        rows.append(dict(
            variant=key[0], N=key[1], BX=key[2], BY=key[3], TX=key[4], TY=key[5],
            SB=key[6], PAD=key[7],
            med_time_ms=f"{np.median(ts)*1000:.4f}",
            std_time_ms=f"{np.std(ts)*1000:.4f}",
            min_gbs=f"{np.min(gbs):.1f}",
            med_gbs=f"{med_gbs:.1f}",
            max_gbs=f"{np.max(gbs):.1f}",
            pct_peak_bw=f"{100*med_gbs/roof['emp_bw']:.1f}",
        ))

    fields = ["variant","N","BX","BY","TX","TY","SB","PAD",
              "med_time_ms","std_time_ms","min_gbs","med_gbs","max_gbs","pct_peak_bw"]
    with open(CSV_AGG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    return rows

# ── Report ──
def report(roof, rows):
    print(f"\n{'='*96}")
    print(f" GPU: {roof['name']}   Emp BW: {roof['emp_bw']:.0f} GB/s   Theo BW: {roof['peak_bw']:.0f} GB/s")
    print(f" Transpose: N={N}  bytes/call = {2*N*N*8/1e9:.2f} GB (read+write)")
    print(f"{'='*96}")

    kern_rows = [r for r in rows if r["variant"] not in LIB_NAMES]
    lib_rows  = [r for r in rows if r["variant"] in LIB_NAMES]

    # Best library row-major as speedup reference
    lib_rm  = [r for r in lib_rows if "blk" not in r["variant"]]
    ref_gbs = float(max(lib_rm, key=lambda x: float(x["med_gbs"]))["med_gbs"]) if lib_rm else None

    # Header
    print(f"\n {'var':<18} {'BX':>3} {'BY':>3} {'TX':>2} {'TY':>2} {'SB':>4} {'P':>1}"
          f"  {'medGB/s':>8} {'minGB/s':>8} {'maxGB/s':>8} {'%pkBW':>6}  {'vs lib':>7}  {'type'}")
    print(f" {'-'*84}")

    # All rows sorted by med_gbs descending — kernels and library interleaved
    for r in sorted(rows, key=lambda x: -float(x["med_gbs"])):
        med    = float(r["med_gbs"])
        vs     = f"{100*med/ref_gbs:.0f}%" if ref_gbs else "—"
        is_lib = r["variant"] in LIB_NAMES
        tag    = " ★ lib" if is_lib else "      "
        print(f" {r['variant']:<18} {r['BX']:>3} {r['BY']:>3} {r['TX']:>2} {r['TY']:>2}"
              f" {r['SB']:>4} {r['PAD']:>1}"
              f"  {med:8.1f} {float(r['min_gbs']):8.1f}"
              f" {float(r['max_gbs']):8.1f}  {r['pct_peak_bw']:>5}%  {vs:>6}  {tag}")

    # Summary lines
    if lib_rows:
        best_lib = max(lib_rows, key=lambda x: float(x["med_gbs"]))
        print(f"\n Best library : {best_lib['variant']} SB={best_lib['SB']}"
              f" -> {best_lib['med_gbs']} GB/s ({best_lib['pct_peak_bw']}% peak)")
    if kern_rows:
        best_kern = max(kern_rows, key=lambda x: float(x["med_gbs"]))
        bk        = float(best_kern["med_gbs"])
        extra     = f"  ({100*bk/ref_gbs:.0f}% of best library)" if ref_gbs else ""
        print(f" Best kernel  : {best_kern['variant']} BX={best_kern['BX']} BY={best_kern['BY']}"
              f" TX={best_kern['TX']} TY={best_kern['TY']} SB={best_kern['SB']}"
              f" -> {best_kern['med_gbs']} GB/s ({best_kern['pct_peak_bw']}% peak){extra}")

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
    sweep()
    if has_lib:
        sweep_lib()

    rows = aggregate(roof)
    if rows:
        report(roof, rows)
    print(f"\nPer-iteration: {CSV_RAW}\nAggregated:    {CSV_AGG}")