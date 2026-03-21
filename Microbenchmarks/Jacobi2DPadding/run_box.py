#!/usr/bin/env python3
"""Jacobi2D GPU sweep: 5pt/box3/box4 stencils × 4 variants × configs, with roofline."""
import subprocess, sys, os, csv
from collections import defaultdict
import numpy as np

import subprocess
# export LD_LIBRARY_PATH=$(spack location -i cuda@12.9)/lib64:$LD_LIBRARY_PATH

# get spack install location for cuda@12.9
cuda_path = subprocess.check_output(
    ["spack", "location", "-i", "cuda@12.9"],
    text=True
).strip()

# prepend to LD_LIBRARY_PATH
lib_path = f"{cuda_path}/lib64"
os.environ["LD_LIBRARY_PATH"] = f"{lib_path}:{os.environ.get('LD_LIBRARY_PATH', '')}"

# (optional) print to verify
print(os.environ["LD_LIBRARY_PATH"])

# ── Spack CUDA env ──
def setup_cuda_env(spec="cuda@12.9"):
    try:
        prefix = subprocess.check_output(["spack","location","-i",spec], text=True).strip()
        ld = os.environ.get("LD_LIBRARY_PATH", "")
        os.environ["LD_LIBRARY_PATH"] = f"{prefix}/lib64:{ld}"
        print(f"  CUDA prefix: {prefix}")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"  [WARN] spack lookup failed ({e}), using system CUDA")

setup_cuda_env()

# ── Config ──
BINARY   = "./box_gpu"
CSV_RAW  = "box_raw.csv"
CSV_AGG  = "box_results.csv"
N        = 4096*4
TSTEPS   = 100

STENCILS = [0, 1, 2]  # 5pt, box3, box4
ST_NAMES = {0: '5pt', 1: "box3", 2: "box4"}
ST_HALO  = {0: 1, 1: 1, 2: 2}
ST_FLOPS = {0: 5, 1: 9, 2: 16}

CONFIGS = [
    (32, 8, 1, 1), (16, 16, 1, 1), (32, 4, 1, 1),
    (32, 8, 2, 1), (16, 16, 2, 1), (32, 4, 2, 1),
    (32, 8, 1, 2), (16, 16, 1, 2),
    (16, 8, 2, 2), (32, 4, 2, 2),
    (32, 4, 4, 1), (16, 8, 4, 1), (32, 4, 1, 4),
    (16, 4, 4, 2), (8, 8, 4, 4),
]
VARIANTS = [0, 1, 2, 3]
PADS     = [1, 2]

RAW_COLS = "stencil,variant,N,tsteps,BX,BY,TX,TY,PAD,iter,time_s,gflops,checksum"

# ── Roofline ──
def get_roofline():
    try:
        import pycuda.driver as cuda
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pycuda", "-q"])
        import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule

    dev = pycuda.autoinit.device
    a = dev.get_attributes()
    sm  = a[cuda.device_attribute.MULTIPROCESSOR_COUNT]
    clk = a[cuda.device_attribute.CLOCK_RATE] / 1000
    mem = a[cuda.device_attribute.MEMORY_CLOCK_RATE] / 1000
    bw  = a[cuda.device_attribute.GLOBAL_MEMORY_BUS_WIDTH]
    cc  = (a[cuda.device_attribute.COMPUTE_CAPABILITY_MAJOR],
           a[cuda.device_attribute.COMPUTE_CAPABILITY_MINOR])

    FP64 = {(7,0):32,(7,5):2,(8,0):32,(8,6):2,(8,7):2,(8,9):2,(9,0):64}
    u = FP64.get(cc, 32 if cc[1]==0 and cc[0]>=7 else 2)
    peak_fp64 = sm * u * 2 * (clk / 1000)
    peak_bw   = 2 * mem * (bw / 8) / 1000

    SZ = 128 * 1024 * 1024
    ag = cuda.mem_alloc(SZ*8); bg = cuda.mem_alloc(SZ*8)
    mod = SourceModule("""
    __global__ void kc(const double*s,double*d,long long n){long long i=blockIdx.x*(long long)blockDim.x+threadIdx.x;if(i<n)d[i]=s[i];}
    __global__ void kd(double a,const double*x,double*y,long long n){long long i=blockIdx.x*(long long)blockDim.x+threadIdx.x;if(i<n)y[i]=a*x[i]+y[i];}
    """)
    kc=mod.get_function("kc"); kd=mod.get_function("kd")
    bl,gr=256,(SZ+255)//256
    def bench(fn, reps=20):
        fn(); cuda.Context.synchronize()
        s,e=cuda.Event(),cuda.Event(); s.record()
        for _ in range(reps): fn()
        e.record(); e.synchronize()
        return s.time_till(e)/reps/1000
    cbw=2*SZ*8/bench(lambda:kc(ag,bg,np.int64(SZ),block=(bl,1,1),grid=(gr,1)))/1e9
    dbw=3*SZ*8/bench(lambda:kd(np.float64(2.),ag,bg,np.int64(SZ),block=(bl,1,1),grid=(gr,1)))/1e9
    ag.free(); bg.free()
    emp=max(cbw,dbw)

    # Per-stencil OI and attainable
    oi, att = {}, {}
    for st in STENCILS:
        H = ST_HALO[st]; S = N + 2*H
        oi[st]  = (float(ST_FLOPS[st]) * N * N) / (2.0 * S * S * 8)
        att[st] = min(peak_fp64, oi[st] * emp)

    return dict(name=dev.name(), cc=cc, sms=sm, peak_fp64=peak_fp64,
                peak_bw=peak_bw, emp_bw=emp, copy_bw=cbw, daxpy_bw=dbw,
                oi=oi, attainable=att)

def compile():
    cmd = f"nvcc -O3 -arch=native -std=c++17 -o {BINARY} box_gpu.cu"
    print(f"Compiling: {cmd}"); subprocess.run(cmd, shell=True, check=True)

def sweep():
    open(CSV_RAW, "w").close()
    for st in STENCILS:
        print(f"\n  ── Stencil: {ST_NAMES[st]} (halo={ST_HALO[st]}, {ST_FLOPS[st]} FLOPs/pt) ──")
        for bx,by,tx,ty in CONFIGS:
            for v in VARIANTS:
                for pad in (PADS if v==2 else [0]):
                    cmd = [BINARY, str(N), str(TSTEPS), str(st), str(v),
                           str(bx), str(by), str(tx), str(ty), str(pad), CSV_RAW]
                    try:
                        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                        if r.returncode != 0:
                            print(f"    FAIL {ST_NAMES[st]} v={v} {bx}x{by}: {r.stderr.strip()[:80]}")
                            continue
                        print(f"    {r.stdout.strip()}")
                    except subprocess.TimeoutExpired:
                        print(f"    TIMEOUT {ST_NAMES[st]} v={v} {bx}x{by}")

def aggregate(roof):
    groups = defaultdict(list)
    with open(CSV_RAW) as f:
        for line in f:
            p = line.strip().split(",")
            if len(p) < 13: continue
            key = tuple(p[:9])  # stencil,variant,N,tsteps,BX,BY,TX,TY,PAD
            groups[key].append((float(p[10]), float(p[11])))

    rows = []
    for key, iters in groups.items():
        st_name = key[0]
        st_id = {v:k for k,v in ST_NAMES.items()}[st_name]
        H = ST_HALO[st_id]; S = N + 2*H
        bpi = 2.0 * S * S * 8  # bytes per iteration

        times = np.array([t for t,_ in iters])
        gfs   = np.array([g for _,g in iters])
        med_t = np.median(times); med_gf = np.median(gfs)
        ebw = bpi / med_t / 1e9

        rows.append(dict(
            stencil=key[0], variant=key[1], N=key[2], tsteps=key[3],
            BX=key[4], BY=key[5], TX=key[6], TY=key[7], PAD=key[8],
            med_time_s=f"{med_t:.6f}", std_time_s=f"{np.std(times):.6f}",
            min_gflops=f"{np.min(gfs):.2f}", med_gflops=f"{med_gf:.2f}",
            max_gflops=f"{np.max(gfs):.2f}", eff_bw_gbs=f"{ebw:.1f}",
            pct_attainable=f"{100*med_gf/roof['attainable'][st_id]:.1f}",
            pct_peak_bw=f"{100*ebw/roof['emp_bw']:.1f}",
        ))

    fields = ["stencil","variant","N","tsteps","BX","BY","TX","TY","PAD",
              "med_time_s","std_time_s","min_gflops","med_gflops","max_gflops",
              "eff_bw_gbs","pct_attainable","pct_peak_bw"]
    with open(CSV_AGG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)
    return rows

def report(roof, rows):
    print(f"\n{'='*88}")
    print(f" GPU: {roof['name']}  CC: {roof['cc'][0]}.{roof['cc'][1]}  SMs: {roof['sms']}")
    print(f" Peak FP64: {roof['peak_fp64']:.0f} GFLOP/s   Emp BW: {roof['emp_bw']:.0f} GB/s")
    for st in STENCILS:
        print(f"   {ST_NAMES[st]:>4}: OI={roof['oi'][st]:.4f} F/B  attainable={roof['attainable'][st]:.1f} GFLOP/s")
    print(f"{'='*88}")

    for st in STENCILS:
        sn = ST_NAMES[st]
        sr = [r for r in rows if r["stencil"] == sn]
        if not sr: continue
        print(f"\n ── {sn} ({ST_FLOPS[st]} FLOPs/pt, halo={ST_HALO[st]}) ──")
        print(f" {'var':<10} {'BX':>3} {'BY':>3} {'TX':>2} {'TY':>2} {'P':>1}"
              f"  {'medGF/s':>8} {'effBW':>6} {'%attn':>6} {'%pkBW':>6}")
        print(f" {'-'*56}")
        for r in sorted(sr, key=lambda x: -float(x["med_gflops"])):
            print(f" {r['variant']:<10} {r['BX']:>3} {r['BY']:>3} {r['TX']:>2} {r['TY']:>2} {r['PAD']:>1}"
                  f"  {float(r['med_gflops']):8.2f} {r['eff_bw_gbs']:>5}  {r['pct_attainable']:>5}% {r['pct_peak_bw']:>5}%")
        best = max(sr, key=lambda x: float(x["med_gflops"]))
        print(f"   Best: {best['variant']} BX={best['BX']} BY={best['BY']} TX={best['TX']} TY={best['TY']}"
              f" -> {best['med_gflops']} GFLOP/s ({best['pct_attainable']}% attn)")

if __name__ == "__main__":
    print("── Roofline ──")
    roof = get_roofline()
    for st in STENCILS:
        print(f"  {ST_NAMES[st]}: OI={roof['oi'][st]:.4f}  attainable={roof['attainable'][st]:.1f} GFLOP/s")
    if not os.path.exists(BINARY) or "--compile" in sys.argv:
        compile()
    print(f"\n── Sweep: N={N} tsteps={TSTEPS} stencils={[ST_NAMES[s] for s in STENCILS]} ──")
    sweep()
    rows = aggregate(roof)
    if rows: report(roof, rows)
    print(f"\nPer-iteration: {CSV_RAW}\nAggregated:    {CSV_AGG}")