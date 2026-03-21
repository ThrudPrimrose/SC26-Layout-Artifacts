#!/usr/bin/env python3
"""Jacobi2D GPU sweep with integrated FP64 roofline analysis."""
import subprocess, sys, os, csv
import numpy as np

# ── Config ──
BINARY = "./jacobi2d_gpu"
CSV_FILE = "jacobi2d_results.csv"
N = 8192
TSTEPS = 100

CONFIGS = [
    # (BX, BY, TX, TY)
    (32, 8, 1, 1), (16, 16, 1, 1), (32, 4, 1, 1),
    (32, 8, 2, 1), (16, 16, 2, 1), (32, 4, 2, 1),
    (32, 8, 1, 2), (16, 16, 1, 2),
    (16, 8, 2, 2), (32, 4, 2, 2),
    (32, 8, 4, 1), (16, 16, 4, 1), (32, 8, 1, 4),
    (16, 8, 4, 2), (32, 8, 4, 4),
]
VARIANTS = [0, 1, 2, 3]  # global, smem, smem_pad, smem_swap
PADS = [1, 2]

# ── Roofline from device properties + empirical BW ──
def get_roofline():
    try:
        import pycuda.driver as cuda
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pycuda", "-q"])
        import pycuda.driver as cuda
    import pycuda.autoinit
    from pycuda.compiler import SourceModule

    dev = pycuda.autoinit.device
    attrs = dev.get_attributes()
    sm_count  = attrs[cuda.device_attribute.MULTIPROCESSOR_COUNT]
    clock_mhz = attrs[cuda.device_attribute.CLOCK_RATE] / 1000
    mem_mhz   = attrs[cuda.device_attribute.MEMORY_CLOCK_RATE] / 1000
    bus_w     = attrs[cuda.device_attribute.GLOBAL_MEMORY_BUS_WIDTH]
    cc = (attrs[cuda.device_attribute.COMPUTE_CAPABILITY_MAJOR],
          attrs[cuda.device_attribute.COMPUTE_CAPABILITY_MINOR])

    FP64_MAP = {(7,0):32,(7,5):2,(8,0):32,(8,6):2,(8,7):2,(8,9):2,(9,0):64}
    fp64u = FP64_MAP.get(cc, 32 if cc[1]==0 and cc[0]>=7 else 2)

    peak_fp64 = sm_count * fp64u * 2 * (clock_mhz / 1000)  # GFLOP/s
    peak_bw   = 2 * mem_mhz * (bus_w / 8) / 1000            # GB/s

    # Empirical BW via copy + daxpy
    SZ = 128 * 1024 * 1024  # 128M doubles = 1GB
    a_gpu = cuda.mem_alloc(SZ * 8)
    b_gpu = cuda.mem_alloc(SZ * 8)
    mod = SourceModule("""
    __global__ void kcopy(const double* s, double* d, long long n) {
        long long i = blockIdx.x*(long long)blockDim.x+threadIdx.x; if(i<n) d[i]=s[i]; }
    __global__ void kdaxpy(double a, const double* x, double* y, long long n) {
        long long i = blockIdx.x*(long long)blockDim.x+threadIdx.x; if(i<n) y[i]=a*x[i]+y[i]; }
    """)
    kcopy = mod.get_function("kcopy")
    kdaxpy = mod.get_function("kdaxpy")
    bl, gr = 256, (SZ + 255) // 256

    def bench(fn, reps=20):
        s, e = cuda.Event(), cuda.Event()
        fn(); cuda.Context.synchronize()  # warmup
        s.record()
        for _ in range(reps): fn()
        e.record(); e.synchronize()
        return s.time_till(e) / reps / 1000

    copy_bw  = 2*SZ*8 / bench(lambda: kcopy(a_gpu,b_gpu,np.int64(SZ),block=(bl,1,1),grid=(gr,1))) / 1e9
    daxpy_bw = 3*SZ*8 / bench(lambda: kdaxpy(np.float64(2.0),a_gpu,b_gpu,np.int64(SZ),block=(bl,1,1),grid=(gr,1))) / 1e9
    emp_bw = max(copy_bw, daxpy_bw)
    a_gpu.free(); b_gpu.free()

    # Jacobi2D OI: read full in[] + write full out[] per iteration
    # bytes_per_iter = 2 * (N+2)^2 * 8,  flops_per_iter = 5 * N^2
    S = N + 2
    oi = (5.0 * N * N) / (2.0 * S * S * 8)
    attainable = min(peak_fp64, oi * emp_bw)

    return dict(name=dev.name(), cc=cc, sms=sm_count,
                peak_fp64=peak_fp64, peak_bw=peak_bw,
                emp_bw=emp_bw, copy_bw=copy_bw, daxpy_bw=daxpy_bw,
                oi=oi, attainable=attainable)

# ── Compile ──
def compile():
    cmd = f"nvcc -O3 -arch=native -o {BINARY} jacobi2d_gpu.cu"
    print(f"Compiling: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

# ── Sweep ──
def sweep(roof):
    # Temp CSV without roofline columns (appended by C++ binary)
    tmp_csv = CSV_FILE + ".tmp"
    open(tmp_csv, "w").close()

    for bx, by, tx, ty in CONFIGS:
        for v in VARIANTS:
            pads = PADS if v == 2 else [0]
            for pad in pads:
                cmd = [BINARY, str(N), str(TSTEPS), str(v),
                       str(bx), str(by), str(tx), str(ty), str(pad), tmp_csv]
                try:
                    r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
                    if r.returncode != 0:
                        print(f"  FAIL v={v} {bx}x{by} TX={tx} TY={ty}: {r.stderr.strip()[:80]}")
                        continue
                    print(f"  {r.stdout.strip()}")
                except subprocess.TimeoutExpired:
                    print(f"  TIMEOUT v={v} {bx}x{by} TX={tx} TY={ty}")

    # Read raw CSV, compute roofline %, write final
    S = N + 2
    bytes_per_run = TSTEPS * 2.0 * S * S * 8
    rows = []
    with open(tmp_csv) as f:
        for line in f:
            parts = line.strip().split(",")
            if len(parts) < 11: continue
            variant, n, ts, bx, by, tx, ty, pad, t_s, gf, ck = parts
            gf_f = float(gf); t_f = float(t_s)
            eff_bw = bytes_per_run / t_f / 1e9
            pct_att = 100 * gf_f / roof["attainable"]
            pct_bw  = 100 * eff_bw / roof["emp_bw"]
            rows.append(dict(variant=variant, N=n, tsteps=ts, BX=bx, BY=by,
                             TX=tx, TY=ty, PAD=pad, time_s=t_s, gflops=gf,
                             eff_bw_gbs=f"{eff_bw:.1f}",
                             pct_attainable=f"{pct_att:.1f}",
                             pct_peak_bw=f"{pct_bw:.1f}", checksum=ck))

    fields = ["variant","N","tsteps","BX","BY","TX","TY","PAD",
              "time_s","gflops","eff_bw_gbs","pct_attainable","pct_peak_bw","checksum"]
    with open(CSV_FILE, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    os.remove(tmp_csv)
    return rows

# ── Report ──
def report(roof, rows):
    print(f"\n{'='*78}")
    print(f" GPU: {roof['name']}   CC: {roof['cc'][0]}.{roof['cc'][1]}   SMs: {roof['sms']}")
    print(f" Peak FP64: {roof['peak_fp64']:.0f} GFLOP/s    Theo BW: {roof['peak_bw']:.0f} GB/s")
    print(f" Empirical BW: {roof['emp_bw']:.0f} GB/s  (copy={roof['copy_bw']:.0f}  daxpy={roof['daxpy_bw']:.0f})")
    print(f" Jacobi2D OI: {roof['oi']:.4f} FLOP/byte   Attainable: {roof['attainable']:.1f} GFLOP/s")
    print(f"{'='*78}")
    hdr = (f" {'variant':<10} {'BX':>3} {'BY':>3} {'TX':>2} {'TY':>2} {'P':>1}"
           f"  {'GFLOP/s':>8} {'eff BW':>7} {'%attn':>6} {'%pkBW':>6} {'time':>8}")
    print(hdr)
    print(f" {'-'*74}")
    for r in sorted(rows, key=lambda x: -float(x["gflops"])):
        print(f" {r['variant']:<10} {r['BX']:>3} {r['BY']:>3} {r['TX']:>2} {r['TY']:>2} {r['PAD']:>1}"
              f"  {float(r['gflops']):8.2f} {r['eff_bw_gbs']:>6}  {r['pct_attainable']:>5}% {r['pct_peak_bw']:>5}% {float(r['time_s']):8.4f}")

    best = max(rows, key=lambda x: float(x["gflops"]))
    print(f"\n Best: {best['variant']} BX={best['BX']} BY={best['BY']} TX={best['TX']} TY={best['TY']}"
          f" PAD={best['PAD']} -> {float(best['gflops']):.2f} GFLOP/s"
          f" ({best['pct_attainable']}% of attainable, {best['pct_peak_bw']}% of peak BW)")

if __name__ == "__main__":
    print("── Roofline analysis ──")
    roof = get_roofline()
    print(f"  Peak FP64: {roof['peak_fp64']:.0f} GFLOP/s  |  Emp BW: {roof['emp_bw']:.0f} GB/s")
    print(f"  Jacobi2D OI={roof['oi']:.4f}  ->  Attainable={roof['attainable']:.1f} GFLOP/s")

    if not os.path.exists(BINARY) or "--compile" in sys.argv:
        compile()
    print(f"\n── Sweep: N={N}  tsteps={TSTEPS} ──")
    rows = sweep(roof)
    if rows:
        report(roof, rows)
    print(f"\nCSV: {CSV_FILE}")