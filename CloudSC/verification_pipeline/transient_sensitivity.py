"""transient_sensitivity.py — Noise-injection sensitivity analysis for CloudSC arrays.

Patches generated C++ so that every write to a target array gets multiplicative
noise:  value *= (1 + eps * uniform(-1,1)).  Arrays that amplify this noise the
most are the ones that need to stay fp64.

IMPORTANT: Run on the fp64 graph (no precision lowering) so that all arrays
are present in their original form (no _lowered variants).

Setup:
    # 1. Generate fp64 code (no precision lowering)
    python cloudsc_cpu_pipeline.py --lowprec fp64
    bash recompile.sh

    # 2. Run the sweep (generates its own baseline, then perturbs each array)
    python transient_sensitivity.py sweep codegen/cloudsc_py.cpp --eps 1e-7 --steps 1

    # 3. Use the ranked table to decide what goes in _LOWERING_EXCLUDE

Single-target usage:
    python transient_sensitivity.py inject codegen/cloudsc_py.cpp --targets zqxn zsolqa --eps 1e-7
    bash recompile.sh && ./cloudsc_cpu_bin 1 --save
    python compare_outputs.py --step 0
    python transient_sensitivity.py strip codegen/cloudsc_py.cpp

List available arrays:
    python transient_sensitivity.py list codegen/cloudsc_py.cpp
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path


def _discover_arrays(file_path: Path, kind: str = "all") -> list[str]:
    """Find array names in generated code.

    Args:
        kind: "transient" — only heap-allocated internal arrays
              "lowered"   — only boundary-cast _lowered arrays
              "param"     — only function parameter arrays
              "scalar"    — only local scalars inside the internal function
              "all"       — everything (arrays + scalars)
    """
    with open(file_path, "r") as f:
        code = f.read()

    # Heap-allocated transients: name = new type DACE_ALIGN(...)[size];
    heap_pattern = re.compile(r"(\w+)\s*=\s*new\s+(?:float|double)\s+DACE_ALIGN\(\d+\)\[")
    heap_names = set(heap_pattern.findall(code))

    # Function parameter arrays: double *__restrict__ name
    # Extract from __program_cloudsc_py_internal signature
    param_pattern = re.compile(r"(?:double|float|int)\s+\*__restrict__\s+(\w+)")
    sig_match = re.search(
        r"void __program_cloudsc_py_internal\((.*?)\)\s*\{",
        code, re.DOTALL,
    )
    if sig_match:
        param_names = set(param_pattern.findall(sig_match.group(1)))
    else:
        param_names = set()
    # Exclude non-arrays
    param_names -= {"handle", "ktype", "ldcum"}

    # Local scalars: double name; (inside the function body, after the opening brace)
    scalar_names = set()
    if sig_match:
        body_start = sig_match.end()
        # Scan local variable declarations at the start of the function body
        scalar_pattern = re.compile(r"^\s+(?:double|float)\s+(\w+)\s*;", re.MULTILINE)
        scalar_names = set(scalar_pattern.findall(code[body_start:]))
        # Exclude DaCe-generated expression temporaries
        _DACE_PREFIXES = ("__", "_", "abs_", "exp_", "min_", "max_", "fmax_",
                          "fmin_", "pow_", "tanh_", "log_", "sqrt_", "neg_", "expr_")
        _DACE_INFIXES = ("_slice", "_plus_", "_minus_", "_times_", "_over_", "_index",
                         "_div_", "_pow_")
        scalar_names = {n for n in scalar_names
                        if not n.startswith(_DACE_PREFIXES)
                        and not any(s in n for s in _DACE_INFIXES)}
        scalar_names -= heap_names  # already covered as arrays

    if kind == "transient":
        return sorted(n for n in heap_names if not n.endswith("_lowered"))
    elif kind == "lowered":
        return sorted(n for n in heap_names if n.endswith("_lowered"))
    elif kind == "param":
        return sorted(param_names)
    elif kind == "scalar":
        return sorted(scalar_names)
    return sorted(heap_names | param_names | scalar_names)


def inject_noise(file_path: Path, targets: list[str], eps: float = 1e-7) -> int:
    """Patch generated C++ to inject noise on target arrays.

    For arrays with write sites (transients, outputs):
      target[idx] = value;  →  target[idx] = transient_noise::perturb(value, eps);

    For read-only arrays (input parameters):
      CopyND<...>::Copy(target + offset, &scalar, stride);
      →  (same) followed by: scalar = transient_noise::perturb(scalar, eps);

    Skips boundary cast static_cast lines.
    Also injects the #include for transient_noise.h.
    """
    with open(file_path, "r") as f:
        code = f.read()

    # Inject include
    include = '#include "transient_noise.h"'
    if include not in code:
        code = re.sub(
            r"(#include\s+<dace/dace\.h>)",
            rf"\1\n{include}",
            code,
            count=1,
        )

    patched_count = 0

    for target in targets:
        # 1a. Patch array writes: target[...] = expr;
        write_pattern = re.compile(
            rf"(\b{re.escape(target)}\[([^\]]+)\])\s*=\s*([^;]+);"
        )

        def write_replacer(m):
            nonlocal patched_count
            lhs = m.group(1)
            value = m.group(3).strip()
            if "static_cast" in value:
                return m.group(0)
            patched_count += 1
            return f"{lhs} = transient_noise::perturb({value}, {eps});"

        code = write_pattern.sub(write_replacer, code)

        # 1b. Patch scalar writes: target = expr; (no index)
        #     Only match plain assignment, not declarations or pointer ops.
        scalar_write_pattern = re.compile(
            rf"(?<=\n)([ \t]+)\b{re.escape(target)}\b\s*=\s*([^;]+);"
        )

        def scalar_write_replacer(m):
            nonlocal patched_count
            indent = m.group(1)
            value = m.group(2).strip()
            if "static_cast" in value or "new " in value:
                return m.group(0)
            patched_count += 1
            return f"{indent}{target} = transient_noise::perturb({value}, {eps});"

        code = scalar_write_pattern.sub(scalar_write_replacer, code)

        # 2. Patch reads via CopyND: ...Copy(target + offset, &scalar, ...);
        #    Add noise to the scalar after the copy.
        read_pattern = re.compile(
            rf"(CopyND<[^>]+>::template\s+Const(?:Dst|Src)<[^>]*>::Copy\(\s*"
            rf"\b{re.escape(target)}\b[^;]*,\s*&(\w+),\s*[^)]*\);)"
        )

        def read_replacer(m):
            nonlocal patched_count
            original = m.group(1)
            scalar = m.group(2)
            patched_count += 1
            return f"{original}\n        {scalar} = transient_noise::perturb({scalar}, {eps});"

        code = read_pattern.sub(read_replacer, code)

    with open(file_path, "w") as f:
        f.write(code)

    print(f"  Injected noise (eps={eps}) into {patched_count} sites for {targets}")
    return patched_count


def strip_noise(file_path: Path):
    """Remove all transient_noise::perturb() calls, restoring original writes."""
    with open(file_path, "r") as f:
        code = f.read()

    code = re.sub(
        r"transient_noise::perturb\(([^,]+),\s*[^)]+\)",
        r"\1",
        code,
    )
    code = code.replace('#include "transient_noise.h"\n', "")

    with open(file_path, "w") as f:
        f.write(code)


def _recompile() -> bool:
    """Run recompile.sh, return True on success."""
    result = subprocess.run(
        ["bash", "recompile.sh"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"    Compile FAILED: {result.stderr[-200:]}")
        return False
    return True


def _run(steps: int = 1) -> bool:
    """Run the binary with --save, return True on success."""
    result = subprocess.run(
        ["./cloudsc_cpu_bin", str(steps), "--save"],
        capture_output=True, text=True,
        timeout=120,
    )
    if result.returncode != 0:
        print(f"    Run FAILED: {result.stderr[-200:]}")
        return False
    return True


def _compute_error_stats(ref_dir: str, test_dir: str, step: int):
    """Compare ref vs test output for given step.

    Finds the worst-affected output field by SNR.
    Returns (field_name, mean_abs_err, stderr_abs_err, field_abs_value, snr).
    """
    import h5py
    import numpy as np

    ref_path = f"{ref_dir}/cpp_output_step_{step}.h5"
    test_path = f"{test_dir}/cpp_output_step_{step}.h5"

    nans = ("", float("nan"), float("nan"), float("nan"), float("nan"))

    if not Path(ref_path).exists() or not Path(test_path).exists():
        return nans

    worst_field = ""
    worst_snr = float("inf")
    worst_mean_abs = 0.0
    worst_stderr_abs = 0.0
    worst_field_abs = 0.0

    with h5py.File(ref_path, "r") as f_ref, h5py.File(test_path, "r") as f_test:
        for key in sorted(set(f_ref.keys()) & set(f_test.keys())):
            ref_data = f_ref[key][()].astype(np.float64).ravel()
            test_data = f_test[key][()].astype(np.float64).ravel()
            if ref_data.shape != test_data.shape:
                continue
            signal = np.mean(np.square(ref_data))
            noise = np.mean(np.square(ref_data - test_data))
            if noise == 0:
                continue
            snr = 10 * np.log10(signal / noise) if signal > 0 else 0.0

            if snr < worst_snr:
                worst_snr = snr
                worst_field = key
                abs_err = np.abs(ref_data - test_data)
                worst_mean_abs = float(np.mean(abs_err))
                worst_stderr_abs = float(np.std(abs_err) / np.sqrt(len(abs_err)))
                worst_field_abs = float(np.mean(np.abs(ref_data)))

    if not worst_field:
        return nans

    return worst_field, worst_mean_abs, worst_stderr_abs, worst_field_abs, worst_snr


def sweep(
    file_path: Path,
    eps: float = 1e-7,
    steps: int = 1,
    samples: int = 10,
    ref_dir: str = "outputs_sweep_ref",
    test_dir: str = "outputs_cpp",
    targets: list[str] | None = None,
):
    """Sweep all (or specified) arrays: inject noise one at a time, measure min SNR.

    First generates an unperturbed C++ baseline in ref_dir, then for each target
    injects noise, recompiles, runs N samples, and averages the error stats.

    Prints a ranked table at the end (worst SNR first = most sensitive).
    """
    import numpy as np

    # Save a clean copy of the source
    backup = file_path.with_suffix(".cpp.bak")
    shutil.copy2(file_path, backup)

    # Generate unperturbed baseline
    print("Generating unperturbed baseline...")
    if not _recompile():
        print("ERROR: baseline compilation failed")
        return []
    if not _run(steps):
        print("ERROR: baseline run failed")
        return []
    ref_path = Path(ref_dir)
    if ref_path.exists():
        shutil.rmtree(ref_path)
    shutil.copytree(test_dir, ref_dir)
    print(f"Baseline saved to {ref_dir}/")

    if targets is None:
        targets = _discover_arrays(file_path)
    print(f"Sweeping {len(targets)} arrays (eps={eps}, steps={steps}, samples={samples})")

    results = []

    for i, name in enumerate(targets):
        print(f"[{i+1}/{len(targets)}] {name}")

        # Restore clean source
        shutil.copy2(backup, file_path)

        # Inject noise for this one target
        n_sites = inject_noise(file_path, [name], eps=eps)
        if n_sites == 0:
            print(f"    Unused (no read/write sites)")
            results.append((name, 0, "", float("nan"), float("nan"), float("nan"), float("inf")))
            continue

        if not _recompile():
            results.append((name, n_sites, "", float("nan"), float("nan"), float("nan"), float("nan")))
            continue

        # Run multiple samples (each run gets different random noise)
        sample_snrs = []
        sample_mean_abs = []
        sample_field_vals = []
        worst_field = ""
        failed = False

        for s in range(samples):
            if not _run(steps):
                failed = True
                break
            field, mean_abs, stderr_abs, field_val, snr = _compute_error_stats(
                ref_dir, test_dir, step=steps - 1)
            if snr != snr:  # nan
                failed = True
                break
            worst_field = field
            sample_snrs.append(snr)
            sample_mean_abs.append(mean_abs)
            sample_field_vals.append(field_val)

        if failed or not sample_snrs:
            results.append((name, n_sites, "", float("nan"), float("nan"), float("nan"), float("nan")))
            continue

        arr_snr = np.array(sample_snrs)
        arr_abs = np.array(sample_mean_abs)
        mean_snr = float(np.mean(arr_snr))
        mean_abs_err = float(np.mean(arr_abs))
        stderr_abs_err = float(np.std(arr_abs) / np.sqrt(len(arr_abs)))
        field_val = float(np.mean(sample_field_vals))

        results.append((name, n_sites, worst_field, mean_abs_err, stderr_abs_err, field_val, mean_snr))
        print(f"    {n_sites} sites, worst = {worst_field} (SNR={mean_snr:.1f} dB, "
              f"|err|={mean_abs_err:.2e} ± {stderr_abs_err:.2e}, |field|={field_val:.2e})")

    # Restore clean source
    shutil.copy2(backup, file_path)
    backup.unlink()

    # Print ranked table (sorted by SNR, worst first)
    results.sort(key=lambda r: r[6] if r[6] == r[6] else float("inf"))

    print("\n" + "=" * 105)
    print(f"Sensitivity ranking (eps={eps}, {steps} steps)")
    print(f"{'Array':<22} {'Sites':>5} {'Worst Field':<18} {'Mean |Err|':>12} {'Stderr':>12} {'|Field|':>12} {'SNR dB':>8}")
    print("-" * 105)
    for name, sites, field, mean_abs, stderr_abs, field_val, snr in results:
        if sites == 0:
            print(f"{name:<22} {'0':>5} {'(unused)':<18} {'---':>12} {'---':>12} {'---':>12} {'---':>8}")
            continue
        def fmt(v):
            return f"{v:.2e}" if v == v else "FAIL"
        snr_s = f"{snr:.1f}" if snr == snr else "FAIL"
        print(f"{name:<22} {sites:>5} {field:<18} {fmt(mean_abs):>12} {fmt(stderr_abs):>12} {fmt(field_val):>12} {snr_s:>8}")
    print("=" * 105)

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Transient sensitivity analysis")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # inject subcommand
    p_inject = sub.add_parser("inject", help="Inject noise into specific targets")
    p_inject.add_argument("file", type=Path)
    p_inject.add_argument("--targets", nargs="+", required=True)
    p_inject.add_argument("--eps", type=float, default=1e-7)

    # strip subcommand
    p_strip = sub.add_parser("strip", help="Remove all injected noise")
    p_strip.add_argument("file", type=Path)

    # sweep subcommand
    p_sweep = sub.add_parser("sweep", help="Sweep arrays and rank by sensitivity")
    p_sweep.add_argument("file", type=Path)
    p_sweep.add_argument("--eps", type=float, default=1e-7)
    p_sweep.add_argument("--steps", type=int, default=1, help="Number of timesteps to run")
    p_sweep.add_argument("--samples", type=int, default=10, help="Number of samples per array")
    p_sweep.add_argument("--targets", nargs="+", default=None, help="Specific targets (default: auto-discover)")
    p_sweep.add_argument("--kind", choices=["all", "transient", "lowered", "param", "scalar"], default="all",
                         help="Which arrays to sweep: transient (internal), lowered (boundary-cast), scalar, or all")
    p_sweep.add_argument("--ref-dir", default="outputs_sweep_ref")
    p_sweep.add_argument("--test-dir", default="outputs_cpp")

    # list subcommand
    p_list = sub.add_parser("list", help="List discoverable arrays")
    p_list.add_argument("file", type=Path)
    p_list.add_argument("--kind", choices=["all", "transient", "lowered", "param", "scalar"], default="all")

    args = parser.parse_args()

    if args.cmd == "inject":
        inject_noise(args.file, args.targets, args.eps)
    elif args.cmd == "strip":
        strip_noise(args.file)
    elif args.cmd == "list":
        names = _discover_arrays(args.file, kind=args.kind)
        print(f"Found {len(names)} arrays (kind={args.kind}):")
        for n in names:
            print(f"  {n}")
    elif args.cmd == "sweep":
        targets = args.targets or _discover_arrays(args.file, kind=args.kind)
        sweep(args.file, eps=args.eps, steps=args.steps,
              samples=args.samples, ref_dir=args.ref_dir,
              test_dir=args.test_dir, targets=targets)
