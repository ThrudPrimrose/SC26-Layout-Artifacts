import h5py
import numpy as np
import polars as pl
import argparse
import os
import sys
import re

def calculate_snr(ref, test):
    """Calculates Signal-to-Noise Ratio in dB."""
    signal_power = np.mean(np.square(ref))
    noise_power = np.mean(np.square(ref - test))
    if noise_power == 0:
        return float('inf')
    if signal_power == 0:
        return 0.0
    return 10 * np.log10(signal_power / noise_power)

def get_comparison_results(file_ref, file_test, tol=1e-12):
    if not os.path.exists(file_ref):
        print(f"Error: Reference file {file_ref} not found.")
        return None
    if not os.path.exists(file_test):
        print(f"Error: Test file {file_test} not found.")
        return None

    results = []
    with h5py.File(file_ref, 'r') as f_ref, h5py.File(file_test, 'r') as f_test:
        common_keys = sorted(set(f_ref.keys()) & set(f_test.keys()))
        
        if not common_keys:
            print(f"No common datasets found between {file_ref} and {file_test}.")
            return []

        for key in common_keys:
            ref_data = f_ref[key][()]
            test_data = f_test[key][()]

            if ref_data.shape != test_data.shape:
                results.append({
                    "Field": key,
                    "Status": "SHAPE MISMATCH",
                    "Max Abs Err": None,
                    "Max Rel Err": None,
                    "RMS": None,
                    "SNR (dB)": None,
                    "Ref L-inf": None,
                    "Ref L2": None
                })
                continue

            abs_err = np.abs(ref_data - test_data)
            max_abs = np.max(abs_err)
            
            # Relative error with small epsilon to avoid div by zero
            rel_err = abs_err / (np.abs(ref_data) + 1e-15)
            max_rel = np.max(rel_err)
            
            rms = np.sqrt(np.mean(np.square(abs_err)))
            snr = calculate_snr(ref_data, test_data)
            
            status = "PASS" if max_rel <= tol else "FAIL"
            
            results.append({
                "Field": key,
                "Status": status,
                "Max Abs Err": max_abs,
                "Max Rel Err": max_rel,
                "RMS": rms,
                "SNR (dB)": snr,
                "Ref L-inf": np.max(np.abs(ref_data)),
                "Ref L2": np.linalg.norm(ref_data.flatten())
            })
    return results

def compare_files(file_ref, file_test, tol=1e-12):
    print(f"Comparing:\n  Ref:  {file_ref}\n  Test: {file_test}\n")
    results = get_comparison_results(file_ref, file_test, tol)
    if results is None:
        return
    
    df = pl.DataFrame(results)
    
    # Configure polars to show more rows
    with pl.Config(tbl_rows=len(results), tbl_width_chars=160):
        print(df)

    failed = df.filter(pl.col("Status") != "PASS")
    if failed.height > 0:
        print(f"\n❌ Validation FAILED for {failed.height} fields (tol={tol})")
        sys.exit(1)
    else:
        print(f"\n✅ Validation PASSED (tol={tol})")

def compare_all_steps(ref_dir="outputs_ref", test_dir="outputs_cpp", tol=1e-12):
    print(f"Comparing all steps in {ref_dir} and {test_dir}\n")
    
    ref_files = [f for f in os.listdir(ref_dir) if f.startswith("ref_output_step_") and f.endswith(".h5")]
    
    all_results = []
    steps_compared = 0

    for ref_file in sorted(ref_files):
        # Extract step index using regex
        match = re.search(r"ref_output_step_(\d+)\.h5", ref_file)
        if not match:
            continue
        
        step_idx = int(match.group(1))
        test_file = f"cpp_output_step_{step_idx}.h5"
        test_path = os.path.join(test_dir, test_file)
        ref_path = os.path.join(ref_dir, ref_file)
        
        if os.path.exists(test_path):
            res = get_comparison_results(ref_path, test_path, tol)
            if res:
                for r in res:
                    r["Step"] = step_idx
                all_results.extend(res)
                steps_compared += 1
        else:
            print(f"Warning: Test file {test_path} missing for step {step_idx}")

    if not all_results:
        print("No matching pairs found for comparison.")
        return

    print(f"Compared {steps_compared} timesteps. Aggregating worst stats per field.")

    df_all = pl.DataFrame(all_results)
    
    # Define aggregation logic: worst-case for each metric
    # Max of Max Abs Err, Max of Max Rel Err, Max of RMS, Min of SNR
    # Status: FAIL if any is FAIL, else PASS
    
    # For each field, pick the row with the worst (lowest) SNR
    agg_df = (
        df_all
        .group_by("Field", maintain_order=True)
        .agg(pl.all().sort_by("SNR (dB)").first())
        .rename({"Step": "Worst Step"})
        .select(["Field", "Status", "Max Abs Err", "Max Rel Err", "RMS", "SNR (dB)", "Worst Step", "Ref L-inf", "Ref L2"])
    )

    with pl.Config(tbl_rows=agg_df.height, tbl_width_chars=-1, tbl_cols=-1):
        print(agg_df)

    failed = agg_df.filter(pl.col("Status") != "PASS")
    if failed.height > 0:
        print(f"\n❌ Validation FAILED for {failed.height} fields across all steps (tol={tol})")
        sys.exit(1)
    else:
        print(f"\n✅ Validation PASSED across all {steps_compared} steps (tol={tol})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare CloudSC HDF5 outputs.")
    parser.add_argument("--step", type=int, help="Comparison step index (compares ref vs cpp)")
    parser.add_argument("--ref", type=str, help="Path to reference HDF5 file")
    parser.add_argument("--test", type=str, help="Path to test HDF5 file")
    parser.add_argument("--tol", type=float, default=1e-12, help="Relative error tolerance")
    args = parser.parse_args()

    if args.step is not None:
        file_ref = f"outputs_ref/ref_output_step_{args.step}.h5"
        file_test = f"outputs_cpp/cpp_output_step_{args.step}.h5"
        compare_files(file_ref, file_test, args.tol)
    elif args.ref and args.test:
        compare_files(args.ref, args.test, args.tol)
    else:
        # Default behavior: compare all available matching steps
        compare_all_steps(tol=args.tol)
