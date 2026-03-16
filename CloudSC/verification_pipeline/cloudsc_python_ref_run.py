import numpy as np
import h5py
import os
import time
import argparse
from types import SimpleNamespace

# Import the reference implementation we just copied
try:
    from cloudsc_python_ref import cloudsc_py
except ImportError:
    print("Error: cloudsc_python_ref.py not found. Please ensure it was copied from the dwarf repo.")
    exit(1)

def load_parameters(path):
    yrecldp = SimpleNamespace()
    ydcst = SimpleNamespace()
    ydthf = SimpleNamespace()
    
    with h5py.File(path, 'r') as f:
        # Map HDF5 keys to the objects based on prefixes
        for k in f.keys():
            if f[k].shape == (1,):
                val = f[k][0]
                if k.startswith('YRECLDP_'):
                    setattr(yrecldp, k[8:].lower(), val)
                elif k.startswith('YDCST_'):
                    setattr(ydcst, k[6:].lower(), val)
                elif k.startswith('YDTHF_'):
                    setattr(ydthf, k[6:].lower(), val)
                elif k in ['RG', 'RD', 'RCPD', 'RETV', 'RLVTT', 'RLSTT', 'RLMLT', 'RTT', 'RV']:
                    setattr(ydcst, k.lower(), val)
                elif k in ['R2ES', 'R3LES', 'R3IES', 'R4LES', 'R4IES', 'R5LES', 'R5IES', 
                           'R5ALVCP', 'R5ALSCP', 'RALVDCP', 'RALSDCP', 'RALFDCP', 
                           'RTWAT', 'RTICE', 'RTICECU', 'RTWAT_RTICE_R', 'RTWAT_RTICECU_R', 
                           'RKOOP1', 'RKOOP2']:
                    setattr(ydthf, k.lower(), val)
    return ydcst, ydthf, yrecldp

def run_ref_simulation():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steps", type=int, default=1)
    parser.add_argument("--input", type=str, default="input.h5")
    parser.add_argument("--save", action="store_true", help="Save outputs to HDF5")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found.")
        return

    print(f"Reading input data from: {args.input}")
    fields = {}
    with h5py.File(args.input, 'r') as f:
        klon = f['KLON'][0]
        klev = f['KLEV'][0]
        ptsphy = f['PTSPHY'][0]
        for k in f.keys():
            if len(f[k].shape) > 0 and f[k].shape != (1,):
                fields[k.lower()] = np.asfortranarray(f[k][()])

    ydcst, ydthf, yrecldp = load_parameters(args.input)
    # ... missing fields init ...
    nclv = 5
    # List of tendency/output arrays to potentially save
    output_arrays = ['tendency_loc_a', 'tendency_loc_t', 'tendency_loc_q', 'tendency_loc_cld',
                     'pcovptot', 'prainfrac_toprfz', 'plude',
                     'pfsqlf', 'pfsqif', 'pfcqnng', 'pfcqlng',
                     'pfsqrf', 'pfsqsf', 'pfcqrng', 'pfcqsng', 'pfsqltur', 'pfsqitur',
                     'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn']
    
    for f_name in output_arrays:
        if f_name not in fields:
            if f_name == 'prainfrac_toprfz': shape = (klon,)
            elif f_name == 'tendency_loc_cld': shape = (nclv, klev, klon)
            elif f_name.startswith('pf'): shape = (klev+1, klon)
            else: shape = (klev, klon)
            fields[f_name] = np.zeros(shape, order='F')

    import inspect
    sig = inspect.signature(cloudsc_py)
    arg_names = list(sig.parameters.keys())
    
    call_args = []
    # (Mapping logic same as before)
    for name in arg_names:
        if name == 'kidia': call_args.append(np.int32(1))
        elif name == 'kfdia': call_args.append(np.int32(klon))
        elif name == 'klon': call_args.append(np.int32(klon))
        elif name == 'klev': call_args.append(np.int32(klev))
        elif name == 'ptsphy': call_args.append(np.float64(ptsphy))
        elif name == 'ydcst': call_args.append(ydcst)
        elif name == 'ydthf': call_args.append(ydthf)
        elif name == 'yrecldp': call_args.append(yrecldp)
        elif name in fields:
            call_args.append(fields[name])
        else:
            raise ValueError(f"Unknown parameter '{name}' in cloudsc_py signature — not in fields or known scalars")

    total_time = 0
    for i in range(args.steps):
        start = time.perf_counter()
        cloudsc_py(*call_args)
        end = time.perf_counter()
        total_time += (end - start)
        
        if args.save:
            os.makedirs("outputs_ref", exist_ok=True)
            out_name = f"outputs_ref/ref_output_step_{i}.h5"
            with h5py.File(out_name, 'w') as f_out:
                for f_name in output_arrays:
                    f_out.create_dataset(f_name.upper(), data=fields[f_name])
            print(f"Saved results to: {os.path.abspath(out_name)}")

if __name__ == "__main__":
    run_ref_simulation()
