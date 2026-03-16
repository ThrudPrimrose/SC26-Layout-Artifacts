import dace
import os
import shutil
import re
from pathlib import Path

from dace.codegen import codegen, compiler
from dace.sdfg import infer_types
from dace.transformation.interstate import LoopToMap


def repl_in_file(file_path: str, src: str, dst: str):
    with open(file_path, "r") as f:
        code = f.read()
    with open(file_path, "w") as f:
        f.write(code.replace(src, dst))


def modify_files_in_directory(directory):
    pattern = re.compile(r"^(\s*)int tmp_struct_symbol")
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if file_path.endswith((".c", ".h", ".cpp", ".cu")):
                modify_file(file_path, pattern)


def modify_file(file_path, pattern):
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    modified = False
    new_lines = []
    for line in lines:
        if pattern.match(line) and "(" not in line and "," not in line and ";" in line:
            line = pattern.sub(r"\1static int tmp_struct_symbol", line)
            modified = True
        new_lines.append(line)

    if modified:
        with open(file_path, "w", encoding="utf-8") as f:
            f.writelines(new_lines)


def flatten_build_folder(build_loc: Path, sdfg_name: str) -> tuple[str, str]:
    """
    Moves necessary files to the parent of build_loc and removes build_loc.
    Matches velocity's flatten_build_folder pattern.
    """
    cpu_path = None
    cpu_ext = "cpp"
    for ext in ["cpp", "cc"]:
        p = build_loc / "src" / "cpu" / f"{sdfg_name}.{ext}"
        if p.exists():
            cpu_path = p
            cpu_ext = ext
            break

    header_path = build_loc / "include" / f"{sdfg_name}.h"
    hash_path = build_loc / "include" / "hash.h"

    target_cpu = build_loc.parent / f"{sdfg_name}.{cpu_ext}"
    target_header = build_loc.parent / f"{sdfg_name}.h"
    target_hash = build_loc.parent / "hash.h"

    if cpu_path and cpu_path.exists():
        shutil.move(str(cpu_path), str(target_cpu))
    else:
        raise FileNotFoundError(
            f"Could not find {sdfg_name}.cpp in {build_loc / 'src' / 'cpu'}"
        )
    if header_path.exists():
        shutil.move(str(header_path), str(target_header))
    if hash_path.exists():
        shutil.move(str(hash_path), str(target_hash))

    # Clean up the original build folder
    shutil.rmtree(build_loc, ignore_errors=True)

    # Patch includes
    repl_in_file(str(target_cpu), "../../include/hash.h", "hash.h")
    repl_in_file(str(target_cpu), f"../../include/{sdfg_name}.h", f"{sdfg_name}.h")

    return str(target_cpu), str(target_header)


# --- Pipeline entry point ---


def main():
    import argparse
    import subprocess
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdfg", type=str, default="../cloudsc_py_unrolled_simplified.sdfgz")
    parser.add_argument("--release", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()

    sdfg = dace.SDFG.from_file(args.sdfg)
    sdfg.name = "cloudsc_py"

    # Set build folder explicitly (like velocity)
    codegen_dir = Path("codegen")
    if codegen_dir.exists():
        shutil.rmtree(codegen_dir)
    codegen_dir.mkdir()
    sdfg.build_folder = str(codegen_dir / sdfg.name)

    print(f"Generating code ({'Release' if args.release else 'Debug'})...")

    # Use velocity's lower-level codegen path
    sdfg.fill_scope_connectors()
    infer_types.infer_connector_types(sdfg)
    infer_types.set_default_schedule_and_storage_types(sdfg, None)
    sdfg.expand_library_nodes()
    infer_types.infer_connector_types(sdfg)
    infer_types.set_default_schedule_and_storage_types(sdfg, None)

    program_objects = codegen.generate_code(sdfg, validate=False)
    compiler.generate_program_folder(sdfg, program_objects, sdfg.build_folder)

    build_loc = Path(sdfg.build_folder)
    modify_files_in_directory(build_loc)

    # Flatten: move files up to codegen/ and clean up build subfolder
    _, header = flatten_build_folder(build_loc, sdfg.name)

    # cloudsc_main.cpp is maintained manually (includes sensitivity mode etc.)

    dace_runtime = Path(dace.__file__).parent / "runtime" / "include"

    # Detect HDF5
    h5_cflags = ""
    h5_libs = ""
    try:
        h5_cflags = subprocess.check_output(["pkg-config", "--cflags", "hdf5"], text=True).strip()
        h5_libs = subprocess.check_output(["pkg-config", "--libs", "hdf5"], text=True).strip()
    except Exception:
        try:
            h5_prefix = subprocess.check_output(["brew", "--prefix", "hdf5"], text=True).strip()
            h5_cflags = f"-I{h5_prefix}/include"
            h5_libs = f"-L{h5_prefix}/lib -lhdf5"
        except Exception:
            pass

    if args.release:
        cpp_flags = "-O3 -std=c++20 -DNDEBUG "
    else:
        cpp_flags = "-O0 -g -std=c++20 -Wall -Wextra "

    cmd = f"c++ {cpp_flags} -Icodegen -Iinclude -I{dace_runtime} {h5_cflags} cloudsc_main.cpp codegen/*.cpp -o cloudsc_cpu_bin -lpthread {h5_libs}"

    with open("recompile.sh", "w") as f:
        f.write(f"#!/bin/bash\nset -e\n{cmd}\n")
    os.chmod("recompile.sh", 0o755)
    print(f"Pipeline ready. Build command updated in recompile.sh (Release={args.release}, HDF5 support included)")


if __name__ == "__main__":
    main()
