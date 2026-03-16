import argparse
import os
import dace
from pathlib import Path
from dace.transformation.passes import GPUKernelLaunchRestructure

def apply_gpu_optimizations(sdfg: dace.SDFG):
    print("Applying GPU optimizations...")
    # Standard GPU transformation
    sdfg.apply_gpu_transformations()
    
    # Set GPU schedule for maps and Global storage for transients
    for state in sdfg.nodes():
        for node in state.nodes():
            if isinstance(node, dace.nodes.MapEntry):
                node.schedule = dace.ScheduleType.GPU_Device
            elif isinstance(node, dace.nodes.AccessNode):
                if node.desc(sdfg).transient:
                    node.desc(sdfg).storage = dace.StorageType.GPU_Global

    # Ensure concurrent streams and C types
    dace.config.Config.set("compiler", "cuda", "max_concurrent_streams", value="10")
    dace.config.Config.set("compiler", "default_data_types", value="C")
    
    # Validation
    sdfg.validate()
    return sdfg

def generate_recompile_script(sdfg: dace.SDFG, main_source: str, output_bin: str):
    # Determine paths based on DaCe's default structure
    build_path = Path(".dacecache") / sdfg.name / "build"
    src_cpu = Path(".dacecache") / sdfg.name / "src" / "cpu"
    src_cuda = Path(".dacecache") / sdfg.name / "src" / "cuda"
    include_path = Path(".dacecache") / sdfg.name / "include"
    
    # Construct the nvcc command
    # We include all generated .cpp and .cu files
    cmd = [
        "nvcc -O3",
        f"-I{include_path}",
        main_source,
        f"{src_cpu}/*.cpp",
        f"{src_cuda}/*.cu",
        f"-o {output_bin}",
        "-lcudart"
    ]
    
    recompile_content = "#!/bin/bash\n"
    recompile_content += " ".join(cmd) + "\n"
    
    with open("recompile.sh", "w") as f:
        f.write(recompile_content)
    os.chmod("recompile.sh", 0o755)
    print("Generated recompile.sh")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdfg", type=str, default="cloudsc_pydace_simplified.sdfgz")
    parser.add_argument("--main", type=str, default="cloudsc_main.cu")
    parser.add_argument("--output", type=str, default="cloudsc_gpu_bin")
    args = parser.parse_args()

    if not os.path.exists(args.sdfg):
        print(f"Error: {args.sdfg} not found.")
        return

    print(f"Loading {args.sdfg}...")
    sdfg = dace.SDFG.from_file(args.sdfg)
    sdfg.name = "cloudsc_py"
    
    # Apply optimizations
    apply_gpu_optimizations(sdfg)
    
    # Generate code only (don't try to build the .so since we might not have the full env)
    print("Generating code...")
    sdfg.generate_code()
    
    # Create the standalone recompile script
    generate_recompile_script(sdfg, args.main, args.output)
    print(f"Pipeline complete. You can now run './recompile.sh' on a GPU machine.")

if __name__ == "__main__":
    main()
