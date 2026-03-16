#!/bin/bash
set -e

# Use our local venv
PYTHON_BIN="./.venv/bin/python"

echo "=== Starting CloudSC GPU Pipeline (Self-Contained) ==="

# 1. Run the Python part: optimization + code generation
echo "[Step 1] Optimizing SDFG and generating code..."
$PYTHON_BIN cloudsc_gpu_pipeline.py

# 2. Inform the user about the generated recompile.sh
echo ""
echo "[Step 2] Build orchestration ready."
echo "A 'recompile.sh' has been generated. Since you are on macOS (no CUDA),"
echo "please run './recompile.sh' on a Linux machine with 'nvcc' to build the final binary."
echo "Note: The 'cloudsc_main.cu' will need to be updated with all array arguments from the SDFG."
echo ""
echo "=== Done ==="
