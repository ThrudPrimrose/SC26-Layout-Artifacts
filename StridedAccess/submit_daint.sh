#!/bin/bash
#SBATCH --job-name=istrid
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --partition=normal
#SBATCH --time=03:00:00
#SBATCH --output=istrid_%j.out
#SBATCH --error=istrid_%j.err

# Load required modules
spack load gcc@15.1/atn6hp2
spack load python@3.15.5/l5hljwm
spack load cuda@12.9

# Set CUDA library path
export LD_LIBRARY_PATH=$(spack location -i cuda@12.9)/lib64:$LD_LIBRARY_PATH

# Print environment info
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "CUDA version:"
nvcc --version
echo ""
echo "GPU info:"
nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv
echo ""

# Set OpenMP threads for CPU parallelism
# Divide CPUs per task
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"

# Launch 1 Python process per GPU automatically using srun
srun --ntasks=$SLURM_NTASKS --gpus-per-task=1 --exact bash -c '
export RANK=$SLURM_PROCID
export LOCAL_RANK=$SLURM_LOCALID
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
export CUDA_ARCH=sm_90a

echo "Starting rank $RANK on GPU $CUDA_VISIBLE_DEVICES"
python3 benchmark_gpu.py
'

# Print completion info
echo ""
echo "Job completed at: $(date)"
