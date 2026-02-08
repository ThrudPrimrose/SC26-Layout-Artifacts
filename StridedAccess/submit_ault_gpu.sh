#!/bin/bash
#SBATCH --job-name=nvstrid
#SBATCH --nodes=1
#SBATCH --ntasks=4                # 4 tasks (one per GPU)
#SBATCH --cpus-per-task=64        # Total 288 CPUs / 4 tasks
#SBATCH --gpus-per-task=1         # 1 GPU per task
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --partition=amda100
#SBATCH --time=04:00:00
#SBATCH --output=istrid_%j_%t.out
#SBATCH --error=istrid_%j_%t.err

# Load required modules
spack load gcc@14.2
spack load python@3.12.9%gcc@14.2
spack load sqlite
spack load cuda@12.9 gcc@14.2

# Set CUDA library path
export LD_LIBRARY_PATH=$(spack location -i cuda)/lib64:$LD_LIBRARY_PATH

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

echo "Starting rank $RANK on GPU $CUDA_VISIBLE_DEVICES"
python3 benchmark_gpu.py
'

# Print completion info
echo ""
echo "Job completed at: $(date)"
