#!/bin/bash
#SBATCH --job-name=elemwise
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-task=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --partition=normal
#SBATCH --time=05:00:00
#SBATCH --output=elemwise_%j.out
#SBATCH --error=elemwise_%j.err

# Load required modules
export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"
pyenv activate dace_py_12

spack load gcc/76jw6nu
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

# Create shared results directory
mkdir -p results

# Launch 1 Python process per GPU
srun --ntasks=$SLURM_NTASKS --gpus-per-task=1 --exact bash -c '
export RANK=$SLURM_PROCID
export LOCAL_RANK=$SLURM_LOCALID
export CUDA_VISIBLE_DEVICES=$LOCAL_RANK
export CUDA_ARCH=sm_90a
echo "Starting rank $RANK on GPU $CUDA_VISIBLE_DEVICES"
python3 benchmark_gpu.py
'

# Merge per-rank CSV files into one
echo ""
echo "Merging results..."
head -1 results/results_rank0.csv > results/results_all.csv
for f in results/results_rank*.csv; do
    tail -n +2 "$f" >> results/results_all.csv
done
echo "Merged CSV: results/results_all.csv"
wc -l results/results_all.csv

echo ""
echo "Job completed at: $(date)"