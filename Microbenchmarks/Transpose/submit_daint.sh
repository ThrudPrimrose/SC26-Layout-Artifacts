#!/bin/bash
#SBATCH --job-name=nbody
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --time=00:30:00
#SBATCH --output=nbody_daint_%j.out
#SBATCH --error=nbody_daint_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=288

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PROC_BIND=true
export OMP_PLACES=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE


echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

spack load gcc/76jw6nu # 14.3
spack loadd cuda@12.9
# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
set -e

source ${SCRATCH}/yakup-dev-env/bin/activate

python run_transpose.py