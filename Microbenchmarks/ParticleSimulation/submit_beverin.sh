#!/bin/bash
#SBATCH --job-name=nbody_cpu
#SBATCH --nodes=1
#SBATCH --partition=mi300
#SBATCH --time=03:00:00
#SBATCH --output=nbody_%j.out
#SBATCH --error=nbody_%j.err

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=192

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=96
export OMP_PROC_BIND=close
export OMP_PLACES=threads

# Strong pinning via Slurm
export SLURM_CPU_BIND=cores

# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
N=8192        # large enough to stress memory (buffer = ~400MB)
STEPS=100      # more time steps
REPS=50       # as requested
VL=8

python3.11 particle_simulation.py \
        --N $N \
        --steps $STEPS \
        --vl $VL \
        --ic random
	--dace

