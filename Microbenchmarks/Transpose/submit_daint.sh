#!/bin/bash
#SBATCH --job-name=transpose_gpu_daint
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --time=02:00:00
#SBATCH --output=transpose_gpu_daint_%j.out
#SBATCH --error=transpose_gpu_daint_%j.err
#SBATCH --ntasks=1
#SBATHC --gpus-per-task=1
#SBATCH --account=g177-1
#SBATCH --cpus-per-task=288
#SBATCH --exclusive

# -------------------------------
# OpenMP configuration
# -------------------------------
export OMP_NUM_THREADS=288
export OMP_PLACES="{0}:72:1,{72}:72:1,{144}:72:1,{216}:72:1"
export OMP_PROC_BIND=close
# Optional: better NUMA behavior
export OMP_DISPLAY_ENV=TRUE


echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

spack load gcc/76jw6nu # 14.3
spack load cutensor
spack load cuda@12.9

# -------------------------------
# Workload parameters (BIG!)
# -------------------------------
set -e
export CUTENSOR_HOME=$(spack location -i cutensor)
source ${SCRATCH}/yakup-dev-env/bin/activate
export CUDA_HOME=$(spack location -i cuda@12.9)
export LD_LIBRARY_PATH=${CUDA_HOME}/lib64:$LD_LIBRARY_PATH


export C_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$SCRATCH/include:$CUTENSOR_HOME/include:$CPLUS_INCLUDE_PATH
export LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LIBRARY_PATH
export LD_LIBRARY_PATH=$SCRATCH/lib:$SCRATCH/lib64:$CUTENSOR_HOME/lib/12:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
export PATH=$SCRATCH/bin:$CUTENSOR_HOME/bin:$PATH

export HPTT_ROOT=$SCRATCH
python run_transpose.py --compile --lib-only