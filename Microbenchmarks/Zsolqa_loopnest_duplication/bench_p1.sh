#!/bin/bash
#SBATCH --job-name=zsolqa
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --exclusive
#SBATCH --time=00:29:59
#SBATCH --output=zsolqa_output.txt
#SBATCH --error=zsolqa_error.txt

export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$(pyenv init -)"
eval "$(pyenv virtualenv-init -)"

pyenv activate dace_py_12

export _DACE_NO_SYNC=1
export LD_LIBRARY_PATH=$(spack location -i cuda@12.9)/lib64:$LD_LIBRARY_PATH
spack load cuda@12.9
spack load gcc/76jw6nu

sleep 1
python zsolqa_loopnest.py --variant=3d-gpu --reps=100 --csv bench_3d_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=split-gpu --reps=100 --csv bench_split_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=3d-reduce-gpu --reps=100 --csv bench_reduce_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=split-reduce-gpu --reps=100 --csv bench_reduce_split_gpu.csv
sleep 1
#python zsolqa_loopnest.py --variant=3d --reps=100 --csv bench_3d.csv
sleep 1
#python zsolqa_loopnest.py --variant=split --reps=100 --csv bench_split.csv
sleep 1
#python zsolqa_loopnest.py --variant=reduce --reps=100 --csv bench_reduce.csv
sleep 1
#python zsolqa_loopnest.py --variant=reduce-split --reps=100 --csv bench_reduce_split.csv

#sleep 1
#ncu --set full --import-source=yes -f -o profile_3d python zsolqa_loopnest.py --variant=3d-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_split python zsolqa_loopnest.py --variant=split-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_reduce python zsolqa_loopnest.py --variant=reduce-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_reduce python zsolqa_loopnest.py --variant=reduce_split-gpu --reps=2 --csv __tr.csv