#!/bin/bash
#SBATCH --job-name=zsolqa
#SBATCH --nodes=1
#SBATCH --partition=normal
#SBATCH --exclusive
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:29:59
#SBATCH --output=zsolqa_output_beverin.txt
#SBATCH --error=zsolqa_error_beverin.txt

export _DACE_NO_SYNC=1
export OMP_NUM_THREADS=96
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=TRUE

echo "Running on $(hostname)"
echo "Threads: $OMP_NUM_THREADS"

sleep 1
python zsolqa_loopnest.py --variant=3d-gpu --reps=100 --csv bench_3d_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=split-gpu --reps=100 --csv bench_split_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=3d-reduce-gpu --reps=100 --csv bench_reduce_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=split-reduce-gpu --reps=100 --csv bench_reduce_split_gpu.csv
sleep 1
python zsolqa_loopnest.py --variant=3d --reps=100 --csv bench_3d.csv
sleep 1
python zsolqa_loopnest.py --variant=split --reps=100 --csv bench_split.csv
sleep 1
python zsolqa_loopnest.py --variant=3d-reduce --reps=100 --csv bench_reduce.csv
sleep 1
python zsolqa_loopnest.py --variant=split-reduce --reps=100 --csv bench_reduce_split.csv

mv *.csv results/beverin

#sleep 1
#ncu --set full --import-source=yes -f -o profile_3d python zsolqa_loopnest.py --variant=3d-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_split python zsolqa_loopnest.py --variant=split-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_reduce python zsolqa_loopnest.py --variant=reduce-gpu --reps=2 --csv __tr.csv
#sleep 1
#ncu --set full --import-source=yes -f -o profile_reduce python zsolqa_loopnest.py --variant=reduce_split-gpu --reps=2 --csv __tr.csv