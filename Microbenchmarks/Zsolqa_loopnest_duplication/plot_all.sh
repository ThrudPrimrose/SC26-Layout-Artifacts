shopt -s nullglob

gpu=( *_gpu.csv )
other=()

for f in results/daint/*.csv; do
    [[ "$f" == *_gpu* ]] && continue
    other+=("$f")
done

python plot_violin.py --csv "${gpu[@]}" --output daint_gpu_results.pdf
python plot_violin.py --csv "${other[@]}" --output daint_cpu_results.pdf