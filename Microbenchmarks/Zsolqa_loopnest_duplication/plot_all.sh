shopt -s nullglob

gpu=( *_gpu.csv )
other=()

for f in *.csv; do
    [[ "$f" == *_gpu* ]] && continue
    other+=("$f")
done

python plot_violin.py --csv "${gpu[@]}" --output gpu_results.pdf
python plot_violin.py --csv "${other[@]}" --output cpu_results.pdf