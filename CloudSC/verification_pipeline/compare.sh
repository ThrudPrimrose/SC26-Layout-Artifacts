python cloudsc_python_ref_run.py --steps 2 --save
python cloudsc_cpu_pipeline.py
./recompile.sh
./cloudsc_cpu_bin 2 --save
python compare_outputs.py