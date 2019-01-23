#!/bin/bash
for item in {1..15}; do	
	cd n_is_$item
	cp ../matrix_pencil_no_plot.slurm .
	sbatch matrix_pencil_no_plot.slurm
	cd ..
done
