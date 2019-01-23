#!/bin/bash -l

#SBATCH --nodes=1 --ntasks-per-node 40
#SBATCH --time=24:00:00
#SBATCH --mail-user=kmylonakis@ucsb.edu
#SBATCH --mail-type=start,end

cd $SLURM_SUBMIT_DIR

module loade MatLab/R2018b

matlab -nodesktop -nojvm -nodisplay -nosplash -r matrix_pencil_no_plot > log.txt </dev/null
