#!/bin/bash

#SBATCH -p parallel
#SBATCH -t 96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --licenses=matlab@matlablm.unige.ch
#SBATCH --licenses=statistics_toolbox@matlablm.unige.ch
#SBATCH --licenses=distrib_computing_toolbox@matlablm.unige.ch
#SBATCH --mem=64000

BASE_MFILE_NAME=runSimFigure3Parallel
MATLAB_MFILE=${BASE_MFILE_NAME}.m

unset DISPLAY

module load matlab/2016b

echo "Starting at $(date)"
echo "Running ${MATLAB_MFILE} on $(hostname)"
srun matlab -nodesktop -nosplash \
-nodisplay -r "simoptions=${SLURM_ARRAY_TASK_ID}; ${BASE_MFILE_NAME}"
echo "Finished at $(date)"



