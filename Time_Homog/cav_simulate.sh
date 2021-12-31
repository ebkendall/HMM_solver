#!/bin/bash

#SBATCH -o Data/sim%a.out
#SBATCH --array=1-50
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --account=jantonelli
#SBATCH --qos=jantonelli-b
#SBATCH --job-name=cavSim
#SBATCH --time=01:00:00
#SBATCH -t 4000
#SBATCH --mem=5gb

export OMP_NUM_THREADS=5

module load R/4.0

R CMD BATCH --no-save cav_simulate.r Data/sim${SLURM_ARRAY_TASK_ID}.Rout
