#!/bin/bash

#SBATCH -o DataOut/a_out/trial%a.out
#SBATCH --array=1-100
#SBATCH --account=jantonelli
#SBATCH --qos=jantonelli-b
#SBATCH --job-name=cavRun
#SBATCH --time=01:00:00
#SBATCH -t 4000
#SBATCH --mem=1gb

module load R/4.0

R CMD BATCH --no-save cav_simulate.r DataOut/a_out/trial${SLURM_ARRAY_TASK_ID}.Rout
