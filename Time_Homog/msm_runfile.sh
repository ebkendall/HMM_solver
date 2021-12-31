#!/bin/bash

#SBATCH -o Model_out/msm/trial%a.out
#SBATCH --array=1-50
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --account=jantonelli
#SBATCH --qos=jantonelli-b
#SBATCH --job-name=homog_msm
#SBATCH --time=01:00:00
#SBATCH -t 4000
#SBATCH --mem=5gb

export OMP_NUM_THREADS=5

module load R/4.0

R CMD BATCH --no-save msm_runfile.r Model_out/msm/trial${SLURM_ARRAY_TASK_ID}.Rout
