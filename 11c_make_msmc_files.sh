#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=msmc_make
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=24G
#SBATCH --array=1-525

module load intel R

input_array=$( head -n${SLURM_ARRAY_TASK_ID} msmc_list.txt | tail -n1 )

Rscript make_MSMC_files.r ${input_array}

