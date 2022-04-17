#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=cut
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-21

output_array=$( head -n${SLURM_ARRAY_TASK_ID} cut_popmap.txt | tail -n1 | cut -f1 )

cut_array=$( head -n${SLURM_ARRAY_TASK_ID} cut_popmap.txt | tail -n1 | cut -f2 )

for i in $( ls *simple.vcf ); do
output_name=${output_array}__${i%.simple.vcf}.vcf

cut $i -f1,2,3,${cut_array} > $output_name

done
