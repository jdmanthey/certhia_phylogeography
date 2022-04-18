#!/bin/bash
#SBATCH --chdir=./
#SBATCH --job-name=msmc
#SBATCH --partition quanah
#SBATCH --nodes=1 --ntasks=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-21

# move to working directory
cd /lustre/scratch/jmanthey/01b_certhia_genomics2/18_msmc

# define variables needed
name_array=$( head -n${SLURM_ARRAY_TASK_ID} msmc_helper.txt | tail -n1 | cut -f1 )

het_array=$( head -n${SLURM_ARRAY_TASK_ID} msmc_helper.txt | tail -n1 | cut -f3 )

output_array=${name_array}__msmc

# change to this individual's directory
cd ${name_array}

# run msmc2 for the full dataset
~/msmc2_linux64bit -o ${output_array} -i 20 -t 2 -m ${het_array} -p 1*2+20*1+1*2+1*3 *txt

# run ten bootstraps
~/multihetsep_bootstrap.py -n 10 -s 1000000 --chunks_per_chromosome 36 --nr_chromosomes 25 --seed 324324 bootstrap *[0-9].txt;


# run each of the ten bootstraps
for i in {1..10}; do
   echo bootstrap_${i};
   # define input files
   boot_input=/lustre/scratch/jmanthey/01b_certhia_genomics2/18_msmc/${name_array}/bootstrap_${i}/*txt;
   # run msmc for this bootstrap
   ~/msmc2_linux64bit -o ${output_array}_bs${i} -i 20 -t 2 -m ${het_array} -p 1*2+20*1+1*2+1*3 ${boot_input};
done

