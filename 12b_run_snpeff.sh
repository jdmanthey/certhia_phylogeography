#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=snpeff1
#SBATCH --nodes=1 --ntasks=2
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-25

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )
input_array=${input_array%.filtered.vcf}

# define main working directory
workdir=/lustre/scratch/jmanthey/01b_certhia_genomics2

# run vcftools with SNP and invariant site output, 50% max missing data, no indels
vcftools --vcf ${workdir}/03_vcf/${input_array}.filtered.vcf --max-missing 0.5 --minGQ 20 --minDP 6 --max-meanDP 50 --max-alleles 2 --mac 1 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/19_snpeff/${input_array}

# run snpEff for this chromosome
java -Xmx16g -jar ~/snpEff/snpEff.jar creeper ${workdir}/19_snpeff/${input_array}.recode.vcf > ${workdir}/19_snpeff/${input_array}.ann.vcf

# remove intergenic, intron, upstream, and downstream regions (i.e., only cds regions)
grep -v "intergenic" ${workdir}/19_snpeff/${input_array}.ann.vcf | grep -v "intron" | grep -v "upstream" | grep -v "downstream" | grep -v "^#" > ${workdir}/19_snpeff/${input_array}.ann2.vcf









