#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter2
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

# define input files from vcf list
vcf_list=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )
input_array=${vcf_list%.g.vcf}

# define main working directory
workdir=/lustre/scratch/jmanthey/01b_certhia_genomics2

# pull out header and add to filtered vcf file
grep "#" ${workdir}/03_vcf/${input_array}.g.vcf > ${workdir}/03_vcf/${input_array}.filtered.vcf

# filter our rows that have low quality filters, genotyped sites with quality less than 20, and null alleles (* in col 4)
grep -v "#" ${workdir}/03_vcf/${input_array}.g.vcf | grep -v "LowQual" | awk '$6 >= 20 || $6 ~ /^\./' | awk '$5 !~ /*/' >> ${workdir}/03_vcf/${input_array}.filtered.vcf

# run vcftools with SNP output, only biallelic SNPs and invariant sites and no outgroup 
vcftools --vcf ${workdir}/03_vcf/${input_array}.filtered.vcf --remove-indv Certhia_familiaris_outgroup_KU92846 --minDP 6 --max-meanDP 50 --max-missing 0.01 --max-alleles 2 --max-maf 0.49 --remove-indels --recode --recode-INFO-all --out ${workdir}/18_msmc/${input_array}

# use bcftools to simplify the vcf files to reduce file size, complexity, and make them easier to work with 
bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' ${workdir}/18_msmc/${input_array}.recode.vcf > ${workdir}/18_msmc/${input_array}.simple.vcf

