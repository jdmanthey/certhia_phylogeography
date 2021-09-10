#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N ca_depth
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=24G

samtools depth -a Certhia_albescens_Chiapas_UWBM113447_final.bam Certhia_albescens_Chiapas_UWBM113461_final.bam \
Certhia_albescens_Chiapas_UWBM113462_final.bam Certhia_albescens_CMexico_UWBM106936_final.bam Certhia_albescens_CMexico_UWBM107890_final.bam \
Certhia_albescens_CMexico_UWBM107897_final.bam Certhia_albescens_Honduras_UWBM103417_final.bam Certhia_albescens_Honduras_UWBM105618_final.bam \
Certhia_albescens_Honduras_UWBM93696_final.bam Certhia_albescens_NuevoLeon_UWBM106768_final.bam Certhia_albescens_NuevoLeon_UWBM111691_final.bam \
Certhia_albescens_NuevoLeon_UWBM111708_final.bam Certhia_americana_California_UWBM112794_final.bam Certhia_americana_California_UWBM112795_final.bam \
Certhia_americana_California_UWBM113198_final.bam Certhia_americana_Utah_UWBM113162_final.bam Certhia_americana_Utah_UWBM113167_final.bam \
Certhia_americana_Utah_UWBM113168_final.bam Certhia_americana_WVirginia_UWBM107061_final.bam Certhia_americana_WVirginia_UWBM112054_final.bam \
Certhia_americana_WVirginia_UWBM112055_final.bam Certhia_familiaris_outgroup_KU92846_final.bam > \
certhia_coverage.txt

# break up the depth files into single column files for each individual (locations dropped)

while read -r name1 number1; do
	number2=$((number1 + 2));
  cut certhia_coverage.txt -f $number2 > ${name1}_depth.txt;
done < popmap.txt
