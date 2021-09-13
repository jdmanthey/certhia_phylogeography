start interactive session
interactive -p quanah

# move to directory
cd /lustre/scratch/jmanthey/01b_certhia_genomics2
cd 13_ds1_admix50kbp

# make one vcf
grep "#" Ca_0001_Tg_2.recode.vcf > total.vcf
for i in $( ls *recode.vcf ); do grep -v "#" $i >> total.vcf; done

# make chromosome map for the vcf
grep -v "#" total.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the combined vcf
vcftools --vcf total.vcf  --plink --chrom-map chrom_map.txt --out total 

# convert  with plink
plink --file total --recode12 --out total2 --noweb

# run admixture 
for K in 1 2 3 4 5 6 7; do admixture --cv total2.ped $K  | tee log_${K}.out; done

# check cv
grep -h CV log_*.out



# run for other thinning setting
cd ../14_ds2_admix100kbp

# make one vcf
grep "#" Ca_0001_Tg_2.recode.vcf > total.vcf
for i in $( ls *recode.vcf ); do grep -v "#" $i >> total.vcf; done

# make chromosome map for the vcf
grep -v "#" total.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom_map.txt

# run vcftools for the combined vcf
vcftools --vcf total.vcf  --plink --chrom-map chrom_map.txt --out total 

# convert  with plink
plink --file total --recode12 --out total2 --noweb

# run admixture 
for K in 1 2 3 4 5 6 7; do admixture --cv total2.ped $K  | tee log_${K}.out; done

# check cv
grep -h CV log_*.out




