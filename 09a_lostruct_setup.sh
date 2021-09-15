interactive -p quanah

cd /lustre/scratch/jmanthey/01b_certhia_genomics2/17_ds5_lostruct

# this process adds an invariant site at position 1 for each chromosome so that the windows may begin at 1
# then gzips the files
for i in $( ls *recode.vcf ); do

	echo $i;
	
	# get the first genotypes line from the file
	grep -v "#" $i | head -n1 > temp2.vcf

	# replace any non-reference alleles with the 0/0
	sed 's/0\/1/0\/0/g' temp2.vcf > temp3.vcf
	sed 's/1\/1/0\/0/g' temp3.vcf > temp4.vcf

	# set the position of the row to 1
	awk '{$2=1 ; OFS="\t"; print ;}' temp4.vcf > temp5.vcf

	# add the header to the new file
	grep "#" $i > ${i%recode.vcf}vcf

	# add the single line to the new file
	cat temp5.vcf >> ${i%recode.vcf}vcf

	# add the rest of the genotype information to the new file
	grep -v "#" $i >> ${i%recode.vcf}vcf

	# convert to bcf and index
	bcftools convert -O b ${i%recode.vcf}vcf > ${i%recode.vcf}bcf
	bcftools index ${i%recode.vcf}bcf

	# remove temp files
	rm temp2.vcf
	rm temp3.vcf
	rm temp4.vcf
	rm temp5.vcf
	rm ${i%recode.vcf}vcf

done

#### copy the gzipped files to a local directory named data to run the next R scripts
