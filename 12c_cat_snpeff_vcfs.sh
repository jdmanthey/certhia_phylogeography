
grep "^#" Ca_0029_Tg_19.ann.vcf > Certhia_cds_snpeff.vcf
for i in $( ls *ann2* ); do grep -v "^#" $i >> Certhia_cds_snpeff.vcf; done

# remove sites with warnings from snpEff
grep -v "WARN" Certhia_cds_snpeff.vcf > Certhia_cds_snpeff2.vcf

# remove the "#" from the header line for reading into R
