
grep "^#" Ca_0029_Tg_19.ann.vcf > Certhia_cds_snpeff.vcf
for i in $( ls *ann2* ); do grep -v "^#" $i >> Certhia_cds_snpeff.vcf; done

