cd /lustre/scratch/jmanthey/01b_certhia_genomics2/20_treemix_50kbp

for i in $( ls *simple.vcf ); do
  cat $i >> treemix_50kbp.vcf;
done

cd ../21_treemix_100kbp/

for i in $( ls *simple.vcf ); do
  cat $i >> treemix_100kbp.vcf;
done
