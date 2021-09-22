for i in $( ls *vcf.gz ); do 

  echo $i

  gunzip -cd $i | grep -v "#" | cut -f 5 >> alt_alleles.txt

done

# count total number of sites
grep "\\." alt_alleles.txt | wc -l
# count number of polymorphic sites
grep -v "\\." alt_alleles.txt | wc -l

