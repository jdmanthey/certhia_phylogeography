for i in {1..21}; do 

dir_array=$( head -n${i} cut_popmap.txt | tail -n1 | cut -f1 )

mkdir $dir_array

done
