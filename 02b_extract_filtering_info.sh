for i in $( ls slurm* ); do
echo $i;
head -n1 $i >> filtering.txt;
grep "Input:" $i >> filtering.txt;
grep "Result:" $i >> filtering.txt;
done
