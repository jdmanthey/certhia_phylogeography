cd 
cd snpEff
mkdir data
cd data
mkdir genomes
mkdir creeper
cd creeper
cp ~/references/06_certhia_reordered.fasta .
mv 06_certhia_reordered.fasta 
mv 06_certhia_reordered.fasta sequences.fa
mv 06_certhia_reordered.all.maker.proteins.fasta protein.fa
mv certhia_genes.maker.gff genes.gff
nano ../../snpEff.config
### add the following lines to the config file at the beginning

data.dir = ~/snpEff/data/

# Certhia
creeper.genome : creeper
creeper.reference : ~/snpEff/data/creeper

# then build the database
java -jar snpEff.jar build -gff3 -v creeper -noCheckCds
