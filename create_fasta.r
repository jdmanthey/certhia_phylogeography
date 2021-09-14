# input args <- input file, popmap
args <- commandArgs(trailingOnly = TRUE)

# add functions for calculations
source("create_fasta_from_vcf.r")

# define minimum number of sites to keep a fasta file 
min_sites <- 10000

# no scientific notation
options(scipen=999)

# read in input file
input_file <- read.table(args[1], sep="\t", stringsAsFactors=F)

# read in populations
populations <- read.table(args[2], stringsAsFactors=F, header=T)

# define output fasta name
output_fasta <- paste(strsplit(args[1], ".simple")[[1]][1], ".fasta", sep="")

# create fasta sequence alignments
create_fasta_from_vcf(input_file, populations, output_fasta, min_sites)
