# create msmc files for each individual, one per chromosome

# input args <- input file (only one arg after script name)
args <- commandArgs(trailingOnly = TRUE)

# no scientific notation
options(scipen=999)

# read in file for this array job
x <- read.table(args[1], stringsAsFactors=F)

# determine chromosome name and output name
chrom_name <- strsplit(strsplit(args[1], "__")[[1]][2], ".vcf")[[1]][1]
output_name <- gsub("__", "/", paste0(substr(args[1], 1, nchar(args[1]) - 4), ".txt"))

# remove sites with missing data 
x <- x[x[,4] != "./.", ]

# replace all phased genotypes to keep it simple
x[,4] <- gsub("\\|", "/", x[,4])

# subset only the variant sites (only biallelic snps kept, so 0/1 genotype)
x_matches <- x[,4] == "0/1"
x_variant <- x[x_matches, ]

# determine number of sites before each match
x_sites <- seq(from=1, to=nrow(x), by=1)
x_sites <- c(0, x_sites[x_matches])
x_sites <- diff(x_sites)

# output
output <- data.frame(
				chrom=as.character(rep(chrom_name, length(x_sites))),
				position=as.numeric(x_variant[,1]),
				sites_called=as.numeric(x_sites),
				genotypes=as.character(paste0(x_variant[,2], x_variant[,3]))
				)

write.table(output, file=output_name, quote=F, row.names=F, col.names=F)



