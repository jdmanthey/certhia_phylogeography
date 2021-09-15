# use a simplified vcf as input for creating a treemix file with snps at least minimum_dist apart

options(scipen=999)

# set up names of individuals
ind_names <- c("Certhia_albescens_CMexico_UWBM106936","Certhia_albescens_CMexico_UWBM10789Certhia_albescens_CMexico_UWBM107897","Certhia_albescens_Chiapas_UWBM113447","Certhia_albescens_Chiapas_UWBM113461","Certhia_albescens_Chiapas_UWBM113462","Certhia_albescens_Honduras_UWBM103417","Certhia_albescens_Honduras_UWBM105618","Certhia_albescens_Honduras_UWBM93696","Certhia_albescens_NuevoLeon_UWBM106768","Certhia_albescens_NuevoLeon_UWBM111691","Certhia_albescens_NuevoLeon_UWBM111708","Certhia_americana_California_UWBM112794","Certhia_americana_California_UWBM112795","Certhia_americana_California_UWBM113198","Certhia_americana_Utah_UWBM113162","Certhia_americana_Utah_UWBM113167","Certhia_americana_Utah_UWBM113168","Certhia_americana_WVirginia_UWBM107061","Certhia_americana_WVirginia_UWBM112054","Certhia_americana_WVirginia_UWBM112055","Certhia_familiaris_outgroup_KU92846")

# set up population names
populations_labels <- sapply(strsplit(ind_names, "_"), "[[", 3)
populations_unique <- unique(populations_labels)
populations_unique <- populations_unique[c(5,6,7,1,4,2,3,8)] # order north to south

# write output 50 kbp 
write(populations_unique, file="certhia_50kbp.treemix", sep=" ", ncolumns=8)

# read in simple vcf
x <- read.table("treemix_50kbp.vcf")

for(a in 1:nrow(x)) {
	allele_freq_list <- c()
	a_rep <- as.character(x[a,5:ncol(x)])
	a_rep <- gsub("\\|", "/", a_rep)
	# loop for each population
	for(b in 1:length(populations_unique)) {
		# get this population's genotypes
		b_rep <- a_rep[populations_labels == populations_unique[b]]
		alt_allele <- length(b_rep[b_rep == "1/1"]) * 2 + length(b_rep[b_rep == "0/1"]) * 1
		ref_allele <- length(b_rep) * 2 - alt_allele
		allele_freq_list <- c(allele_freq_list, paste(ref_allele, alt_allele, sep=","))
	}
	write(allele_freq_list, file="certhia_50kbp.treemix", sep=" ", ncolumns=8, append=T)
}


# write output 100 kbp
write(populations_unique, file="certhia_100kbp.treemix", sep=" ", ncolumns=8)

# read in simple vcf
x <- read.table("treemix_100kbp.vcf")

for(a in 1:nrow(x)) {
	allele_freq_list <- c()
	a_rep <- as.character(x[a,5:ncol(x)])
	a_rep <- gsub("\\|", "/", a_rep)
	# loop for each population
	for(b in 1:length(populations_unique)) {
		# get this population's genotypes
		b_rep <- a_rep[populations_labels == populations_unique[b]]
		alt_allele <- length(b_rep[b_rep == "1/1"]) * 2 + length(b_rep[b_rep == "0/1"]) * 1
		ref_allele <- length(b_rep) * 2 - alt_allele
		allele_freq_list <- c(allele_freq_list, paste(ref_allele, alt_allele, sep=","))
	}
	write(allele_freq_list, file="certhia_100kbp.treemix", sep=" ", ncolumns=8, append=T)
}



