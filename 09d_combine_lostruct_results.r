options(scipen=999)

# for 50 kbp
x_pca_files <- list.files("lostruct_certhia_50kbp", pattern="*pca.csv")
x_mds <- read.table("lostruct_certhia_50kbp/mds_coords.csv", sep=",", header=T, stringsAsFactors=F)
x_mds <- na.omit(x_mds)
# loop to read in and combine pca and regions files
total_pca <- c()
for(a in 1:length(x_pca_files)) {
	a_pca <- read.table(paste0("lostruct_certhia_50kbp/", x_pca_files[a]), sep=",", header=T, stringsAsFactors=F)
	a_regions <- read.table(paste0("lostruct_certhia_50kbp/", strsplit(x_pca_files[a], ".pca")[[1]][1], ".regions.csv"), sep=",", header=T, stringsAsFactors=F)
	
	total_pca <- rbind(total_pca, cbind(a_regions, a_pca[,4:ncol(a_pca)]))	
}

# combine mds and pca files
x_mds <- cbind(x_mds, total_pca)

write.table(x_mds, file="b_lostruct_50kbp.txt", sep="\t", quote=F, row.names=F)

# for 100 kbp
x_pca_files <- list.files("lostruct_certhia_100kbp", pattern="*pca.csv")
x_mds <- read.table("lostruct_certhia_100kbp/mds_coords.csv", sep=",", header=T, stringsAsFactors=F)
x_mds <- na.omit(x_mds)
# loop to read in and combine pca and regions files
total_pca <- c()
for(a in 1:length(x_pca_files)) {
	a_pca <- read.table(paste0("lostruct_certhia_100kbp/", x_pca_files[a]), sep=",", header=T, stringsAsFactors=F)
	a_regions <- read.table(paste0("lostruct_certhia_100kbp/", strsplit(x_pca_files[a], ".pca")[[1]][1], ".regions.csv"), sep=",", header=T, stringsAsFactors=F)
	
	total_pca <- rbind(total_pca, cbind(a_regions, a_pca[,4:ncol(a_pca)]))	
}

# combine mds and pca files
x_mds <- cbind(x_mds, total_pca)

write.table(x_mds, file="b_lostruct_100kbp.txt", sep="\t", quote=F, row.names=F)
