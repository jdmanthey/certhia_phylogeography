# for 50 kbp directory

options(scipen=999)

# list all the files in the trees directory
x_files <- list.files("windows", full.names=T)

# find the chromosome, start, and end for each tree
x_names <- list.files("windows")
x_chrom <- sapply(strsplit(sapply(strsplit(x_names, "RAxML_bipartitions."), "[[", 2), "__"), "[[", 1)
x_start <- sapply(strsplit(x_names, "__"), "[[", 2)
x_end <- sapply(strsplit(sapply(strsplit(x_names, "__"), "[[", 3), ".tre"), "[[", 1)

# write tree info
write.table(cbind(x_chrom, x_start, x_end), file="certhia_50kbp_tree_info.txt", sep="\t", quote=F, row.names=F, col.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="certhia_50kbp.trees", ncolumns=1)


# for 100 kbp directory

options(scipen=999)

# list all the files in the trees directory
x_files <- list.files("windows", full.names=T)

# find the chromosome, start, and end for each tree
x_names <- list.files("windows")
x_chrom <- sapply(strsplit(sapply(strsplit(x_names, "RAxML_bipartitions."), "[[", 2), "__"), "[[", 1)
x_start <- sapply(strsplit(x_names, "__"), "[[", 2)
x_end <- sapply(strsplit(sapply(strsplit(x_names, "__"), "[[", 3), ".tre"), "[[", 1)

# write tree info
write.table(cbind(x_chrom, x_start, x_end), file="certhia_100kbp_tree_info.txt", sep="\t", quote=F, row.names=F, col.names=F)

# trees into one file
tree_list <- list()
for(a in 1:length(x_files)) {
	tree_list[[a]] <- scan(x_files[a], what="character")
}
tree_list <- unlist(tree_list)
write(tree_list, file="certhia_100kbp.trees", ncolumns=1)
