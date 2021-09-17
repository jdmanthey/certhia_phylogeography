# have all phylogenies in one file
# have a file with the chromosome and window for each tree in the correct order
# have a popmap with the individual names and groupings wished to test
# have the outgroup labeled once in the popmap as "outgroup"

library(ape)
library(phytools)

options(scipen=999)

# read in trees, info, and popmap for 50kbp trees
x_trees <- read.tree("certhia_50kbp.trees")
x_info <- read.table("certhia_50kbp_tree_info.txt", sep="\t", stringsAsFactors=F)
x_popmap <- read.table("gsi_popmap.txt", sep="\t", stringsAsFactors=F)
x_output <- "b_gsi_50kbp_output.txt"

# define outgroup
outgroup <- x_popmap[x_popmap[,2] == "outgroup", 1]

# remove outgroup from popmap
x_popmap <- x_popmap[x_popmap[,2] != "outgroup", ]

# write initial line of output
write(c("chr", "start", "end", "pop", "gsi"), file=x_output, sep="\t", ncolumns=5)

# loop for each tree
for(a in 1:length(x_trees)) {
	# select the tree
	x <- x_trees[a][[1]]
	# reroot
	x <- midpoint.root(x)
	x <- root(x, outgroup, resolve.root=T)
	
	# define the groups
	groups <- unique(x_popmap[,2])
	
	# loop for each group of interest
	for(g in 1:length(groups)) {
		# define the group of interest
		group_of_interest <- x_popmap[x_popmap[,2] == groups[g],1]
		
		# calculate MRCA of group of interest
		group_mrca <- getMRCA(x, tip = group_of_interest)
		
		# calculate number of ndoes to reach common ancestor for each individual in group
		# and then remove redundant nodes
		nodes_needed <- c()
		for(b in 1:length(group_of_interest)) {
			# what is the number of this individual?
			b_number <- match(group_of_interest[b], x$tip.label)
				
			# loop throup edge table until reaching MRCA
			while_loop <- 0
			while(while_loop != group_mrca) {
				# add node
				nodes_needed <- c(nodes_needed, x$edge[x$edge[,2]==b_number,1])
		
				# change new number to that node and update the while_loop object to the node
				b_number <- x$edge[x$edge[,2]==b_number,1]
				while_loop <- b_number
			}
			# add the MRCA to the nodes_needed object
			nodes_needed <- c(nodes_needed, group_mrca)
	
			# only keep unique nodes
			nodes_needed <- unique(nodes_needed)
		}

		# calculate gsi

		# gs 
		gs <- (length(group_of_interest) - 1) / length(nodes_needed)

		# max gs = 1
		max_gs <- 1

		# min gs = minimum number of nodes to connect all individuals (n - 1) / total number of nodes
		min_gs <- (length(group_of_interest) - 1) / length(unique(x$edge[,1]))

		# equation 4 of Cummings et al. 2008 (GSI paper)
		gsi <- (gs  - min_gs) / (max_gs - min_gs)

		# write output
		write(c(x_info[a,1], x_info[a,2], x_info[a,3], groups[g], gsi), file=x_output, sep="\t", ncolumns=5, append=T)
	}
}










# read in trees, info, and popmap for 100kbp trees
x_trees <- read.tree("certhia_100kbp.trees")
x_info <- read.table("certhia_100kbp_tree_info.txt", sep="\t", stringsAsFactors=F)
x_popmap <- read.table("gsi_popmap.txt", sep="\t", stringsAsFactors=F)
x_output <- "b_gsi_100kbp_output.txt"

# define outgroup
outgroup <- x_popmap[x_popmap[,2] == "outgroup", 1]

# remove outgroup from popmap
x_popmap <- x_popmap[x_popmap[,2] != "outgroup", ]

# write initial line of output
write(c("chr", "start", "end", "pop", "gsi"), file=x_output, sep="\t", ncolumns=5)

# loop for each tree
for(a in 1:length(x_trees)) {
	# select the tree
	x <- x_trees[a][[1]]
	# reroot
	x <- midpoint.root(x)
	x <- root(x, outgroup, resolve.root=T)
	
	# define the groups
	groups <- unique(x_popmap[,2])
	
	# loop for each group of interest
	for(g in 1:length(groups)) {
		# define the group of interest
		group_of_interest <- x_popmap[x_popmap[,2] == groups[g],1]
		
		# calculate MRCA of group of interest
		group_mrca <- getMRCA(x, tip = group_of_interest)
		
		# calculate number of ndoes to reach common ancestor for each individual in group
		# and then remove redundant nodes
		nodes_needed <- c()
		for(b in 1:length(group_of_interest)) {
			# what is the number of this individual?
			b_number <- match(group_of_interest[b], x$tip.label)
				
			# loop throup edge table until reaching MRCA
			while_loop <- 0
			while(while_loop != group_mrca) {
				# add node
				nodes_needed <- c(nodes_needed, x$edge[x$edge[,2]==b_number,1])
		
				# change new number to that node and update the while_loop object to the node
				b_number <- x$edge[x$edge[,2]==b_number,1]
				while_loop <- b_number
			}
			# add the MRCA to the nodes_needed object
			nodes_needed <- c(nodes_needed, group_mrca)
	
			# only keep unique nodes
			nodes_needed <- unique(nodes_needed)
		}

		# calculate gsi

		# gs 
		gs <- (length(group_of_interest) - 1) / length(nodes_needed)

		# max gs = 1
		max_gs <- 1

		# min gs = minimum number of nodes to connect all individuals (n - 1) / total number of nodes
		min_gs <- (length(group_of_interest) - 1) / length(unique(x$edge[,1]))

		# equation 4 of Cummings et al. 2008 (GSI paper)
		gsi <- (gs  - min_gs) / (max_gs - min_gs)

		# write output
		write(c(x_info[a,1], x_info[a,2], x_info[a,3], groups[g], gsi), file=x_output, sep="\t", ncolumns=5, append=T)
	}
}
