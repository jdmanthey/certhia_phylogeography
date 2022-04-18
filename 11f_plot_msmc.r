library(MetBrewer)
# list all final files (to get names of all files)
x_files <- list.files(pattern="*__msmc.final.txt")
x_files_base <- sapply(strsplit(x_files, ".final.txt"), "[[", 1)
x_names <- x_files_base
x_names <- gsub("Chiapas", "C.America.North", x_names)
x_names <- gsub("CMexico", "C.Mexico", x_names)
x_names <- gsub("Honduras", "C.America.South", x_names)
x_names <- gsub("NuevoLeon", "S.M.Oriental", x_names)
x_names <- gsub("California", "Pacific", x_names)
x_names <- gsub("Utah", "Rocky.Mtns", x_names)
x_names <- gsub("WVirginia", "E.N.America", x_names)

# define parameters
gen <- 2
mu <- 2.506e-9 * gen # needs to be per generation
min_age <- 000
max_age <- 500000
plotting_age_range <- 1000 # years for x axis


# loop through each output directory
output <- list()
output_bootstraps <- list()
for(a in 1:length(x_files_base)) {
	
	# identify main output file
	a_main <- paste(x_files_base[a], ".final.txt", sep="")	
	# identify bootstrap outputs
	a_bootstraps <- paste(x_files_base[a], "_bs", seq(from=1, to=10, by=1), ".final.txt", sep="")

	# read in main file
	a_rep <- read.table(a_main, sep="\t", header=T)
	# rearrange main file for plotting lines
	for(d in 1:nrow(a_rep)) {
		if(d == 1) {
			a_rep2 <- rbind(c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		} else {
			a_rep2 <- rbind(a_rep2, c(a_rep[d,2], a_rep[d,4]), c(a_rep[d,3], a_rep[d,4]))
		}
	}
	a_rep <- a_rep2
	# scale by mutation rate
	a_rep[,1] <- a_rep[,1] / mu #don't multiply by generation because mu is in change per year
	a_rep[,2] <- (1 / a_rep[,2]) / (2 * mu)
	# remove very young and old time frames prone to error
	a_rep <- a_rep[a_rep[,1] >= min_age & a_rep[,1] <= max_age,]
	# scale by plotting age range and pop size range
	a_rep <- a_rep / plotting_age_range
	# add to output list
	output[[a]] <- a_rep
	
	# output for each bootstrap
	output_bootstraps[[a]] <- list()
	for(b in 1:length(a_bootstraps)) {
		b_rep <- read.table(a_bootstraps[b], sep="\t", header=T)
		# rearrange main file for plotting lines
		for(d in 1:nrow(b_rep)) {
			if(d == 1) {
				b_rep2 <- rbind(c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			} else {
				b_rep2 <- rbind(b_rep2, c(b_rep[d,2], b_rep[d,4]), c(b_rep[d,3], b_rep[d,4]))
			}
		}
		b_rep <- b_rep2
		# scale by mutation rate
		b_rep[,1] <- b_rep[,1] / mu #don't multiply by generation because mu is in change per year
		b_rep[,2] <- (1 / b_rep[,2]) / (2 * mu)
		# remove very young and old time frames prone to error
		b_rep <- b_rep[b_rep[,1] >= min_age & b_rep[,1] <= max_age,]
		# scale by plotting age range and pop size range
		b_rep <- b_rep / plotting_age_range
		# add to output list
		output_bootstraps[[a]][[b]] <- b_rep		
	}
}

# set plot items
a_col <- "darkgreen"
par(mfrow=c(4,6))
par(mar=c(1,1,2,0.2))
for(a in 1:length(output)) {
	plot_name1 <- sapply(strsplit(x_names[a], "_"), "[[", 3)
	plot_name2 <- sapply(strsplit(x_names[a], "_"), "[[", 4)
	plot_name <- paste(plot_name1, plot_name2)
	plot(c(-1,1), xlim=c(15, 200), ylim=c(0,1500), pch=19, cex=0.01, log="x", xlab="", ylab="", main="", xaxt="n", yaxt="n")
	title(main=plot_name, adj=0,line=0.5, cex.main=0.7)
	axis(side=2, at=c(0, 200, 400, 600, 800, 1000, 1200, 1400), labels=FALSE)
	axis(side=1, at=c(10, 20, 50, 100, 200), labels=F)		
	# plot bootstraps
	for(b in 1:length(output_bootstraps[[1]])) {
		lines(output_bootstraps[[a]][[b]][,1], output_bootstraps[[a]][[b]][,2], col=a_col, lwd=0.3)
	}
	lines(output[[a]][,1], output[[a]][,2], col=a_col, lwd=3)
}


# harmonic mean of pop sizes from most recent to 200k years ago
harmonic_pop <- c()
for(a in 1:length(output)) {
	out_rep <- output[[a]]
	
	# define time series
	time_series <- seq(from=as.integer(out_rep[2,1])+1, to=200, by=1)
	# time series pops
	time_pops <- c()
	for(b in 1:length(time_series)) {
		time_pops <- c(time_pops, out_rep[time_series[b] < out_rep[,1],][1,2])
	}
	# harmonic mean of this individual
	harm_rep <- length(time_pops) / sum((1 / time_pops))
	
	# add to output element
	harmonic_pop <- c(harmonic_pop, harm_rep)
}

harmonic_pop <- data.frame(id=as.character(x_files_base), harmonic_pop=as.numeric(harmonic_pop))
harmonic_pop

# plot groups of three from same population together
# define plot groupings and order (groups of 3)
plot_groups <- list()
plot_groups[[1]] <- c(13, 14, 15)
plot_groups[[2]] <- c(16, 17, 18)
plot_groups[[3]] <- c(19, 20, 21)
plot_groups[[4]] <- c(4, 5, 6)
plot_groups[[5]] <- c(10, 11, 12)
plot_groups[[6]] <- c(1, 2, 3)
plot_groups[[7]] <- c(7, 8, 9)
plot_group_names <- c("Pacific", "Rocky Mountains", "Eastern North America", "Central Mexico", "Sierra Madre Oriental", "Central America North", "Central America South")
# set plot items
a_cols <- met.brewer(name="Egypt", n=4, type="discrete")[c(1,2,4)]
a_cols2 <- c(rgb(221,81,41,max=255,alpha=200), rgb(15,123,162,max=255,alpha=200), rgb(250,178,85,max=255,alpha=200))
par(mfrow=c(1,7))
par(mar=c(1,1,2,0.2))
for(a in 1:length(plot_groups)) {
	plot_name1 <- sapply(strsplit(x_names[a], "_"), "[[", 3)
	plot_name2 <- sapply(strsplit(x_names[a], "_"), "[[", 4)
	plot_name <- plot_group_names[a]
	plot(c(-1,1), xlim=c(15, 200), ylim=c(0,1500), pch=19, cex=0.01, log="x", xlab="", ylab="", main="", xaxt="n", yaxt="n")
	title(main=plot_name, adj=0,line=0.5, cex.main=0.9)
	axis(side=2, at=c(0, 200, 400, 600, 800, 1000, 1200, 1400), labels=FALSE)
	axis(side=1, at=c(10, 20, 50, 100, 200), labels=F)
	# loop for all three
	for(d in 1:3) {
		d_rep <- plot_groups[[a]][d]
		# plot bootstraps
		for(b in 1:length(output_bootstraps[[1]])) {
			lines(output_bootstraps[[d_rep]][[b]][,1], output_bootstraps[[d_rep]][[b]][,2], col=a_cols[d], lwd=0.15)
		}
	}
	d_rep <- plot_groups[[a]][1]
	lines(output[[d_rep]][,1], output[[d_rep]][,2], col=a_cols2[1], lwd=3)
	d_rep <- plot_groups[[a]][2]
	lines(output[[d_rep]][,1], output[[d_rep]][,2], col=a_cols2[2], lwd=2.5)
	d_rep <- plot_groups[[a]][3]
	lines(output[[d_rep]][,1], output[[d_rep]][,2], col=a_cols2[3], lwd=2)
	
	
}











