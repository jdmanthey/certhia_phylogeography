
x_files <- list.files(pattern="*txt")
sample_names <- substr(x_files, 1, nchar(x_files)-10)

coverage_max_to_check <- 50
plot_max <- coverage_max_to_check + 2

output1 <- list()
output2 <- list()
mean_coverage <- c()
for(a in 1:length(x_files)) {
	a_rep <- scan(x_files[a])
	a_output <- c()
	for(b in 1:plot_max) {
		if(b == plot_max) {
			a_output <- c(a_output, length(a_rep[a_rep >= (b - 1)]))
		} else {
			a_output <- c(a_output, length(a_rep[a_rep == (b - 1)]))
		}
	}
	a_output2 <- a_output / sum(a_output)
	output1[[a]] <- a_output
	output2[[a]] <- a_output2
	mean_coverage <- c(mean_coverage, mean(a_rep))
}

rm(a_rep)
save.image("certhia2_coverage.RData")







load("certhia2_coverage.RData")

par(mfrow=c(4,6))
for(a in 1:length(output2)) {
	plot(0:51, output2[[a]], pch=19, cex=0.1, xlab="Coverage", ylab="Proportion C. americana genome", main=sample_names[a], ylim=c(0,0.14))
	poly.plot <- rbind(cbind(0:51, output2[[a]]), c(51, 0), c(0,0))
	polygon(poly.plot, col="gray")
	abline(v=mean_coverage[a], col="red")
}
