# popmap = individual base names of fastq files, one line per individual

# make sure reference is indexed with bwa and samtools before use, and use CreateSequenceDictionary in GATK 
	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/01b_certhia_genomics2"
	directory_name <- "11_genotype_scripts_certhia"
	reference_genome_location <- "/home/jmanthey/references/06_certhia_reordered.fasta"
	cluster <- "quanah"
	output_name <- "certhia_genotype"
	popmap <- "popmap.txt"
	individuals <- read.table(popmap, sep="\t")
	faidx <- read.table("06_certhia_reordered.fasta.fai", stringsAsFactors=F)
	
	min_scaffold_size <- 1000000
	max_genotype_job_size <- 10000000
	max_individual_genotype_job_size <- 100000000
	

	# make directories
	dir.create(directory_name)
	dir.create(paste(directory_name, "/01_gatk_split", sep=""))
	dir.create(paste(directory_name, "/02b_gatk_database", sep=""))
	dir.create(paste(directory_name, "/03b_group_genotype_database", sep=""))

	# subset the index
	faidx <- faidx[faidx[,2] >= min_scaffold_size, ]

	# finds scaffolds too big to genotype at once
	faidx_keep <- faidx[faidx[,2] < max_individual_genotype_job_size,1:2]
	faidx_change <- faidx[faidx[,2] >= max_individual_genotype_job_size,1:2]
	# paste the interval to use to each of the faidx objects
	faidx_keep <- cbind(faidx_keep, faidx_keep[,1], rep(1, nrow(faidx_keep)), faidx_keep[,2])
	faidx_change <- cbind(faidx_change, rep("x", nrow(faidx_change)))
	new_faidx_change <- c()
	for(a in 1:nrow(faidx_change)) {
		a_rep <- faidx_change[a,]
		a_breaks <- floor(as.numeric(a_rep[1,2]) / 2)
		a_break1 <- c(a_rep[1,1], a_breaks, paste(a_rep[1,1], ":1-", a_breaks, sep=""), 1, a_breaks)
		a_break2 <- c(paste(a_rep[1,1], "b", sep=""), as.numeric(a_rep[1,2]) - a_breaks, paste(a_rep[1,1], ":", a_breaks + 1, "-", as.numeric(a_rep[1,2]), sep=""), a_breaks + 1, as.numeric(a_rep[1,2]))
		new_faidx_change <- rbind(new_faidx_change, a_break1, a_break2)
	}
	colnames(faidx_keep) <- c("id", "length", "interval", "start", "end")
	colnames(new_faidx_change) <- c("id", "length", "interval", "start", "end")
	faidx <- rbind(faidx_keep, new_faidx_change)
	faidx[,3] <- as.character(faidx[,3])
	faidx[,1] <- as.character(faidx[,1])
	faidx[,2] <- as.numeric(faidx[,2])
	faidx[,4] <- as.numeric(faidx[,4])
	faidx[,5] <- as.numeric(faidx[,5])
	faidx <- na.omit(faidx)
	
	# write the helper files for the 1st genotyping step
	for(a in 1:nrow(individuals)) {
		if(a == 1) {
			helper1 <- cbind(faidx[,1], as.character(faidx[,3]), rep(individuals[a,1], nrow(faidx)))
		} else {
			helper1 <- rbind(helper1, cbind(faidx[,1], as.character(faidx[,3]), rep(individuals[a,1], nrow(faidx))))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper1[,2], file=paste(directory_name, "/01_gatk_split", "/helper1.txt", sep=""), ncolumns=1)
	# write the individual to genotype for each job
	write(helper1[,3], file=paste(directory_name, "/01_gatk_split", "/helper2.txt", sep=""), ncolumns=1)
	# write the chromosome name for each job
	write(helper1[,1], file=paste(directory_name, "/01_gatk_split", "/helper1b.txt", sep=""), ncolumns=1)
	
	
	# step 1
	# genotype all individuals using GATK, one array job using the above two helper files

	a.script <- paste(directory_name, "/01_gatk_split/step1_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "step1", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=8", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", nrow(helper1), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} helper1.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("ind_array=$( head -n${SLURM_ARRAY_TASK_ID} helper2.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper1b.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	#gatk 4.0
	a_name <- paste(project_directory, "/01_bam_files/", "${ind_array}", "_final.bam", sep="")
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx64g" HaplotypeCaller -R ', reference_genome_location, " -I ", a_name, " -ERC GVCF -O ", project_directory, "/02_vcf/", "${name_array}", "._${ind_array}_.g.vcf", " --QUIET --intervals ", "${chr_array}", sep="")
	write(gatk_command, file=a.script, append=T)
	
	
	
	
	# write the helper files for the 2nd genotyping step
	# write the chromosome/scaffold to database for each job
	write(faidx[,1], file=paste(directory_name, "/02b_gatk_database", "/helper3.txt", sep=""), ncolumns=1)
	write(faidx[,3], file=paste(directory_name, "/02b_gatk_database", "/helper3b.txt", sep=""), ncolumns=1)

	
	
	# step 2
	# create genotyping database for each of the chromosomes
	a.script <- paste(directory_name, "/02b_gatk_database/step2_array.sh", sep="")

	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "step2", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=8", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", nrow(faidx), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper3.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} helper3b.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
		
	#make list of all vcfs to database
	for(b in 1:nrow(individuals)) {
		if(b == 1) {
			vcf_total <- paste("-V ", project_directory, "/02_vcf/", "${name_array}", "._", individuals[b,1], "_.g.vcf", sep="")
		} else {
			vcf_total <- paste(vcf_total, " -V ", project_directory, "/02_vcf/", "${name_array}", "._", individuals[b,1], "_.g.vcf", sep="")
		}
	}
	
	#gatk 4.0
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx60g" GenomicsDBImport ', vcf_total, " --genomicsdb-workspace-path ", project_directory, "/02_vcf/", "${name_array}", " -L ", "${interval_array}", sep="")
				
	write(gatk_command, file=a.script, append=T)
		




	# write the helper files for the 3rd genotyping step
	helper <- c()
	job_suffixes <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")
	for(a in 1:nrow(faidx)) {
		a_rep <- faidx[a,]
		job_number <- ceiling(a_rep[1,2] / max_genotype_job_size)
		job_number <- job_suffixes[1:job_number]
		
		if(length(job_number) > 1) {
			# define interval to group genotype if more than one job
			interval_start <- a_rep[1,4]
			interval_end <- a_rep[1,4] + max_genotype_job_size - 1
			for(b in 1:length(job_number)) {
				helper <- rbind(helper, c(faidx[a,1], 
					paste(strsplit(faidx[a,3], ":")[[1]][1], ":", interval_start, "-", interval_end, sep=""),
					paste(faidx[a,1], "__", job_number[b], sep="")
					))
					interval_start <- interval_start + max_genotype_job_size
				if((b+1) != length(job_number)) {
					interval_end <- interval_end + max_genotype_job_size
				} else {
					interval_end <- faidx[a,5]
				}
			}
		} else {
			helper <- rbind(helper, c(faidx[a,1], 
				paste(faidx[a,1], ":", 1, "-", faidx[a,2], sep=""),
				paste(faidx[a,1], sep="")
				))
		}
	}
	# write the chromosome/scaffold to genotype for each job
	write(helper[,1], file=paste(directory_name, "/03b_group_genotype_database", "/helper4.txt", sep=""), ncolumns=1)
	# write the interval range to genotype for each job
	write(helper[,2], file=paste(directory_name, "/03b_group_genotype_database", "/helper5.txt", sep=""), ncolumns=1)
	# write the output base name for each job
	write(helper[,3], file=paste(directory_name, "/03b_group_genotype_database", "/helper6.txt", sep=""), ncolumns=1)
	
	
	


	

	# step 3
	# group genotyping for each interval
	a.script <- paste(directory_name, "/03b_group_genotype_database/step3_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "step3", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=10", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", nrow(helper), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel java", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("chr_array=$( head -n${SLURM_ARRAY_TASK_ID} helper4.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("interval_array=$( head -n${SLURM_ARRAY_TASK_ID} helper5.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("name_array=$( head -n${SLURM_ARRAY_TASK_ID} helper6.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)

	#gatk 4.0
	gatk_command <- paste('/lustre/work/jmanthey/gatk-4.1.0.0/gatk --java-options "-Xmx80g" GenotypeGVCFs -R ', reference_genome_location, " -V gendb://", project_directory, "/02_vcf/", "${chr_array}", " --include-non-variant-sites -O ", project_directory, "/03_vcf/", "${name_array}", ".g.vcf", " -L ", "${interval_array}", sep="")
	write(gatk_command, file=a.script, append=T)
