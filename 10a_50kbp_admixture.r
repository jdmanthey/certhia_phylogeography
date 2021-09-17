	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/01b_certhia_genomics2/18_ds6_stats50kbp"
	directory_name <- "stat_50kbp"
	cluster <- "quanah"
	max_number_jobs <- 400
	
	# read in reference index
	# filtered to only include genotyped chromosomes
	ref_index <- read.table("06_certhia_reordered.fasta.fai", stringsAsFactors=F)
	
	# define window size
	window_size <- 50000
	
	# make directories
	dir.create(directory_name)
	
	# define intervals and write to helper files
	stat_helper1 <- list()
	stat_helper2 <- list()
	stat_helper3 <- list()
	counter <- 1
	for(a in 1:nrow(ref_index)) {
		a_start <- 1
		a_end <- a_start + window_size - 1
		a_max <- ref_index[a,2]
		a_windows <- ceiling((a_max - a_start) / window_size)
		a_chromosome <- ref_index[a,1]
		
		# loop for defining helper info for each window
		for(b in 1:a_windows) {
			if(b == a_windows) {
				a_end <- a_max
			}
			stat_helper1[[counter]] <- a_chromosome
			stat_helper2[[counter]] <- a_start
			stat_helper3[[counter]] <- a_end

			a_start <- a_start + window_size
			a_end <- a_end + window_size
			counter <- counter + 1
		}
	}
	stat_helper1 <- unlist(stat_helper1)
	stat_helper2 <- unlist(stat_helper2)
	stat_helper3 <- unlist(stat_helper3)
	
	# calculate number of array jobs
	if(length(stat_helper3) > max_number_jobs) {
		n_jobs_per_array <- ceiling(length(stat_helper3) / max_number_jobs)
		n_array_jobs <- ceiling(length(stat_helper3) / n_jobs_per_array)
	} else {
		n_array_jobs <- length(stat_helper3)
		n_jobs_per_array <- 1
	}
	
	stat_helper1 <- c(stat_helper1, rep("x", n_jobs_per_array - length(stat_helper1) %% n_jobs_per_array))
	stat_helper2 <- c(stat_helper2, rep(1, n_jobs_per_array - length(stat_helper2) %% n_jobs_per_array))
	stat_helper3 <- c(stat_helper3, rep(1, n_jobs_per_array - length(stat_helper3) %% n_jobs_per_array))
	length(stat_helper3)
	write(stat_helper1, file=paste(directory_name, "/stat_helper_chrom.txt", sep=""), ncolumns=1)
	write(stat_helper2, file=paste(directory_name, "/stat_helper_start.txt", sep=""), ncolumns=1)
	write(stat_helper3, file=paste(directory_name, "/stat_helper_end.txt", sep=""), ncolumns=1)

	# write the array script
	a.script <- paste(directory_name, "/stat50kbp_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "stat", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=1", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", n_array_jobs, sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("module load intel R", file=a.script, append=T)
	write("", file=a.script, append=T)

	write("# Set the number of runs that each SLURM task should do", file=a.script, append=T)
	write(paste("PER_TASK=", n_jobs_per_array, sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("# Calculate the starting and ending values for this task based", file=a.script, append=T)
	write("# on the SLURM task and the number of runs per task.", file=a.script, append=T)
	write("START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))", file=a.script, append=T)
	write("END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("# Print the task and run range", file=a.script, append=T)
	write("echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM", file=a.script, append=T)
	write("", file=a.script, append=T)

	write("# Run the loop of runs for this task.", file=a.script, append=T)	
	write("for (( run=$START_NUM; run<=$END_NUM; run++ )); do", file=a.script, append=T)
	write("\techo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("\tchrom_array=$( head -n${run} stat_helper_chrom.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("\tstart_array=$( head -n${run} stat_helper_start.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("\tend_array=$( head -n${run} stat_helper_end.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	# add header to output file
	header <- paste('\tgunzip -cd ', project_directory, '/Ca_0029_Tg_19.recode.vcf.gz | grep "#" > ', project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep="")
	write(header, file=a.script, append=T)
	write("", file=a.script, append=T)
	
	#tabix command
	tabix_command <- paste("\ttabix ", project_directory, "/${chrom_array}.recode.vcf.gz ${chrom_array}:${start_array}-${end_array} >> ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep="")
	write(tabix_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	
	# bcftools command
	bcf_tools_command <- paste("\tbcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf > ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf", sep="")
	write(bcf_tools_command, file=a.script, append=T)
	write("", file=a.script, append=T)

	# ADMIXTURE commands
	# make chromosome map
	admixture_commands <- paste0('\tgrep -v "#" ', project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf | cut -f 1 | uniq | awk '{print $0\"\\t\"$0}' > ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}_chrom_map.txt")
	write(admixture_commands, file=a.script, append=T)
	write("", file=a.script, append=T)
	# convert to plink
	admixture_commands <- paste0("\tvcftools --vcf ",project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf --plink --chrom-map ",project_directory, "/windows/${chrom_array}__${start_array}__${end_array}_chrom_map.txt --out ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}")
	write(admixture_commands, file=a.script, append=T)
	write("", file=a.script, append=T)
	# reformat in plink
	admixture_commands <- paste0("\tplink --file ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array} --recode12 --out ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b --noweb")
	write(admixture_commands, file=a.script, append=T)
	write("", file=a.script, append=T)
	# run admixture
	admixture_commands <- paste0("\tfor K in 2 7; do admixture -C 0.000001 --cv ",project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b.ped $K; done")
	write(admixture_commands, file=a.script, append=T)
	write("", file=a.script, append=T)
		
	# remove unnecessary files at end
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.recode.vcf", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.simple.vcf", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.ped", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}.map", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}_chrom_map.txt", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b.ped", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b.nosex", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b.map", sep=""), file=a.script, append=T)
	write(paste("\trm ", project_directory, "/windows/${chrom_array}__${start_array}__${end_array}__b.log", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	# finish
	write("done", file=a.script, append=T)
	
