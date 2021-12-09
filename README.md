## certhia_phylogeography

Steps to analyze resequencing data for Certhia americana samples from 7 populations. 

    01_setup.sh -> make directories for organizing all files throughout analyses
    02_align.sh -> filters, trims, and aligns data to reference genome
                    -> requires input basenames.txt to run 
                    -> all fastq files should be gzipped and in the formats: basename_R1.fastq.gz
                                                                             basename_R2.fastq.gz
    02b_extract_filtering_info.sh -> pulls out filtering info from bbduk slurm.out files
    03a_coverage.sh -> uses samtools to measure alignment coverage from bam files
    03b_plot_coverage.r -> uses R to plot the output from 03a_coverage.sh
    04_create_genotype_scripts.r -> R script that makes slurm array scripts and helper text files for genotyping in GATK
                    -> requires popmap.txt and 06_certhia_reordered.fasta.fai files for job creation
    05_concatenate_vcf_files.sh + 05b_concatenate_vcf_files.sh -> combine all vcfs from the same chromosome into single vcf files
    06_filter_vcf.sh -> filter the whole-chromosome vcf files for downstream analyses
    07_whole_genome_admixture.sh -> run whole-genome admixture analyses with thinned SNP datasets
    08a_phylo_stats_50kbp.r + 08b_phylo_stats_100kbp.r -> R scripts to build SLURM array jobs to estimate phylogenies and popgen summary statistics in sliding windows
                    ->subsequent phylogeny functions in SLURM jobs need the create_fasta_from_vcf.r and create_fasta.r R scripts
                    ->subsequent popgen functions in SLURM jobs need the calculate_windows.r and window_stat_calculations.r R scripts
    08c_combine_trees.r -> combine all raxml output .tre files into a single .trees file
    08d_count_sites_in_tree_files.sh -> count invariable and variable sites per vcf
    08d_species_trees.sh -> estimate species trees from gene trees
    08e_gsi.r -> calculate GSI from each of the gene trees
    09 sh and r scripts -> run LOSTRUCT
    10a_50kbp_admixture.r + 10b_100kbp_admixture.r -> R scripts that make SLURM array scripts for sliding window ADMIXTURE analyses
    
    
    
    
    
    
    
    
    
    
    
    
    
