#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 1 

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=25G

#Save standard error and out to files:
#$ -e stderr_eQTl_01.file
#$ -o stdout_eQTL_01.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript 01_eQTL_subset_files.R \
	P140343-Results_FinalReport_clean_SNPs_autosome_individuals.A-transpose.matrixQTL.geno \
	BEST-D_phenotype_file_final.tsv  
