#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 10

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=20G

#Save standard error and out to files:
#$ -e stderr_qfam.file
#$ -o stdout_qfam.file

#Run the job - program must have full path:
/ifs/apps/bio/plink-1.07/plink-1.07-x86_64/plink \
--noweb --bfile all_awk_loose_processed_filtered_MAF_PASS_IDs \
--allow-no-sex --pheno metadata_final_pseudomonas.tsv --missing-phenotype -99999 \
--pheno-name $pheno_name $cmd --out $outfile 
