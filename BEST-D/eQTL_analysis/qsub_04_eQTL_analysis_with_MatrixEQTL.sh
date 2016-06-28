#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 2

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=20G

#Save standard error and out to files:
#$ -e stderr_eQTl_04.file
#$ -o stdout_eQTL_04.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript /ifs/devel/antoniob/projects/BEST-D/04_eQTL_analysis_with_MatrixEQTL.R $geno $expr $covar $snp_pos $probe_pos
