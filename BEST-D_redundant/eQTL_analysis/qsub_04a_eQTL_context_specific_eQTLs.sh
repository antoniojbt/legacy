#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 2

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=10G

#Save standard error and out to files:
#$ -e stderr_eQTl_context_specific.file
#$ -o stdout_eQTL_context_specific.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript /ifs/devel/antoniob/projects/BEST-D/04a_eQTL_context_specific_eQTLs.R \
$EQTL_1 \
$EQTL_2 \
$GENES_INTEREST \
$SNPS_INTEREST \
$total_SNPs_tested \
$total_probes_tested \
$total_SNPprobe_pairs_tested
