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
#$ -e stderr_eQTl_02.file
#$ -o stdout_eQTL_02.file

#Run the job - program must have full path:
/ifs/apps/apps/R-3.1.3/bin/Rscript \
	02_eQTL_order_and_match_vers3.R \
	genotype_data_all_treated_baseline.tsv \
	GEx_baseline_4000_and_2000.tsv \
	principal_components_normalised_filtered_PC20.tsv
	#biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt \
	#snp146Common_MatrixEQTL_snp_pos.txt

/ifs/apps/apps/R-3.1.3/bin/Rscript \
        02_eQTL_order_and_match_vers3.R \
        genotype_data_all_treated_final.tsv \
        GEx_treated_4000_and_2000.tsv \
        principal_components_normalised_filtered_PC20.tsv
