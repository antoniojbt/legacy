#!/bin/bash
qsub -v pheno_name=TM_new -v cmd='--qfam --perm' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_TM qsub_qfam.sh
qsub -v pheno_name=IP_new -v cmd='--qfam --perm' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_IP qsub_qfam.sh
qsub -v pheno_name=AT_new -v cmd='--qfam --perm' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_AT qsub_qfam.sh
qsub -v pheno_name=CI_new -v cmd='--qfam --perm' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_CI qsub_qfam.sh
qsub -v pheno_name=CO_new -v cmd='--qfam --perm' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_CO qsub_qfam.sh
