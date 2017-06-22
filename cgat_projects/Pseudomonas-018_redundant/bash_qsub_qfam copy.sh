#!/bin/bash
qsub -v mpheno=43 -v cmd='--qfam-total perm emp-se' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_TM qsub_qfam.sh
qsub -v mpheno=44 -v cmd='--qfam-total perm emp-se' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_IP qsub_qfam.sh
qsub -v mpheno=45 -v cmd='--qfam-total perm emp-se' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_AT qsub_qfam.sh
qsub -v mpheno=46 -v cmd='--qfam-total perm emp-se' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_CI qsub_qfam.sh
qsub -v mpheno=47 -v cmd='--qfam-total perm emp-se' -v outfile=all_awk_loose_processed_filtered_MAF_PASS_IDs_CO qsub_qfam.sh
