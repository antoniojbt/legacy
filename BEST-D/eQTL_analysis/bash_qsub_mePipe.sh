#!/bin/bash

set -e

qsub -v GENO=cut_genotype_data_all_treated_baseline.tsv_matched.tsv \
	-v EXPR=cut_GEx_baseline_4000_and_2000.tsv_matched.tsv \
	-v SNP_POS=snp146Common_MatrixEQTL_snp_pos.txt \
	-v PROBE_POS=biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt \
	-v threshold="1e-5" \
	-v cisThreshold="0.001" \
	qsub_mePipe.sh


qsub -v GENO=cut_genotype_data_all_treated_final.tsv_matched.tsv \
	-v EXPR=cut_GEx_treated_4000_and_2000.tsv_matched.tsv \
	-v SNP_POS=snp146Common_MatrixEQTL_snp_pos.txt \
	-v PROBE_POS=biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt \
	-v threshold="1e-05" \
        -v cisThreshold="0.001" \
	qsub_mePipe.sh

