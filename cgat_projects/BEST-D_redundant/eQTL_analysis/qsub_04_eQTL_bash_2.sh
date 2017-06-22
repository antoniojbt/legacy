#!/bin/bash
qsub \
-v geno=cut_genotype_data_all_treated_baseline.tsv_matched.tsv \
-v expr=cut_GEx_baseline_4000_and_2000.tsv_matched.tsv \
-v covar=cut_principal_components_GEx_baseline_matched.tsv \
-v snp_pos=snp146Common_MatrixEQTL_snp_pos.txt \
-v probe_pos=biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt \
qsub_04_eQTL_analysis_with_MatrixEQTL.sh

qsub \
-v geno=cut_genotype_data_all_treated_final.tsv_matched.tsv \
-v expr=cut_GEx_treated_4000_and_2000.tsv_matched.tsv \
-v covar=cut_principal_components_GEx_treated__matched.tsv \
-v snp_pos=snp146Common_MatrixEQTL_snp_pos.txt \
-v probe_pos=biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt \
qsub_04_eQTL_analysis_with_MatrixEQTL.sh

