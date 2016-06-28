#!/bin/bash

qsub \
-v EQTL_1=cis_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv \
-v EQTL_2=cis_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_final.tsv_matched.tsv \
-v GENES_INTEREST=VD_genes.txt \
-v SNPS_INTEREST=VD_SNPs_GWAS_list.txt \
qsub_04a_eQTL_context_specific_eQTLs.sh

qsub \
-v EQTL_1=trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv \
-v EQTL_2=trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_final.tsv_matched.tsv \
-v GENES_INTEREST=VD_genes.txt \
-v SNPS_INTEREST=VD_SNPs_GWAS_list.txt \
qsub_04a_eQTL_context_specific_eQTLs.sh

