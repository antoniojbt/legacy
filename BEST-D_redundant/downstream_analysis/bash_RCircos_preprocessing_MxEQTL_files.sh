#!/bin/bash
Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/RCircos_plotting_2/RCircos_preprocessing_MxEQTL_files.R"

set -e

#eQTL files from MxEQTL, large p-value cis:
$Rscript $R_script \
cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p1.0_1MB.cis \
SNP \
gene

#eQTL files from MxEQTL, large p-value cis:
$Rscript $R_script \
cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p1.0_1MB.cis \
SNP \
gene

#eQTL files from MxEQTL trans:
$Rscript $R_script \
trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv  \
SNP \
gene

#eQTL files from MxEQTL trans:
$Rscript $R_script \
trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_final.tsv_matched.tsv  \
SNP \
gene