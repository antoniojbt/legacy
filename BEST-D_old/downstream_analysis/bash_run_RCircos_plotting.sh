#!/bin/bash

#reQTL files with run_RCircos_plotting_reQTLs.R
/Library/Frameworks/R.framework/Resources/bin/Rscript /Users/antoniob/Desktop/scripts_to_upload/run_RCircos_plotting_reQTLs.R \
biomart_annot_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_cis.txt.txt \
biomart_annot_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_trans.txt.txt \
biomart_annot_full_topTable_pairing_all_treated.txt.txt \
biomart_annot_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_cis.txt.txt \
5 \
1 \
chrs_to_exclude_RCircos.txt

mXEQTL files with run_RCircos_plotting_MxEQTL_files.R
/Library/Frameworks/R.framework/Resources/bin/Rscript /Users/antoniob/Desktop/scripts_to_upload/run_RCircos_plotting_MxEQTL_files.R \
biomart_annot_cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p5_1MB.cis.txt \
biomart_annot_cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p8_1MB.trans.txt \
biomart_annot_full_topTable_pairing_all_treated.txt.txt \
biomart_annot_cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p5_1MB.cis.txt \
5 \
1 \
chrs_to_exclude_RCircos.txt


/Library/Frameworks/R.framework/Resources/bin/Rscript /Users/antoniob/Desktop/scripts_to_upload/run_RCircos_plotting_MxEQTL_files.R \
biomart_annot_cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p5_1MB.cis.txt \
biomart_annot_cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p8_1MB.trans.txt \
biomart_annot_full_topTable_pairing_all_treated.txt.txt \
biomart_annot_cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p5_1MB.cis.txt \
5 \
1 \
chrs_to_exclude_RCircos.txt