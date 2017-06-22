#!/bin/bash

# Renaming files links:

ln -s cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p1_1e+06.cis 2000+4000-baseline-1.eQTL_cis

ln -s cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p1_1e+06.cis 2000+4000-12months-1.eQTL_cis 

ln -s trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv 2000+4000-baseline-1.eQTL_trans

ln -s trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_final.tsv_matched.tsv 2000+4000-12months-1.eQTL_trans

ln -s full_topTable_pairing_all_treated.txt 2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx

ln -s reQTLs_FDR5_cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p1_1e+06.cis.txt 2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis

ln -s reQTLs_FDR5_trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt 2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans

#ln -s /ifs/projects/proj043/analysis.dir/files_renamed.dir/files_renamed/2000* .
