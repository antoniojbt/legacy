#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/Venn_two_sets.R"

set -e

#############
# cis v diff exp
hits_file1="2000+4000-baseline-1.eQTL_cis"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="gene"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############


#############
# cis v diff exp
hits_file1="2000+4000-12months-1.eQTL_cis"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="gene"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############


#############
# trans v diff exp
hits_file1="2000+4000-baseline-1.eQTL_trans"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="gene"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############

#############
# trans v diff exp
hits_file1="2000+4000-12months-1.eQTL_trans"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="gene"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############


#############
# reQTLs v diff exp
hits_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="Probe_ID"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR.final"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############

#############
# reQTLs v diff exp
hits_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.diffGEx"
column_file1="Probe_ID"
column_file2="probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR.final"
col_adj_pvalue2="adj.P.Val"
background_file="background_probeID_BESTD_expressed.txt"
adj_pvalue_file2="0.20"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
$adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file $adj_pvalue_file2
#############