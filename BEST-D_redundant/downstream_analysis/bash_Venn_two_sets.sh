#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/Venn_two_sets.R"

set -e


#############
# cis v cis genes
hits_file1="2000+4000-baseline-1.eQTL_cis"
hits_file2="2000+4000-12months-1.eQTL_cis"
column_file1="gene"
column_file2="gene"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="FDR"
background_file="background_probeID_BESTD_expressed.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v cis snps
column_file1="SNP"
column_file2="SNP"
background_file="background_SNPs_MatrixEQTL.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############


#############
# trans v trans genes
hits_file1="2000+4000-baseline-1.eQTL_trans"
hits_file2="2000+4000-12months-1.eQTL_trans"
column_file1="gene"
column_file2="gene"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="FDR"
background_file="background_probeID_BESTD_expressed.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# trans v trans snps
column_file1="SNP"
column_file2="SNP"
background_file="background_SNPs_MatrixEQTL.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans genes baseline
hits_file1="2000+4000-baseline-1.eQTL_cis"
hits_file2="2000+4000-baseline-1.eQTL_trans"
column_file1="gene"
column_file2="gene"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="FDR"
background_file="background_probeID_BESTD_expressed.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans SNPs baseline
column_file1="SNP"
column_file2="SNP"
background_file="background_SNPs_MatrixEQTL.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans genes final
hits_file1="2000+4000-12months-1.eQTL_cis"
hits_file2="2000+4000-12months-1.eQTL_trans"
column_file1="gene"
column_file2="gene"
adj_pvalue="0.05"
col_adj_pvalue="FDR"
col_adj_pvalue2="FDR"
background_file="background_probeID_BESTD_expressed.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans SNPs final
column_file1="SNP"
column_file2="SNP"
background_file="background_SNPs_MatrixEQTL.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans genes reQTLs
hits_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis"
hits_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans"
column_file1="Probe_ID"
column_file2="Probe_ID"
adj_pvalue="0.05"
col_adj_pvalue="FDR.final"
col_adj_pvalue2="FDR.final"
background_file="background_probeID_BESTD_expressed.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############

#############
# cis v trans SNPs reQTLs
column_file1="SNP"
column_file2="SNP"
background_file="background_SNPs_MatrixEQTL.txt"
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############


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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
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
$Rscript $R_script $hits_file1 $column_file1 $hits_file2 $column_file2 \
  $adj_pvalue $col_adj_pvalue $col_adj_pvalue2 $background_file
#############