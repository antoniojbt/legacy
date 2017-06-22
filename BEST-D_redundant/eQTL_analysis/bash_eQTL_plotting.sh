#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="eQTL_plotting.R"
eQTL_file="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis"

set -e

###########
# Plots at baseline:
geno_file="cut_genotype_data_all_treated_baseline.tsv_matched.tsv"
gex_file="cut_GEx_baseline_4000_and_2000.tsv_matched.tsv"
PC_file='PCs_to_adjust_for_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt'
PCs_to_correct='40'

SNP="rs6435757"
PROBE="ILMN_1807425"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs371671"
PROBE="ILMN_1758633"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs11231141"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7124057"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs12225213"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2509982"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2926806"
PROBE="ILMN_1787345"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10203838"
PROBE="ILMN_1807425"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10108662"
PROBE="ILMN_1656310"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs1044158"
PROBE="ILMN_1659523"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

###########
# Plots at 12 months:
geno_file="cut_genotype_data_all_treated_final.tsv_matched.tsv"
gex_file="cut_GEx_treated_4000_and_2000.tsv_matched.tsv"
PC_file='PCs_to_adjust_for_cut_genotype_data_all_treated_final.tsv_matched.tsv.txt'
PCs_to_correct='35'

SNP="rs6435757"
PROBE="ILMN_1807425"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs371671"
PROBE="ILMN_1758633"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs11231141"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7124057"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs12225213"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2509982"
PROBE="ILMN_1796968"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2926806"
PROBE="ILMN_1787345"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10203838"
PROBE="ILMN_1807425"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10108662"
PROBE="ILMN_1656310"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs1044158"
PROBE="ILMN_1659523"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct
