#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/eQTL_plotting.R"
eQTL_file="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans"

set -e

###########
# Plots at baseline:
geno_file="cut_genotype_data_all_treated_baseline.tsv_matched.tsv"
gex_file="cut_GEx_baseline_4000_and_2000.tsv_matched.tsv"
PC_file='PCs_to_adjust_for_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt'
PCs_to_correct='35'

SNP="rs2249742"
PROBE="ILMN_2088124"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7036709"
PROBE="ILMN_1678799"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10469365"
PROBE="ILMN_1759341"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2395943"
PROBE="ILMN_1736238"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7615782"
PROBE="ILMN_1793724"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs978458"
PROBE="ILMN_1715569"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7204270"
PROBE="ILMN_1774901"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs3019776"
PROBE="ILMN_1671818"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs35120848"
PROBE="ILMN_1717261"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs13192471"
PROBE="ILMN_1717261"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2017409"
PROBE="ILMN_1739885"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct
###########

###########
# Same plots at 12 months:
geno_file="cut_genotype_data_all_treated_final.tsv_matched.tsv"
gex_file="cut_GEx_treated_4000_and_2000.tsv_matched.tsv"
PC_file='PCs_to_adjust_for_cut_genotype_data_all_treated_final.tsv_matched.tsv.txt'
PCs_to_correct='25'

SNP="rs2249742"
PROBE="ILMN_2088124"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7036709"
PROBE="ILMN_1678799"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10469365"
PROBE="ILMN_1759341"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2395943"
PROBE="ILMN_1736238"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7615782"
PROBE="ILMN_1793724"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs978458"
PROBE="ILMN_1715569"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs7204270"
PROBE="ILMN_1774901"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs3019776"
PROBE="ILMN_1671818"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs35120848"
PROBE="ILMN_1717261"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs13192471"
PROBE="ILMN_1717261"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs2017409"
PROBE="ILMN_1739885"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct