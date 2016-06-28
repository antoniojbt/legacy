#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Documents/github.dir/cgat_projects/projects/BEST-D/eQTL_plotting.R"

set -e

###########
# Plots at baseline for cis:
eQTL_file="2000+4000-baseline-1.eQTL_cis"
geno_file="cut_genotype_data_all_treated_baseline.tsv_matched.tsv"
gex_file="cut_GEx_baseline_4000_and_2000.tsv_matched.tsv"
PC_file='PCs_to_adjust_for_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt'
PCs_to_correct='40'

SNP="rs7143764"
PROBE="ILMN_1798177"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs11150882"
PROBE="ILMN_1707137"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs104664"
PROBE="ILMN_1809147"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs4360063"
PROBE="ILMN_1743145"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs10044354"
PROBE="ILMN_1743145"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs3747240"
PROBE="ILMN_1809147"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs6006992"
PROBE="ILMN_1809147"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs3827393"
PROBE="ILMN_1809147"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs2551949"
PROBE="ILMN_3250243"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs2256941"
PROBE="ILMN_3250243"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

###########
# Plots for trans:
eQTL_file="2000+4000-baseline-1.eQTL_trans"
PCs_to_correct='35'


SNP="rs3761268"
PROBE="ILMN_1742442"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10876864"
PROBE="ILMN_1678522"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs10876864"
PROBE="ILMN_1750636"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs10876864"
PROBE="ILMN_1726647"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs17334797"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct


SNP="rs17426064"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs17763596"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs12185233"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs12185268"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct

SNP="rs12373123"
PROBE="ILMN_2393693"
$Rscript $R_script $eQTL_file $geno_file $gex_file $SNP $PROBE $PC_file $PCs_to_correct
