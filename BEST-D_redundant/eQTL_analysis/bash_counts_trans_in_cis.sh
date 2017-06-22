#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/count_trans_in_cis.R"

set -e

total_SNPs_tested="477422"
total_probes_tested="14972"
total_SNPprobe_pairs_tested="5755203"
pvalue="0.05"
pvalue_col="FDR"

eQTL_file1='2000+4000-baseline-1.eQTL_trans'
eQTL_file2='2000+4000-baseline-1.eQTL_cis'
$Rscript $R_script $eQTL_file1  $eQTL_file2 $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested $pvalue $pvalue_col

eQTL_file1="2000+4000-12months-1.eQTL_trans"
eQTL_file2="2000+4000-12months-1.eQTL_cis"
$Rscript $R_script $eQTL_file1  $eQTL_file2 $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested $pvalue $pvalue_col

# eQTL_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans"
# eQTL_file2="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis"
# pvalue_col="FDR.final"
# $Rscript $R_script $eQTL_file1  $eQTL_file2 $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested $pvalue $pvalue_col

#####
# Run concatenation:
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/eQTL_counting_concat.R"

$Rscript $R_script