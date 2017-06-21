#!/bin/bash

Rscript="/Library/Frameworks/R.framework/Resources/bin/Rscript"
R_script="/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/eQTL_counting.R"

set -e

column_adj_pvalue="FDR"
pvalue="0.05"
total_SNPs_tested="477422"
total_probes_tested="14972"
total_SNPprobe_pairs_tested="5755203"

eQTL_file1="2000+4000-12months-1.eQTL_cis"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested

eQTL_file1="2000+4000-12months-1.eQTL_trans"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested

eQTL_file1="2000+4000-baseline-1.eQTL_cis"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested

eQTL_file1="2000+4000-baseline-1.eQTL_trans"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested

eQTL_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested

eQTL_file1="2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans"
$Rscript $R_script $eQTL_file1 $column_adj_pvalue $pvalue $total_SNPs_tested $total_probes_tested $total_SNPprobe_pairs_tested