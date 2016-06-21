#!/bin/bash

qsub \
-v EQTL_1=2000+4000-baseline-1.eQTL_cis \
-v EQTL_2=2000+4000-12months-1.eQTL_cis \
-v GENES_INTEREST=VD_genes.txt \
-v SNPS_INTEREST=VD_SNPs_GWAS_list.txt \
-v total_SNPs_tested=477422 \
-v total_probes_tested=14972 \
-v total_SNPprobe_pairs_tested=5755203 \
qsub_04a_eQTL_context_specific_eQTLs.sh

qsub \
-v EQTL_1=2000+4000-baseline-1.eQTL_trans \
-v EQTL_2=2000+4000-12months-1.eQTL_trans \
-v GENES_INTEREST=VD_genes.txt \
-v SNPS_INTEREST=VD_SNPs_GWAS_list.txt \
-v total_SNPs_tested=477422 \
-v total_probes_tested=14972 \
-v total_SNPprobe_pairs_tested=5755203 \
qsub_04a_eQTL_context_specific_eQTLs.sh

