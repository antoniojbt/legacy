#!/bin/bash
qsub -v geno=genotype_data_placebo_baseline.tsv -v expr=GEx_baseline_placebo.tsv -v covar=covar_PCAs_t_placebo_baseline.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R		
qsub -v geno=genotype_data_placebo_final.tsv -v expr=GEx_finalVisit_placebo.tsv -v covar=covar_PCAs_t_placebo_final.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R
qsub -v geno=genotype_data_2000_baseline.tsv -v expr=GEx_baseline_2000.tsv -v covar=covar_PCAs_t_2000_baseline.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R
qsub -v geno=genotype_data_2000_final.tsv -v expr=GEx_finalVisit_2000.tsv -v covar=covar_PCAs_t_2000_final.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R
qsub -v geno=genotype_data_4000_baseline.tsv -v expr=GEx_baseline_4000.tsv -v covar=covar_PCAs_t_4000_baseline.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R
qsub -v geno=genotype_data_4000_final.tsv -v expr=GEx_finalVisit_4000.tsv -v covar=covar_PCAs_t_4000_final.tsv qsub_04_eQTL_analysis_with_MatrixEQTL.R
