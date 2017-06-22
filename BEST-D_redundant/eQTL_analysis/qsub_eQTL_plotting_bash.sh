!#/bin/bash
qsub -v eQTL_file=genotype_data_placebo_baseline.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_placebo_baseline.tsv -v expr=subset_baseline_placebo.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
qsub -v eQTL_file=genotype_data_placebo_final.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_placebo_final.tsv -v expr=subset_finalVisit_placebo.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
qsub -v eQTL_file=genotype_data_2000_baseline.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_2000_baseline.tsv -v expr=subset_baseline_2000.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
qsub -v eQTL_file=genotype_data_2000_final.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_2000_final.tsv -v expr=subset_finalVisit_2000.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
qsub -v eQTL_file=genotype_data_4000_baseline.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_4000_baseline.tsv -v expr=subset_baseline_4000.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
qsub -v eQTL_file=genotype_data_4000_final.tsv_MatrixEQTL_loose_1MB.cis -v geno=genotype_data_4000_final.tsv -v expr=subset_finalVisit_4000.tsv -v snp=rs7143764 -v probe=ILMN_1798177 qsub_eQTL_plotting.sh 
