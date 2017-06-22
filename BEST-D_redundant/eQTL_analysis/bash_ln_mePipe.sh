#!/bin/bash
# Scripts and files needed to process files and run context specific eQTL analysis with mePipe R package script

# Get scripts:
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/eQTL_analysis/run_mePipe.R .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/eQTL_analysis/bash_qsub_mePipe.sh .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/eQTL_analysis/qsub_mePipe.sh .

# Get files from eQTL results:
#ln -s /ifs/projects/proj043/analysis.dir/files_renamed.dir/2000* .
ln -s ../cut* .

# Get files for genes of interest, SNPs of interest, and SNP and probe locations:
ln -s /ifs/projects/proj043/analysis.dir/snp146Common_MatrixEQTL_snp_pos.txt .
ln -s /ifs/projects/proj043/analysis.dir/biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt .
#ln -s /ifs/projects/proj043/analysis.dir/VD_genes.txt .
#ln -s /ifs/projects/proj043/analysis.dir/VD_SNPs_GWAS_list.txt .

