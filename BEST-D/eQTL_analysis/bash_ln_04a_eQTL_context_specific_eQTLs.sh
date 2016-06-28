#!/bin/bash
# Scripts and files needed to process files and run context specific eQTL analysis with R script
# Run script 04a_eQTL_xxx 


# Get scripts:
ln -s /ifs/devel/antoniob/projects/BEST-D/04a_eQTL_context_specific_eQTLs.R .
ln -s /ifs/devel/antoniob/projects/BEST-D/qsub_04a_eQTL_context_specific_eQTLs.sh .
cp /ifs/devel/antoniob/projects/BEST-D/bash_qsub_04a_eQTL_context_specific_eQTLs.sh . 

# Get files from eQTL results:
#ln -s /ifs/projects/proj043/analysis.dir/mePipe_runs.dir/cis* .
#ln -s /ifs/projects/proj043/analysis.dir/mePipe_runs.dir/trans* .
ln -s /ifs/projects/proj043/analysis.dir/files_renamed.dir/2000* .

# Get files for genes of interest, SNPs of interest, and SNP and probe locations:
ln -s /ifs/projects/proj043/analysis.dir/snp146Common_MatrixEQTL_snp_pos.txt .
ln -s /ifs/projects/proj043/analysis.dir/biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt .
ln -s /ifs/projects/proj043/analysis.dir/VD_genes.txt .
ln -s /ifs/projects/proj043/analysis.dir/VD_SNPs_GWAS_list.txt .

