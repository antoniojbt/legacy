#!/bin/bash
# Scripts and files needed to process files and run eQTL analysis with MatrixEQTL
# Run script 00_eQTL_xxx then 01, etc. after this file. 

# Get scripts:
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/eQTL_analysis/*eQTL* .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/moveme.R .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R .

# Get files from gene expression and genotype QC results, update depending on recent analysis for these files:
ln -s /ifs/projects/proj043/analysis.dir/gene_expression_3.dir/GEx_* .
ln -s /ifs/projects/proj043/analysis.dir/gene_expression_3.dir/normalised_filtered*.tab .
ln -s /ifs/projects/proj043/analysis.dir/gene_expression_3.dir/principal_components_normalised_filtered.tsv .
ln -s /ifs/projects/proj043/analysis.dir/gene_expression_3.dir/membership_file_cleaned_all.tab .
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.bed .
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.bim .
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.fam .

# Get files for covariates and SNP and probe locations:
ln -s /ifs/projects/proj043/analysis.dir/BEST-D_phenotype_file_final.tsv .
ln -s /ifs/projects/proj043/analysis.dir/BEST-D_lab_kit_IDs_for_array_QC.tsv .
ln -s /ifs/projects/proj043/analysis.dir/03_form_lab_kits.csv.tr .
ln -s /ifs/projects/proj043/analysis.dir/snp146Common_MatrixEQTL_snp_pos.txt .
ln -s /ifs/projects/proj043/analysis.dir/biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt .


# Remove files not needed (TO DO: clean up):
rm -rf 03*QTL* 04a_eQTL_context_specific_eQTLs* 05_eQTL_analysis_with_MatrixEQTL_multi_tissue.R 
rm -rf 06_eQTL_analysis_with_MatrixEQTL_multi_tissue.R 07_eQTL_analysis_with_MatrixEQTL_multi_tissue.R 
rm -rf bash_ln_04a_eQTL_context_specific_eQTLs* bash_qsub_04a_eQTL_context_specific_eQTLs* eQTL_analysis_Peter_course_notes.R 
rm -rf qsub_04a_eQTL_context_specific_eQTLs*
rm -rf eQTL_notes_to_delete.R 04_eQTL_analysis_with_MatrixEQTL_no_covar.R 01_eQTL_process_and_clean_legacy.R
rm -rf 02_eQTL_order_and_match_vers2.R
rm -rf 02_eQTL_order_and_match.R 

# Use only the top 20 PCs:
cat principal_components_normalised_filtered.tsv | cut -f1-21 > principal_components_normalised_filtered_PC20.tsv
