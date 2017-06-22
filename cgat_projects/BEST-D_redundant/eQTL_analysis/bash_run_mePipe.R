#!/bin/bash

# Get files needed:
ln -s /ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/cut_principal_components_GEx_* .
ln -s /ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/cut_GEx_* .
ln -s /ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/cut_genotype_data_all_treated_*.tsv .
ln -s /ifs/projects/proj043/analysis.dir/snp146Common_MatrixEQTL_snp_pos.txt .
ln -s /ifs/projects/proj043/analysis.dir/biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt .

# Get scripts needed:
ln -s /ifs/devel/antoniob/projects/BEST-D/run_mePipe.R .
ln -s /ifs/devel/antoniob/projects/BEST-D/qsub_mePipe.sh .
ln -s /ifs/devel/antoniob/projects/BEST-D/bash_qsub_mePipe.sh .
