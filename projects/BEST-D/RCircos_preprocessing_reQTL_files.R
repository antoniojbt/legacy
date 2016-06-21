##############
# Antonio Berlanga
# Takes MatrixEQTL output and adds chr coordinates from Biomart
# Adds genomic coordinates for link plot in RCircos and re-orders so that 
# chr, start, end, chr2, start2, end2 is output
# based on SNP location followed by Illumina probe location
# 18 May 2016
##############

####################
##Set working directory and file locations and names of required inputs:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/scripts to upload/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_illumina_probes_chr_mismatch_check.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_RCircos_preprocessing','.RData', sep='')
R_session_saved_image
##############

##############
# Load libraries:
library(biomaRt)
library(data.table)
library(plyr)
##############

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Illumina probes with chr locations:
Mx_eQTL_file <- as.character(args[1])
# Mx_eQTL_file <- 'cis_tx_fdr5_reQTLs_annot_all_Tx_joint_cis.txt'
# Mx_eQTL_file <- 'cis_tx_fdr5_reQTLs_annot_all_Tx_joint_trans.txt'
# Mx_eQTL_file <- 'cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p5_1MB.cis'
Mx_eQTL_file <- 'cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p8_1MB.trans'
 
# Mx_eQTL_file <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p5_1MB.cis'
# Mx_eQTL_file <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p8_1MB.trans'

column_to_use <- as.character(args[2])
# column_to_use <- 'Symbol' # use 'Symbol' for GO files
####################


####################
# Read data:
Mx_eQTL_data <- read.csv(Mx_eQTL_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, 
                         strip.white = TRUE, na.strings = c('', ' ', 'NA', '-', 'na'))
class(Mx_eQTL_data)
str(Mx_eQTL_data)
head(Mx_eQTL_data)
dim(Mx_eQTL_data)
####################

####################
# Run biomaRt to get annotations:
# See: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf
# Biomart has hg38 as default
# SNPS are from Ensembl release 84, Variation data are from dbSNP build 146
# http://www.ensembl.info/blog/2016/01/06/whats-coming-in-ensembl-release-84/

# Use SNP position as first set of coordinates:
# Check what's available:
listMarts()
# Check what datasets are included in the biomart database:
mart <- useMart(biomart = 'ENSEMBL_MART_SNP')
listDatasets(mart)
# Set the source (database) and dataset to use:
snpmart <- useMart(biomart = 'ENSEMBL_MART_SNP', dataset = 'hsapiens_snp')
# Check what information is available:
filters <- listFilters(snpmart)
head(filters)
tail(filters)
# View(filters)

attributes_biomart <- listAttributes(snpmart)
# View(attributes_biomart)
class(listAttributes(snpmart))

# Search by keyword:
listAttributes(snpmart)[1:5, 1:3]
listAttributes(snpmart)[which(grepl('chr', listAttributes(snpmart)[, 1], ignore.case = TRUE)), ]

# Obtain values from the database and dataset selected:
snp_ids <- Mx_eQTL_data[, c('SNP')]
head(snp_ids)
length(snp_ids)
snp_attributes <- c('refsnp_id', 'chr_name', 'chrom_start', 'chrom_end')
snp_locations <- getBM(attributes = snp_attributes, filters = 'snp_filter', 
                       values = snp_ids, mart = snpmart) 

# Check results:
head(snp_locations)
tail(snp_locations)
dim(snp_locations)
length(unique(snp_locations$refsnp_id))

# Order by chr, then chr start:
snp_locations <- snp_locations[order(snp_locations$chr_name, snp_locations$chrom_start), ]
head(snp_locations)
count(snp_locations$chr_name)

# Get rid of patches, haplotypes, etc:
snp_locations$chr_name <- as.integer(snp_locations$chr_name)
count(snp_locations$chr_name)
dim(snp_locations[which(!is.na(snp_locations$chr_name)), ])
snp_locations <- snp_locations[which(!is.na(snp_locations$chr_name)), ]
dim(snp_locations)
head(snp_locations)
##############

##############
# Get coordinates using Illumina IDs present in the MatrixEQTL file:
# Set the source (database) and dataset to use:
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check what information is available
filters <- listFilters(ensembl)
head(filters)
# View(listAttributes(ensembl))
class(listAttributes(ensembl))
# Search by keyword:
listAttributes(ensembl)[1:5, 1:3]
listAttributes(ensembl)[which(grepl('illumina', listAttributes(ensembl)[, 1], ignore.case = TRUE)), ]
# Get the values:
probe_IDs <- Mx_eQTL_data[, c('gene')]
head(probe_IDs)
length(probe_IDs)
length(unique(probe_IDs))

probe_attributes <- c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id', 'entrezgene',
                      'chromosome_name', 'start_position', 'end_position')
probe_locations <- getBM(attributes = probe_attributes, filters = 'illumina_humanht_12_v4', 
                         values = probe_IDs, mart = ensembl)

# Check results:
head(probe_locations)
tail(probe_locations)
dim(probe_locations)
length(Mx_eQTL_data$Probe_ID)
length(unique(Mx_eQTL_data$Probe_ID))

# Some Illumina probes map to more than one location on Biomart:
dim(probe_locations)
length(unique(probe_locations[, 1]))
length(which(duplicated(probe_locations[, 1])))
length(which(!duplicated(probe_locations[, 1])))

# Removing without preference:
probe_locations <- probe_locations[which(!duplicated(probe_locations[, 1])), ]
dim(probe_locations)
####################

####################
# Merge MatrixEQTL file with new annotations (SNP and probe locations):
# Change column names to allow merge for SNPs:
head(Mx_eQTL_data)
colnames(Mx_eQTL_data)
colnames(Mx_eQTL_data)[2] <- 'Probe_ID'
head(snp_locations)
colnames(snp_locations)[1] <- 'SNP'
colnames(snp_locations)[2] <- 'SNP_chr_name'
colnames(snp_locations)[3] <- 'SNP_chrom_start'
colnames(snp_locations)[4] <- 'SNP_chrom_end'
colnames(snp_locations)

# Merge with SNP annotations:
Mx_eQTL_data_annot_snp <- merge(Mx_eQTL_data, snp_locations)#, by = 'SNP')
dim(Mx_eQTL_data)
dim(Mx_eQTL_data_annot_snp)
head(Mx_eQTL_data_annot_snp)
length(which(is.na(Mx_eQTL_data_annot_snp$SNP)))

# Change column names to allow merge for probe IDs:
head(probe_locations)
colnames(probe_locations)
colnames(probe_locations)[1] <- 'Probe_ID'
colnames(probe_locations)[2] <- 'hgnc_symbol_biomart'
colnames(probe_locations)[3] <- 'ensembl_gene_id_biomart'
colnames(probe_locations)[4] <- 'entrezgene_biomart'
colnames(probe_locations)[5] <- 'Probe_ID_chromosome_name'
colnames(probe_locations)[6] <- 'Probe_ID_start_position'
colnames(probe_locations)[7] <- 'Probe_ID_end_position'
colnames(probe_locations)

# Merge with SNP annotations:
Mx_eQTL_data_annot_final <- merge(Mx_eQTL_data_annot_snp, probe_locations)#, by = 'Probe_ID')
dim(Mx_eQTL_data)
dim(Mx_eQTL_data_annot_snp)
dim(Mx_eQTL_data_annot_final)
head(Mx_eQTL_data_annot_final)
tail(Mx_eQTL_data_annot_final)
length(which(is.na(Mx_eQTL_data_annot_final$SNP)))
length(which(is.na(Mx_eQTL_data_annot_final$Probe_ID)))
colnames(Mx_eQTL_data_annot_final)
####################

####################
# TO DO: biomaRt doesn't find about a third of the trans probe IDs in some cases (eg reQTL files), 
# so lost for circos plots for now. 
# Doesn't affect other results, only figures... Try illuminaHumanv4.db again or Illumina's original file?
####################


####################
# Write results to disk:
write.table(Mx_eQTL_data_annot_final, sprintf('biomart_annot_%s.txt', Mx_eQTL_file), 
            row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)

# Generate file for RCircos:
# reQTL files have more columns (baseline results and 12 months results):
file_name <- sprintf('biomart_annot_%s.txt', Mx_eQTL_file)
cmd_run <- sprintf('cat %s | cut -f13-15,19-21 > RCircos_%s', file_name, file_name)
system(cmd_run)
####################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run the script for xxx.
####################
