##############
# Antonio Berlanga
# Takes Limma's toptable output and adds chr coordinates from Biomart
# based on Illumina probe location
# 19 May 2016
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
R_session_saved_image <- paste('R_session_saved_image_RCircos_preprocessing_toptable','.RData', sep='')
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
# Mx_eQTL_file <- 'full_topTable_pairing_all_treated.txt'

column_to_use <- as.character(args[2])
# column_to_use <- 'probe_ID' # use 'Symbol' for GO files
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
probe_IDs <- Mx_eQTL_data[, c(column_to_use)]
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
length(Mx_eQTL_data[, c(column_to_use)])
length(unique(Mx_eQTL_data[, c(column_to_use)]))

# Some Illumina probes map to more than one location on Biomart:
dim(probe_locations)
length(unique(probe_locations[, 1]))
length(which(duplicated(probe_locations[, 1])))
length(which(!duplicated(probe_locations[, 1])))

# Removing without preference:
probe_locations <- probe_locations[which(!duplicated(probe_locations[, 1])), ]
dim(probe_locations)
head(probe_locations)
####################

####################
# Change column names to allow merge for probe IDs:
head(probe_locations)
colnames(probe_locations)
colnames(probe_locations)[1] <- 'probe_ID'
colnames(probe_locations)[2] <- 'hgnc_symbol_biomart'
colnames(probe_locations)[3] <- 'ensembl_gene_id_biomart'
colnames(probe_locations)[4] <- 'entrezgene_biomart'
colnames(probe_locations)[5] <- 'Probe_ID_chromosome_name'
colnames(probe_locations)[6] <- 'Probe_ID_start_position'
colnames(probe_locations)[7] <- 'Probe_ID_end_position'
colnames(probe_locations)

# Merge with SNP annotations:
Mx_eQTL_data_annot <- merge(Mx_eQTL_data, probe_locations)
dim(Mx_eQTL_data)
dim(Mx_eQTL_data_annot)
head(Mx_eQTL_data_annot)
tail(Mx_eQTL_data_annot)
length(which(is.na(Mx_eQTL_data_annot[, c(column_to_use)])))
colnames(Mx_eQTL_data_annot)
####################


####################
# Write results to disk:
write.table(Mx_eQTL_data_annot, sprintf('biomart_annot_%s.txt', Mx_eQTL_file), 
            row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)

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
