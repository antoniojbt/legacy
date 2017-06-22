##############
# Antonio Berlanga
# Takes GO annotation file and adds chr coordinates from Biomart to symbol column
# Adds genomic coordinates
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
R_session_saved_image <- paste('R_session_saved_image_RCircos_preprocessing_GO_file','.RData', sep='')
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
# Mx_eQTL_file <- 'GO_VD_proteins.txt'

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

# Get coordinates using Symbols provided by GO files, e.g.:
# http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0033280#term=annotation

# Set the source (database) and dataset to use:
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check what information is available
filters <- listFilters(ensembl)
head(filters)
# View(filters)
filters[1:5, ]
filters[which(grepl('hgnc', filters[, 1], ignore.case = TRUE)), ]

# View(listAttributes(ensembl))
class(listAttributes(ensembl))
# Search by keyword:
listAttributes(ensembl)[1:5, 1:3]
listAttributes(ensembl)[which(grepl('hgnc', listAttributes(ensembl)[, 1], ignore.case = TRUE)), ]
# Get the values:
probe_IDs <- Mx_eQTL_data[, c(column_to_use)]
head(probe_IDs)
length(probe_IDs)
length(unique(probe_IDs))

probe_attributes <- c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene',
                      'chromosome_name', 'start_position', 'end_position')
probe_locations <- getBM(attributes = probe_attributes, filters = 'hgnc_symbol', 
                         values = probe_IDs, mart = ensembl)

# Check results:
head(probe_locations)
tail(probe_locations)
dim(probe_locations)
length(Mx_eQTL_data[, c(column_to_use)])
length(unique(Mx_eQTL_data[, c(column_to_use)]))
count(Mx_eQTL_data[, c(column_to_use)])
count(is.na(probe_locations$hgnc_symbol))
####################


####################
# Write results to disk:
write.table(probe_locations, sprintf('biomart_annot_%s.txt', Mx_eQTL_file), 
            row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)

# Generate file for RCircos:
probe_locations_rcircos <- probe_locations[, c(4, 5, 6, 1)]
head(probe_locations_rcircos)

write.table(probe_locations_rcircos, sprintf('RCircos_biomart_annot_%s.txt', Mx_eQTL_file), 
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
