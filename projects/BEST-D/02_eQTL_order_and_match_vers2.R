#############################################
# eQTL file processing script
# Antonio J Berlanga-Taylor
# 28 July 2015

# Input: loads .RData file saved image from 01_eQTL_xxx script (which processed genotype and phenotype 
# file to produce subsets needed). 

# TO DO: This script is OK for small files but not large genotype files (or saved in multiple chr files):
# Concatenate first with GATK, VCFtools, etc.
# Use plink to remove samples eg: plink2 --keep filename or --remove filename
# https://www.cog-genomics.org/plink2/filter#cluster

# Plink2 can re-order:
# https://www.cog-genomics.org/plink2/data#indiv_sort
# plink --bfile [original fileset] --indiv-sort f [file describing new order] --make-bed --out [new prefix]

# Order and match PC, pheno and GEx files, save the IDs and order to file and pass to plink for file/subsets

# Steps:
# Read one file from 01 script - 
# Order each according to unique ID
# Match against each other, 

#############################################

#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir')
# setwd('Desktop/BEST_D_03_MAR.DIR/eqtl_files.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_order_match",Sys.Date(),".txt", sep=""), open = 'a')
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
load('R_session_saved_image_processed_files.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
# R_session_saved_image <- paste('R_session_saved_image_order_and_match_2','.RData', sep='')
R_session_saved_image <- paste('R_session_saved_image_order_and_match','.RData', sep='')
R_session_saved_image
#############################################


#############################################
## Load packages:
library(data.table)

# Get script with functions needed:
source('functions_for_MatrixeQTL.R')
#############################################


#############################################
# Run with command line arguments:
options(echo=TRUE) # to see commands in output file. TO DO: check how it works with sink() above.
args <- commandArgs(trailingOnly = TRUE)
print(args)
#############################################


#############################################
#  Set variables:
# TO DO: pass to configuration file
input_file <- as.chracter(args[1])
# GEx_file <- 'normalised_filtered_expression_values.tab'
ordered_complete_cases_IDs_file <- as.character(args[2])
# ordered_complete_cases_IDs_file <- ''

#############################################


#############################################
# Read data in:
input_data <- fread(input_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
# setkey(input_data, ) # optional, could add as parameter
ordered_complete_cases_IDs_data <- read.csv(ordered_complete_cases_IDs_file, sep = '\t', header = TRUE)
#############################################

#############################################
# Transpose subsetted files for MatrixeQTL, use data.table's transpose:
input_data <- transpose(genotype_data)

# Order data:
input_data <- input_data[, order(names(input_data))]

# Remove samples based on complete cases list:
input_data <- input_data[, which(colnames(input_data) %in% ordered_complete_cases_IDs_data)]
identical(colnames(input_data), ordered_complete_cases_IDs_data)


#############################################


#############################################
# Write to disk:
write.table(covar_PCAs_t_placebo_baseline, 'covar_PCAs_t_placebo_baseline.tsv', 
            sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)



# All columns must be numeric,convert first column to row names;
rownames(GEx_baseline_placebo) <- GEx_baseline_placebo[, 1]
GEx_baseline_placebo[, 1] <- NULL
GEx_baseline_placebo[1:5, 1:5]

rownames(genotype_data_placebo_baseline) <- genotype_data_placebo_baseline[, 1]
genotype_data_placebo_baseline[, 1] <- NULL
genotype_data_placebo_baseline[1:5, 1:5]

# TO DO/check: Floating points can create massive matrices?, round off expression data:
GEx_baseline_placebo <- round(GEx_baseline_placebo, 4)
GEx_baseline_placebo[1:5, 1:5]

# Write files to disk as MatrixEQTL doesn't seem to slice data if loaded directly from an object.
write.table(GEx_baseline_placebo,'GEx_baseline_placebo.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(genotype_data_placebo_baseline, 'genotype_data_placebo_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:

rm(moveme, object_sizes)

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()

q()

# Next: run the script for xxx.
#############################
