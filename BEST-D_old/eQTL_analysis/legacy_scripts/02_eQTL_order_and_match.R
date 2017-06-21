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

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""), open = 'a')
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
phenotype_file <- as.character(args[1])
# phenotype_file <- 'BEST-D_phenotype_file_final.tsv'
GEx_file <- as.chracter(args[2])
# GEx_file <- 'normalised_filtered_expression_values.tab'
genotype_file <- as.character(args[3]) 
# genotype_file <- 'xxxgenotype_data.geno' # From 00_eqtl_scriptxxx and if subsetted then from 01_eqtl_scriptxxx
principal_components_file <- as.character(args[4]) 
# principal_components_file <- 'principal_components_normalised_filtered_PC20.tsv'
#############################################


#############################################
# Read data in:
phenotype_data <- read.csv(phenotype_file, sep = '\t', header = TRUE)

GEx_data <- fread(GEx_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE) #check if first cell is not empty
setkey(GEx_data, ) # data.table doesn't take column numbers

genotype_data <- fread(genotype_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE) #If straight from plink then space separated not tab
setkey(genotype_data, SNP) # data.table doesn't take column numbers

principal_components_data <- read.csv(principal_components_file, sep = '\t', header = TRUE)
#############################################


#############################################
# Phenotype/phenotype file needs to be transposed for MatrixeQTL (so that individuals are columns not rows):
# Get column number that will serve as column names:
column_ID <- which(colnames(phenotype_data) == 'kit_id_randomisation')
column_ID
# Run function to transpose and name columns (in functions file sourced above):
phenotype_data_transposed <- transpose_file(phenotype_data, column_ID)
dim(phenotype_data)
dim(phenotype_data_transposed)
head(phenotype_data_transposed)
tail(phenotype_data_transposed)
#############################################

#############################################
# Plink outputs FID and IID as FID_IID, so new headers are needed to match the gene expression headers:
# Not needed if genotype file was subsetted (script 01_eqtl_xxx):
# genotype_data <- split_plink_IDs(genotype_data)
# colnames(genotype_data)

# Transpose subsetted files for MatrixeQTL, use data.table's transpose:
genotype_data <- transpose(genotype_data)

GEx_data <- transpose(GEx_data)
#############################################


#############################################
## Match each set of files (number, name, order)
# Get subset files from expression analysis (from normalisation script output, written to files):

# subset_baseline_placebo <- readr::read_tsv('subset_baseline_placebo.tab') # This causes some issues downstream with classes and MatrixEQTL reading the file.
GEx_baseline_placebo <- read.csv('GEx_baseline_placebo.tsv', sep = '\t', check.names = FALSE)
head(GEx_baseline_placebo)
dim(GEx_baseline_placebo)
class(GEx_baseline_placebo)

# GEx_baseline_2000 <- readr::read_tsv('GEx_baseline_2000.tsv')
GEx_baseline_2000 <- read.csv('GEx_baseline_2000.tsv', sep = '\t', check.names = FALSE)
# GEx_baseline_4000 <- readr::read_tsv('GEx_baseline_4000.tsv')
GEx_baseline_4000 <- read.csv('GEx_baseline_4000.tsv', sep = '\t', check.names = FALSE)
# GEx_finalVisit_placebo <- readr::read_tsv('GEx_finalVisit_placebo.tsv')
GEx_finalVisit_placebo <- read.csv('GEx_finalVisit_placebo.tsv', sep = '\t', check.names = FALSE)
# GEx_finalVisit_2000 <- readr::read_tsv('GEx_finalVisit_2000.tsv')
GEx_finalVisit_2000 <- read.csv('GEx_finalVisit_2000.tsv', sep = '\t', check.names = FALSE)
# GEx_finalVisit_4000 <- readr::read_tsv('GEx_finalVisit_4000.tsv')
GEx_finalVisit_4000 <- read.csv('GEx_finalVisit_4000.tsv', sep = '\t', check.names = FALSE)

# TO DO: I should create lists with all files and run with function to avoid repetition...
# GEx_baseline_placebo[1:5, 1:5]
# GEx_list <- list(GEx_baseline_placebo, GEx_finalVisit_placebo,
#                     GEx_baseline_2000, GEx_finalVisit_2000, 
#                     GEx_baseline_4000, GEx_finalVisit_4000)
 
# sapply(GEx_list, dim)

# Convert data.tables to data.frames:
genotype_data_placebo_baseline <- data.frame(genotype_data_placebo_baseline, check.names = FALSE)
dim(genotype_data_placebo_baseline)
class(genotype_data_placebo_baseline)
genotype_data_placebo_baseline[1:5, 1:5]

genotype_data_2000_baseline <- data.frame(genotype_data_2000_baseline, check.names = FALSE)
genotype_data_4000_baseline <- data.frame(genotype_data_4000_baseline, check.names = FALSE)
genotype_data_placebo_final  <- data.frame(genotype_data_placebo_final, check.names = FALSE)
genotype_data_2000_final <- data.frame(genotype_data_2000_final, check.names = FALSE)
genotype_data_4000_final <- data.frame(genotype_data_4000_final, check.names = FALSE)


# Genotype, expression and covariate files should be with the right sample names, 
# transposed and GExted:
genotype_data_placebo_baseline[1:5, 1:5]
names(genotype_data_placebo_baseline)[1]
class(genotype_data_placebo_baseline)
# tables() # gives a summary of what data.table exist

names(GEx_baseline_placebo)[1] <- 'FID'
names(GEx_baseline_2000)[1] <- 'FID' 
names(GEx_baseline_4000)[1] <- 'FID'
names(GEx_finalVisit_placebo)[1] <- 'FID'
names(GEx_finalVisit_2000)[1] <- 'FID'
names(GEx_finalVisit_4000)[1] <- 'FID'

# First column in data.table is called 'variable', so naming expression files 
# the same: 
# Change setnames(x,old,new) to avoid copying everything.
names(genotype_data_placebo_baseline)[1] <- 'FID'
names(genotype_data_2000_baseline)[1] <- 'FID'
names(genotype_data_4000_baseline)[1] <- 'FID'
names(genotype_data_placebo_final)[1] <- 'FID'
names(genotype_data_2000_final)[1] <- 'FID'
names(genotype_data_4000_final)[1] <- 'FID'

# Order files:
GEx_baseline_placebo <- GEx_baseline_placebo[, order(names(GEx_baseline_placebo))]
GEx_baseline_placebo[1:5, 1:5]
GEx_baseline_2000 <- GEx_baseline_2000[, order(names(GEx_baseline_2000))]
GEx_baseline_4000 <- GEx_baseline_4000[, order(names(GEx_baseline_4000))]
GEx_finalVisit_placebo <- GEx_finalVisit_placebo[, order(names(GEx_finalVisit_placebo))]
GEx_finalVisit_2000 <- GEx_finalVisit_2000[, order(names(GEx_finalVisit_2000))]
GEx_finalVisit_4000 <- GEx_finalVisit_4000[, order(names(GEx_finalVisit_4000))]

genotype_data_placebo_baseline <- genotype_data_placebo_baseline[, order(names(genotype_data_placebo_baseline))]
genotype_data_placebo_baseline[1:5, 1:5]
genotype_data_2000_baseline <- genotype_data_2000_baseline[, order(names(genotype_data_2000_baseline))]
genotype_data_4000_baseline <- genotype_data_4000_baseline[, order(names(genotype_data_4000_baseline))]
genotype_data_placebo_final  <- genotype_data_placebo_final[, order(names(genotype_data_placebo_final))]
genotype_data_2000_final <- genotype_data_2000_final[, order(names(genotype_data_2000_final))]
genotype_data_4000_final <- genotype_data_4000_final[, order(names(genotype_data_4000_final))]
genotype_data_4000_final[1:5, 1:5]

# Move last column (FID) to first:
#source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
source('moveme.R')
GEx_baseline_placebo <- GEx_baseline_placebo[, moveme(names(GEx_baseline_placebo), 'FID first')]
GEx_baseline_2000 <- GEx_baseline_2000[, moveme(names(GEx_baseline_2000), 'FID first')]
GEx_baseline_4000 <- GEx_baseline_4000[, moveme(names(GEx_baseline_4000), 'FID first')]
GEx_finalVisit_placebo <- GEx_finalVisit_placebo[, moveme(names(GEx_finalVisit_placebo), 'FID first')]
GEx_finalVisit_2000 <- GEx_finalVisit_2000[, moveme(names(GEx_finalVisit_2000), 'FID first')]
GEx_finalVisit_4000 <- GEx_finalVisit_4000[, moveme(names(GEx_finalVisit_4000), 'FID first')]

dim(GEx_finalVisit_4000)
GEx_finalVisit_4000[1:5, 1:5]

genotype_data_placebo_baseline <- genotype_data_placebo_baseline[, moveme(names(genotype_data_placebo_baseline), 'FID first')]
genotype_data_2000_baseline <- genotype_data_2000_baseline[, moveme(names(genotype_data_2000_baseline), 'FID first')]
genotype_data_4000_baseline <- genotype_data_4000_baseline[, moveme(names(genotype_data_4000_baseline), 'FID first')]
genotype_data_placebo_final  <- genotype_data_placebo_final[, moveme(names(genotype_data_placebo_final), 'FID first')]
genotype_data_2000_final <- genotype_data_2000_final[, moveme(names(genotype_data_2000_final), 'FID first')]
genotype_data_4000_final <- genotype_data_4000_final[, moveme(names(genotype_data_4000_final), 'FID first')]

dim(genotype_data_4000_final)
genotype_data_4000_final[1:5, 1:5]

# Now matching files to have the same individuals (and same order) for genotype 
# and expression files:
dim(GEx_baseline_placebo)
dim(genotype_data_placebo_baseline)
length(which((colnames(GEx_baseline_placebo) %in% colnames(genotype_data_placebo_baseline))))
GEx_baseline_placebo <- GEx_baseline_placebo[, which(colnames(GEx_baseline_placebo) %in% colnames(genotype_data_placebo_baseline))]
genotype_data_placebo_baseline <- genotype_data_placebo_baseline[, which(colnames(genotype_data_placebo_baseline) %in% colnames(GEx_baseline_placebo))]
identical(colnames(GEx_baseline_placebo), colnames(genotype_data_placebo_baseline))
dim(GEx_baseline_placebo)
dim(genotype_data_placebo_baseline)
GEx_baseline_placebo[1:5, 1:5]
genotype_data_placebo_baseline[1:5, 1:5]

GEx_baseline_2000 <- GEx_baseline_2000[, which(colnames(GEx_baseline_2000) %in% colnames(genotype_data_2000_baseline))]
genotype_data_2000_baseline <- genotype_data_2000_baseline[, which(colnames(genotype_data_2000_baseline) %in% colnames(GEx_baseline_2000))]
identical(colnames(GEx_baseline_2000), colnames(genotype_data_2000_baseline))

GEx_baseline_4000 <- GEx_baseline_4000[, which(colnames(GEx_baseline_4000) %in% colnames(genotype_data_4000_baseline))]
genotype_data_4000_baseline <- genotype_data_4000_baseline[, which(colnames(genotype_data_4000_baseline) %in% colnames(GEx_baseline_4000))]
identical(colnames(GEx_baseline_4000), colnames(genotype_data_4000_baseline))

GEx_finalVisit_placebo <- GEx_finalVisit_placebo[, which(colnames(GEx_finalVisit_placebo) %in% colnames(genotype_data_placebo_final))]
genotype_data_placebo_final <- genotype_data_placebo_final[, which(colnames(genotype_data_placebo_final) %in% colnames(GEx_finalVisit_placebo))]
identical(colnames(GEx_finalVisit_placebo), colnames(genotype_data_placebo_final))

GEx_finalVisit_2000 <- GEx_finalVisit_2000[, which(colnames(GEx_finalVisit_2000) %in% colnames(genotype_data_2000_final))]
genotype_data_2000_final <- genotype_data_2000_final[, which(colnames(genotype_data_2000_final) %in% colnames(GEx_finalVisit_2000))]
identical(colnames(GEx_finalVisit_2000), colnames(genotype_data_2000_final))

GEx_finalVisit_4000 <- GEx_finalVisit_4000[, which(colnames(GEx_finalVisit_4000) %in% colnames(genotype_data_4000_final))]
genotype_data_4000_final <- genotype_data_4000_final[, which(colnames(genotype_data_4000_final) %in% colnames(GEx_finalVisit_4000))]
identical(colnames(GEx_finalVisit_4000), colnames(genotype_data_4000_final))

# TO DO: clean up, separate covar PCA to other covars, etc.
# covar PCA matching:
covar_PCAs <- read.csv('principal_components_normalised_filtered_PC20.tsv', header = TRUE, 
                       sep = '\t')
head(covar_PCAs)
dim(covar_PCAs)
covar_PCAs$kit_id <- as.character(covar_PCAs$kit_id)
str(covar_PCAs)
object.size(covar_PCAs)

covar_PCAs_t <- data.table(covar_PCAs, check.names = FALSE)
covar_PCAs_t <- dcast.data.table(melt(covar_PCAs_t, id.vars = 1), variable ~ kit_id)
covar_PCAs_t <- data.frame(covar_PCAs_t, check.names = FALSE)
head(covar_PCAs_t)
str(covar_PCAs_t)
dim(covar_PCAs_t)
object.size(covar_PCAs_t)

covar_PCAs_t <- covar_PCAs_t[, order(names(covar_PCAs_t))]
head(covar_PCAs_t)
covar_PCAs_t[1:5, 1:5]

covar_PCAs_t_placebo_baseline <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(genotype_data_placebo_baseline))]
dim(covar_PCAs_t_placebo_baseline)
dim(genotype_data_placebo_baseline) # Plus one column from FID
covar_PCAs_t_placebo_baseline[1:5, 1:5]
genotype_data_placebo_baseline[1:5, 1:5]

identical(colnames(covar_PCAs_t_placebo_baseline), colnames(genotype_data_placebo_baseline)[-1])
identical(colnames(covar_PCAs_t_placebo_baseline), colnames(GEx_baseline_placebo)[-1])

covar_PCAs_t_2000_baseline <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(genotype_data_2000_baseline))]
identical(colnames(covar_PCAs_t_2000_baseline), colnames(genotype_data_2000_baseline)[-1])

covar_PCAs_t_4000_baseline <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(genotype_data_4000_baseline))]
identical(colnames(covar_PCAs_t_4000_baseline), colnames(genotype_data_4000_baseline)[-1])

covar_PCAs_t_placebo_final <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(GEx_finalVisit_placebo))]
identical(colnames(covar_PCAs_t_placebo_final), colnames(genotype_data_placebo_final)[-1])

covar_PCAs_t_2000_final <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(genotype_data_2000_final))]
identical(colnames(covar_PCAs_t_2000_final), colnames(genotype_data_2000_final)[-1])

covar_PCAs_t_4000_final <- covar_PCAs_t[, which(colnames(covar_PCAs_t) %in% colnames(genotype_data_4000_final))]
identical(colnames(covar_PCAs_t_4000_final), colnames(genotype_data_4000_final)[-1])

# Write to disk:
write.table(covar_PCAs_t_placebo_baseline, 'covar_PCAs_t_placebo_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(covar_PCAs_t_2000_baseline, 'covar_PCAs_t_2000_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(covar_PCAs_t_4000_baseline, 'covar_PCAs_t_4000_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(covar_PCAs_t_placebo_final, 'covar_PCAs_t_placebo_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(covar_PCAs_t_2000_final, 'covar_PCAs_t_2000_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(covar_PCAs_t_4000_final, 'covar_PCAs_t_4000_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)


# Order and match covariates file:
dim(covar_transposed)
covar_transposed[1:5, 1:5]

covar_transposed <- covar_transposed[, order(names(covar_transposed))]
covar_transposed[1:5, 1:5]
dim(covar_transposed)

covar_transposed_placebo_baseline <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(genotype_data_placebo_baseline))]
dim(covar_transposed_placebo_baseline)
dim(genotype_data_placebo_baseline)
covar_transposed_placebo_baseline[1:5, 1:5]
genotype_data_placebo_baseline[1:5, 1:5]

identical(colnames(covar_transposed_placebo_baseline), colnames(genotype_data_placebo_baseline)[-1])
identical(colnames(covar_transposed_placebo_baseline), colnames(GEx_baseline_placebo)[-1])

covar_transposed_2000_baseline <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(genotype_data_2000_baseline))]
identical(colnames(covar_transposed_2000_baseline), colnames(genotype_data_2000_baseline)[-1])

covar_transposed_4000_baseline <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(genotype_data_4000_baseline))]
identical(colnames(covar_transposed_4000_baseline), colnames(genotype_data_4000_baseline)[-1])

# TO DO: covar headers, FID and IID are from baseline kits and doesn't have final visit kits.
# covar_transposed_placebo_final <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(GEx_finalVisit_placebo))]
# covar_transposed_placebo_final <- covar_transposed[, which(covar_transposed['kit_id_finalVisit', ] %in% colnames(GEx_finalVisit_placebo))]
# identical(colnames(covar_transposed_placebo_final), colnames(genotype_data_placebo_final)[-1])
# dim(covar_transposed_placebo_final)
# which(colnames(covar_transposed) == '120005003') # Final visit
# which(colnames(covar_transposed) == '120000035') # Randomisation
# 
# dim(genotype_data_placebo_final)
# dim(covar_transposed)
# genotype_data_placebo_final[1:5, 1:5]
# covar_transposed[1:5, 1:5]
# 
# covar_transposed_2000_final <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(genotype_data_2000_final))]
# identical(colnames(covar_transposed_2000_final), colnames(genotype_data_2000_final)[-1])
# 
# covar_transposed_4000_final <- covar_transposed[, which(colnames(covar_transposed) %in% colnames(genotype_data_4000_final))]
# identical(colnames(covar_transposed_4000_final), colnames(genotype_data_4000_final)[-1])

# All columns must be numeric,convert first column to row names;
rownames(GEx_baseline_placebo) <- GEx_baseline_placebo[, 1]
GEx_baseline_placebo[, 1] <- NULL
GEx_baseline_placebo[1:5, 1:5]

rownames(GEx_baseline_2000) <- GEx_baseline_2000[, 1]
GEx_baseline_2000[, 1] <- NULL

rownames(GEx_baseline_4000) <- GEx_baseline_4000[, 1]
GEx_baseline_4000[, 1] <- NULL

rownames(GEx_finalVisit_placebo) <- GEx_finalVisit_placebo[, 1]
GEx_finalVisit_placebo[, 1] <- NULL

rownames(GEx_finalVisit_2000) <- GEx_finalVisit_2000[, 1]
GEx_finalVisit_2000[, 1] <- NULL

rownames(GEx_finalVisit_4000) <- GEx_finalVisit_4000[, 1]
GEx_finalVisit_4000[, 1] <- NULL

rownames(genotype_data_placebo_baseline) <- genotype_data_placebo_baseline[, 1]
genotype_data_placebo_baseline[, 1] <- NULL
genotype_data_placebo_baseline[1:5, 1:5]

rownames(genotype_data_2000_baseline) <- genotype_data_2000_baseline[, 1]
genotype_data_2000_baseline[, 1] <- NULL
genotype_data_2000_baseline[1:5, 1:5]

rownames(genotype_data_4000_baseline) <- genotype_data_4000_baseline[, 1]
genotype_data_4000_baseline[, 1] <- NULL

rownames(genotype_data_placebo_final) <- genotype_data_placebo_final[, 1]
genotype_data_placebo_final[, 1] <- NULL

rownames(genotype_data_2000_final) <- genotype_data_2000_final[, 1]
genotype_data_2000_final[, 1] <- NULL

rownames(genotype_data_4000_final) <- genotype_data_4000_final[, 1]
genotype_data_4000_final[, 1] <- NULL

# TO DO/check: Floating points can create massive matrices?, round off expression data:
GEx_baseline_placebo <- round(GEx_baseline_placebo, 4)
GEx_baseline_placebo[1:5, 1:5]
GEx_baseline_2000 <- round(GEx_baseline_2000, 4)
GEx_baseline_4000 <- round(GEx_baseline_4000, 4)
GEx_finalVisit_placebo <- round(GEx_finalVisit_placebo, 4)
GEx_finalVisit_2000 <- round(GEx_finalVisit_2000, 4)
GEx_finalVisit_4000 <- round(GEx_finalVisit_4000, 4)

# Write files to disk as MatrixEQTL doesn't seem to slice data if loaded directly from an object.
write.table(GEx_baseline_placebo,'GEx_baseline_placebo.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(GEx_baseline_2000,'GEx_baseline_2000.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(GEx_baseline_4000,'GEx_baseline_4000.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(GEx_finalVisit_placebo,'GEx_finalVisit_placebo.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(GEx_finalVisit_2000,'GEx_finalVisit_2000.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(GEx_finalVisit_4000,'GEx_finalVisit_4000.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(genotype_data_placebo_baseline, 'genotype_data_placebo_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(genotype_data_2000_baseline, 'genotype_data_2000_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(genotype_data_4000_baseline, 'genotype_data_4000_baseline.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(genotype_data_placebo_final, 'genotype_data_placebo_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(genotype_data_2000_final, 'genotype_data_2000_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(genotype_data_4000_final, 'genotype_data_4000_final.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
#############################################


#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])

rm(moveme, object_sizes)

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()

q()

# Next: run the script for xxx.
#############################
