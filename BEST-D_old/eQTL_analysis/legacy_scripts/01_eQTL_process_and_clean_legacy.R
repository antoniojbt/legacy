#############################################
# eQTL file processing script
# Antonio J Berlanga-Taylor
# 28 July 2015

# Current input: quality controlled gene expression and genotyping data, annotation files
# for SNPs (positions) and probes (positions and associated annotations)

# Outputs: various plots and tables for eQTL associations.
# See libraries and packages required below

# Inputs are:
# Genotype file
# Gene expression file
# Covariates file

# e.g.
# ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport_clean-individuals_and_SNPs* .
# ln -s /ifs/projects/proj043/analysis.dir/BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt .
# ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport.map .
# ln -s /ifs/projects/proj043/analysis.dir/gene_expression.dir/normalised_filtered_all.tab .
# ln -s /ifs/projects/proj043/analysis.dir/gene_expression.dir/membership_file_cleaned_all.tab .
# mv P140343-Results_FinalReport.map P140343-Results_FinalReport_clean-individuals_and_SNPs.map

# These must be in the format that Matrix eQTL needs (SNPs/genes/covariates in rows, individuals in columns), \
# see below for file formatting and conversions.
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html

#############################################


#############################################
# Check Peter H. files and workflow:
# https://registry.hub.docker.com/u/humburg/eqtl-intro/dockerfile/

# Other sources, summaries, slides, etc. to check:
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
# https://github.com/jknightlab
# http://link.springer.com/protocol/10.1007%2F978-1-61779-785-9_14


# Genotype files need to be processed for use in Matrix eQTL. First convert plink formats to oxford format:
# https://www.cog-genomics.org/plink2/data#recode
# Run the /ifs/devel/antoniob/projects/BEST-D/00_eQTL_genotype_process.sh script to do this.
# Convert the covariate file to the format specified (individual IDs as columns, phenotypes as rows). This is done in 
# the next script (02_eQTL_xxx).

# The gene expression file should already be in the right format after limma processing (as probes though, more processing \
# for genes is necessary).

#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir')
# setwd('/Users/antoniob/Desktop/eQTL_27_Aug/')

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
# load('R_session_saved_image_processed_files.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_processed_files','.RData', sep='')
R_session_saved_image

#############################################


#############################################
## Update packages if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite('data.table')

# Had a segmentation fault with devtools:
# https://github.com/hadley/devtools/issues/653
# remove.packages('RCurl')
# remove.packages('curl')
# biocLite('RCurl')
# biocLite('curl')
# install.packages("devtools", dependencies = TRUE)
# library(reshape)

## Download Peter's scripts for processing scripts for MatrixQTL:
# devtools::install_github('humburg/Rsge')
# devtools::install_github('jknightlab/mePipe')
# devtools::install_github('hadley/readr')


#and load them:
packages <-c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse', 'Rsge', 'mePipe', 
             'readr', 'reshape2', 'data.table')
lapply(packages, require, character.only = TRUE)
sessionInfo()
#############################################


#############################################
# Run with command line arguments:
options(echo=TRUE) # to see commands in output file. TO DO: check how it works with sink() above.
args <- commandArgs(trailingOnly = TRUE)
print(args)


# Read data in:
# TO DO: pass to configuration file
genotype_data <- as.character(args[1]) 
#'P140343-Results_FinalReport_clean-individuals_and_SNPs.A-transpose.matrixQTL.geno' 
expression_data <- as.character(args[2]) 
#'normalised_filtered_expression_values_all.tab'
covariates_data <- as.character(args[3])
#'BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt'

# TO DO: pass this as a parameter
# Read file with master IDs for both kits and patients:
lab_kit_IDs <- read.csv('03_form_lab_kits.csv.tr')
head(lab_kit_IDs)
dim(lab_kit_IDs)

# genotype_data <- ('P140343-Results_FinalReport_clean-individuals_and_SNPs.A-transpose.matrixQTL.geno')
# expression_data <-('normalised_filtered_expression_values_all.tab')
# covariates_data <- ('BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt')

geno <- readr::read_tsv(genotype_data)
expr <- readr::read_tsv(expression_data)
covar <- read.csv(covariates_data, sep = '\t', header = TRUE)

# Also read in membership file will all kit and patient IDs (needed to cross match IDs as genotype file only has baseline kit IDs for example but
# MatrixQTL needs matching sample names and order for pairs of genotype/phenotype files):
membership_file_cleaned_all <- readr::read_tsv('membership_file_cleaned_all.tab')
head(membership_file_cleaned_all)
# Missing first column name:
colnames(membership_file_cleaned_all)[1] <- 'FID'
head(membership_file_cleaned_all)
length(rownames(membership_file_cleaned_all))
str(membership_file_cleaned_all)
class(membership_file_cleaned_all)

# Inspect data:
# View(geno)

str(geno)
head(geno, 2)
colnames(geno)
geno[1:5, 1:5]

str(expr)
head(expr, 2)
colnames(expr)
expr[1:5, 1:5]

str(covar)
class(covar)
head(covar)
covar[1:5, 1:5]
covar[nrow(covar), ]
dim(covar)

# Phenotype/covariates file needs to be transposed for MatrixeQTL (so that individuals are columns not rows):

covar_names <- covar[, 1] # Keep first column row names to use as column names for transposed file.
covar_names
covar_transposed <- as.data.frame(t(covar))
colnames(covar_transposed) <- covar_names

str(covar_transposed)
class(covar_transposed)
head(covar_transposed)
covar[1:5,1:5]
covar_transposed[1:5,1:5]
covar_transposed[, ncol(covar_transposed)]
dim(covar_transposed)
dim(covar)

# Some sanity checks:
length(which(covar_transposed[4, ] == 'Randomisation')) # Should match the number of samples
length(which(covar_transposed[4, ] != 'Randomisation'))

#############################################


#############################################
## Process data and subset files so that there are matching (number, label and order) genotype and expression pairs for each arm.

# Plink outputs FID and IID as FID_IID, so new headers are needed to match the gene expression headers:
geno_new_colnames <- strsplit(colnames(geno), split='_', fixed=TRUE) #[1,]
# View(geno_new_colnames)
geno_new_colnames[[2]][[1]]
geno_new_colnames <- matrix(unique(unlist(geno_new_colnames), ncol=1, byrow=T))

length(geno_new_colnames)
length(geno) # Lengths should match

geno2 <- geno
geno2[1:5,1:5]
names(geno2) <- geno_new_colnames
geno2[1:5,1:5]
dim(geno2)

## Genotype file has kit IDs, not patient IDs, and they correspond to baseline.
# Add patient IDs along with final visit kit IDs so that sub-setting is easier afterwards:

# Geno2 file (after splitting strings for 'kitID_kitID' names from Plink transposing) has baseline kit IDs, not patient IDs.
# Expression files have kit IDs.
# Covariates file has baseline IDs only:
# membership_file_cleaned_all (from gene expression analysis) has all IDs for samples which passed gene expression QC.

# TO DO: clean up, expression isn't used hin this script?
# Expression and membership sample number should match as these files come from the gene expression QC.
# There's one extra column in the expression file with the probe ID).
dim(expr)
dim(membership_file_cleaned_all)
head(membership_file_cleaned_all)
head(expr)

dim(geno2)
dim(covar_transposed)

geno2[1:5, 1:5]
expr[1:5, 1:5]
covar_transposed[1:5, 1:5]
membership_file_cleaned_all[1:5, 1:5]

# Check master file with kit and patient IDs to match against genotype file: 
head(lab_kit_IDs)
tail(lab_kit_IDs)
dim(lab_kit_IDs)
lab_kit_IDs$pt_id <- as.character(lab_kit_IDs$pt_id)
lab_kit_IDs$kit_id <- as.character(lab_kit_IDs$kit_id)
str(lab_kit_IDs)
summary(lab_kit_IDs)

# Keep only samples which correspond to baseline and final visit. This corresponds to the full set which would have been available for
# gene expression and genotyping (check this is correct):
lab_kit_IDs_baseline_final <- lab_kit_IDs[which(lab_kit_IDs[, 3] != 'Six Month'), ]
lab_kit_IDs_baseline_final <- lab_kit_IDs_baseline_final[which(lab_kit_IDs_baseline_final[, 3] != 'One Month'), ]
head(lab_kit_IDs_baseline_final)
tail(lab_kit_IDs_baseline_final)
dim(lab_kit_IDs_baseline_final)
summary(lab_kit_IDs_baseline_final)

# Cross with QC'd genotype file (which has baseline kit IDs):
geno2_headers <- as.data.frame(names(geno2)[-1])
names(geno2_headers)[1] <- 'kit_id'
head(geno2_headers)
genotype_IDs <- merge(lab_kit_IDs_baseline_final, geno2_headers, by='kit_id', all.y = TRUE)
head(genotype_IDs)
tail(genotype_IDs)
dim(genotype_IDs)
summary(genotype_IDs)

# Add patient IDs to this:
genotype_IDs_double_kit <- merge(lab_kit_IDs_baseline_final, genotype_IDs, by='pt_id', all.y = TRUE)
head(genotype_IDs_double_kit)
tail(genotype_IDs_double_kit)
dim(genotype_IDs_double_kit)
str(genotype_IDs_double_kit)
summary(genotype_IDs_double_kit)

# Re-label and delete repeated columns for clarity:
genotype_IDs_double_kit_clean <- genotype_IDs_double_kit[, 1:5]
head(genotype_IDs_double_kit_clean)
dim(genotype_IDs_double_kit_clean)
names(genotype_IDs_double_kit_clean)
names(genotype_IDs_double_kit_clean)[2] <- 'arm'
names(genotype_IDs_double_kit_clean)[3] <- 'visit_type'
names(genotype_IDs_double_kit_clean)[4] <- 'kit_id_final_visit'
names(genotype_IDs_double_kit_clean)[5] <- 'kit_id_randomisation'
head(genotype_IDs_double_kit_clean)
dim(genotype_IDs_double_kit_clean)
genotype_IDs_double_kit_clean_ordered <- genotype_IDs_double_kit_clean[order(genotype_IDs_double_kit_clean[, 1], genotype_IDs_double_kit_clean[, 3]), ]
head(genotype_IDs_double_kit_clean_ordered)
tail(genotype_IDs_double_kit_clean_ordered)

# Keep unique patient IDs with both kit IDs as columns:
genotype_IDs_double_kit_clean_unique <- genotype_IDs_double_kit_clean_ordered[!duplicated(genotype_IDs_double_kit_clean_ordered[, 1]), ] 
head(genotype_IDs_double_kit_clean_unique)
tail(genotype_IDs_double_kit_clean_unique)
dim(genotype_IDs_double_kit_clean_unique)
write.table(genotype_IDs_double_kit_clean_unique, 'genotype_IDs_double_kit_clean_unique.tsv', sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
#############################################


#############################################
## Generate file subsets according to arm and visit type.
# Add label names (columns) to genotype file so as to be able to subset based on visit type afterwards:
# Easier to transpose genotype data and merge than to merge row-wise:
geno2_by_columns <- t(geno2)
geno2_by_columns[1:5, 1:5]
dim(geno2_by_columns)
colnames(geno2_by_columns) <- geno2_by_columns[1, ]
geno2_by_columns_df <- cbind(kit_id = rownames(geno2_by_columns), geno2_by_columns)
geno2_by_columns_df[1:5, 1:5]
geno2_by_columns_df <- as.data.frame(geno2_by_columns_df[-1, ])
rownames(geno2_by_columns_df) <- NULL
names(geno2_by_columns_df)[1] <- 'kit_id_randomisation'
geno2_by_columns_df[1:5, 1:5]
dim(geno2_by_columns_df)
# write.table(geno2_by_columns_df, 'geno2_by_columns.tsv', quote = FALSE, row.names = TRUE, col.names = FALSE, sep = '\t')

# TO DO: use tmp, free some spaces:
rm(geno, geno2, geno2_by_columns)
gc()

# Merge SNP data with patient and kit IDs:
geno2_by_columns_df_pt_IDs <- merge(genotype_IDs_double_kit_clean_unique, geno2_by_columns_df, by='kit_id_randomisation')
geno2_by_columns_df_pt_IDs[1:5, 1:10]
dim(geno2_by_columns_df_pt_IDs)


# Subset files:
length(which(geno2_by_columns_df_pt_IDs[, 3] == 2)) # Placebo
length(which(geno2_by_columns_df_pt_IDs[, 3] == 1)) # 2000 U
length(which(geno2_by_columns_df_pt_IDs[, 3] == 0)) # 4000 U
length(which(geno2_by_columns_df_pt_IDs[, 4] == 'Randomisation'))

columns_to_delete <- c(2, 3, 4, 5)
genotype_data_placebo_baseline <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 2), -columns_to_delete]
genotype_data_placebo_baseline[1:5, 1:10]
dim(genotype_data_placebo_baseline)
class(genotype_data_placebo_baseline)

genotype_data_2000_baseline <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 1), -columns_to_delete]
genotype_data_4000_baseline <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 0), -columns_to_delete]

columns_to_delete <- c(1, 2, 3, 4)
genotype_data_placebo_final <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 2), -columns_to_delete]
genotype_data_placebo_final[1:5, 1:5]
dim(genotype_data_placebo_final)
genotype_data_2000_final <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 1), -columns_to_delete]
genotype_data_4000_final <- geno2_by_columns_df_pt_IDs[which(geno2_by_columns_df_pt_IDs[, 3] == 0), -columns_to_delete]

# Tranpose subsetted files for MatrixeQTL, use reshape2 and data.table:

genotype_data_placebo_baseline <- data.table(genotype_data_placebo_baseline)
genotype_data_2000_baseline <- data.table(genotype_data_2000_baseline)
genotype_data_4000_baseline <- data.table(genotype_data_4000_baseline)

genotype_data_placebo_final <- data.table(genotype_data_placebo_final)
genotype_data_2000_final <- data.table(genotype_data_2000_final)
genotype_data_4000_final <- data.table(genotype_data_4000_final)

genotype_data_placebo_baseline <- dcast.data.table(melt(genotype_data_placebo_baseline, id.vars = 1), variable ~ kit_id_randomisation)
genotype_data_2000_baseline <- dcast.data.table(melt(genotype_data_2000_baseline, id.vars = 1), variable ~ kit_id_randomisation)
genotype_data_4000_baseline <- dcast.data.table(melt(genotype_data_4000_baseline, id.vars = 1), variable ~ kit_id_randomisation)

genotype_data_placebo_final <- dcast.data.table(melt(genotype_data_placebo_final, id.vars = 1), variable ~ kit_id_final_visit)
genotype_data_2000_final <- dcast.data.table(melt(genotype_data_2000_final, id.vars = 1), variable ~ kit_id_final_visit)
genotype_data_4000_final <- dcast.data.table(melt(genotype_data_4000_final, id.vars = 1), variable ~ kit_id_final_visit)

# TO DO: These should be functions...
# transpose_genos <- function(input_file) {
#   input_file_t <- dcast.data.table(melt(input_file, id.vars = 1), variable ~ kit_id_randomisation)
#   file_name <- paste(input_file, '_t_df.txt', sep = '')
#   assign(x = input_file_t_df, value = file_name, pos = .GlobalEnv)
#   #  result <- input_file_t_df
# }
 

## TO DO: get basic counts for number of samples per group:
# print('Samples that passed genotyping QC and gene expression QC for final visit:')
# length(which(final_w_GEx[, 16] == 'FinalVisit'))
#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
ls()
object_sizes <- sapply(ls(), function(x) object.size(get(x)))
as.matrix(rev(sort(object_sizes))[1:10])

# Delete:
rm(columns_to_delete, covar, covar_names, geno2_by_columns_df, geno2_by_columns_df_pt_IDs, 
   geno2_headers, geno_new_colnames, genotype_IDs, genotype_IDs_double_kit, 
   genotype_IDs_double_kit_clean, genotype_IDs_double_kit_clean_ordered, genotype_data, 
   lab_kit_IDs_baseline_final, object_sizes, output_file, packages, covariates_data,
   expression_data, expr)

sessionInfo()

# # Keep:
#   covar_transposed, genotype_IDs_double_kit_clean_unique, 
# genotype_data_2000_baseline
# genotype_data_2000_final
# genotype_data_4000_baseline
# genotype_data_4000_final
# genotype_data_placebo_baseline
# genotype_data_placebo_final
# membership_file_cleaned_all
# lab_kit_IDs
# expr

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

# objects_to_save <- (c('', '', '', ''))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script 02_eQTL_xxx for matching files between genotypes, covariates and expression.
#############################################
