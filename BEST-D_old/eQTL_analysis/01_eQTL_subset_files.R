#############################################
# eQTL file processing script
# Antonio J Berlanga-Taylor
# 28 July 2015

# The files are processed to generate subsets of interest for later analysis with  MatrixeQTL
# which needs (SNPs/genes/phenotype in rows, individuals in columns), see:
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html

# Inputs needed are:
# QCd plink processed genotype file and phenotype (phenotype) file.

# Outputs: files individuals in columns and variables as rows (required by MatrixeQTL).
# Outputs haven't been ordered or matched to each other (ie geno to expression to phenotype). This is done in 02_eQTL xxx file. 

# Notes:
# Genotype files need to be processed for use in Matrix eQTL. First convert plink formats to oxford format:
# Run the /ifs/devel/antoniob/projects/BEST-D/00_eQTL_genotype_process.sh script first to do this.
# Convert the phenotype file to the format specified (individual IDs as columns, phenotypes as rows). 
# This is done in the next script (02_eQTL_xxx).

# The gene expression file should already be in the right format after limma processing (as probes though, more processing \
# for genes is necessary).
#############################################


#############################################
# Examples and other scripts related to MatrixeQTL to check:
# https://registry.hub.docker.com/u/humburg/eqtl-intro/dockerfile/
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
# https://github.com/jknightlab
# http://link.springer.com/protocol/10.1007%2F978-1-61779-785-9_14
#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/gex_FC_tests/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_", Sys.Date(), ".txt", sep=""))
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

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_eQTL_subset_files','.RData', sep='')
R_session_saved_image
#############################################


#############################################
## Load packages:
library(data.table)

# Get script with functions needed:
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/moveme.R')
source('/ifs/devel/antoniob/projects/BEST-D/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
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
genotype_file <- as.character(args[1]) 
# genotype_file <- '~/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/P140343-Results_FinalReport_clean_SNPs_autosome.A-transpose.matrixQTL.geno'
phenotype_file <- as.character(args[2])
# phenotype_file <- '~/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/BEST-D_phenotype_file_final.tsv'
#############################################

#############################################
# Read data in:
genotype_data <- fread(genotype_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setkey(genotype_data, SNP) # data.table doesn't take column numbers

phenotype_data <- read.csv(phenotype_file, sep = '\t', header = TRUE)

## Inspect data:
# View(genotype_data)
class(genotype_data)
head(sapply(genotype_data, class))
str(genotype_data)
head(genotype_data, 2)
head(colnames(genotype_data))
head(rownames(genotype_data))
genotype_data[1:5, 1:5, with = F]
tables()

str(phenotype_data)
class(phenotype_data)
head(phenotype_data)
tail(phenotype_data)
phenotype_data[1:5, 1:5]
phenotype_data[nrow(phenotype_data), ]
dim(phenotype_data)

#############################################
## Process data and subset files for comparisons of interest.

# Plink outputs FID and IID as FID_IID, so new headers are needed to match the gene expression headers:
genotype_data <- split_plink_IDs(genotype_data)
colnames(genotype_data)
# TO DO: run with plink instead: http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#updatefam
# plink2: https://www.cog-genomics.org/plink2/data

## Generate file subsets according to arm and visit type.
## Genotype file has kit IDs, not patient IDs, and they correspond to baseline.
# Baseline and final visit IDs for kits are in phenotype data file.

# genotype_data file (after splitting strings for 'kitID_kitID' names from Plink transposing) has baseline kit IDs, not patient IDs.
# Expression files have kit IDs (both as it contains both pre and post samples).

# Get IDs to subset genotype file on:
placebo_baseline_samples_IDs <- get_subset_IDs(phenotype_data, 'arm', 2, 'kit_id_randomisation') # placebo
Tx2000_baseline_samples_IDs <- get_subset_IDs(phenotype_data, 'arm', 1, 'kit_id_randomisation') # 2000 UI
Tx4000_baseline_samples_IDs <- get_subset_IDs(phenotype_data, 'arm', 0, 'kit_id_randomisation') # 4000 UI
all_treated_baseline_IDs <- phenotype_data[which(phenotype_data$arm == 0 | phenotype_data$arm == 1), ]
all_treated_baseline_IDs <- all_treated_baseline_IDs$kit_id_randomisation
head(all_treated_baseline_IDs)

# Final visit samples for genotypes need the baseline kit but with samples removed that have NAs at kit_id final visit:
placebo_final_samples_IDs <- get_subset_IDs_two_vars(phenotype_data, 'arm', 2, 'kit_id_finalVisit', 'NA', 'kit_id_randomisation')
Tx2000_final_samples_IDs <- get_subset_IDs_two_vars(phenotype_data, 'arm', 1, 'kit_id_finalVisit', 'NA', 'kit_id_randomisation')
Tx4000_final_samples_IDs <- get_subset_IDs_two_vars(phenotype_data, 'arm', 0, 'kit_id_finalVisit', 'NA', 'kit_id_randomisation')
all_treated_final_IDs <- phenotype_data[which(phenotype_data$arm == 0 | phenotype_data$arm == 1), ]
all_treated_final_IDs <- all_treated_final_IDs[which(all_treated_final_IDs$kit_id_finalVisit != 'NA'), ]
all_treated_final_IDs <- all_treated_final_IDs$kit_id_randomisation
head(all_treated_final_IDs)

# Subset using data.table with function from script sourced above:
# TO DO: save IDs from subset above (to disk) and pass these to plink2 --remove/--keep instead.
# eg plink2 --keep filename or --remove filename
# https://www.cog-genomics.org/plink2/filter#cluster
genotype_data_placebo_baseline <- get_subset_dt(genotype_data, placebo_baseline_samples_IDs, 1)
genotype_data_2000_baseline <- get_subset_dt(genotype_data, Tx2000_baseline_samples_IDs, 1)
genotype_data_4000_baseline <- get_subset_dt(genotype_data, Tx4000_baseline_samples_IDs, 1)
genotype_data_all_treated_baseline <- get_subset_dt(genotype_data, all_treated_baseline_IDs, 1)
genotype_data_all_treated_baseline[1:5, 1:5, with = F]

genotype_data_placebo_final <- get_subset_dt(genotype_data, placebo_final_samples_IDs, 1)
genotype_data_2000_final <- get_subset_dt(genotype_data, Tx2000_final_samples_IDs, 1)
genotype_data_4000_final <- get_subset_dt(genotype_data, Tx4000_final_samples_IDs, 1)
genotype_data_all_treated_final <- get_subset_dt(genotype_data, all_treated_final_IDs, 1)
genotype_data_all_treated_final[1:5, 1:5, with = F]
################################

################################
# Final visit genotypes need final kit IDs to match GEx and covar PC files:
phenotype_data[1:5, 1:5]
master_IDs <- phenotype_data[, c('kit_id_randomisation', 'kit_id_finalVisit', 'pt_id')]
master_IDs$kit_id_randomisation <- as.character(master_IDs$kit_id_randomisation)
master_IDs$kit_id_finalVisit <- as.character(master_IDs$kit_id_finalVisit)
str(master_IDs)
class(master_IDs)
head(master_IDs)
# names(master_IDs)[2] <- 'FID'
names(master_IDs)
# Get colnames from genotype file and match to final kit:
cols_geno <- as.data.frame(colnames(genotype_data_all_treated_final)[-1]) # Minus 'SNP' column name
names(cols_geno)[1] <- 'kit_id_randomisation'
head(cols_geno)

# Merge files:
cols_geno_all_IDs <- merge(cols_geno, master_IDs)
cols_geno_all_IDs$kit_id_randomisation <- as.character(cols_geno_all_IDs$kit_id_randomisation)
str(cols_geno_all_IDs)
cols_geno_all_IDs <- rbind(cols_geno_all_IDs, 'SNP') # Add row for binding with geno data
cols_geno_all_IDs <- cols_geno_all_IDs[order(cols_geno_all_IDs$kit_id_randomisation), ]
# vi(cols_geno_all_IDs)
dim(cols_geno_all_IDs)
head(cols_geno_all_IDs)
tail(cols_geno_all_IDs)

# Genotype file needs to be ordered first, otherwise rbind doesn't match:
genotype_data_all_treated_final <- as.data.frame(genotype_data_all_treated_final)
genotype_data_all_treated_final <- genotype_data_all_treated_final[, order(colnames(genotype_data_all_treated_final))]
head(genotype_data_all_treated_final)
# Add colnames as rows to genotype file:
genotype_data_all_treated_final <- rbind(genotype_data_all_treated_final, as.list(cols_geno_all_IDs$kit_id_randomisation))
genotype_data_all_treated_final <- rbind(genotype_data_all_treated_final, as.list(cols_geno_all_IDs$kit_id_finalVisit))
head(genotype_data_all_treated_final)
tail(genotype_data_all_treated_final)
identical(colnames(genotype_data_all_treated_final), as.character(genotype_data_all_treated_final[497137, ]))
colnames(genotype_data_all_treated_final) <- genotype_data_all_treated_final[497138, ]
colnames(genotype_data_all_treated_final)
genotype_data_all_treated_final <- genotype_data_all_treated_final[-c(497137, 497138), ]
colnames(genotype_data_all_treated_final)[ncol(genotype_data_all_treated_final)] <- 'FID'
genotype_data_all_treated_final <- genotype_data_all_treated_final[, moveme(names(genotype_data_all_treated_final), 'FID first')]
head(genotype_data_all_treated_final)
################################

################################
# TO DO: write to file so that the next script can pick one file at a time (for BEST-D and 018, 
# for others transfer to cmd with plink after getting IDs to subset from).
write.table(genotype_data_all_treated_baseline, paste(deparse(substitute(genotype_data_all_treated_baseline)), '.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)
write.table(genotype_data_all_treated_final, paste(deparse(substitute(genotype_data_all_treated_final)), '.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)

write.table(all_treated_baseline_IDs, paste(deparse(substitute(all_treated_baseline_IDs)), '.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)
write.table(all_treated_final_IDs, paste(deparse(substitute(all_treated_final_IDs)), '.tsv', sep = ''), 
            sep='\t', quote = FALSE, col.names = NA)

# TO DO: function for naming on the fly and writing to disk works but get var name just passes literally
# var1 <- 1
# deparse(x)
# deparse(substitute(var1))
# substitute(var1)
# deparse(var1)
# eval(substitute(var1))
# quote(var1)
# quote(x)
# enquote(var1)
# 
# test_file <- genotype_data_placebo_baseline[1:10, 1:10, with = F]
# var_name <- get_var_name(test_file)
# var_name
# test_file
# 
# write.table(test_file, paste(var_name, '.txt', sep = ''), sep='\t', 
#             quote = FALSE, col.names = NA)
# 
# get_var_name(genotype_data_placebo_baseline[1:10, 1:10, with = F])
# get_var_name(xxx)
# write_DT_tsv(genotype_data_placebo_baseline[1:10, 1:10, with = F])
################################

################################
## TO DO: get basic counts for number of samples per group:
# print('Samples that passed genotyping QC and gene expression QC for final visit:')
# length(which(final_w_GEx[, 16] == 'FinalVisit'))
#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:

rm("all_treated_baseline", 
"all_treated_final",            
"column_ID",                
"genotype_data",
"genotype_file",
"get_subset_dt",
"get_subset_IDs",
"get_subset_IDs_two_vars",
"output_file",
"phenotype_data",
"phenotype_file",
"placebo_baseline_samples",
"placebo_final_samples",
"split_plink_IDs",
"transpose_file",
"Tx2000_baseline_samples",
"Tx2000_final_samples",
"Tx4000_baseline_samples",
"Tx4000_final_samples")

# Print session info:
sessionInfo()

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run the script 02_eQTL_xxx for matching files between genotypes, phenotype and expression.
#############################################
