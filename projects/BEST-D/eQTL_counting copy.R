########################
# Script for counting SNPs and probes/genes following MatrixEQTL analysis, not for reQTLs
# Antonio J Berlanga-Taylor
# 30 Sept 2015
########################


#############################################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/')

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
#load('R_session_saved_image_order_and_match.RData', verbose=T)
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_eQTL_counts','.RData', sep='')
R_session_saved_image
####################


########################
# Load packages:
library(data.table)
########################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# MatrixeQTL files:
eQTL_file1 <- as.character(args[1])
# eQTL_file1 <- '2000+4000-baseline-1.eQTL_cis'
# eQTL_file1 <- '2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis'

column_adj_pvalue <- as.character(args[2])
# column_adj_pvalue <- 'FDR'
# column_adj_pvalue <- 'FDR.final'

pvalue <- as.numeric(args[3])
# pvalue <- 0.05

total_SNPs_tested <- as.numeric(args[4])
# total_SNPs_tested <- 477422

total_probes_tested <- as.numeric(args[5])
# total_probes_tested <- 14972

total_SNPprobe_pairs_tested <- as.numeric(args[6])
# total_SNPprobe_pairs_tested <- 5755203

print(args)
####################


####################
# Read data:
eQTL_data1 <- fread(eQTL_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(eQTL_data1)
eQTL_data1
setkey(eQTL_data1)
########################

####################
# Switch names for reQTLs:
# TO DO: in reQTLs switch '.final' to the .x and .y outputs and print file names at top with #.
# This works OK because MxEQTL files only have 6 columns, otherwise it will change column names as only testing 
# for NA.

if (colnames(eQTL_data1)[1] == 'Probe_ID') {
  setnames(eQTL_data1, "Probe_ID", "gene")
}

if (is.na(colnames(eQTL_data1)[7])) {
  column_adj_pvalue <- column_adj_pvalue
} else {
  setnames(eQTL_data1, "beta.final", "beta")  
}

if (is.na(colnames(eQTL_data1)[8])) {
  column_adj_pvalue <- column_adj_pvalue
} else {
  setnames(eQTL_data1, "t-stat.final", "t-stat")
}

if (is.na(colnames(eQTL_data1)[9])) {
  column_adj_pvalue <- column_adj_pvalue
} else {
  setnames(eQTL_data1, "p-value.final", "p-value")
}

if (is.na(colnames(eQTL_data1)[10])) {
  column_adj_pvalue <- column_adj_pvalue
} else {
  setnames(eQTL_data1, "FDR.final", "FDR")
}

colnames(eQTL_data1)
eQTL_data1
####################

####################
# Basic counts:
dim(eQTL_data1)
colnames(eQTL_data1)

# Sanity check (no NAs should be present unless it's an reQTL file):
dim(eQTL_data1)
length(which(complete.cases(eQTL_data1)))

# Concatenate SNP-probe pairs in new column, sanity check (each row should correspond to a unique SNP-probe
# pair combination):
# Runs but takes a long time and not needed.
# eQTL_data1 <- within(eQTL_data1,  pairs <- paste(SNP, gene, sep="-"))
# eQTL_data1
# dim(eQTL_data1)
# dim(eQTL_data1[, 'pairs', with = F])
# length(which(!duplicated(eQTL_data1[, 'pairs', with = F])))
# identical(dim(eQTL_data1)[1], length(which(!duplicated(eQTL_data1[, 'pairs', with = F]))))

#######
# Full count:
# Get counts and summaries for each column:
snp_probe_pairs <- nrow(eQTL_data1)
snp_probe_pairs
dim(eQTL_data1) # Equal to SNP-probe pairs (which is equal to each row that MxEQTL outputs)

snp_count_full <- length(which(!duplicated(eQTL_data1[, 'SNP', with = F])))
snp_count_full

probe_count_full <- length(which(!duplicated(eQTL_data1[, 'gene', with = F])))
probe_count_full

beta_summary_full <- t(summary(eQTL_data1[, "beta", with = F]))
beta_summary_full

tstat_summary_full <- t(summary(eQTL_data1[, "t-stat", with = F]))
tstat_summary_full

pvalue_summary_full<- t(summary(eQTL_data1[, "p-value", with = F]))
pvalue_summary_full

min_pvalue_full<- min(eQTL_data1[, "p-value", with = F])
min_pvalue_full

FDR_summary_full<- t(summary(eQTL_data1[, "FDR", with = F]))
FDR_summary_full

min_FDR_full<- min(eQTL_data1[, "FDR", with = F])
min_FDR_full


#######
# Count at significant threshold:
# Get subset which is significant:
eQTL_data1_FDR <- eQTL_data1[which(eQTL_data1[, column_adj_pvalue, with = F] < pvalue), ]
eQTL_data1_FDR

# Get counts and summaries for each column for FDR cut-off:
snp_probe_pairs_FDR <- nrow(eQTL_data1_FDR)
snp_probe_pairs_FDR
dim(eQTL_data1) # Equal to SNP-probe pairs (which is equal to each row that MxEQTL outputs)

snp_count_FDR <- length(which(!duplicated(eQTL_data1_FDR[, 'SNP', with = F])))
snp_count_FDR

probe_count_FDR <- length(which(!duplicated(eQTL_data1_FDR[, 'gene', with = F])))
probe_count_FDR

beta_summary_FDR <- t(summary(eQTL_data1_FDR[, "beta", with = F]))
beta_summary_FDR

tstat_summary_FDR <- t(summary(eQTL_data1_FDR[, "t-stat", with = F]))
tstat_summary_FDR

pvalue_summary_FDR <- t(summary(eQTL_data1_FDR[, "p-value", with = F]))
pvalue_summary_FDR

min_pvalue_FDR <- min(eQTL_data1_FDR[, "p-value", with = F])
min_pvalue_FDR

FDR_summary_FDR <- t(summary(eQTL_data1_FDR[, "FDR", with = F]))
FDR_summary_FDR

min_FDR_FDR <- min(eQTL_data1_FDR[, "FDR", with = F])
min_FDR_FDR
###########
####################

####################
# Save results of counts above to file:
# TO DO: Warning that percentages are extracted from file read but may not be all the tested 
# elements (eg if MxEQTL wasn't run with full p-value).

options(digits = 2)

cat(file = sprintf('counts_%s.txt', eQTL_file1), 
    paste('File: '), '\t', paste(eQTL_file1), '\n',
    # Full count:
    paste('Number passed as total SNPs tested:'), '\t', total_SNPs_tested, '\n', 
    paste('Number passed as total probes tested:'), '\t', total_probes_tested, '\n',
    paste('Number passed as total SNP-probe pairs tested:'), '\t', total_SNPprobe_pairs_tested, '\n',
    paste('Total number of SNP-probe pairs (in file):'), '\t', snp_probe_pairs, '\n',
    paste('Total number of SNPs (in file):'), '\t', snp_count_full, '\n',
    paste('Total number of probes (in file):'), '\t', probe_count_full, '\n',
    #  sprintf('unique_genes_%s:', eQTL), '\t', length(unique(eQTL_data1$gene)), '\n',
    paste('Summary stats beta values:'), '\t', beta_summary_full, '\n', 
    paste('Summary stats t-stat values:'), '\t', tstat_summary_full, '\n', 
    paste('Summary stats p-values:'), '\t', pvalue_summary_full, '\n', 
    paste('Lowest p-value:'), '\t', min_pvalue_full, '\n', 
    paste('Summary stats FDR values:'), '\t', FDR_summary_full, '\n', 
    paste('Lowest FDR value:'), '\t', min_FDR_full, '\n', 
    #Counts under adj.pvalue specified:
    sprintf('Total number of SNP-probe pairs under FDR %s%%:', pvalue*100), '\t', nrow(eQTL_data1_FDR), '\n', 
    sprintf('Total number of SNPs under FDR %s%%:', pvalue*100), '\t', snp_count_FDR, '\n', 
    sprintf('Total number of probes under FDR %s%%:', pvalue*100), '\t', probe_count_FDR, '\n', 
    #  sprintf('unique_genes_%s:', eQTL), '\t', length(unique(eQTL_data1$gene)), '\n',
    sprintf('Summary stats beta values under FDR %s%%:', pvalue*100), '\t', beta_summary_FDR, '\n', 
    sprintf('Summary stats t-stat values under FDR %s%%:', pvalue*100), '\t', tstat_summary_FDR, '\n', 
    sprintf('Summary stats p-values under FDR %s%%:', pvalue*100), '\t', pvalue_summary_FDR, '\n', 
    sprintf('Summary stats FDR values under FDR %s%%:', pvalue*100), '\t', FDR_summary_FDR, '\n', 
    # Additional numbers:
    paste('Total number of SNP-probe pairs under FDR 10%:'), '\t', length(which(eQTL_data1[, 'FDR', with = F] < 0.10)), '\n',
    paste('Total number of SNP-probe pairs under FDR 20%:'), '\t', length(which(eQTL_data1[, 'FDR', with = F] < 0.20)), '\n',
    sprintf('Percentage of SNP-probe pairs under FDR %s%% (of total tested):', pvalue*100), '\t', (nrow(eQTL_data1_FDR) / total_SNPprobe_pairs_tested) * 100, '\n',
    sprintf('Percentage of SNPs under FDR %s%% (of total unique):', pvalue*100), '\t', (snp_count_FDR / total_SNPs_tested) * 100, '\n',
    sprintf('Percentage of probes under FDR %s%% (of total unique):', pvalue*100), '\t', (probe_count_FDR / total_probes_tested) * 100, '\n',
append = FALSE)

########################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run XGR, IPA, etc, cross with GWAS, ENCODE, etc.
####################