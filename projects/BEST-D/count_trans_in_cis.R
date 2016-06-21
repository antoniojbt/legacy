#############################################
# eQTL counts of proportions of trans in cis
# Antonio J Berlanga-Taylor
# 2 June 2016
# Input: eQTL files from MatrixeQTL
# Outputs: file with counts and IDs of SNPs, probes and SNP-probe pairs shared/not shared
#############################################


#############################################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_trans_in_cis_counts",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output"))

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
R_session_saved_image <- paste('R_session_saved_image_eQTL_trans_in_cis_counts','.RData', sep='')
R_session_saved_image
####################


####################
# Load packages:
library(data.table)
library(illuminaHumanv4.db)
# source('/ifs/devel/antoniob/projects/BEST-D/functions_for_MatrixeQTL.R')
source('/Users/antoniob/Desktop/Downloads_to_delete/ifs_scripts_backup.dir/projects/BEST-D/functions_for_MatrixeQTL.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# MatrixeQTL files:
eQTL_file1 <- as.character(args[1])
# eQTL_file1 <- '2000+4000-baseline-1.eQTL_trans'
# eQTL_file1 <- '2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_trans'

eQTL_file2 <- as.character(args[2])
# eQTL_file2 <- '2000+4000-baseline-1.eQTL_cis'
# eQTL_file2 <- '2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis'

total_SNPs_tested <- as.numeric(args[3])
# total_SNPs_tested <- 477422

total_probes_tested <- as.numeric(args[4])
# total_probes_tested <- 14972

total_SNPprobe_pairs_tested <- as.numeric(args[5])
# total_SNPprobe_pairs_tested <- 5755203

pvalue <- as.numeric(args[6])
# pvalue <- 0.05

pvalue_col <- as.character(args[7])
# pvalue_col <- 'FDR'

print(args)
####################

####################
# Set up file naming for outputs:
# to create: 2000+4000-baseline-1.trans_in_cis
eQTL_file1_base <- strsplit(eQTL_file1, '[.]')
eQTL_file1_base <- eQTL_file1_base[[1]][1]
eQTL_file1_base
eQTL_file1_ext <- strsplit(eQTL_file1, '_')
eQTL_file1_ext <- eQTL_file1_ext[[1]][2]
eQTL_file1_ext

eQTL_file2_base <- strsplit(eQTL_file2, '[.]')
eQTL_file2_base <- eQTL_file2_base[[1]][1]
eQTL_file2_base
eQTL_file2_ext <- strsplit(eQTL_file2, '_')
eQTL_file2_ext <- eQTL_file2_ext[[1]][2]
eQTL_file2_ext

output_file_name <- sprintf('%s_in_%s_%s.eQTL', eQTL_file1_ext, eQTL_file2_ext, eQTL_file1_base)
output_file_name
####################


####################
# Read data:
eQTL_data1 <- fread(eQTL_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data1, 'gene', 'Probe_ID')
dim(eQTL_data1)
colnames(eQTL_data1)
setkey(eQTL_data1, SNP, Probe_ID)

eQTL_data2 <- fread(eQTL_file2, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data2, 'gene', 'Probe_ID')
dim(eQTL_data2)
colnames(eQTL_data2)
setkey(eQTL_data2, SNP, Probe_ID)

tables()

eQTL_data1
eQTL_data2
####################

####################
# Subset data to pvalue given:
eQTL_data1_FDR <- eQTL_data1[which(eQTL_data1[, pvalue_col, with = F] < pvalue), ]
eQTL_data1_FDR
dim(eQTL_data1_FDR)
eQTL_data2_FDR <- eQTL_data2[which(eQTL_data2[, pvalue_col, with = F] < pvalue), ]
eQTL_data2_FDR
dim(eQTL_data2_FDR)
####################

####################
# Get proportions of trans which are cis: pass pairs of files (cis/trans (baseline, 12months, etc.))

#######
# SNP-probe pairs in each file and shared:
# Number of SNP-probe pairs in each file, sanity check (no NAs should be present):
dim(eQTL_data1_FDR)
length(which(complete.cases(eQTL_data1_FDR)))
dim(eQTL_data2_FDR)
length(which(complete.cases(eQTL_data2_FDR)))

# Concatenate SNP-probe pairs in new column:
# TO DO: this errors with reQTL files: 'Can't assign to the same column twice in the same query (duplicates detected).'
eQTL_data1_FDR <- within(eQTL_data1_FDR,  pairs <- paste(SNP, Probe_ID, sep="-"))
eQTL_data1_FDR
colnames(eQTL_data1_FDR)

# eQTL_data1_FDR[which(eQTL_data1_FDR[, Probe_ID] == 'ILMN_1793287'), ]

eQTL_data2_FDR <- within(eQTL_data2_FDR,  pairs <- paste(SNP, Probe_ID, sep="-"))
eQTL_data2_FDR

# Shared pairs:
shared_pairs <- getShared(eQTL_data1_FDR, eQTL_data2_FDR, 'pairs', pvalue_col)
shared_pairs
dim(shared_pairs)
class(shared_pairs)

# Pairs which only appear in one file:
pairs_eQTL_data1_FDR <- getPrivate(eQTL_data1_FDR, eQTL_data2_FDR, 'pairs', pvalue_col)
head(pairs_eQTL_data1_FDR)
dim(pairs_eQTL_data1_FDR)

pairs_eQTL_data2_FDR <- getPrivate(eQTL_data2_FDR, eQTL_data1_FDR, 'pairs', pvalue_col)
head(pairs_eQTL_data2_FDR)
dim(pairs_eQTL_data2_FDR)

# Sanity, the maximum number of private pairs would be the FDR cut or private plus shared (use non-duplicated values):
length(unique(eQTL_data1_FDR$pairs)) == length(unique(pairs_eQTL_data1_FDR$pairs))
length(unique(eQTL_data1_FDR$pairs)) == length(unique(pairs_eQTL_data1_FDR$pairs)) + nrow(shared_pairs)
length(unique(eQTL_data2_FDR$pairs)) == length(unique(pairs_eQTL_data2_FDR$pairs))
length(unique(eQTL_data2_FDR$pairs)) == length(unique(pairs_eQTL_data2_FDR$pairs)) + nrow(shared_pairs)
# Sanity, should be zero:
length(which(as.character(pairs_eQTL_data1_FDR$pairs) %in% as.character(pairs_eQTL_data2_FDR$pairs)))

#######
# Shared SNPs (returns SNP-probe pairs where shared SNPs appear) :
shared_SNPs <- getShared(eQTL_data1_FDR, eQTL_data2_FDR, 'SNP', pvalue_col)
head(shared_SNPs)
dim(shared_SNPs)
length(unique(shared_SNPs$SNP))

# SNPs private to each group, returns SNP-probe pairs in data1 appear:
SNPs_eQTL_data1_FDR <- getPrivate(eQTL_data1_FDR, eQTL_data2_FDR, 'SNP', pvalue_col)
head(SNPs_eQTL_data1_FDR)
dim(SNPs_eQTL_data1_FDR)
length(unique(SNPs_eQTL_data1_FDR$SNP))

SNPs_eQTL_data2_FDR <- getPrivate(eQTL_data2_FDR, eQTL_data1_FDR, 'SNP', pvalue_col)
head(SNPs_eQTL_data2_FDR)
dim(SNPs_eQTL_data2_FDR)
length(unique(SNPs_eQTL_data2_FDR$SNP))
length(which(complete.cases(SNPs_eQTL_data2_FDR)))
length(which(complete.cases(shared_SNPs)))

# Sanity, the maximum number of private SNPs would be the FDR cut or private plus shared (use non-duplicated values):
length(unique(eQTL_data1_FDR$SNP)) == length(unique(SNPs_eQTL_data1_FDR$SNP))
length(unique(eQTL_data1_FDR$SNP)) == length(unique(SNPs_eQTL_data1_FDR$SNP)) + length(unique(shared_SNPs$SNP))
length(unique(eQTL_data2_FDR$SNP)) == length(unique(eQTL_data2_FDR$SNP))
length(unique(eQTL_data2_FDR$SNP)) == length(unique(SNPs_eQTL_data2_FDR$SNP)) +length(unique(shared_SNPs$SNP))
# Sanity, should be zero:
length(which(as.character(SNPs_eQTL_data1_FDR$SNP) %in% as.character(SNPs_eQTL_data2_FDR$SNP)))

#######
# Shared probes:
# Shared SNPs (returns SNP-probe pairs where shared SNPs appear) :
shared_probes <- getShared(eQTL_data1_FDR, eQTL_data2_FDR, 'Probe_ID', pvalue_col)
head(shared_probes)
dim(shared_probes)
length(unique(shared_probes$Probe_ID))

# Probes private to each group, returns SNP-probe pairs in data1 appear:
probes_eQTL_data1_FDR <- getPrivate(eQTL_data1_FDR, eQTL_data2_FDR, 'Probe_ID', pvalue_col)
head(probes_eQTL_data1_FDR)
dim(probes_eQTL_data1_FDR)
length(unique(probes_eQTL_data1_FDR$Probe_ID))

probes_eQTL_data2_FDR <- getPrivate(eQTL_data2_FDR, eQTL_data1_FDR, 'Probe_ID', pvalue_col)
head(probes_eQTL_data2_FDR)
dim(probes_eQTL_data2_FDR)
length(unique(probes_eQTL_data2_FDR$Probe_ID))
length(which(complete.cases(probes_eQTL_data2_FDR)))
length(which(complete.cases(shared_probes)))

# Sanity, the maximum number of private probes would be the FDR cut or private plus shared (use non-duplicated values):
length(unique(eQTL_data1_FDR$Probe_ID)) == length(unique(probes_eQTL_data1_FDR$Probe_ID))
length(unique(eQTL_data1_FDR$Probe_ID)) == length(unique(probes_eQTL_data1_FDR$Probe_ID)) + length(unique(shared_probes$Probe_ID))
length(unique(eQTL_data2_FDR$Probe_ID)) == length(unique(eQTL_data2_FDR$Probe_ID))
length(unique(eQTL_data2_FDR$Probe_ID)) == length(unique(probes_eQTL_data2_FDR$Probe_ID)) +length(unique(shared_probes$Probe_ID))
# Sanity, should be zero:
length(which(as.character(probes_eQTL_data1_FDR$Probe_ID) %in% as.character(probes_eQTL_data2_FDR$Probe_ID)))
####################

####################
# Add annotations:
shared_pairs <- get_illumina_annot(shared_pairs, 2)
shared_pairs
class(shared_pairs)

pairs_eQTL_data1_FDR  <- get_illumina_annot(pairs_eQTL_data1_FDR, 2)
pairs_eQTL_data1_FDR
pairs_eQTL_data2_FDR <- get_illumina_annot(pairs_eQTL_data2_FDR, 2)

shared_SNPs <- get_illumina_annot(shared_SNPs, 2)
SNPs_eQTL_data1_FDR <- get_illumina_annot(SNPs_eQTL_data1_FDR, 2)
SNPs_eQTL_data2_FDR <- get_illumina_annot(SNPs_eQTL_data2_FDR, 2)

shared_probes <- get_illumina_annot(shared_probes, 2)
probes_eQTL_data1_FDR <- get_illumina_annot(probes_eQTL_data1_FDR, 2)
probes_eQTL_data2_FDR <- get_illumina_annot(probes_eQTL_data2_FDR, 2)

####################

####################
# Save to file private and shared SNPs, probes and pairs:
shared_pairs <- shared_pairs[order(shared_pairs[, pvalue_col]), ]
shared_pairs
write.table(shared_pairs, sprintf('shared_pairs_%s', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

pairs_eQTL_data1_FDR <- pairs_eQTL_data1_FDR[order(pairs_eQTL_data1_FDR[, pvalue_col]), ]
pairs_eQTL_data1_FDR
write.table(pairs_eQTL_data1_FDR, sprintf('private_pairs_%s', eQTL_file1, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

pairs_eQTL_data2_FDR <- pairs_eQTL_data2_FDR[order(pairs_eQTL_data2_FDR[, pvalue_col]), ]
write.table(pairs_eQTL_data2_FDR, sprintf('private_pairs_%s', eQTL_file2, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

shared_SNPs <- shared_SNPs[order(shared_SNPs[, pvalue_col]), ]
write.table(shared_SNPs, sprintf('shared_SNPs_%s', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

SNPs_eQTL_data1_FDR <- SNPs_eQTL_data1_FDR[order(SNPs_eQTL_data1_FDR[, pvalue_col]), ]
write.table(SNPs_eQTL_data1_FDR, sprintf('private_SNPs_%s', eQTL_file1, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

SNPs_eQTL_data2_FDR <- SNPs_eQTL_data2_FDR[order(SNPs_eQTL_data2_FDR[, pvalue_col]), ]
write.table(pairs_eQTL_data2_FDR, sprintf('private_SNPs_%s', eQTL_file2, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

shared_probes <- shared_probes[order(shared_probes[, pvalue_col]), ]
write.table(shared_probes, sprintf('shared_probes_%s', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

probes_eQTL_data1_FDR <- probes_eQTL_data1_FDR[order(probes_eQTL_data1_FDR[, pvalue_col]), ]
write.table(probes_eQTL_data1_FDR, sprintf('private_probes_%s', eQTL_file1, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

probes_eQTL_data2_FDR <- probes_eQTL_data2_FDR[order(probes_eQTL_data2_FDR[, pvalue_col]), ]
write.table(probes_eQTL_data2_FDR, sprintf('private_probes_%s', eQTL_file2, output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
####################


####################
# Save results of counts above to file:
# TO DO: Warning that percentages are extracted from file read but may not be all the tested 
# elements (eg if MxEQTL wasn't run with full p-value).

options(digits = 2)
cat(file = sprintf('counts_%s.txt', output_file_name),
    paste('File1 is: '), '\t', paste(eQTL_file1), '\n',
    paste('File2 is: '), '\t', paste(eQTL_file2), '\n',
    paste('Number passed as total SNPs tested: '), '\t', total_SNPs_tested, '\n', 
    paste('Number passed as total probes tested: '), '\t', total_probes_tested, '\n',
    paste('Number passed as total SNP-probe pairs tested: '), '\t', total_SNPprobe_pairs_tested, '\n',
    # Numbers for individual files (also processed with eQTL_counting.R though):
    paste('Number of SNP-probe pairs in file1'), '\t', nrow(eQTL_data1), '\n', 
    paste('Number of SNP-probe pairs in file2'), '\t', nrow(eQTL_data2), '\n', 
    paste('Number of unique SNPs in file1'), '\t', length(unique(eQTL_data1$SNP)), '\n',
    paste('Number of unique SNPs in file2'), '\t', length(unique(eQTL_data2$SNP)), '\n',
    #  paste('unique genes %s', eQTL_file1), '\t', length(unique(eQTL_data1$gene)), '\n',
    #  paste('unique genes %s', eQTL_file2), '\t', length(unique(eQTL_data2$gene)), '\n',
    paste('Number of unique probes in file1'), '\t', length(unique(eQTL_data1$Probe_ID)), '\n',
    paste('Number of unique probes in file2'), '\t', length(unique(eQTL_data2$Probe_ID)), '\n',
    # Summary stats:
    paste('Summary stats of beta values in file1:'), '\t', t(summary(eQTL_data1[, "beta", with = F])), '\n',
    paste('Summary stats of beta values in file2:'), '\t', t(summary(eQTL_data2[, "beta", with = F])), '\n',
    paste('Summary stats of t-stat values in file1:'), '\t', t(summary(eQTL_data1[, "t-stat", with = F])), '\n',
    paste('Summary stats of t-stat values in file2:'), '\t', t(summary(eQTL_data2[, "t-stat", with = F])), '\n',
    paste('Summary stats of p-values file1:'), '\t', t(summary(eQTL_data1[, "p-value", with = F])), '\n',
    paste('Summary stats of p-values file2:'), '\t', t(summary(eQTL_data2[, "p-value", with = F])), '\n',
    paste('Lowest p-value in file1:'), '\t', min(eQTL_data1[, "p-value", with = F], na.rm = TRUE), '\n',
    paste('Lowest p-value in file2:'), '\t', min(eQTL_data2[, "p-value", with = F], na.rm = TRUE), '\n',
    paste('Summary stats of FDR values in file1:'), '\t', t(summary(eQTL_data1[, "FDR", with = F])), '\n',
    paste('Summary stats of FDR values in file2:'), '\t', t(summary(eQTL_data2[, "FDR", with = F])), '\n',
    paste('Lowest FDR value in file1:'), '\t', min(eQTL_data1[, "FDR", with = F], na.rm = TRUE), '\n',
    paste('Lowest FDR value in file2:'), '\t', min(eQTL_data2[, "FDR", with = F], na.rm = TRUE), '\n',
    #Counts under adj.pvalue specified:
    sprintf('Total number of SNP-probe pairs in file 1 under FDR %s%%:', pvalue*100), '\t', nrow(eQTL_data1_FDR), '\n', 
    sprintf('Total number of SNP-probe pairs in file 2 under FDR %s%%:', pvalue*100), '\t', nrow(eQTL_data2_FDR), '\n', 
    sprintf('Total number of (unique) SNPs in file 1 under FDR %s%%:', pvalue*100), '\t', length(unique(eQTL_data1_FDR$SNP)), '\n', 
    sprintf('Total number of (unique) SNPs in file 2 under FDR %s%%:', pvalue*100), '\t', length(unique(eQTL_data2_FDR$SNP)), '\n', 
    sprintf('Total number of (unique) probes in file 1 under FDR %s%%:', pvalue*100), '\t', length(unique(eQTL_data1_FDR$Probe_ID)), '\n', 
    sprintf('Total number of (unique) probes in file 2 under FDR %s%%:', pvalue*100), '\t', length(unique(eQTL_data2_FDR$Probe_ID)), '\n', 
    # Shared elements at given pvalue:
    sprintf('File1 SNPs found in file2 (both at FDR < %s%%):', pvalue*100), '\t', length(which(unique(eQTL_data1_FDR[, SNP]) %in% unique(eQTL_data2_FDR[, SNP]))), '\n',
    paste('Percentage of file1 SNPs found in file2 (of total SNPs tested):'), '\t', (length(which(unique(eQTL_data1_FDR[, SNP]) %in% unique(eQTL_data2_FDR[, SNP]))) / total_SNPs_tested) * 100, '\n',
    paste('Percentage of shared SNPs (trans which are also cis of trans total):'), '\t', (length(which(unique(eQTL_data1_FDR[, SNP]) %in% unique(eQTL_data2_FDR[, SNP]))) / length(unique(eQTL_data1_FDR$SNP))) * 100, '\n',
    sprintf('File1 probes found in file2 (both at FDR < %s%%):', pvalue*100), '\t', length(which(unique(eQTL_data1_FDR[, Probe_ID]) %in% unique(eQTL_data2_FDR[, Probe_ID]))), '\n', 
    paste('Percentage of file1 probes found in file2  (of total probes tested):'), '\t', (length(which(unique(eQTL_data1_FDR[, Probe_ID]) %in% unique(eQTL_data2_FDR[, Probe_ID]))) / total_probes_tested) * 100, '\n',
    paste('Percentage of shared probes (trans which are also cis of trans total):'), '\t', (length(which(unique(eQTL_data1_FDR[, Probe_ID]) %in% unique(eQTL_data2_FDR[, Probe_ID]))) / length(unique(eQTL_data1_FDR$Probe_ID))) * 100, '\n',
    sprintf('File1 SNP-probe pairs found in file2 (shared pairs) both under FDR %s%%):', pvalue*100), '\t', nrow(shared_pairs), '\n', 
    paste('Percentage of file1 SNP-probe pairs found in file2 (shared pairs) both under FDR 5% (of total tested):'), '\t', (nrow(shared_pairs) / total_SNPprobe_pairs_tested) * 100, '\n', 
    sprintf('SNP-probe pairs private to file1 FDR < %s%%):', pvalue*100), '\t', nrow(pairs_eQTL_data1_FDR), '\n', 
    paste('Percentage of SNP-probe pairs private to file1 FDR < 5% (of total tested):'), '\t', (nrow(pairs_eQTL_data1_FDR) / total_SNPprobe_pairs_tested) * 100, '\n', 
    sprintf('SNP-probe pairs private to file2 FDR < %s%%):', pvalue*100), '\t', nrow(pairs_eQTL_data2_FDR), '\n',
    sprintf('Percentage of SNP-probe pairs private to file2 FDR < %s%%):', pvalue*100), '\t', (nrow(pairs_eQTL_data2_FDR) / total_SNPprobe_pairs_tested) * 100, '\n',
    append = FALSE)
####################


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