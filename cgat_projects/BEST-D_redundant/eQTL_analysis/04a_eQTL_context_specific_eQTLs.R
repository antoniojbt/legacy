#############################################
# eQTL analysis of context/treatment specific eSNPs
# Antonio J Berlanga-Taylor
# March 2016

# Input: eQTL files from MatrixeQTL, genes of interest

# Outputs: various plots and tables for eQTL associations.
# See libraries and packages required below

#############################################

#############################################
## Questions

# 
#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/functional_annotation.dir/')

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
R_session_saved_image <- paste('R_session_saved_image_eQTL_responseQTLs','.RData', sep='')
R_session_saved_image
####################


####################
# Load packages:
library(data.table)
library(illuminaHumanv4.db)
source('/ifs/devel/antoniob/projects/BEST-D/functions_for_MatrixeQTL.R')
# source('/Users/antoniob/Desktop/Downloads_to_delete/ifs_scripts_backup.dir/projects/BEST-D//functions_for_MatrixeQTL.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# MatrixeQTL files:
eQTL_file1 <- as.character(args[1])
# eQTL_file1 <- '2000+4000-baseline-1.eQTL_trans'

eQTL_file2 <- as.character(args[2])
# eQTL_file2 <- '2000+4000-12months-1.eQTL_trans'

# Genes of interest (eg VD gene list):
genes_of_interest_file <- as.character(args[3])
# genes_of_interest_file <- 'VD_genes.txt'
  
# List of SNPs in LD or of interest: 
# VD GWAS SNP list LD r^1:
SNPs_of_interest_file <- as.character(args[4])
# SNPs_of_interest_file <- 'VD_SNPs_GWAS_list.txt'

total_SNPs_tested <- as.numeric(args[5])
# total_SNPs_tested <- 477422

total_probes_tested <- as.numeric(args[6])
# total_probes_tested <- 14972

total_SNPprobe_pairs_tested <- as.numeric(args[7])
# total_SNPprobe_pairs_tested <- 5755203
####################

####################
# Set up file naming for outputs:
# to create: 2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis
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

output_file_name <- sprintf('%s-VS-%s.reQTL_%s', eQTL_file2_base, eQTL_file1_base, eQTL_file1_ext)
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
dim(eQTL_data2)

tables()

# eQTL_data1 <- as.data.frame(eQTL_data1)
# eQTL_data2 <- as.data.frame(eQTL_data2)
eQTL_data1
eQTL_data2

SNPs_of_interest_data <- fread(SNPs_of_interest_file, header = F, stringsAsFactors = FALSE)
class(SNPs_of_interest_data)
head(SNPs_of_interest_data)
genes_of_interest_data <- fread(genes_of_interest_file, header = F, stringsAsFactors = FALSE)
# genes_of_interest_data <- rbind('VDR', 'GC', 'CYP2R1', 'VDBP', 'CYP27B1', 'CYP27A1', 'DHCR7', 'NADSYN1')
head(genes_of_interest_data)
####################


####################
# TO DO: this works but takes a long time
# Get annotations:
# eQTL_data1 <- get_illumina_annot(eQTL_data1, 2)
# head(eQTL_data1)
# tail(eQTL_data1)
# dim(eQTL_data1)
# 
# eQTL_data2 <- get_illumina_annot(eQTL_data2, 2)
# head(eQTL_data2)
# tail(eQTL_data2)
# dim(eQTL_data2)
####################


####################
# Basic counts:
dim(eQTL_data1)
dim(eQTL_data2)

length(which(eQTL_data1$FDR < 0.05))
length(which(eQTL_data2$FDR < 0.05))
range(eQTL_data1$beta)
range(eQTL_data2$beta)

length(unique(eQTL_data1$SNP))
length(unique(eQTL_data2$SNP))
# Sanity check against number passed in args:
# TO DO: if else raise warning:
identical(length(unique(eQTL_data1$SNP)), as.integer(total_SNPs_tested))
identical(length(unique(eQTL_data2$SNP)), as.integer(total_SNPs_tested))

#length(unique(eQTL_data1$gene))
#length(unique(eQTL_data2$gene))

length(unique(eQTL_data1$Probe_ID))
length(unique(eQTL_data2$Probe_ID))
# Sanity check against number passed in args:
# TO DO: if else raise warning:
identical(length(unique(eQTL_data1$Probe_ID)), as.integer(total_probes_tested))
identical(length(unique(eQTL_data2$Probe_ID)), as.integer(total_probes_tested))

length(which(eQTL_data1[, SNP] %in% eQTL_data2[, SNP]))
length(which(eQTL_data1[, Probe_ID] %in% eQTL_data2[, Probe_ID]))

# Get values under FDR 5%:
eQTL_data1_FDR5 <- eQTL_data1[which(eQTL_data1$FDR < 0.05), ]
dim(eQTL_data1_FDR5)
eQTL_data2_FDR5 <- eQTL_data2[which(eQTL_data2$FDR < 0.05), ]
dim(eQTL_data2_FDR5)
####################

####################
# SNP-probe pairs in each file and shared:
# Number of SNP-probe pairs in each file, sanity check (no NAs should be present):
dim(eQTL_data1_FDR5)
length(which(complete.cases(eQTL_data1_FDR5)))
dim(eQTL_data2_FDR5)
length(which(complete.cases(eQTL_data2_FDR5)))
# Concatenate SNP-probe pairs in new column:
eQTL_data1_new_FDR5 <- within(eQTL_data1_FDR5,  pairs <- paste(SNP, Probe_ID, sep="-"))
eQTL_data1_FDR5
eQTL_data1_new_FDR5

eQTL_data2_new_FDR5 <- within(eQTL_data2_FDR5,  pairs <- paste(SNP, Probe_ID, sep="-"))
eQTL_data2_FDR5
eQTL_data2_new_FDR5

# Check how many pairs appear in both sets:
dim(eQTL_data1_FDR5)
dim(eQTL_data2_FDR5)
shared_pairs <- length(which(as.character(eQTL_data1_new_FDR5$pairs) %in% as.character(eQTL_data2_new_FDR5$pairs)))
shared_pairs
# SNP-probe pairs unique to each group:
eQTL_data1_unique_pairs <- eQTL_data1_new_FDR5[which(!as.character(eQTL_data1_new_FDR5$pairs) %in% as.character(eQTL_data2_new_FDR5$pairs)), ]
eQTL_data2_unique_pairs <- eQTL_data2_new_FDR5[which(!as.character(eQTL_data2_new_FDR5$pairs) %in% as.character(eQTL_data1_new_FDR5$pairs)), ]
eQTL_data1_unique_pairs <- eQTL_data1_unique_pairs[order(eQTL_data1_unique_pairs$FDR), ]
eQTL_data2_unique_pairs <- eQTL_data2_unique_pairs[order(eQTL_data2_unique_pairs$FDR), ]
eQTL_data1_unique_pairs
eQTL_data2_unique_pairs
dim(eQTL_data1_unique_pairs)
dim(eQTL_data2_unique_pairs)
dim(eQTL_data1_FDR5)
dim(eQTL_data2_FDR5)
shared_pairs
# Sanity check, unique and total pairs should be TRUE:
identical(shared_pairs + nrow(eQTL_data1_unique_pairs), nrow(eQTL_data1_FDR5))
identical(shared_pairs + nrow(eQTL_data2_unique_pairs), nrow(eQTL_data2_FDR5))
length(which(as.character(eQTL_data1_unique_pairs$pairs) %in% as.character(eQTL_data2_unique_pairs$pairs)))
####################


####################
# Treatment (group 2) only SNPs (those with FDR < 0.05 in group 2 and NA or FDR > 0.5 in group 1):
# So 'context-specific':
# Use data.table join based on 2nd group. x[y] is an inner join (keeps those which match from both tables):
# https://rstudio-pubs-static.s3.amazonaws.com/52230_5ae0d25125b544caab32f75f0360e775.html

# eQTL_data_joint <- eQTL_data2[eQTL_data1]
# eQTL_data_joint
# View(eQTL_data_joint)


# Get context specific eSNPs:
# Full outer join will keep everything and fill in with NAs:
all_eQTL_data_joint <- merge(eQTL_data1, eQTL_data2, all = TRUE)
all_eQTL_data_joint #'.x' are baseline and '.y' are treated
dim(all_eQTL_data_joint)
length(which(all_eQTL_data_joint$FDR.y < 0.05)) # Number of SNP-probe pairs under FDR 5%
length(which(is.na(all_eQTL_data_joint$FDR.x))) # Number of SNP-probe pairs which are NA for baseline
range(all_eQTL_data_joint$FDR.x, na.rm = TRUE)
range(all_eQTL_data_joint$FDR.y, na.rm = TRUE)
range(all_eQTL_data_joint$`p-value.x`, na.rm = TRUE)
range(all_eQTL_data_joint$`p-value.y`, na.rm = TRUE)
length(which(all_eQTL_data_joint$`p-value.x` < 10e-5))
length(which(all_eQTL_data_joint$`p-value.y` < 10e-5))

# Keep only SNP-probe pairs which are context specific (FDR < 5% in group 2 and absent or FDR > 50% in group 2):
eQTL_tx_fdr5 <- all_eQTL_data_joint[which(all_eQTL_data_joint$FDR.y < 0.05), ]
eQTL_tx_fdr5
dim(eQTL_tx_fdr5)
# Order by significance:
eQTL_tx_fdr5 <- eQTL_tx_fdr5[order(eQTL_tx_fdr5$FDR.y, eQTL_tx_fdr5$FDR.x), ] # Ordered by FDR final, then baseline
eQTL_tx_fdr5

length(which(eQTL_tx_fdr5$FDR.x > 0.50))
length(which(eQTL_tx_fdr5$FDR.x > 0.50 | is.na(eQTL_tx_fdr5$FDR.x)))
eQTL_tx_fdr5_reQTLs <- eQTL_tx_fdr5[which(eQTL_tx_fdr5$FDR.x > 0.50 | is.na(eQTL_tx_fdr5$FDR.x)), ]

eQTL_tx_fdr5_reQTLs <- eQTL_tx_fdr5_reQTLs[order(eQTL_tx_fdr5_reQTLs$FDR.y), ]
eQTL_tx_fdr5_reQTLs
names(eQTL_tx_fdr5_reQTLs)[7] <- 'beta.final'
names(eQTL_tx_fdr5_reQTLs)[8] <- 't-stat.final'
names(eQTL_tx_fdr5_reQTLs)[9] <- 'p-value.final'
names(eQTL_tx_fdr5_reQTLs)[10] <- 'FDR.final'
names(eQTL_tx_fdr5_reQTLs)
eQTL_tx_fdr5_reQTLs
####################

####################
# Genes of interest in eSNPs:
# Get annotations:
eQTL_tx_fdr5_reQTLs_annot <- get_illumina_annot(eQTL_tx_fdr5_reQTLs, 2)
class(eQTL_tx_fdr5_reQTLs_annot)
eQTL_tx_fdr5_reQTLs_annot
dim(eQTL_tx_fdr5_reQTLs_annot)
# View(eQTL_tx_fdr5_reQTLs_annot)
eQTL_tx_fdr5_reQTLs_annot <- eQTL_tx_fdr5_reQTLs_annot[order(eQTL_tx_fdr5_reQTLs_annot$FDR.final), ]
eQTL_tx_fdr5_reQTLs_annot
####################

####################
# Save results of counts above to file:
options(digits = 2)
cat(file = sprintf('counts_%s.txt', output_file_name),
    paste('File1 is: '), '\t', paste(eQTL_file1), '\n',
    paste('File2 is: '), '\t', paste(eQTL_file2), '\n',
    paste('Number passed as total SNPs tested: '), '\t', total_SNPs_tested, '\n', 
    paste('Number passed as total probes tested: '), '\t', total_probes_tested, '\n',
    paste('Number passed as total SNP-probe pairs tested: '), '\t', total_SNPprobe_pairs_tested, '\n',
    # Look mainly at shared numbers only, numbers for individual files are processed with eQTL_counting.R
    paste('Number of SNP-probe pairs in file1:'), '\t', nrow(eQTL_data1), '\n', 
    paste('Number of SNP-probe pairs in file2:'), '\t', nrow(eQTL_data2), '\n', 
    paste('Number of SNP-probe pairs at FDR <5% in file1:'), '\t', length(which(eQTL_data1$FDR < 0.05)), '\n',
    paste('Number of SNP-probe pairs at FDR <5% in file2:'), '\t', length(which(eQTL_data2$FDR < 0.05)), '\n',
    paste('Number of unique SNPs in file1 under FDR 5%:'), '\t', length(unique(eQTL_data1_FDR5$SNP)), '\n',
    paste('Number of unique SNPs in file2 under FDR 5%:'), '\t', length(unique(eQTL_data2_FDR5$SNP)), '\n',
    #  paste('unique genes %s', eQTL_file1), '\t', length(unique(eQTL_data1$gene)), '\n',
    #  paste('unique genes %s', eQTL_file2), '\t', length(unique(eQTL_data2$gene)), '\n',
    paste('Number of unique probes in file1 under FDR 5%:'), '\t', length(unique(eQTL_data1_FDR5$Probe_ID)), '\n',
    paste('Number of unique probes in file2 under FDR 5%:'), '\t', length(unique(eQTL_data2_FDR5$Probe_ID)), '\n',
    # Stats for reQTL files:
    paste('Summary stats of beta values in reQTLs (file1):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "beta.x", with = F])), '\n',
    paste('Summary stats of beta values in reQTLs (file2):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "beta.final", with = F])), '\n',    
    paste('Summary stats of t-stat values in reQTLs (file1):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "t-stat.x", with = F])), '\n',
    paste('Summary stats of t-stat values in reQTLs (file2):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "t-stat.final", with = F])), '\n',
    paste('Summary stats of p-values in reQTLs (file1):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "p-value.x", with = F])), '\n',
    paste('Summary stats of p-values values in reQTLs (file2):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "p-value.final", with = F])), '\n',
    paste('Lowest p-value in reQTLs (file1):'), '\t', min(eQTL_tx_fdr5_reQTLs_annot[, "p-value.x", with = F], na.rm = TRUE), '\n',
    paste('Lowest p-value in reQTLs (file2):'), '\t', min(eQTL_tx_fdr5_reQTLs_annot[, "p-value.final", with = F], na.rm = TRUE), '\n',
    paste('Summary stats of FDR values in reQTLs (file1):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "FDR.x", with = F])), '\n',
    paste('Summary stats of FDR values in reQTLs (file2):'), '\t', t(summary(eQTL_tx_fdr5_reQTLs_annot[, "FDR.final", with = F])), '\n',
    paste('Lowest FDR value in reQTLs (file1):'), '\t', min(eQTL_tx_fdr5_reQTLs_annot[, "FDR.x", with = F], na.rm = TRUE), '\n',
    paste('Lowest FDR value in reQTLs (file2):'), '\t', min(eQTL_tx_fdr5_reQTLs_annot[, "FDR.final", with = F], na.rm = TRUE), '\n',
    # Shared elements:
    paste('File1 SNPs found in file2 (both at FDR < 5%):'), '\t', length(which(eQTL_data1_FDR5[, SNP] %in% eQTL_data2_FDR5[, SNP])), '\n',
    paste('Percentage of file1 SNPs found in file2 (of total SNPs tested):'), '\t', (length(which(eQTL_data1_FDR5[, SNP] %in% eQTL_data2_FDR5[, SNP])) / total_SNPs_tested) * 100 , '\n',
    paste('File1 probes found in file2 (both at FDR < 5%):'), '\t', length(which(unique(eQTL_data1_FDR5[, Probe_ID]) %in% unique(eQTL_data2_FDR5[, Probe_ID]))), '\n', 
    paste('Percentage of file1 probes found in file2  (of total probes tested):'), '\t', (length(which(unique(eQTL_data1_FDR5[, Probe_ID]) %in% unique(eQTL_data2_FDR5[, Probe_ID]))) / total_probes_tested) * 100, '\n',
    paste('File1 SNP-probe pairs found in file2 (shared pairs) both under FDR 5%':), '\t', shared_pairs, '\n', 
    paste('Percentage of file1 SNP-probe pairs found in file2 (shared pairs) both under FDR 5% (of total tested):'), '\t', (shared_pairs / total_SNPprobe_pairs_tested) * 100, '\n', 
    paste('SNP-probe pairs unique to file1 FDR < 5%:'), '\t', nrow(eQTL_data1_unique_pairs), '\n', 
    paste('Percentage of SNP-probe pairs unique to file1 FDR < 5% (of total tested):'), '\t', (nrow(eQTL_data1_unique_pairs) / total_SNPprobe_pairs_tested) * 100, '\n', 
    paste('SNP-probe pairs unique to file2 FDR < 5%:'), '\t', nrow(eQTL_data2_unique_pairs), '\n',
    paste('Percentage of SNP-probe pairs unique to file2 FDR < 5%:'), '\t', (nrow(eQTL_data2_unique_pairs) / total_SNPprobe_pairs_tested) * 100, '\n',
    paste('SNP-probe pairs FDR >50% in file1 in reQTLs (equals SNP-probe pairs unless pvalue <1 was run):'), '\t', length(which(eQTL_tx_fdr5_reQTLs_annot$FDR.x > 0.50)), '\n',
    paste('SNP-probe pairs absent in file1 and present in file2:'), '\t', length(which(is.na(eQTL_tx_fdr5_reQTLs_annot$FDR.x))), '\n',
    append = FALSE)
####################


####################
# Some possible hits for genes of interest:
# Save results to file:
save_output <- capture.output(
  paste("eQTL_data1[which(grepl(pattern = 'CYP', x = eQTL_data1$gene)), ]"),
  eQTL_data1[which(grepl(pattern = 'CYP', x = eQTL_data1$gene)), ],
  paste("eQTL_data2[which(grepl(pattern = 'CYP', x = eQTL_data2$gene)), ]"),
  eQTL_data2[which(grepl(pattern = 'CYP', x = eQTL_data2$gene)), ],
  paste("eQTL_data1[which(grepl(pattern = 'VDR', x = eQTL_data1$gene)), ]"),
  eQTL_data1[which(grepl(pattern = 'VDR', x = eQTL_data1$gene)), ],
  paste("eQTL_data2[which(grepl(pattern = 'IL', x = eQTL_data2$gene)), ]"),
  eQTL_data2[which(grepl(pattern = 'IL', x = eQTL_data2$gene)), ],
  paste("eQTL_data1[which(grepl(pattern = 'VDR', x = eQTL_data1$gene)), ]"),
  eQTL_data1[which(grepl(pattern = 'VDR', x = eQTL_data1$gene)), ],
  paste("eQTL_data2[which(grepl(pattern = 'IL', x = eQTL_data2$gene)), ]"),
  eQTL_data2[which(grepl(pattern = 'IL', x = eQTL_data2$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'VD', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'VD', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'IL', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'IL', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'CYP', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'CYP', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'IFN', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'IFN', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'COX', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'COX', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'TOR', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'TOR', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'PIP', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'PIP', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ],
  paste("eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'TSC', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]"),
  eQTL_tx_fdr5_reQTLs_annot[which(grepl(pattern = 'TSC', x = eQTL_tx_fdr5_reQTLs_annot$gene)), ]
)
save_output
# Save results to file:
write.table(save_output, sprintf('select_SNPs_genes_of_interest_%s.txt', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, 
            row.names = FALSE)
####################

####################
# Gene hits against list provided:
#genes_of_interest_data_group1 <- eQTL_data1[which(eQTL_data1[, gene] %in% genes_of_interest_data[, V1]), ]
#genes_of_interest_data_group1

#genes_of_interest_data_group2 <- eQTL_data2[which(eQTL_data2[, gene] %in% genes_of_interest_data[, V1]), ]
#genes_of_interest_data_group2

genes_of_interest_hits_joint <- eQTL_tx_fdr5_reQTLs_annot[which(genes_of_interest_data[, V1] %in% eQTL_tx_fdr5_reQTLs_annot[, gene]), ]
genes_of_interest_hits_joint
genes_of_interest_hits_joint[, 'SNP', with = F]

# Wang et al GWAS, rs2060793 at chr11:14915310 in CYP2R1, these all locate on chr11, prob same haplotype
# No pubmed result COPB1 and vitamin D; PMID:26780889 

# SNPs of interest in eSNPs against list provided:
SNPs_of_interest_hits_file1 <- eQTL_data1[which(eQTL_data1[, SNP] %in% SNPs_of_interest_data[, V1]), ]
SNPs_of_interest_hits_file2 <- eQTL_data2[which(eQTL_data2[, SNP] %in% SNPs_of_interest_data[, V1]), ]
SNPs_of_interest_hits_file1
SNPs_of_interest_hits_file2

SNPs_of_interest_hits_spec <- eQTL_tx_fdr5_reQTLs[which(eQTL_tx_fdr5_reQTLs[, SNP] %in% SNPs_of_interest_data[, V1]), ]
SNPs_of_interest_hits_spec

# Save results to file:

#write.table(genes_of_interest_data_group1, sprintf('%s.txt', deparse(substitute(genes_of_interest_data_group1))),
#            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

#write.table(genes_of_interest_data_group2, sprintf('%s.txt', deparse(substitute(genes_of_interest_data_group2))),
#            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(genes_of_interest_hits_joint, sprintf('%s_%s.txt', deparse(substitute(genes_of_interest_hits_joint)), output_file_name),
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(SNPs_of_interest_hits_file1, sprintf('%s_%s.txt', deparse(substitute(SNPs_of_interest_hits_file1)), output_file_name),
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(SNPs_of_interest_hits_file2, sprintf('%s_%s.txt', deparse(substitute(SNPs_of_interest_hits_file2)), output_file_name),
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

write.table(SNPs_of_interest_hits_spec, sprintf('%s_%s.txt', deparse(substitute(SNPs_of_interest_hits_spec)), output_file_name),
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)

####################


####################
# Write results to disk for further processing:
options(digits = 2)
eQTL_tx_fdr5_reQTLs_annot
output_file_name
write.table(eQTL_tx_fdr5_reQTLs_annot, output_file_name, 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
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

