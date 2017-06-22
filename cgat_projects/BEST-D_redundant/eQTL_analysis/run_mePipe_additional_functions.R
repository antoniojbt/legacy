#################
# Script to run mePipe additional options such as multi-peak, LD, etc.
# Antonio J Berlanga-Taylor
# 20 May 2016
# Requires output from mePipe (from run_mePipe.R script)
# Outputs plots and files to disk
#################

#################
# See:
# https://github.com/jknightlab/mePipe
#################


#################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/mePipe_runs.dir')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/mePipe_results')

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

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_mePipe_additional','.RData', sep='')
R_session_saved_image
####################

#################
# Load libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("MatrixEQTL", "trio", "optparse", "XML", "snow"))
# library(devtools)
# biocLite('optparse')
# # devtools::install_github("humburg/Rsge")
# system("wget https://github.com/jknightlab/mePipe/releases/download/mePipe_v1.3.5/mePipe_1.3.5.tar.gz")
# install.packages("mePipe_1.3.5.tar.gz", repos=NULL)
library(mePipe)
library(ggplot2)
library(data.table)
# library(MatrixEQTL)
source('run_mePipe_functions.R')
#################

#################
# Switch Rsge off:
sge.options(sge.use.cluster = FALSE) #sge.user.options=opt$sgeoptions)
# Also check the following as errors in cgat150
# make_option("--sgeoptions", default="-S /bin/bash -V",
# help="Options to pass to qsub (ignored unless `--cluster` is used. [default: %default])"),
# mePipe::
mePipe::getOptions()

#################


####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# File with maximum number of eQTL results based on covariate loop are soft linked 
# and start with cis/trans _covSelect...

eQTL_file <- as.character(args[1])
# eQTL_file <- 'cis_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv'

# trans_file <- as.character(args[2])
# trans_file <- 'trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv'

geno_file <- as.character(args[2])
# geno_file <- 'cut_genotype_data_all_treated_baseline.tsv_matched.tsv'

# Gene expression levels:
Gex_file <- as.character(args[3])
# Gex_file <- 'cut_GEx_baseline_4000_and_2000.tsv_matched.tsv'

covars_file <- as.character(args[4])
# covars_file <- 'PCs_to_adjust_for_cut_genotype_data_all_treated_baseline.tsv_matched.tsv.txt'

snp_pos = as.character(args[4])
# snp_pos = 'snp146Common_MatrixEQTL_snp_pos.txt'
##################

##################
# Read files:
eQTL_data <- fread(eQTL_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
eQTL_data
dim(eQTL_data)
# mePipe script asks for 'snps' column (instead of 'SNP' as output by MatrixEQTL):
setnames(eQTL_data, 'SNP', 'snps')
eQTL_data

geno_data <- fread(geno_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
geno_data
dim(geno_data)

Gex_data <- fread(Gex_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
Gex_data
dim(Gex_data)

snp_pos_data <- fread(snp_pos, sep = ' ', header = FALSE, stringsAsFactors = FALSE)
snp_pos_data
dim(snp_pos_data)
# mePipe LD blocks functions require columns to be 'chrom' and 'pos', 'pos' being rsID:
# setnames(snp_pos_data, 'V1', 'rsID')
setnames(snp_pos_data, 'V1', 'snp')
setnames(snp_pos_data, 'V2', 'chrom')
# Cut third column:
snp_pos_data[, V3:=NULL]
# Re-order columns:
# setcolorder(snp_pos_data, c("pos", "chrom"))
snp_pos_data
# Convert to data.frame only:
# snp_pos_data <- as.data.frame(snp_pos_data)
# head(snp_pos_data)

# mePipe getLDblocks() requires rownames to be the SNP positions:
# rownames(snp_pos_data) <- snp_pos_data$pos
# head(rownames(snp_pos_data))
# head(snp_pos_data)
#################


#################
# Obtain model summaries for Matrix-eQTL model
# Errors:
# See:
# https://github.com/jknightlab/mePipe/blob/master/R/utils.R 
#################

#################
# R options for reading files are required for 
# getEffectSize()
# but give errors when passing file name only.

# Original code is in:
# https://github.com/jknightlab/mePipe/blob/master/R/effectSize.R
# https://github.com/jknightlab/mePipe/blob/master/R/utils.R
# https://github.com/jknightlab/mePipe/blob/master/exec/runMatrixEQTL.R

# So loading files directly using MatrixEQTL function:
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(Gex_file);

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile(geno_file)

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values; This is from the plink encoding.
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covars_file);

# snp_pos_data = SlicedData$new();
# snp_pos_data$fileDelimiter = "\t";      # the TAB character
# snp_pos_data$fileOmitCharacters = "NA"; # denote missing values; This is from the plink encoding.
# snp_pos_data$fileSkipRows = 1;          # one row of column labels
# snp_pos_data$fileSkipColumns = 1;       # one column of row labels
# snp_pos_data$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
# snp_pos_data$LoadFile(snp_pos);
#################

#################
# Compute effect size for all SNPs with significant associations:
# ?getEffectSize
effect_sizes <- getEffectSize(hits = as.data.frame(eQTL_data), expression = gene, genotype = snps, covariate = cvrt, 
                              minFDR = 0.05)
# is(gene, "SlicedData")
# is(snps, "SlicedData")
# is(cvrt, "SlicedData")

# TO DO: errors:
# Error in geno[[as.character(hits$snps[i])]] : 
#   attempt to select less than one element in get1index
#################

#################
# Identify region of high LD around a given SNP and summarise eQTLs by LD block:
# ?getLDblocks
# Errors, seems circular, pbs with pos variable column names...
LD_blocks <- getLDblocks(eqtls = as.data.frame(eQTL_data), genotype = snps, pos = as.data.frame(snp_pos_data))
LD_blocks <- getLDblocks(eqtls = eQTL_data, genotype = geno_file, pos = snp_pos_data)

eqtls <- eQTL_data
pos <- snp_pos_data
head(rownames(pos))
head(which(tolower(colnames(pos)) == "snp" | tolower(colnames(pos)) == "id"))

idx <- which(tolower(colnames(pos)) == "snp" | tolower(colnames(pos)) == "id")
if(length(idx) == 1){
  rownames(pos) <- as.character(pos[,idx])
  
  chroms <- unique(as.character(pos[rownames(pos) %in% eqtls$snps, "chrom"]))
  chroms


# make_option(c("--maxDist"), default=1e5L,
#             help="Maximum distance allowed between two adjacent SNPs within an LD block. [default: %default]"),

# make_option(c("--maxSNPs"), default=200L,
#             help="Maximum number of SNPs to consider for each LD block. [default: %default]"),

# make_option(c("--ldFDR"), default=0.05,
#             help="Maximum FDR of eQTLs to be included in list of SNPs for each block. Only blocks with at least one SNP significant at this level will be reported. [default: %default]"),

# make_option(c("--ldPvalR2"),
#             help="Maximum p-value for correlation between two eSNPs for that should be considered part of the same signal by `ldPairs`. Use this option to test for significant correlation between SNPs instead of using the fixed R^2 threshold. Note that this may group SNPs with relatively low LD together, especially if the samplesize is large."),

# make_option(c("--ldR2"), default=0.8,
#             help="Cut-off for R^2 between two eSNPs. All SNPs that have a higher R^2 with a peak SNP will be considered part of the same signal. [default: %default]"),

# make_option(c("--ldOnly"), action="store_true", default=FALSE,
#             help="Compute LD blocks for existing eQTL results. This assumes that previous results can be loaded from the file implied by '--output'. Implies '--ldBlocks'"),
#################

#################
# # getLDpairs
# make_option(c("--ldPairs"), action="store_true", default=FALSE,
#             help="Flag indicating whether eQTLs should be grouped based on pairwise measures of LD. This differs from `ldBlocks` in that all SNPs that have R^2 > `ldR2` with the peak SNP will be grouped together and no attempt is made to infer local LD patterns. [default: %default]"),
#################

#################
# # getMultiPeak
# make_option(c("--multiPeak"), action="store_true", default=FALSE,
#             help="Flag indicating whether an attempt should be made to resolve multiple independent eSNPs for the same gene via multiple regression. [default: %default]"),

# make_option(c("--multiPvalue"), default=1e-6,
#             help="P-value threshold to use for associations between secondary eSNPs and gene expression measurements. [default: %default]"),

# make_option(c("--multiGene"), meta="GENE",
#             help="Restrict analysis to this gene (for analyses that support this feature)"),

#################

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

# Next run:
####################