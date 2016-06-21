#############################################
# eQTL analysis with MatrixQTL
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

# These must be in the format that Matrix eQTL needs (SNPs/genes/covariates in rows, individuals in columns), \
# see below for file formatting and conversions.
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
# Run R script for BEST-D data first.

#############################################

#############################################
## eQTL analysis comparisons for BEST-D:

# 1) Diff. exp. changes between groups, then eQTL (maybe better for dose effects).
# 2) eQTL in each group separately, then compare eQTLs for presence/absence.
# 3) Difference in difference estimators
# 4) 

#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')

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
R_session_saved_image <- paste('R_session_saved_image_MatrixQTL_analysis','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
# args[1] = 'genotype_data_placebo_baseline.tsv'
# args[2] = 'subset_baseline_placebo.tsv'
# args[3] = 'covar_PCAs_t_placebo_baseline.tsv'

#############################################


#############################################
## Update packages if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse'))
## Download Peter's scripts for processing scripts for MatrixQTL:

# devtools::install_github('humburg/Rsge')
# devtools::install_github('jknightlab/mePipe')
# devtools::install_github('hadley/readr')

#and load them:
packages <-c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse', 'Rsge', 'mePipe', 
             'readr', 'reshape2', 'data.table')
lapply(packages, require, character.only = TRUE)
#############################################


#############################################
## Start actual analysis once files are read, processed and subsets created.

# Calculate and assign minor allele frequency labels:
# Check how alleles were encoded. With plink --recode A-tranpose, A1 alleles are counted as default:
# https://www.cog-genomics.org/plink2/data#recode
# A1 alleles should correspond to the minor allele:
# http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml


# summary(geno2_matched[,2:15])
# maf <- colMeans(geno2_matched[,-1], na.rm=TRUE)/2
# which(maf > 0.50)
# maf2 <- pmin(maf, 1-maf)
# maf2[1:10]
# maf[1:10]
#############################################

#############################################
# Run MatrixEQTL:
# See mini-tutorials:
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html
# https://adairama.wordpress.com/2013/10/25/278/

## Get specific SNPs:
# SNP_file_name <- genotype_data_4000_final
# SNPs_of_interest <- c('rs2060793', 'rs1993116', 'rs3829251')
# SNPs_of_interest
# # rownames(genotype_data_placebo_final)
# index_SNPs <- which(rownames(SNP_file_name) %in% SNPs_of_interest)


# Set up MatrixEQTL
useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file = as.character(args[1])
# SNP_file = 'genotype_data_4000_baseline.tsv'
SNP_file
expression_file = as.character(args[2])
# expression_file = 'subset_baseline_4000.tsv'
expression_file

# Specify name of condition (if running multiple files/groups/tissues/etc.):
tissue <- SNP_file
tissue

# A separate file may be provided with extra covariates. 
#covariates_file = as.character(args[3])
#covariates_file
# covariates_file = 'covar_PCAs_t_4000_baseline.tsv'

# In case of no covariates set the variable covariates_file to character().
covariates_file = character()

# Load annotations for SNP location, probe annotation and probe location:
# snpspos = as.character(args[4])
# genepos = as.character(args[5])

snpspos = 'SNP_genomic_locations.txt'
genepos = 'illumina_probes_genomic_locations_loose.txt'
# snpspos = 'test_snp142_MatrixEQTL_snp_pos.txt'
# genepos = 'test_illumina_probes_genomic_locations.txt'
snpspos
genepos

snpPos <- fread(snpspos, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
head(snpPos)
dim(snpPos)

probePos <- fread(genepos, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
head(probePos)
dim(probePos)
# probeAnnot <- readr::read_tsv('')

# Switch classes as these are read by data.table and cause problems later:
snpPos <- as.data.frame(snpPos)
class(snpPos)

probePos <- as.data.frame(probePos)

# Set output file names:
output_file_name = paste(SNP_file, '_MatrixEQTL_loose_p8_1MB_no_covar.trans', sep = '')
output_file_name
#'./4000_baseline_MatrixEQTL.trans'
output_file_name.cis = paste(SNP_file, '_MatrixEQTL_loose_p5_1MB_no_covar.cis', sep = '')
output_file_name.cis
#'./4000_baseline_MatrixEQTL.cis'


# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
# Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset 
# may cause excessively large output files.

pvOutputThreshold = 1e-8
pvOutputThreshold.cis = 1e-5

# Define the covariance matrix for the error term. 
# This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

errorCovariance = numeric()

# Load the files:
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file)
# snps$CreateFromMatrix(data.matrix(SNP_file)); # Data need to be passed as a matrix


gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file);
# gene$CreateFromMatrix(data.matrix(expression_file));

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "-99999"; # denote missing values; This is from the plink encoding.
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
#cvrt$LoadFile(covariates_file);
cvrt$CreateFromMatrix(data.matrix(covariates_file));

# Check files:
str(cvrt)
str(snps)
str(gene)


# Finally, the main Matrix eQTL function is called:

# Set parameters for analysis for distance window and p-value threshold (1 Mb, 10^-3 (cis), 10^-5 (trans):
# Consider an error covariance matrix if necessary (correlated variables or errors)
# ?Matrix_eQTL_main

me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  output_file_name.cis = output_file_name.cis,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = 'qqplot',
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE,
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpPos,
  genepos = probePos,
  cisDist = 1e6) # TO DO: set cis distance as parameter, 1e5 = 100kb

# unlink(output_file_name)
# unlink(output_file_name.cis)

show(me$all$eqtls)
# qplot(me$all$eqtls)
# plot(me)

#run_me <- assign(me, paste('MeQTL_results_', snps, sep = '' ))
#run_me


## Save degrees of freedom in order to be able to run multi-tissue Matrix EQTL:
cat(file = 'df.txt', tissue, "\t", me$param$dfFull, '\n', append = TRUE)


# Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. 
# In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans. 
# Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.


# TO DO: run with trailing args so as to input files (covar, expr, geno) from the command line and submit as batch.
# Run all subsets:

# expression_file = subset_final_placebo
# output_file='./MatrixEQTL.subset_baseline_placebo.trans'
# output_file.cis='./MatrixEQTL.subset_baseline_placebo.cis'
# gene = SlicedData$new();
# gene$LoadFile('subset_baseline_placebo.tab')
# me
 
# # TO DO: clean up / test:
# # Inspect results:
# results <- readr::read_tsv('MatrixEQTL.trans')
# head(me$all$eqtls)
# plot(me)
# qqnorm(me$all$eqtls[, 4])
# head(results$cis$eqtls)
# nrow($cis$eqtls)
# head($cis$ntests)
# 
# head($trans$eqtls)
# nrow($trans$eqtls)
# head($trans$ntests)
# 
# # Include a covariate matrix in the analysis (use PCs generated previously, how to specify how many?):
# covar_MatrixQTL <- SlicedData$new()
# covar_MatrixQTL$CreateFromMatrix(t(pc[,1:10]))
# covar_MatrixQTL
# eQTL.pc10 <- Matrix_eQTL_main(snps, genes, cvrt=covar_MatrixQTL, 
#                               output_file_name = './eQTL.pc10.trans', 
#                               output_file_name.cis = './eQTL.pc10.cis', 
#                               pvOutputThreshold.cis = 1e-3, snpspos = as.data.frame(snpPos), 
#                               genepos = as.data.frame(probePos))
# 
# # Inspect results:
# head(eQTL.pc10$cis$eqtls)
# nrow(eQTL.pc10$cis$eqtls)
# head(eQTL.pc10$cis$ntests)
# 
# head(eQTL.pc10$trans$eqtls)
# nrow(eQTL.pc10$trans$eqtls)
# head(eQTL.pc10$trans$ntests)


#############################################
## Basic counting:
# TO DO:
# file_snp <- 'MatrixEQTL.trans'
# 
# read_file <- read.table(file_snp, header=T, sep='\t')
# 
# head(read_file)
# 
# unique_snps <- unique(read_file$SNP)
# head(unique_snps)
# length(unique_snps)
# 
# unique_genes <- unique(read_file$gene)
# head(unique_genes)
# length(unique_genes)
# 
# merged_gene <- merge(unique_genes, read_file)
# head(merged_gene)
# length(merged_gene$SNP)
#############################################


#############################################
## Plot results and provide summaries:
# Check circos plots for overall view
# Plot cis and trans
# Descriptives on p-values and effect sizes
# Boxplots for SNPs/genes of interest
# 
#############################################


#############################################
## Interpreting eQTL results
# GWAS catalogue
# Blueprint
# ENCODE
# Visualise in context in UCSC genome browser
# Consider imputing SNPs in area of interest
# Check LD structure in area of interest
# Check disease, evolution, context specificity, etc.
# Consider specific mechanisms of function of variants (eg TF and microRNA binding sites)
# Check super-enhancer, enhancer, promoter, etc.  areas for SNPs of interest
# Master regulatory loci


#############################################


#############################################
## Consider further downstream analysis
# Co-expression networks - WGCNA
# Pathway analysis - GO and IPA
# Compare to other tissues, conditions, datasets:
# http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/

#############################################

# TO DO:
# Bayesian and non-parametric eQTL analysis:
# http://www.bioconductor.org/packages/release/bioc/vignettes/iBMQ/inst/doc/iBMQ.pdf
# See also Panama:
# http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002330
# http://pmbio.github.io/envGPLVM/


#############################


#############################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])


#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run the script for xxx.
#############################
