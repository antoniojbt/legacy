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
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/gex_FC_tests/')

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
#############################################


#############################################
## Load libraries:
packages <-c('MatrixEQTL', 'data.table')
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
# useModel = modelLINEAR_CROSS # Should work better for Tx specific eQTLs
#useModel = modelANOVA

SNP_file = as.character(args[1])
# SNP_file = 'cut_genotype_data_all_treated_baseline.tsv_matched.tsv'
# SNP_file = 'genotype_2000+4000.tsv'
SNP_file

expression_file = as.character(args[2])
# expression_file = 'cut_GEx_baseline_4000_and_2000.tsv_matched.tsv'
# expression_file = 'cut_gex_FC_GEx_treated_4000_and_2000.tsv_over_GEx_baseline_4000_and_2000.tsv_matched.tsv'
# expression_file = 'FC_GEx_4000+2000.tsv'
expression_file

# Specify name of condition (if running multiple files/groups/tissues/etc.):
tissue <- SNP_file
tissue

# A separate file may be provided with extra covariates. 
covariates_file = as.character(args[3])
# covariates_file = 'cut_principal_components_GEx_baseline_matched.tsv'
# In case of no covariates set the variable covariates_file_name to character():
# covariates_file = character()
covariates_file

# Load annotations for SNP location, probe annotation and probe location:
snpspos = as.character(args[4])
genepos = as.character(args[5])
# snpspos = '~/Desktop/BEST_D.DIR/mac_runs_to_upload/annotation_external/snp146Common_MatrixEQTL_snp_pos.txt'
# snpspos = 'snp146Common_MatrixEQTL_snp_pos.txt'
# genepos = '~/Desktop/BEST_D.DIR/mac_runs_to_upload/annotation_internal/biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt'
# genepos = 'biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt'
snpspos
genepos

cisDist = as.numeric(args[6])
# cisDist = 1e6

# The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. 
# Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset 
# may cause excessively large output files.
pvOutputThreshold = as.numeric(args[7])
# pvOutputThreshold = 1e-8
pvOutputThreshold.cis = as.numeric(args[8])
# pvOutputThreshold.cis = 1e-5
#############################################


#############################################
# Set output file names:
expression_file_split <- strsplit(expression_file, '[.]')[[1]][1]
output_file_name = sprintf('%s_%s_%s_MxEQTL.trans', expression_file_split, pvOutputThreshold, cisDist)
output_file_name
#'./4000_baseline_MatrixEQTL.trans'
output_file_name.cis = sprintf('%s_%s_%s_MxEQTL.cis', expression_file_split, pvOutputThreshold.cis, cisDist)
output_file_name.cis
#############################################


#############################################
# Read files:
snpPos <- fread(snpspos, sep = ' ', header = TRUE, stringsAsFactors = FALSE)
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
#############################################


#############################################
# Set-up MatrixeQTL:
# Define the covariance matrix for the error term. 
# Consider an error covariance matrix if necessary (correlated variables or errors)
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
snps
# snps$CreateFromMatrix(data.matrix(SNP_file)); # Data need to be passed as a matrix

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile(expression_file);
gene
# gene$CreateFromMatrix(data.matrix(expression_file));

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values; This is from the plink encoding.
cvrt$fileSkipRows = 0;          # one row of column labels
cvrt$fileSkipColumns = 0;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file);
# cvrt$CreateFromMatrix(data.matrix(covariates_file)); # use this if covariates_file = character()

# Check files:
str(cvrt)
str(snps)
str(gene)

dim(cvrt)
dim(snps)
dim(gene)
#############################################

#############################################
# Call the main Matrix eQTL function
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  output_file_name.cis = output_file_name.cis,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = 'qqplot',
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE,
  pvOutputThreshold = pvOutputThreshold,
  pvOutputThreshold.cis = pvOutputThreshold.cis,
  snpspos = snpPos,
  genepos = probePos,
  cisDist = cisDist)

# unlink(output_file_name)
# unlink(output_file_name.cis)

# Inspect results:
# Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object me. 
# In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans. 
# Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.
# show(me$all$eqtls)
# head(me$all$eqtls)
# plot(me)
# qqnorm(me$all$eqtls[, 4])
# head(results$cis$eqtls)
# nrow($cis$eqtls)
# head($cis$ntests)
# head($trans$eqtls)
# nrow($trans$eqtls)
# head($trans$ntests)
#############################################

#############################################
# Plot:
png(paste('qqplot_', SNP_file, '.png', sep = ''))
plot(me)
dev.off()
#############################################

#############################################
## Save degrees of freedom in order to be able to run multi-tissue Matrix EQTL:
cat(file = 'df_anova.txt', tissue, "\t", me$param$dfFull, '\n', append = TRUE)
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


#############################################
# TO DO:
# Bayesian and non-parametric eQTL analysis:
# http://www.bioconductor.org/packages/release/bioc/vignettes/iBMQ/inst/doc/iBMQ.pdf
# See also Panama:
# http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002330
# http://pmbio.github.io/envGPLVM/
#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))
#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for downstream analyses, plotting, counting, etc.
#############################