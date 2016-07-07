####################
# Script to obtain ratio between two files of gene expression values
# 05 June 2016
####################

####################
# See:

####################

####################
options(echo = TRUE)
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/gex_FC_tests/')

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
# load('R_session_saved_image_XGR.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_gex_FC_ratios','.RData', sep='')
R_session_saved_image
####################

####################
# Libraries:
library(data.table)
library(gvlma)
source('functions_for_MatrixeQTL.R')
source('moveme.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene expression files:
gex_file1 <- as.character(args[1])
gex_file1 <- 'GEx_baseline_4000_and_2000.tsv'

gex_file2 <- as.character(args[2])
gex_file2 <- 'GEx_treated_4000_and_2000.tsv'

pheno_file <- as.character(args[3])
pheno_file <- 'membership_file_cleaned_all.tab'

print(args)
####################


####################
# Read data:
gex_1 <- fread(gex_file1, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gex_1[1:5, 1:5, with = F]
dim(gex_1)

gex_2 <- fread(gex_file2, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
gex_2[1:5, 1:5, with = F]
dim(gex_2)

pheno_data <- fread(pheno_file, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
pheno_data[1:5, 1:5, with = F]
dim(pheno_data)
class(pheno_data)
colnames(pheno_data)
setkey(pheno_data)
####################

####################
# Match files in order to have the same individuals and order
# Needs matching kit IDs to patient IDs

# Transpose expression files to allow merging later (pheno data is in the right format):
gex_1[1:5, 1:5, with = F]
gex_1_t <- transpose_file(gex_1, 1)
setkey(gex_1_t)
gex_1_t[1:5, 'rownames', with = F]
class(gex_1_t)
gex_1_t[1:5, 1:5, with = F]
gex_1_t[1:5, c((ncol(gex_1_t)-5):ncol(gex_1_t)), with = F]
dim(gex_1_t)
dim(gex_1)

gex_2[1:5, 1:5, with = F]
gex_2_t <- transpose_file(gex_2, 1)
setkey(gex_2_t)
gex_2_t[1:5, 'rownames', with = F]
class(gex_2_t)
gex_2_t[1:5, 1:5, with = F]
gex_2_t[1:5, c((ncol(gex_2_t)-5):ncol(gex_2_t)), with = F]

length(which(colnames(gex_1) %in% colnames(pheno_IDs))) # First column is probe IDs
length(which(colnames(gex_2) %in% colnames(pheno_IDs))) # First column is probe IDs


# Merge files with pheno data to add pt IDs
# Rename columns to match each other:
pheno_data[1:5, 1:5, with = F]
setnames(pheno_data, 'V1', 'kit_id')
setnames(gex_1_t, 'rownames', 'kit_id')
setnames(gex_2_t, 'rownames', 'kit_id')
# Convert column classes:
gex_1_t <- gex_1_t[, kit_id:=as.character(kit_id)]
gex_2_t <- gex_2_t[, kit_id:=as.character(kit_id)]
pheno_data <- pheno_data[, kit_id:=as.character(kit_id)]
tables()
str(gex_1_t[, 'kit_id', with = F])
str(gex_2_t[, 'kit_id', with = F])
str(pheno_data[, 'kit_id', with = F])
# Changing column types seems to mess the keys:
setkey(gex_1_t, 'kit_id')
setkey(gex_2_t, 'kit_id')
setkey(pheno_data, 'kit_id')

# Check matches:
length(which(pheno_data[['kit_id']] %in% gex_1_t[['kit_id']])) # Must be tested as vectors
length(which(pheno_data[['kit_id']] %in% gex_2_t[['kit_id']])) # Must be tested as vectors

# Merge to get patient IDs with kit IDs for each expression file:
gex_1_pt_IDs <- gex_1_t[pheno_data, nomatch = 0]
dim(gex_1_pt_IDs)
gex_1_pt_IDs[1:5, 1:5, with = F]
gex_1_pt_IDs[1:5, c((ncol(gex_1_pt_IDs)-5):ncol(gex_1_pt_IDs)), with = F]

gex_2_pt_IDs <- gex_2_t[pheno_data, nomatch = 0]
dim(gex_2_pt_IDs)
gex_2_pt_IDs[1:5, 1:5, with = F]
gex_2_pt_IDs[1:5, c((ncol(gex_2_pt_IDs)-5):ncol(gex_2_pt_IDs)), with = F]
####################

####################
# Order and match gene expression files by patient ID:
# Change column classes:
gex_1_pt_IDs <- gex_1_pt_IDs[, pt_id:=as.character(pt_id)]
gex_2_pt_IDs <- gex_2_pt_IDs[, pt_id:=as.character(pt_id)]
pheno_data <- pheno_data[, pt_id:=as.character(pt_id)]
str(gex_1_pt_IDs[, 'pt_id', with = F])
str(gex_2_pt_IDs[, 'pt_id', with = F])
str(pheno_data[, 'pt_id', with = F])
tables()
# Set keys:
setkey(gex_1_pt_IDs, 'pt_id')
setkey(gex_2_pt_IDs, 'pt_id')
setkey(pheno_data, 'pt_id')

dim(gex_1_pt_IDs)
dim(gex_2_pt_IDs)

# Keep only files which appear in both:
gex_1_keep<- which(gex_1_pt_IDs[['pt_id']] %in% gex_2_pt_IDs[['pt_id']])
length(gex_1_keep)
gex_1_matched <- gex_1_pt_IDs[gex_1_keep, , ]
dim(gex_1_matched)
gex_1_matched[1:5, 1:5, with = F]
gex_1_matched[1:5, c((ncol(gex_1_matched) - 5):ncol(gex_1_matched)), with = F]

gex_2_keep<- which(gex_2_pt_IDs[['pt_id']] %in% gex_1_pt_IDs[['pt_id']])
length(gex_2_keep)
gex_2_matched <- gex_2_pt_IDs[gex_2_keep, , ]
dim(gex_2_matched)
gex_2_matched[1:5, 1:5, with = F]
gex_2_matched[1:5, c((ncol(gex_2_matched)-5):ncol(gex_2_matched)), with = F]

# Check all columns and order match:
# TO DO: stop if false:
identical(colnames(gex_1_matched), colnames(gex_2_matched))

# Keep only gene expression values:
gex_1_matched <- gex_1_matched[, -c((ncol(gex_1_matched) - (ncol(pheno_data) - 1)):ncol(gex_1_matched)), with = F]
gex_1_matched[1:5, 1:5, with = F]
gex_1_matched[1:5, c((ncol(gex_1_matched) - 5):ncol(gex_1_matched)), with = F]

gex_2_matched <- gex_2_matched[, -c((ncol(gex_2_matched) - (ncol(pheno_data) - 1)):ncol(gex_2_matched)), with = F]
gex_2_matched[1:5, 1:5, with = F]
gex_2_matched[1:5, c((ncol(gex_2_matched) - 5):ncol(gex_2_matched)), with = F]

dim(gex_1_matched)
dim(gex_2_matched)
dim(gex_1)
dim(gex_2)

# Check all columns and order match:
# TO DO: stop if false:
identical(colnames(gex_1_matched), colnames(gex_2_matched))
####################

####################
# Get PCs (from max MatrixeQTL loop) for each file:
# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca_gex_1_matched <- prcomp(gex_1_matched, center=TRUE, scale=TRUE)
pca_gex_2_matched <- prcomp(gex_2_matched, center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc_gex_1_matched <- data.frame(round(pca_gex_1_matched$x, 2))
pc_gex_2_matched <- data.frame(round(pca_gex_2_matched$x, 2))
# pc_gex_1_matched$sample_id <- rownames(pc_gex_1_matched) 
# pc_gex_2_matched$sample_id <- rownames(pc_gex_2_matched) 
dim(pc_gex_1_matched)
dim(pc_gex_2_matched)

# pc_gex_1_matched <- pc_gex_1_matched[, moveme(names(pc_gex_1_matched), 'sample_id first')]
# pc_gex_2_matched <- pc_gex_2_matched[, moveme(names(pc_gex_2_matched), 'sample_id first')]

names(pc_gex_1_matched)[1:10]
names(pc_gex_2_matched)[1:10]
class(pc_gex_1_matched)
class(pc_gex_2_matched)

# First column names are probe IDs and should match:
# identical(pc_gex_1_matched[, 1], pc_gex_2_matched[, 1]) # If adding rownames as extra column
identical(rownames(pc_gex_1_matched), rownames(pc_gex_2_matched))

# Explore dimensions and plot first 10 or so components:
dim(pc_gex_1_matched)
dim(pc_gex_2_matched)
str(pca_gex_1_matched)
pc_gex_1_matched[1:5, 1:5]
pc_gex_2_matched[1:5, 1:5]

summary(pca_gex_1_matched)
summary(pca_gex_2_matched)

# Plot PCA results:
plot_PCA <- ('plot_PCA_gex1_v_gex2.png')
png(plot_PCA, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_gex_1_matched, main = sprintf('Normalised expression values - %s', gex_file1))
plot(pca_gex_2_matched, main = sprintf('Normalised expression values - %s', gex_file2))
# Scatterplot of PC1 and PC2:
biplot(pca_gex_1_matched, main = sprintf('PC1 v PC2 - %s', gex_file1))
biplot(pca_gex_2_matched, main = sprintf('PC1 v PC2 - %s', gex_file2))
par(mfrow=c(1,1))
dev.off()

# Regress out PCs from gene expression values for each file
# MatrixeQTL loop PCs to regress:
# 25 baseline?
# 35 final?
PCs_to_adjust_1 <- 25
PCs_to_adjust_2 <- 35

pc_gex_1_matched_to_adjust <- pc_gex_1_matched[1:nrow(pc_gex_1_matched), 1:PCs_to_adjust_1]
pc_gex_2_matched_to_adjust <- pc_gex_2_matched[1:nrow(pc_gex_2_matched), 1:PCs_to_adjust_2]

pc_gex_1_matched_to_adjust[1:5, ]
pc_gex_2_matched_to_adjust[1:5, ]

dim(pc_gex_1_matched_to_adjust)
dim(pc_gex_2_matched_to_adjust)

gex_1_matched[1:5, 1:5, with = F]
gex_2_matched[1:5, 1:5, with = F]

dim(gex_1_matched)
dim(gex_2_matched)

class(pc_gex_1_matched)
class(gex_1_matched)

# Create dataframe with data (expression value plus PCs)
gex_1_matched_lm <- data.frame(cbind(gex_1_matched[, 1, with = F], pc_gex_1_matched_to_adjust))
# gex_2_matched_lm <- data.frame(cbind(gex_2_matched[, 1, with = F], pc_gex_2_matched_to_adjust))
gex_1_matched_lm[1:5, 1:5]
colnames(gex_1_matched_lm)[1] <- 'probe'

# Linear regression corrected for PCs:
pass_formula <- as.formula(sprintf('%s ~ .', 'probe'))
pass_formula
lm_gex_1_matched <- lm(formula = pass_formula, data = gex_1_matched_lm)
# lm_gex_2_matched <- lm(formula = pass_formula, data = gex_2_matched_lm)
summary.lm(lm_gex_1_matched)

# Test assumptions:
gvmodel <- gvlma(lm_gex_1_matched)
summary(gvmodel)

# Model comparison:
# AIC(lm_1, lm_2)
# anova(lm_1, lm_2)

# Get results:
coefficients(lm_gex_1_matched) # lm_fit_PCs$coefficients
str(lm_gex_1_matched)
# Get coefficients (estimate, t-stat, p-value) from lm:
lm_summary <- summary.lm(lm_gex_1_matched)
lm_summary$coefficients
lm_summary_coef <- as.data.frame(lm_summary$coefficients)
str(lm_summary_coef)

# Correct expression value:
gex_1_matched_lm[1:5, 1:5]
gex_corrected <- gex_1_matched_lm[, 'probe'] - rowSums(coef(lm_gex_1_matched)[c(-1)] * gex_1_matched_lm[, 2:ncol(gex_1_matched_lm)])
# gex_corrected_2 <- gex_2_matched_lm[, 'probe'] - rowSums(coef(lm_gex_2_matched)[c(-1)] * gex_2_matched_lm[, 2:ncol(gex_2_matched_lm)])
gex_corrected

# TO DO: Do the same for 16,700 probes and get ratios afterwards
gex_corrected_2 / gex_corrected
 
####################

####################
# Get ratio of file2 over file1:
gex_FC <- gex_2_matched / gex_1_matched
dim(gex_FC)
class(gex_FC)
gex_FC[1:5, 1:5, with = F]

gex_FC[, list = mean(gex_FC)]

####################

####################
# Write to disk:
gex_file1
gex_file2

write.table(gex_FC, sprintf('gex_FC_%s_over_%s', gex_file2, gex_file1), sep='\t', quote = FALSE, col.names = NA)
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

# Next: run the script for xxx.
####################
