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
# load('R_session_saved_image_gex_FC_ratios.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_gex_FC_ratios','.RData', sep='')
R_session_saved_image
####################

####################
# Libraries:
library(data.table)
library(ggplot2)
library(gridExtra)
library(gvlma)
library(biglm)

# Get script with functions needed:
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/moveme.R')
source('/ifs/devel/antoniob/projects/BEST-D/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/eQTL_analysis/functions_for_MatrixeQTL.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene expression files:
gex_file1 <- as.character(args[1])
# gex_file1 <- 'GEx_baseline_4000_and_2000.tsv'

gex_file2 <- as.character(args[2])
# gex_file2 <- 'GEx_treated_4000_and_2000.tsv'

pheno_file <- as.character(args[3])
# pheno_file <- 'membership_file_cleaned_all.tab'

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

# length(which(colnames(gex_1) %in% colnames(pheno_data))) # First column is probe IDs
# length(which(colnames(gex_2) %in% colnames(pheno_data))) # First column is probe IDs
length(which(colnames(gex_1) %in% pheno_data[['V1']])) # First column is probe IDs
length(which(colnames(gex_2) %in% pheno_data[['V1']])) # First column is probe IDs

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

# Check all columns and order match, pt_ids will as they are unique but kit_ids won't as these
# specify baseline and 12 months in BEST-D:
# TO DO: stop if false:
identical(colnames(gex_1_matched), colnames(gex_2_matched))
identical(gex_1_matched[['pt_id']], gex_2_matched[['pt_id']])
identical(gex_1_matched[['kit_id']], gex_2_matched[['kit_id']]) # Will be false

## Keep only gene expression values
# kit_ids are the sample IDs, these are kept as they were set as keys for data.table
# Keep sample IDs pt_ids as sanity though
# Delete the last rows that were added from the IDs and pheno file:
gex_1_matched[, c((ncol(gex_1_matched) - (ncol(pheno_data) - 3)):ncol(gex_1_matched)), with = F]
gex_1_matched <- gex_1_matched[, -c((ncol(gex_1_matched) - (ncol(pheno_data) - 3)):ncol(gex_1_matched)), with = F]
gex_1_matched[1:5, 1:5, with = F]
gex_1_matched[1:5, c((ncol(gex_1_matched) - 5):ncol(gex_1_matched)), with = F]
colnames(gex_1_matched)[1:5]
colnames(gex_1_matched)[ncol(gex_1_matched)]

gex_2_matched <- gex_2_matched[, -c((ncol(gex_2_matched) - (ncol(pheno_data) - 3)):ncol(gex_2_matched)), with = F]
gex_2_matched[1:5, 1:5, with = F]
gex_2_matched[1:5, c((ncol(gex_2_matched) - 5):ncol(gex_2_matched)), with = F]
colnames(gex_2_matched)[1:5]
colnames(gex_2_matched)[ncol(gex_2_matched)]

dim(gex_1_matched)
dim(gex_2_matched)
dim(gex_1)
dim(gex_2)

gex_1_matched[1:5, c((ncol(gex_1_matched) - 5):ncol(gex_1_matched)), with = F]
gex_2_matched[1:5, c((ncol(gex_2_matched) - 5):ncol(gex_2_matched)), with = F]

# Check all columns and order match:
# TO DO: stop if false:
identical(colnames(gex_1_matched), colnames(gex_2_matched))
identical(gex_1_matched[['pt_id']], gex_2_matched[['pt_id']])
####################

####################
# Get PCs (from max MatrixeQTL loop) for each file:
# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
gex_1_matched[, c('pt_id', 'kit_id'), with = F]
gex_2_matched[, c('pt_id', 'kit_id'), with = F]
pca_gex_1_matched <- prcomp(gex_1_matched[, -c('pt_id', 'kit_id'), with = F], center=TRUE, scale=TRUE)
pca_gex_2_matched <- prcomp(gex_2_matched[, -c('pt_id', 'kit_id'), with = F], center=TRUE, scale=TRUE)

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
names(pc_gex_1_matched)[ncol(pc_gex_1_matched)]
names(pc_gex_2_matched)[1:10]
names(pc_gex_2_matched)[ncol(pc_gex_2_matched)]
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

sum_pc_gex_1_matched <- summary(pca_gex_1_matched)
sum_pc_gex_2_matched <- summary(pca_gex_2_matched)

sum_pc_gex_1_matched$importance[, 1:10]
sum_pc_gex_2_matched$importance[, 1:10]
sum_pc_gex_1_matched_df <- as.data.frame(sum_pc_gex_1_matched$importance)
sum_pc_gex_2_matched_df <- as.data.frame(sum_pc_gex_2_matched$importance)

sum_pc_gex_1_matched_df <- t(sum_pc_gex_1_matched_df)
sum_pc_gex_2_matched_df <- t(sum_pc_gex_2_matched_df)

sum_pc_gex_1_matched_df <- as.data.frame(sum_pc_gex_1_matched_df)
sum_pc_gex_2_matched_df <- as.data.frame(sum_pc_gex_2_matched_df)

sum_pc_gex_1_matched_df$percent_var <- round(100 * (sum_pc_gex_1_matched_df$`Proportion of Variance`), 3)
sum_pc_gex_2_matched_df$percent_var <- round(100 * (sum_pc_gex_2_matched_df$`Proportion of Variance`), 3)

sum_pc_gex_1_matched_df$PC <- factor(row.names(sum_pc_gex_1_matched_df), levels = row.names(sum_pc_gex_1_matched_df),
                        labels = row.names(sum_pc_gex_1_matched_df))
sum_pc_gex_2_matched_df$PC <- factor(row.names(sum_pc_gex_2_matched_df), levels = row.names(sum_pc_gex_2_matched_df),
                        labels = row.names(sum_pc_gex_2_matched_df))

head(sum_pc_gex_1_matched_df)
tail(sum_pc_gex_1_matched_df)
sum_pc_gex_1_matched_df[1:20, ]
names(sum_pc_gex_1_matched_df)
str(sum_pc_gex_1_matched_df)

head(sum_pc_gex_2_matched_df)
tail(sum_pc_gex_2_matched_df)
sum_pc_gex_2_matched_df[1:20, ]
names(sum_pc_gex_2_matched_df)
str(sum_pc_gex_2_matched_df)
####################

####################
# Plot PCA results:
png('plot_PCA_gex1_v_gex2.png', width = 12, height = 12, units = 'in', res = 300)
# Plot proportion of variance of first x PCs:
p1 <- ggplot(sum_pc_gex_1_matched_df, aes(y = percent_var, x = PC)) + 
  geom_bar(stat = 'identity') +
  theme_classic() +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  ggtitle(gex_file1) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
p2 <- ggplot(sum_pc_gex_2_matched_df, aes(y = percent_var, x = PC)) + 
  geom_bar(stat = 'identity') +
  theme_classic() +
  labs(x = 'Principal component', y = 'Proportion of variance (%)') +
  ggtitle(gex_file2) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
# Scatterplot of PC1 and PC2:
p3 <- ggplot(data=pc_gex_1_matched, aes(x=PC1, y=PC2)) + 
  theme_bw() + geom_point(alpha = 0.8) +
  ggtitle(gex_file1)
p4 <- ggplot(data=pc_gex_2_matched, aes(x=PC1, y=PC2)) + 
  theme_bw() + geom_point(alpha = 0.8) +
  ggtitle(gex_file2)
grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()
####################

####################
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
####################

####################
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
####################

####################
# Do the same for 16,700 probes and get ratios
# Linear regression corrected for PCs:

# Function to run linear regression on each probe:
get_lm <- function(DT, PC_DT, col_name) {
  # Create dataframe with data (expression value plus PCs)
  DT2 <- data.frame(cbind(DT[, as.character(col_name), with = F], PC_DT))
  # Linear regression corrected for PCs:
  col_names <- names(DT2)
  pass_formula <- as.formula(sprintf('%s ~ %s', col_name, paste(col_names[!col_names %in% col_name], collapse = " + ")))
  biglm_fit <- biglm(formula = pass_formula, data = DT2)
  # summary(biglm_fit)
  # Get coefficients from lm:
  # biglm_fit$qr$thetab # These are the coefficients
  # Correct expression value:
  DT2_corrected <- DT2[, col_name] - rowSums(biglm_fit$qr$thetab[c(-1)] * DT2[, 2:ncol(DT2)])
  # DT2_corrected
  return(DT2_corrected)
}
####################

####################
# Run regressions for file gex 1:
# Get all probe names minus sample ID column:
colnames(gex_1_matched)[c(ncol(gex_1_matched))]
dim(gex_1_matched)
gex_1_matched[, c('pt_id', 'kit_id'), with = F]
gex_2_matched[, c('pt_id', 'kit_id'), with = F]

cols <- colnames(gex_1_matched)[1:c(ncol(gex_1_matched)-2)]
length(cols) # Should be two less than the file
dim(gex_1_matched)
head(cols)
tail(cols)
# Create data.table with sample IDs:
# DT_corrected <- data.table(c(1:nrow(gex_1_matched)), c('dummy')) # Dummy table
DT_corrected <- gex_1_matched[, c('pt_id', 'kit_id'), with = F]
DT_corrected

for(i in cols) {
  corrected <- get_lm(gex_1_matched, pc_gex_1_matched_to_adjust, i)
  DT_corrected[, paste0(i) := corrected]
}
DT_corrected[1:5, 1:5, with = F]
dim(DT_corrected)
DT_corrected_1 <- DT_corrected
DT_corrected_1[1:5, 1:5, with = F]

# TO DO: Should use lapply, eg:
# DT_corrected <- gex_1_matched[, paste0(i) := lapply(.SD, get_lm, args, args),]
# http://stackoverflow.com/questions/25443658/r-add-new-columns-to-a-data-table-containing-many-variables

# # Delete dummy columns:
# DT_corrected_1[, c('V1', 'V2') := NULL]
# DT_corrected_1[1:5, 1:5, with = F]
# DT_corrected_1[, ncol(DT_corrected_1), with = F]
# dim(DT_corrected_1)
####################

####################
# Repeat for gex_2 file:
# Get all probe names minus sample ID column:
colnames(gex_2_matched)[c(ncol(gex_2_matched))]
dim(gex_2_matched)
gex_2_matched[, c('pt_id', 'kit_id'), with = F]

cols <- colnames(gex_2_matched)[1:c(ncol(gex_2_matched)-2)]
length(cols) # Should be one less than the file
dim(gex_2_matched)
head(cols)
tail(cols)

# Create data.table with sample IDs:
DT_corrected <- gex_2_matched[, c('pt_id', 'kit_id'), with = F]
DT_corrected

for(i in cols) {
  corrected <- get_lm(gex_2_matched, pc_gex_2_matched_to_adjust, i)
  DT_corrected[, paste0(i) := corrected]
}
DT_corrected[1:5, 1:5, with = F]
dim(DT_corrected)
DT_corrected_2 <- DT_corrected
DT_corrected_2[1:5, 1:5, with = F]

# # Delete dummy columns:
# DT_corrected_2[, c('V1', 'V2') := NULL]
# DT_corrected_2[1:5, 1:5, with = F]
# dim(DT_corrected_2)

# TO DO: Test assumptions (summary?):
# gvmodel <- gvlma(biglm_gex_1_matched)
# summary(gvmodel)
####################

####################
# Get ratio of file2 over file1
# Check columns and rownames match:
DT_corrected_1[1:5, 1:5, with = F]
DT_corrected_2[1:5, 1:5, with = F]
identical(colnames(DT_corrected_1), colnames(DT_corrected_2))
identical(DT_corrected_1[['pt_id']], DT_corrected_2[['pt_id']])

# Keep IDs, pt_id plus baseline kit_id:
sample_IDs <- DT_corrected_1[, c('pt_id', 'kit_id'), with = F]
sample_IDs

# Get ratios:
gex_FC <- DT_corrected_2[, -c('pt_id', 'kit_id'), with = F] / DT_corrected_1[, -c('pt_id', 'kit_id'), with = F]
gex_FC <- cbind(sample_IDs, gex_FC)
dim(gex_FC)
class(gex_FC)
gex_FC[1:5, 1:5, with = F]

# TO DO: check these ratios are OK:
# range(as.list(gex_FC[, -c('pt_id', 'kit_id'), .(all_means = lapply(.SD, mean)), with = F]))
####################


####################
# Prepare for MatrixeQTL, samples as columns, values as rows.
# gex_FC <- fread('gex_FC_GEx_treated_4000_and_2000.tsv_over_GEx_baseline_4000_and_2000.tsv', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
gex_FC[1:5, 1:5, with = F]
str(gex_FC)
gex_FC_t <- transpose_file(gex_FC, 2) # data.table transpose() turns factors into character
str(gex_FC_t)
class(gex_FC_t)
gex_FC_t <- as.data.frame(gex_FC_t, stringsAsFactors = FALSE)
gex_FC_t <- gex_FC_t[, moveme(names(gex_FC_t), 'rownames first')]
# Delete first row with pt_id:
gex_FC_t <- gex_FC_t[-1, ]
head(gex_FC_t)[1:5]
str(gex_FC_t)

# Convert all columns to numeric:
for(i in c(2:ncol(gex_FC_t))) {
  gex_FC_t[, i] <- as.numeric(as.character(gex_FC_t[, i]))
  }
str(gex_FC_t)
head(gex_FC_t)[1:5]

range(colMeans(gex_FC_t[, -1]))
####################


####################
# Write to disk:
gex_file1
gex_file2

write.table(gex_FC_t, sprintf('gex_FC_%s_over_%s', gex_file2, gex_file1), sep='\t', quote = FALSE, col.names = NA)
####################

####################
# The end:
# Remove objects that are not necessary to save:
#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx.
####################