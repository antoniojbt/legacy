#############################
# To be run after 02 normalisation of array data
# Antonio J Berlanga-Taylor
# 7 Nov 2016
# BEST-D project differential expression - VD deficient comparisons, all samples
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/')
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/tables_and_plots_for_draft/final_draft_BEST-D/tables/all_GEx_tables/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_VD_deficient_all",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is very project specific. Check ways of making count comparisons.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
load('R_session_saved_image_pheno_file_check.RData', verbose=T)
# load('R_session_saved_image_diff_expression.RData', verbose=T) # Has subsetted objects for array data.
# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_VD_deficient_all', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')
#install.packages("statmod")


library(limma)
library(ggplot2)
library(ellipse)
library(Hmisc)
library(splines)
library(plyr)
library(statmod)
library(illuminaHumanv4.db)

# Get additional functions needed:
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/microarray_analysis/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/moveme.R')
#############################

#############################
# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

#Check dimensions between annotation file with meta-data (must have the same number of rows, otherwise
#errors downstream):
#TO DO: Raise error if not the same.

dim(membership_file_cleaned)
str(membership_file_cleaned)
head(membership_file_cleaned)
tail(membership_file_cleaned)
head(normalised_filtered)
dim(normalised_filtered)
str(normalised_filtered)
dim(normalised_filtered_annotated)

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))


# Load full phenotype data for covariates adjustment:
phenotype_data <- read.csv('BEST-D_phenotype_file_final.tsv', sep = '\t', 
                           header = TRUE, na.string = c('-99999', "", " ", "NA"))
dim(phenotype_data)
length(which(complete.cases(phenotype_data)))
#View(phenotype_data)
head(phenotype_data)
tail(phenotype_data)
summary(phenotype_data)
str(phenotype_data)
class(phenotype_data)
names(phenotype_data)
summary(phenotype_data[, c('kit_id_randomisation', 'vitd0', 'vitd12')])
#############################

#############################
# Subset by 50 nmol/L for all samples for vitD0 values:
phenotype_data$vitd0_less_50 <- ifelse(phenotype_data$vitd0 < 50, 1, 0)
summary(as.factor(phenotype_data$vitd0_less_50))
# View(phenotype_data[, c('vitd0', 'arm')])

# Get subsets 50 nmol:
summary(phenotype_data$vitd0_less_50 == 1)
VD_def_50_all <- phenotype_data[which(phenotype_data$vitd0_less_50 == 1), ]
dim(VD_def_50_all)
head(VD_def_50_all)
head(VD_def_50_all$kit_id_randomisation)

# Subset kit_id array file 50 nmol:
dim(membership_file_cleaned)
head(membership_file_cleaned)
membership_file_cleaned_vd_50_def_all <- membership_file_cleaned
membership_file_cleaned_vd_50_def_all$kit_id <- row.names(membership_file_cleaned_vd_50_def_all)
head(membership_file_cleaned_vd_50_def_all)

membership_file_cleaned_vd_50_def_all <- merge(membership_file_cleaned_vd_50_def_all, VD_def_50_all)
dim(membership_file_cleaned_vd_50_def_all)
head(membership_file_cleaned_vd_50_def_all)
row.names(membership_file_cleaned_vd_50_def_all) <- membership_file_cleaned_vd_50_def_all$kit_id
membership_file_cleaned_vd_50_def_all <- membership_file_cleaned_vd_50_def_all[order(membership_file_cleaned_vd_50_def_all$kit_id), ]
names(membership_file_cleaned_vd_50_def_all)
head(membership_file_cleaned_vd_50_def_all)
head(membership_file_cleaned_vd_50_def_all[, c('vitd0', 'vitd12', 'pt_id', 'visit_type', 'vitd0_less_50', 'vitd0_less_25')])
dim(membership_file_cleaned)
dim(membership_file_cleaned_vd_50_def_all)
count(membership_file_cleaned_vd_50_def_all$group_membership)

# Subset array file:
dim(normalised_filtered)
normalised_filtered_VD_50_def_all <- normalised_filtered[, which(colnames(normalised_filtered) %in% row.names(membership_file_cleaned_vd_50_def_all))]
dim(normalised_filtered_VD_50_def_all)
head(normalised_filtered_VD_50_def_all)

# Check rows and columns match in pheno and array files:
identical(row.names(membership_file_cleaned_vd_50_def_all), colnames(normalised_filtered_VD_50_def_all))
#############################


#############################
# Compare joint 2000+4000 vs baseline and placebo with pairing for those with baseline values < 50 nmol/L.
# Pairing and treatment:
pairing_joint_vd_def_50 <- factor(membership_file_cleaned_vd_50_def_all$pt_id)
head(pairing_joint_vd_def_50)
length(pairing_joint_vd_def_50)

# Get treatment groups:
treatment_joint_vd_def_50 <- factor(membership_file_cleaned_vd_50_def_all$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
dim(membership_file_cleaned_vd_50_def_all)
length(treatment_joint_vd_def_50)
head(treatment_joint_vd_def_50)
summary(treatment_joint_vd_def_50)

#Define design:
design_all_treated_pairs_vd_def_50 <- model.matrix(~pairing_joint_vd_def_50+treatment_joint_vd_def_50)
head(design_all_treated_pairs_vd_def_50)[1:5, 1:5]
tail(design_all_treated_pairs_vd_def_50)[, -1]
dim(design_all_treated_pairs_vd_def_50)
colnames(design_all_treated_pairs_vd_def_50)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated_vd_def_50 <- lmFit(normalised_filtered_VD_50_def_all, design_all_treated_pairs_vd_def_50)
fit_all_treated_2_vd_def_50 <- eBayes(fit_all_treated_vd_def_50)
dim(fit_all_treated_2_vd_def_50)
str(fit_all_treated_2_vd_def_50)
summary(fit_all_treated_2_vd_def_50)
names(fit_all_treated_2_vd_def_50)
head(fit_all_treated_2_vd_def_50$coefficients)#[1:5, 1:5]
head(fit_all_treated_2_vd_def_50$lods)#[1:5, 1:5]
head(fit_all_treated_2_vd_def_50$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2_vd_def_50$coefficients
colnames(fit_all_treated_2_vd_def_50)

# Get results:
topTable(fit_all_treated_2_vd_def_50, adjust='BH')
topTable(fit_all_treated_2_vd_def_50, coef="treatment_joint_vd_def_50treated_2000+4000", adjust='BH')
topTable(fit_all_treated_2_vd_def_50, coef="treatment_joint_vd_def_50treated_placebo", adjust='BH')

topTable_pairing_joint_treated_vd_def_50 <- topTable(fit_all_treated_2_vd_def_50, coef="treatment_joint_vd_def_50treated_2000+4000", adjust='BH', n=Inf)
topTable_pairing_joint_placebo_vd_def_50 <- topTable(fit_all_treated_2_vd_def_50, coef="treatment_joint_vd_def_50treated_placebo", adjust='BH', n=Inf)
class(topTable_pairing_joint_placebo_vd_def_50)
head(topTable_pairing_joint_treated_vd_def_50)
dim(topTable_pairing_joint_treated_vd_def_50)
head(topTable_pairing_joint_placebo_vd_def_50)
dim(topTable_pairing_joint_placebo_vd_def_50)

# Basic counts
count(topTable_pairing_joint_treated_vd_def_50$adj.P.Val < 0.05)
count(topTable_pairing_joint_placebo_vd_def_50$adj.P.Val < 0.05)

count(topTable_pairing_joint_treated_vd_def_50$adj.P.Val < 0.10)
count(topTable_pairing_joint_placebo_vd_def_50$adj.P.Val < 0.10)

count(topTable_pairing_joint_treated_vd_def_50$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated_vd_def_50$logFC) > 1.1)
count(topTable_pairing_joint_treated_vd_def_50$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated_vd_def_50$logFC) < 0.9)

# Plot overlaps:
fit_all_treated_2_vd_def_50$coef[1:3, 160:ncol(fit_all_treated_2_vd_def_50$coef)]
fit_all_treated_2_vd_def_50_subset_coef <- fit_all_treated_2_vd_def_50[, 160:161]
results_by_pairing_all_treated_def_50 <- decideTests(fit_all_treated_2_vd_def_50_subset_coef)
head(results_by_pairing_all_treated_def_50)
str(results_by_pairing_all_treated_def_50)
heatDiagram(results_by_pairing_all_treated_def_50, fit_all_treated_2_vd_def_50_subset_coef)#, primary = 'treatmenttreated_4000')
vennDiagram(results_by_pairing_all_treated_def_50)

# Get correlation:
genas(fit_all_treated_2_vd_def_50, coef = c(160, 161), subset = 'logFC', plot = TRUE)

# Get annotations for topTable:
topTable_pairing_joint_placebo_vd_def_50 <- get_gene_symbols(topTable_pairing_joint_placebo_vd_def_50)
topTable_pairing_joint_treated_vd_def_50 <- get_gene_symbols(topTable_pairing_joint_treated_vd_def_50)
head(topTable_pairing_joint_treated_vd_def_50)
dim(topTable_pairing_joint_treated_vd_def_50)
head(topTable_pairing_joint_placebo_vd_def_50)
dim(topTable_pairing_joint_placebo_vd_def_50)
# View(topTable_pairing_joint_treated)

# Interpretation: There are very few significant differences for paired tests (before vs after) when only looking at those
# deficient at baseline (vitd0 < 50 nmol/L) when joining treatment groups (placebo and 2000 + 4000).

# Write results to disk:
# write.table(x=topTable_pairing_joint_treated_vd_def_50, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
#             file='full_topTable_pairing_all_treated_vd_def_50.txt')
# 
# write.table(x=topTable_pairing_joint_placebo_vd_def_50, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
#             file='full_topTable_pairing_all_placebo_vd_def_50.txt')
# Write to disk:
write_topTable("treatment_joint_vd_def_50treated_2000+4000", fit_all_treated_2_vd_def_50)
write_topTable("treatment_joint_vd_def_50treated_placebo", fit_all_treated_2_vd_def_50)
#############################


#############################
# Get subsets for paired analyses of all samples:
# Subset by 25 nmol/L for all samples for vitD0 values:
phenotype_data$vitd0_less_25 <- ifelse(phenotype_data$vitd0 < 25, 1, 0)
summary(as.factor(phenotype_data$vitd0_less_25))

# Get subsets 25 nmol:
summary(phenotype_data$vitd0_less_25 == 1)
VD_def_25_all <- phenotype_data[which(phenotype_data$vitd0_less_25 == 1), ]
dim(VD_def_25_all)
head(VD_def_25_all)
head(VD_def_25_all$kit_id_randomisation)

# Subset kit_id array file 25 nmol:
dim(membership_file_cleaned)
head(membership_file_cleaned)
membership_file_cleaned_vd_25_def_all <- membership_file_cleaned
membership_file_cleaned_vd_25_def_all$kit_id <- row.names(membership_file_cleaned_vd_25_def_all)
dim(membership_file_cleaned_vd_25_def_all)
head(membership_file_cleaned_vd_25_def_all)

membership_file_cleaned_vd_25_def_all <- merge(membership_file_cleaned_vd_25_def_all, VD_def_25_all)
dim(membership_file_cleaned_vd_25_def_all)
head(membership_file_cleaned_vd_25_def_all)
row.names(membership_file_cleaned_vd_25_def_all) <- membership_file_cleaned_vd_25_def_all$kit_id
membership_file_cleaned_vd_25_def_all <- membership_file_cleaned_vd_25_def_all[order(membership_file_cleaned_vd_25_def_all$kit_id), ]
names(membership_file_cleaned_vd_25_def_all)
head(membership_file_cleaned_vd_25_def_all)
head(membership_file_cleaned_vd_25_def_all[, c('vitd0', 'vitd12', 'pt_id', 'visit_type', 'vitd0_less_50', 'vitd0_less_25')])
dim(membership_file_cleaned)
dim(membership_file_cleaned_vd_25_def_all)

# Subset array file:
dim(normalised_filtered)
normalised_filtered_VD_25_def_all <- normalised_filtered[, which(colnames(normalised_filtered) %in% row.names(membership_file_cleaned_vd_25_def_all))]
dim(normalised_filtered_VD_25_def_all)
head(normalised_filtered_VD_25_def_all)

# Check rows and columns match in pheno and array files:
identical(row.names(membership_file_cleaned_vd_25_def_all), colnames(normalised_filtered_VD_25_def_all))
#############################


#############################
# Compare joint 2000+4000 vs baseline and placebo with pairing for those with baseline values < 25 nmol/L.
# Pairing and treatment:
pairing_joint_vd_def_25 <- factor(membership_file_cleaned_vd_25_def_all$pt_id)
head(pairing_joint_vd_def_25)
length(pairing_joint_vd_def_25)

# Get treatment groups:
treatment_joint_vd_def_25 <- factor(membership_file_cleaned_vd_25_def_all$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
dim(membership_file_cleaned_vd_25_def_all)
length(treatment_joint_vd_def_25)
head(treatment_joint_vd_def_25)
summary(treatment_joint_vd_def_25)

#Define design:
design_all_treated_pairs_vd_def_25 <- model.matrix(~pairing_joint_vd_def_25+treatment_joint_vd_def_25)
head(design_all_treated_pairs_vd_def_25)[1:5, 1:5]
tail(design_all_treated_pairs_vd_def_25)[, -1]
dim(design_all_treated_pairs_vd_def_25)
colnames(design_all_treated_pairs_vd_def_25)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated_vd_def_25 <- lmFit(normalised_filtered_VD_25_def_all, design_all_treated_pairs_vd_def_25)
fit_all_treated_2_vd_def_25 <- eBayes(fit_all_treated_vd_def_25)
dim(fit_all_treated_2_vd_def_25)
str(fit_all_treated_2_vd_def_25)
summary(fit_all_treated_2_vd_def_25)
names(fit_all_treated_2_vd_def_25)
head(fit_all_treated_2_vd_def_25$coefficients)#[1:5, 1:5]
head(fit_all_treated_2_vd_def_25$lods)#[1:5, 1:5]
head(fit_all_treated_2_vd_def_25$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2_vd_def_25$coefficients
colnames(fit_all_treated_2_vd_def_25)

# Get results:
topTable(fit_all_treated_2_vd_def_25, adjust='BH')
topTable(fit_all_treated_2_vd_def_25, coef="treatment_joint_vd_def_25treated_2000+4000", adjust='BH')
topTable(fit_all_treated_2_vd_def_25, coef="treatment_joint_vd_def_25treated_placebo", adjust='BH')

topTable_pairing_joint_treated_vd_def_25 <- topTable(fit_all_treated_2_vd_def_25, coef="treatment_joint_vd_def_25treated_2000+4000", adjust='BH', n=Inf)
topTable_pairing_joint_placebo_vd_def_25 <- topTable(fit_all_treated_2_vd_def_25, coef="treatment_joint_vd_def_25treated_placebo", adjust='BH', n=Inf)
class(topTable_pairing_joint_placebo_vd_def_25)
head(topTable_pairing_joint_treated_vd_def_25)
dim(topTable_pairing_joint_treated_vd_def_25)
head(topTable_pairing_joint_placebo_vd_def_25)
dim(topTable_pairing_joint_placebo_vd_def_25)

# Basic counts
count(topTable_pairing_joint_treated_vd_def_25$adj.P.Val < 0.05)
count(topTable_pairing_joint_placebo_vd_def_25$adj.P.Val < 0.05)

count(topTable_pairing_joint_treated_vd_def_25$adj.P.Val < 0.10)
count(topTable_pairing_joint_placebo_vd_def_25$adj.P.Val < 0.10)

count(topTable_pairing_joint_treated_vd_def_25$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated_vd_def_25$logFC) > 1.1)
count(topTable_pairing_joint_treated_vd_def_25$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated_vd_def_25$logFC) < 0.9)

# Plot overlaps:
fit_all_treated_2_vd_def_25$coef[1:3, 14:ncol(fit_all_treated_2_vd_def_25$coef)]
fit_all_treated_2_vd_def_25_subset_coef <- fit_all_treated_2_vd_def_25[, 14:15]
results_by_pairing_all_treated_def_25 <- decideTests(fit_all_treated_2_vd_def_25_subset_coef)
head(results_by_pairing_all_treated_def_25)
str(results_by_pairing_all_treated_def_25)
heatDiagram(results_by_pairing_all_treated_def_25, fit_all_treated_2_vd_def_25_subset_coef)#, primary = 'treatmenttreated_4000')
vennDiagram(results_by_pairing_all_treated_def_25)

# Get correlation:
genas(fit_all_treated_2_vd_def_25, coef = c(14, 15), subset = 'logFC', plot = TRUE)

# Get annotations for topTable:
topTable_pairing_joint_placebo_vd_def_25 <- get_gene_symbols(topTable_pairing_joint_placebo_vd_def_25)
topTable_pairing_joint_treated_vd_def_25 <- get_gene_symbols(topTable_pairing_joint_treated_vd_def_25)
head(topTable_pairing_joint_treated_vd_def_25)
dim(topTable_pairing_joint_treated_vd_def_25)
head(topTable_pairing_joint_placebo_vd_def_25)
dim(topTable_pairing_joint_placebo_vd_def_25)
# View(topTable_pairing_joint_treated)

# Interpretation: There are very few significant differences for paired tests (before vs after) when only looking at those
# deficient at baseline (vitd0 < 25 nmol/L) when joining treatment groups (placebo and 2000 + 4000).

# Write results to disk:
# write.table(x=topTable_pairing_joint_treated_vd_def_25, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
#             file='full_topTable_pairing_all_treated_vd_def_25.txt')
# 
# write.table(x=topTable_pairing_joint_placebo_vd_def_25, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
#             file='full_topTable_pairing_all_placebo_vd_def_25.txt')
# Write to disk:
write_topTable("treatment_joint_vd_def_25treated_2000+4000", fit_all_treated_2_vd_def_25)
write_topTable("treatment_joint_vd_def_25treated_placebo", fit_all_treated_2_vd_def_25)
#############################

#############################
# Compare those which increased by the median or more (vitd12 minus vitd0)
# Subset all samples on vitD delta values:
phenotype_data$delta <- phenotype_data$vitd12 - phenotype_data$vitd0
summary(phenotype_data$delta)
dim(phenotype_data)

phenotype_data$delta_bin <- ifelse(phenotype_data$delta >= 44.79, 1, 0)
summary(as.factor(phenotype_data$delta_bin))

# Get subsets with high delta:
summary(phenotype_data$delta_bin == 1)
VD_delta <- phenotype_data[which(phenotype_data$delta_bin == 1), ]
dim(VD_delta)
head(VD_delta)
head(VD_delta$kit_id_randomisation)

# Subset kit_id array file:
dim(membership_file_cleaned)
head(membership_file_cleaned)
membership_file_cleaned_vd_delta <- membership_file_cleaned
membership_file_cleaned_vd_delta$kit_id <- row.names(membership_file_cleaned_vd_delta)
head(membership_file_cleaned_vd_delta)

membership_file_cleaned_vd_delta <- merge(membership_file_cleaned_vd_delta, VD_delta)
dim(membership_file_cleaned_vd_delta)
head(membership_file_cleaned_vd_delta)
row.names(membership_file_cleaned_vd_delta) <- membership_file_cleaned_vd_delta$kit_id
membership_file_cleaned_vd_delta <- membership_file_cleaned_vd_delta[order(membership_file_cleaned_vd_delta$kit_id), ]
names(membership_file_cleaned_vd_delta)
head(membership_file_cleaned_vd_delta)
head(membership_file_cleaned_vd_delta[, c('vitd0', 'vitd12', 'pt_id', 'visit_type', 'vitd0_less_50', 'vitd0_less_25')])
dim(membership_file_cleaned)
dim(membership_file_cleaned_vd_delta)

# Subset array file:
dim(normalised_filtered)
normalised_filtered_VD_delta <- normalised_filtered[, which(colnames(normalised_filtered) %in% row.names(membership_file_cleaned_vd_delta))]
dim(normalised_filtered_VD_delta)
head(normalised_filtered_VD_delta)

# Check rows and columns match in pheno and array files:
identical(row.names(membership_file_cleaned_vd_delta), colnames(normalised_filtered_VD_delta))
#############################

#############################
# Compare joint 2000+4000 vs baseline and placebo with pairing for those with VD delta values > median.
# Pairing and treatment:
pairing_joint_vd_delta <- factor(membership_file_cleaned_vd_delta$pt_id)
head(pairing_joint_vd_delta)
length(pairing_joint_vd_delta)

# Get treatment groups, exclude placebo as these have no delta with > median:
count(membership_file_cleaned_vd_delta$two_group_Tx)
count(membership_file_cleaned_vd_delta$delta_bin)
treatment_joint_vd_delta <- factor(membership_file_cleaned_vd_delta$two_group_Tx, levels = c('untreated', 'treated_2000+4000'))
dim(membership_file_cleaned_vd_delta)
length(treatment_joint_vd_delta)
head(treatment_joint_vd_delta)
summary(treatment_joint_vd_delta)

#Define design:
design_all_treated_pairs_vd_delta <- model.matrix(~pairing_joint_vd_delta+treatment_joint_vd_delta)
head(design_all_treated_pairs_vd_delta)[1:5, 1:5]
tail(design_all_treated_pairs_vd_delta)[, -1]
dim(design_all_treated_pairs_vd_delta)
colnames(design_all_treated_pairs_vd_delta)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated_vd_delta <- lmFit(normalised_filtered_VD_delta, design_all_treated_pairs_vd_delta)
fit_all_treated_2_vd_delta <- eBayes(fit_all_treated_vd_delta)
dim(fit_all_treated_2_vd_delta)
str(fit_all_treated_2_vd_delta)
summary(fit_all_treated_2_vd_delta)
names(fit_all_treated_2_vd_delta)
head(fit_all_treated_2_vd_delta$coefficients)#[1:5, 1:5]
head(fit_all_treated_2_vd_delta$lods)#[1:5, 1:5]
head(fit_all_treated_2_vd_delta$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2_vd_delta$coefficients
colnames(fit_all_treated_2_vd_delta)

# Get results:
topTable(fit_all_treated_2_vd_delta, adjust='BH')
topTable(fit_all_treated_2_vd_delta, coef="treatment_joint_vd_deltatreated_2000+4000", adjust='BH')

topTable_pairing_joint_treated_vd_delta <- topTable(fit_all_treated_2_vd_delta, coef="treatment_joint_vd_deltatreated_2000+4000", adjust='BH', n=Inf)
class(topTable_pairing_joint_treated_vd_delta)
head(topTable_pairing_joint_treated_vd_delta)
dim(topTable_pairing_joint_treated_vd_delta)

# Basic counts
count(topTable_pairing_joint_treated_vd_delta$adj.P.Val < 0.05)
count(topTable_pairing_joint_treated_vd_delta$adj.P.Val < 0.10)
count(topTable_pairing_joint_treated_vd_delta$adj.P.Val < 0.05 & 2^(topTable_pairing_joint_treated_vd_delta$logFC) > 1.1)

# Plot overlaps:
fit_all_treated_2_vd_delta$coef[1:3, 140:ncol(fit_all_treated_2_vd_delta$coef)]
fit_all_treated_2_vd_delta_subset_coef <- fit_all_treated_2_vd_delta[, 144]
results_by_pairing_all_treated_delta <- decideTests(fit_all_treated_2_vd_delta_subset_coef)
head(results_by_pairing_all_treated_delta)
str(results_by_pairing_all_treated_delta)
heatDiagram(results_by_pairing_all_treated_delta, fit_all_treated_2_vd_delta_subset_coef)#, primary = 'treatmenttreated_4000')
vennDiagram(results_by_pairing_all_treated_delta)

# Get annotations for topTable:
topTable_pairing_joint_treated_vd_delta <- get_gene_symbols(topTable_pairing_joint_treated_vd_delta)
head(topTable_pairing_joint_treated_vd_delta)
dim(topTable_pairing_joint_treated_vd_delta)
# View(topTable_pairing_joint_treated_vd_delta)

# Interpretation: There are very few significant differences for paired tests (before vs after) when only looking at those
# whose vitamin D levels increased by more than the median (44 nmol/L) when joining treatment groups (2000 + 4000).
# There weren't any samples from placebo which had this, the median VD change was 2.58:
summary(phenotype_data[which(phenotype_data$arm == 2), c('delta')])

# Write results to disk:
write_topTable("treatment_joint_vd_deltatreated_2000+4000", fit_all_treated_2_vd_delta)
#############################

#############################
# Compare individual arms vs baseline and placebo with contrasts for those with VD delta values > median.
# Pairing and treatment:
# pairing_arms_vd_delta <- factor(membership_file_cleaned_vd_delta$pt_id)
# head(pairing_arms_vd_delta)
# length(pairing_arms_vd_delta)

# Get treatment groups, exclude placebo as these have no delta with > median:
count(membership_file_cleaned_vd_delta$group_membership)
count(membership_file_cleaned_vd_delta$delta_bin)
treatment_arms_vd_delta <- factor(membership_file_cleaned_vd_delta$group_membership, 
                                  levels = c('baseline_2000', 'baseline_4000', 'final_2000', 'final_4000'))
dim(membership_file_cleaned_vd_delta)
length(treatment_arms_vd_delta)
head(treatment_arms_vd_delta)
summary(treatment_arms_vd_delta)

#Define design:
design_arms_pairs_vd_delta <- model.matrix(~0+treatment_arms_vd_delta)
head(design_arms_pairs_vd_delta)
tail(design_arms_pairs_vd_delta)
dim(design_arms_pairs_vd_delta)
colnames(design_arms_pairs_vd_delta)
colSums(design_arms_pairs_vd_delta)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_arms_vd_delta <- lmFit(normalised_filtered_VD_delta, design_arms_pairs_vd_delta)

# Set contrasts:
cont.matrix <- makeContrasts(delta_f4000vsb4000=treatment_arms_vd_deltafinal_4000-treatment_arms_vd_deltabaseline_4000, 
                             delta_f2000vsb2000=treatment_arms_vd_deltafinal_2000-treatment_arms_vd_deltabaseline_2000,
                             levels=design_arms_pairs_vd_delta)
cont.matrix

# Obtain differentially expressed genes based on contrasted factors:
fit_arms_2_vd_delta <- contrasts.fit(fit_arms_vd_delta, cont.matrix)
fit_arms_2_vd_delta <- eBayes(fit_arms_2_vd_delta)
names(fit_arms_2_vd_delta)
colnames(fit_arms_2_vd_delta)

# Get results:
topTable(fit_arms_2_vd_delta, adjust='BH')
topTable(fit_arms_2_vd_delta, coef="delta_f4000vsb4000", adjust='BH')
topTable(fit_arms_2_vd_delta, coef="delta_f2000vsb2000", adjust='BH')

# Interpretation: There are no significant differences for tests (before vs after) when only looking at those
# whose vitamin D levels increased by more than the median (44 nmol/L) within each treatment arm.

# Write results to disk:
# write_topTable("treatment_arms_vd_deltatreated_2000+4000", fit_arms_2_vd_delta)
#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')
sessionInfo()

q()
# Next: run script for higher level analyses.
#############################