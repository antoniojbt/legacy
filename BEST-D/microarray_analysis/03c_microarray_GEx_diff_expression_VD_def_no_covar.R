#############################
# To be run after 02 normalisation of array data
# Antonio J Berlanga-Taylor
# 25 Feb 2016
# BEST-D project differential expression - VD deficient comparisons
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_VD_deficient",".txt", sep=""), open='a')
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
R_session_saved_image <- paste('R_session_saved_image_diff_expression_VD_deficient', '.RData', sep='')

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

# TO DO:
# Subset phenotype data so that it only contains data from those which have GEx array data:
head(membership_file_cleaned)
count(membership_file_cleaned$group_membership)
baseline_placebo <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'baseline_placebo'), ]
baseline_2000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'baseline_2000'), ]
baseline_4000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'baseline_4000'), ]
final_placebo <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'final_placebo'), ]
final_2000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'final_2000'), ]
final_4000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'final_4000'), ]

baseline_4000_and_2000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'baseline_2000' |
                                                          membership_file_cleaned$group_membership == 'baseline_4000'), ]
final_4000_and_2000 <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'final_2000' |
                                                          membership_file_cleaned$group_membership == 'final_4000'), ]

subsets_list <- list(baseline_placebo, baseline_2000, baseline_4000, 
                     final_placebo, final_2000, final_4000,
                     baseline_4000_and_2000, final_4000_and_2000)
sapply(subsets_list, dim)
head(baseline_4000_and_2000)

pheno_baseline_keep <- which(as.character(phenotype_data$kit_id_randomisation) %in% row.names(baseline_4000_and_2000))
length(pheno_baseline_keep)
dim(phenotype_data)
# summary(phenotype_data[, c('kit_id_randomisation', 'kit_id_finalVisit')])
phenotype_data_array_baseline <- phenotype_data[pheno_baseline_keep, ]
phenotype_data_array_baseline <- phenotype_data_array_baseline[order(phenotype_data_array_baseline$kit_id_randomisation), ]
dim(phenotype_data_array_baseline)
head(phenotype_data_array_baseline$kit_id_randomisation)
summary(phenotype_data_array_baseline$kit_id_randomisation)

# Process phenotype file for correlations, vitd0:
phenotype_data_array_baseline$kit_id_randomisation
summary(phenotype_data_array_baseline[, c('kit_id_randomisation', 'vitd0', 'albumin0', 'creatinine0')])
vitd0 <- phenotype_data_array_baseline[, c('kit_id_randomisation', 'vitd0', 'albumin0', 'creatinine0')]
dim(vitd0)

# Order and transpose it so that it matches array data order by kit id:
row.names(vitd0) <- vitd0$kit_id_randomisation
head(vitd0)
vitd0 <- vitd0[order(vitd0$kit_id_randomisation), ]
vitd0_t <- t(vitd0[order(vitd0$kit_id_randomisation), ]) # Transposed
head(vitd0)
head(vitd0_t)

# Check match with array data from baseline samples:
dim(vitd0)
dim(array_baseline_4000_and_2000)
head(colnames(array_baseline_4000_and_2000))
length(which(as.character(phenotype_data_array_baseline$kit_id_randomisation) %in% row.names(vitd0)))
length(which(colnames(array_baseline_4000_and_2000) %in% row.names(vitd0)))
identical(row.names(vitd0), colnames(array_baseline_4000_and_2000))
#############################




#########################
## Compare low vs high baseline VD:
# Experimental design matrix specification
# Subset by 50 nmol/L of baseline VD values:
phenotype_data_array_baseline$vitd0_less_50 <- ifelse(phenotype_data_array_baseline$vitd0 < 50, 1, 0)
summary(as.factor(phenotype_data_array_baseline$vitd0_less_50))

#Define design:
group_vitd0_less_50 <- factor(phenotype_data_array_baseline$vitd0_less_50)
str(group_vitd0_less_50)
count(group_vitd0_less_50)
head(group_vitd0_less_50)

design_vitd0_less_50 <- model.matrix(~group_vitd0_less_50)
head(design_vitd0_less_50)
dim(design_vitd0_less_50)
dim(array_baseline_4000_and_2000)

#Run linear model:
fit_def_vitd0 <- lmFit(array_baseline_4000_and_2000, design_vitd0_less_50)
fit_def_vitd0_2 <- eBayes(fit_def_vitd0)

#Get results and plot:
topTable(fit_def_vitd0_2, adjust="BH")

#volcanoplot(fit2)
results_vitd0_def_50 <- decideTests(fit_def_vitd0_2) 
results_vitd0_def_50
vennDiagram(results_vitd0_def_50)

# Interpretation: No significant differences vitd0 < 50 vs > 50 in GEx at baseline.
# Write to disk:
# write_topTable("treatment_jointtreated_2000+4000", fit_all_treated_2)
###############

#########################
## Compare first lowest vs highest baseline VD:
# Experimental design matrix specification
# Subset by deficient vs high (< 25 vs > 75 nmol/L) of baseline VD values:
phenotype_data_array_baseline$vitd0_less_25 <- ifelse(phenotype_data_array_baseline$vitd0 < 25, 1, 
                                                      ifelse(phenotype_data_array_baseline$vitd0 > 75, 0,
                                                             NA))
summary(as.factor(phenotype_data_array_baseline$vitd0_less_25))

# Remove from array data those not present in vitd0 < 25:
to_remove <- which(is.na(phenotype_data_array_baseline$vitd0_less_25))
to_remove <- phenotype_data_array_baseline[to_remove, 'kit_id_randomisation']
head(to_remove)
array_baseline_4000_and_2000_vitd0_less_25 <- array_baseline_4000_and_2000[, -which(colnames(array_baseline_4000_and_2000) %in% to_remove)]
dim(array_baseline_4000_and_2000_vitd0_less_25)

# Check match between rows and columns for phenotype and array file:
identical(colnames(array_baseline_4000_and_2000_vitd0_less_25), 
          as.character(phenotype_data_array_baseline$kit_id_randomisation[
            which(!is.na(phenotype_data_array_baseline$vitd0_less_25))]))

#Define design:
group_vitd0_less_25 <- factor(phenotype_data_array_baseline$vitd0_less_25)
str(group_vitd0_less_25)
count(group_vitd0_less_25)
head(group_vitd0_less_25)

design_vitd0_less_25 <- model.matrix(~group_vitd0_less_25)
head(design_vitd0_less_25)
dim(design_vitd0_less_25)
dim(array_baseline_4000_and_2000_vitd0_less_25)


#Run linear model:
fit_def_vitd0_25 <- lmFit(array_baseline_4000_and_2000_vitd0_less_25, design_vitd0_less_25)
fit_def_vitd0_25_2 <- eBayes(fit_def_vitd0_25)

#Get results and plot:
topTable(fit_def_vitd0_25_2, adjust="BH")

#volcanoplot(fit2)
results_vitd0_def_25 <- decideTests(fit_def_vitd0_25_2) 
results_vitd0_def_25
vennDiagram(results_vitd0_def_25)

# Interpretation: No significant differences vitd0 < 25 vs > 75 in GEx at baseline.
# Write to disk:
# write_topTable("treatment_jointtreated_2000+4000", fit_all_treated_2)
###############

########################
# Compare joint 2000+4000 vs baseline and placebo with pairing for those with baseline values < 50 nmol/L.

# Get subsets:
summary(phenotype_data_array_baseline$vitd0_less_50 == 1)
VD_def_50_baseline <- phenotype_data_array_baseline[which(phenotype_data_array_baseline$vitd0_less_50 == 1), ]
dim(VD_def_50_baseline)
head(VD_def_50_baseline)
head(VD_def_50_baseline$kit_id_randomisation)

# Subset kit_id array file:
head(membership_file_cleaned)
membership_file_cleaned_vd_50_def <- membership_file_cleaned
membership_file_cleaned_vd_50_def$kit_id <- row.names(membership_file_cleaned_vd_50_def)
head(membership_file_cleaned_vd_50_def)

membership_file_cleaned_vd_50_def <- merge(membership_file_cleaned_vd_50_def, VD_def_50_baseline)
dim(membership_file_cleaned_vd_50_def)
head(membership_file_cleaned_vd_50_def)
row.names(membership_file_cleaned_vd_50_def) <- membership_file_cleaned_vd_50_def$kit_id
membership_file_cleaned_vd_50_def <- membership_file_cleaned_vd_50_def[order(membership_file_cleaned_vd_50_def$kit_id), ]
head(membership_file_cleaned_vd_50_def)
head(membership_file_cleaned_vd_50_def[, c('vitd0', 'vitd12', 'pt_id', 'visit_type', 'vitd0_less_50', 'vitd0_less_25')])
dim(membership_file_cleaned)
dim(membership_file_cleaned_vd_50_def)

# Subset array file:
dim(normalised_filtered)
normalised_filtered_VD_50_def <- normalised_filtered[, which(colnames(normalised_filtered) %in% row.names(membership_file_cleaned_vd_50_def))]
dim(normalised_filtered_VD_50_def)
head(normalised_filtered_VD_50_def)

# Check rows and columns match in pheno and array files:
identical(row.names(membership_file_cleaned_vd_50_def), colnames(normalised_filtered_VD_50_def))

# Pairing and treatment:
pairing_joint_vd_def_50 <- factor(membership_file_cleaned_vd_50_def$pt_id)
head(pairing_joint_vd_def_50)
length(pairing_joint_vd_def_50)

# Get treatment groups:
treatment_joint_vd_def_50 <- factor(membership_file_cleaned_vd_50_def$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
head(treatment_joint_vd_def_50)
summary(treatment_joint_vd_def_50)

#Define design:
design_all_treated_pairs_vd_def_50 <- model.matrix(~pairing_joint_vd_def_50+treatment_joint_vd_def_50)
head(design_all_treated_pairs_vd_def_50)[1:5, 1:5]
tail(design_all_treated_pairs_vd_def_50)[, -1]
dim(design_all_treated_pairs_vd_def_50)
colnames(design_all_treated_pairs_vd_def_50)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated_vd_def_50 <- lmFit(normalised_filtered_VD_50_def, design_all_treated_pairs_vd_def_50)
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
genas(fit_all_treated_2, coef = c(160, 161), subset = 'logFC', plot = TRUE)

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
write.table(x=topTable_pairing_joint_treated_vd_def_50, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_treated_vd_def_50.txt')

write.table(x=topTable_pairing_joint_placebo_vd_def_50, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE,
            file='full_topTable_pairing_all_placebo_vd_def_50.txt')
# Write to disk:
# write_topTable("treatment_jointtreated_2000+4000", fit_all_treated_2)
#########################

########################
# Compare joint 2000+4000 vs baseline and placebo with pairing for those with baseline values < 25 nmol/L.

# Get subsets:
summary(phenotype_data_array_baseline$vitd0_less_25 == 1)
VD_def_25_baseline <- phenotype_data_array_baseline[which(phenotype_data_array_baseline$vitd0_less_25 == 1), ]
dim(VD_def_25_baseline)
head(VD_def_25_baseline)
head(VD_def_25_baseline$kit_id_randomisation)

# Subset kit_id array file:
head(membership_file_cleaned)
membership_file_cleaned_vd_25_def <- membership_file_cleaned
membership_file_cleaned_vd_25_def$kit_id <- row.names(membership_file_cleaned_vd_25_def)
head(membership_file_cleaned_vd_25_def)

membership_file_cleaned_vd_25_def <- merge(membership_file_cleaned_vd_25_def, VD_def_25_baseline)
dim(membership_file_cleaned_vd_25_def)
head(membership_file_cleaned_vd_25_def)
row.names(membership_file_cleaned_vd_25_def) <- membership_file_cleaned_vd_25_def$kit_id
membership_file_cleaned_vd_25_def <- membership_file_cleaned_vd_25_def[order(membership_file_cleaned_vd_25_def$kit_id), ]
head(membership_file_cleaned_vd_25_def)
head(membership_file_cleaned_vd_25_def[, c('vitd0', 'vitd12', 'pt_id', 'visit_type', 'vitd0_less_50', 'vitd0_less_25')])
dim(membership_file_cleaned)
dim(membership_file_cleaned_vd_25_def)

# Subset array file:
dim(normalised_filtered)
normalised_filtered_vd_25_def <- normalised_filtered[, which(colnames(normalised_filtered) %in% row.names(membership_file_cleaned_vd_25_def))]
dim(normalised_filtered_vd_25_def)
head(normalised_filtered_vd_25_def)

# Check rows and columns match in pheno and array files:
identical(row.names(membership_file_cleaned_vd_25_def), colnames(normalised_filtered_vd_25_def))

# Pairing and treatment:
pairing_joint_vd_def_25 <- factor(membership_file_cleaned_vd_25_def$pt_id)
head(pairing_joint_vd_def_25)
length(pairing_joint_vd_def_25)

# Get treatment groups:
treatment_joint_vd_def_25 <- factor(membership_file_cleaned_vd_25_def$all_treated, levels = c('untreated', 'treated_2000+4000', 'treated_placebo'))
head(treatment_joint_vd_def_25)
summary(treatment_joint_vd_def_25)

#Define design:
design_all_treated_pairs_vd_def_25 <- model.matrix(~pairing_joint_vd_def_25+treatment_joint_vd_def_25)
head(design_all_treated_pairs_vd_def_25)[1:5, 1:5]
tail(design_all_treated_pairs_vd_def_25)[, -1]
dim(design_all_treated_pairs_vd_def_25)
colnames(design_all_treated_pairs_vd_def_25)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated_vd_def_25 <- lmFit(normalised_filtered_vd_25_def, design_all_treated_pairs_vd_def_25)
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
genas(fit_all_treated_2, coef = c(14, 15), subset = 'logFC', plot = TRUE)

# Get annotations for topTable:
topTable_pairing_joint_placebo_vd_def_25 <- get_gene_symbols(topTable_pairing_joint_placebo_vd_def_25)
topTable_pairing_joint_treated_vd_def_25 <- get_gene_symbols(topTable_pairing_joint_treated_vd_def_25)
head(topTable_pairing_joint_treated_vd_def_25)
dim(topTable_pairing_joint_treated_vd_def_25)
head(topTable_pairing_joint_placebo_vd_def_25)
dim(topTable_pairing_joint_placebo_vd_def_25)
# View(topTable_pairing_joint_treated)

# Interpretation: There are no significant differences for paired tests (before vs after) when only looking at those
# severely deficient at baseline (vitd0 < 25 nmol/L) when joining treatment groups (placebo and 2000 + 4000).
# Write to disk:
# write_topTable("treatment_jointtreated_2000+4000", fit_all_treated_2)
#########################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################