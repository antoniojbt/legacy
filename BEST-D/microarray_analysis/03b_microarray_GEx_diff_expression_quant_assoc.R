#############################
# To be run after 02 normalisation of array data
# Antonio J Berlanga-Taylor
# 25 Feb 2016
# BEST-D project differential expression - quantitative trait tests
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/results_1/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_3",".txt", sep=""), open='a')
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
#load('R_session_saved_image_pheno_file_check.RData', verbose=T)
load('R_session_saved_image_diff_expression.RData', verbose=T) # Has subsetted objects for array data.
# load('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/R_session_saved_image_diff_expression_3.RData', verbose=T)
# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_3', '.RData', sep='')

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
source('/Users/antoniob/Desktop/BEST_D_12_FEB.DIR/scripts/gene_expression_functions.R')
#source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
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

# Subset phenotype data so that it only contains data from those which have GEx array data:
colnames(array_baseline_4000_and_2000)
dim(array_baseline_4000_and_2000)
pheno_baseline_keep <- which(as.character(phenotype_data$kit_id_randomisation) %in% colnames(array_baseline_4000_and_2000))
length(pheno_baseline_keep)
dim(phenotype_data)
summary(phenotype_data[, c('kit_id_randomisation', 'kit_id_finalVisit')])
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
## Test for correlation between vitamin D and log-expression simple correlations:
# Analyse one group at a time first.

# Subsets of array data:
dim(array_treated_4000_and_2000) # These are 12 month samples for arms 0 and 1 (4000 and 2000)
dim(array_baseline_4000_and_2000) # These are all baseline (including baseline placebo, hence ~90 more)

head(array_treated_4000_and_2000)[1:5, 1:5]
head(array_baseline_4000_and_2000)[1:5, 1:5]
head(phenotype_data)[1:5, 1:5]
class(array_treated_4000_and_2000)
class(phenotype_data)

# TO DO: check VD-sepsis corr method:
# vitd0 correlation to gene expression levels:
dim(vitd0_t)
dim(as.matrix(array_baseline_4000_and_2000))
vitd0_corr_baseline_logExp <- rcorr(vitd0_t[2, ], as.matrix(array_baseline_4000_and_2000))
names(vitd0_corr_baseline_logExp)
dim(vitd0_corr_baseline_logExp$r)
head(vitd0_corr_baseline_logExp$r)
head(vitd0_corr_baseline_logExp$P)
# heatmap(vitd0_corr_baseline_logExp$r)
# heatmap(vitd0_corr_baseline_logExp$P)

# vitd0 correlation with albumin0 and creatinine0:
row.names(vitd0_t)
rcorr(vitd0_t[2, ], vitd0_t[3, ])
rcorr(vitd0_t[2, ], vitd0_t[4, ])

# albumin0 correlation to gene expression levels:
albumin0_corr_baseline_logExp <- rcorr(vitd0_t[3, ], as.matrix(array_baseline_4000_and_2000))
dim(albumin0_corr_baseline_logExp$r)
head(albumin0_corr_baseline_logExp$r)
head(albumin0_corr_baseline_logExp$P)
# heatmap(albumin0_corr_baseline_logExp$r)
# heatmap(albumin0_corr_baseline_logExp$P)

# Albumin0 correlation to gene expression levels:
creatinine0_corr_baseline_logExp <- rcorr(vitd0_t[3, ], as.matrix(array_baseline_4000_and_2000))
dim(creatinine0_corr_baseline_logExp$r)
head(creatinine0_corr_baseline_logExp$r)
head(creatinine0_corr_baseline_logExp$P)
# heatmap(creatinine0_corr_baseline_logExp$r)
# heatmap(creatinine0_corr_baseline_logExp$P)
######################### 


######################### 
## Test for correlation between vitamin D and log-expression adjusted for covariates:

# Adjust for cofactors:
# Pass covariates for adjustment here, consider spline package for non-linear if needed. Cut down to these variables
# looking at correlations (association script for BEST-D). Some variables still correlate (age and gender particularly).
# TO DO: I haven't assessed non-linear relationships. check spline package.
# Variables adjusted for in association analysis (plink checks co-linearity an excludes some so not final list):
# incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, 
# bmi0, calendar_age_ra, season_randomisation_2, male
# vitd0, 'vitd6', 'vitd12'

# Check rows and columns between pheno file and array file:
identical(colnames(array_baseline_4000_and_2000), as.character(phenotype_data_array_baseline$kit_id_randomisation))

incident_fracture <- factor(phenotype_data_array_baseline$incident_fracture)
levels(incident_fracture)
count(incident_fracture)

incident_resp_infection <- factor(phenotype_data_array_baseline$incident_resp_infection)
diabetes <- factor(phenotype_data_array_baseline$diabetes)
heart_disease <- factor(phenotype_data_array_baseline$heart_disease)
copd_asthma <- factor(phenotype_data_array_baseline$copd_asthma)
basemed_vitamind <- factor(phenotype_data_array_baseline$basemed_vitamind)
currsmoker <- factor(phenotype_data_array_baseline$currsmoker)
bmi0 <- phenotype_data_array_baseline$bmi0
calendar_age_ra <- phenotype_data_array_baseline$calendar_age_ra
season_randomisation_2 <- factor(phenotype_data_array_baseline$season_randomisation_2)
male <- factor(phenotype_data_array_baseline$male)
count(male)

vitd12 <- phenotype_data_array_baseline$vitd12
vitd0 <- phenotype_data_array_baseline$vitd0
summary(vitd12)
summary(vitd0)

covariates_list <- cbind(vitd12, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, 
                         currsmoker, bmi0, calendar_age_ra, season_randomisation_2, male)

summary(covariates_list)
count(covariates_list)
covariates_list_df <- as.data.frame(covariates_list)
head(covariates_list_df)
which(!complete.cases(covariates_list_df))

# Others for testing:
creatinine0 <- phenotype_data_array_baseline$creatinine0
albumin0 <- phenotype_data_array_baseline$albumin0

# TO DO: haven't adjusted for cell count. Should be at least partially accounted for in paired analysis and randomisation. Check if
# I can pass top PCs for example.

#Define design:
# R command for limma for covariates would be: 
# Testing ony baseline here, so no groupings based on treatment:
design_quant_assoc <- model.matrix(~vitd0 +
                                          incident_fracture +
                                          incident_resp_infection +
                                          diabetes +
                                          heart_disease +
                                          copd_asthma +
                                          basemed_vitamind +
                                          currsmoker +
                                          bmi0 +
                                          calendar_age_ra +
                                          season_randomisation_2 +
                                          male
                                   )

design_quant_assoc
dim(design_quant_assoc)

#Run linear model on quantitative trait:
fit_vitd0 <- lmFit(array_baseline_4000_and_2000, design_quant_assoc)
fit_vitd0_2 <- eBayes(fit_vitd0)
names(fit_vitd0_2)

#Get results and plot:
topTable(fit_vitd0_2, adjust="BH")
topTable(fit_vitd0_2, adjust="BH", coef = 'vitd0')
topTable(fit_vitd0_2, adjust="BH", coef = 'bmi0')
topTable(fit_vitd0_2, adjust="BH", coef = 'male1')
topTable(fit_vitd0_2, adjust="BH", coef = 'calendar_age_ra')
topTable(fit_vitd0_2, adjust="BH", coef = 'season_randomisation_23')
topTable(fit_vitd0_2, adjust="BH", coef = 'season_randomisation_24')
topTable(fit_vitd0_2, adjust="BH", coef = 'currsmoker1')
topTable(fit_vitd0_2, adjust="BH", coef = 'basemed_vitamind1')
topTable(fit_vitd0_2, adjust="BH", coef = 'heart_disease1')
topTable(fit_vitd0_2, adjust="BH", coef = 'diabetes1')
topTable(fit_vitd0_2, adjust="BH", coef = 'incident_resp_infection1')
topTable(fit_vitd0_2, adjust="BH", coef = 'incident_fracture1')
topTable(fit_vitd0_2, adjust="BH", coef = 'copd_asthma1')

#volcanoplot(fit2)
results_vitd0_quant_adj <- decideTests(fit_vitd0_2) 
results_vitd0_quant_adj
str(results_vitd0_quant_adj)
colnames(results_vitd0_quant_adj)
vennDiagram(results_vitd0_quant_adj[, 2:6])
vennDiagram(results_vitd0_quant_adj[, 7:11])
vennDiagram(results_vitd0_quant_adj[, 11:14])
vennDiagram(results_vitd0_quant_adj[, c(2, 7:10)])
# Interpretation: No significant genes that correlate between baseline vitamin D and baseline Gex values.

#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################
