#############################
# To be run after 02 normalisation of array data
# Antonio J Berlanga-Taylor
# 25 Feb 2016
# BEST-D project differential expression - VD deficient comparisons, baseline only for <50 vs >50 
# and <25 vs >75 nmol/L
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_VD_deficient_baseline",".txt", sep=""), open='a')
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
R_session_saved_image <- paste('R_session_saved_image_diff_expression_VD_deficient_baseline', '.RData', sep='')

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
all_baseline <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'baseline_placebo' | 
                                                membership_file_cleaned$group_membership == 'baseline_2000' |
                                                membership_file_cleaned$group_membership == 'baseline_4000'), ]
all_final <- membership_file_cleaned[which(membership_file_cleaned$group_membership == 'final_placebo' | 
                                                membership_file_cleaned$group_membership == 'final_2000' |
                                                membership_file_cleaned$group_membership == 'final_4000'), ]

subsets_list <- list(baseline_placebo, baseline_2000, baseline_4000, 
                     final_placebo, final_2000, final_4000,
                     baseline_4000_and_2000, final_4000_and_2000,
                     all_baseline, all_final)
sapply(subsets_list, dim)

pheno_baseline_keep <- which(as.character(phenotype_data$kit_id_randomisation) %in% row.names(all_baseline))
length(pheno_baseline_keep)
dim(phenotype_data)

phenotype_data_array_baseline <- phenotype_data[pheno_baseline_keep, ]
phenotype_data_array_baseline <- phenotype_data_array_baseline[order(phenotype_data_array_baseline$kit_id_randomisation), ]
row.names(phenotype_data_array_baseline) <- phenotype_data_array_baseline$kit_id_randomisation
dim(phenotype_data_array_baseline)
head(phenotype_data_array_baseline$kit_id_randomisation)

# Subset gene expression normalised data to keep only those with baseline data:
dim(all_baseline)
length(which(colnames(normalised_filtered) %in% rownames(all_baseline)))
GEx_data_all_baseline <- normalised_filtered[, which(colnames(normalised_filtered) %in% rownames(all_baseline))]
dim(GEx_data_all_baseline)
identical(colnames(GEx_data_all_baseline), rownames(all_baseline))

summary(phenotype_data_array_baseline[, c('kit_id_randomisation', 'vitd0', 'vitd12')])

# Check match with array data from baseline samples:
dim(phenotype_data_array_baseline)
identical(row.names(phenotype_data_array_baseline), rownames(all_baseline))
identical(row.names(phenotype_data_array_baseline), colnames(GEx_data_all_baseline))
#############################


#########################
## Compare low vs high baseline VD:
# Experimental design matrix specification
# Subset by 50 nmol/L of baseline VD values:
phenotype_data_array_baseline$vitd0_less_50 <- ifelse(phenotype_data_array_baseline$vitd0 < 50, 1, 0)
summary(as.factor(phenotype_data_array_baseline$vitd0_less_50))
summary(phenotype_data$vitd0)
summary(phenotype_data_array_baseline$vitd0)

#Define design:
group_vitd0_less_50 <- factor(phenotype_data_array_baseline$vitd0_less_50)
str(group_vitd0_less_50)
count(group_vitd0_less_50)
head(group_vitd0_less_50)

design_vitd0_less_50 <- model.matrix(~group_vitd0_less_50)
head(design_vitd0_less_50)
dim(design_vitd0_less_50)
dim(GEx_data_all_baseline)

#Run linear model:
fit_def_vitd0 <- lmFit(GEx_data_all_baseline, design_vitd0_less_50)
fit_def_vitd0_2 <- eBayes(fit_def_vitd0)

#Get results and plot:
topTable(fit_def_vitd0_2, adjust="BH")

#volcanoplot(fit2)
results_vitd0_def_50 <- decideTests(fit_def_vitd0_2) 
results_vitd0_def_50
vennDiagram(results_vitd0_def_50)

# Interpretation: No significant differences vitd0 < 50 vs > 50 in GEx at baseline.
# Write to disk:
write.table(topTable(fit_def_vitd0_2, adjust="BH", number = Inf), sep='\t', quote = FALSE, 
            col.names = NA, row.names = TRUE, 
            file=paste0('topTable_', 'all_baseline_', 'VDdef50', '.txt'))
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
GEx_data_all_baseline_vitd0_less_25 <- GEx_data_all_baseline[, -which(colnames(GEx_data_all_baseline) %in% to_remove)]
dim(GEx_data_all_baseline_vitd0_less_25)

# Check match between rows and columns for phenotype and array file:
identical(colnames(GEx_data_all_baseline_vitd0_less_25), 
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
dim(GEx_data_all_baseline_vitd0_less_25)


#Run linear model:
fit_def_vitd0_25 <- lmFit(GEx_data_all_baseline_vitd0_less_25, design_vitd0_less_25)
fit_def_vitd0_25_2 <- eBayes(fit_def_vitd0_25)

#Get results and plot:
topTable(fit_def_vitd0_25_2, adjust="BH")

#volcanoplot(fit2)
results_vitd0_def_25 <- decideTests(fit_def_vitd0_25_2) 
results_vitd0_def_25
vennDiagram(results_vitd0_def_25)

# Interpretation: No significant differences vitd0 < 25 vs > 75 in GEx at baseline.
# Write to disk:
write.table(topTable(fit_def_vitd0_2, adjust="BH", number = Inf), sep='\t', quote = FALSE, 
            col.names = NA, row.names = TRUE, 
            file=paste0('topTable_', 'all_baseline_', 'VDdef25', '.txt'))
###############

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################