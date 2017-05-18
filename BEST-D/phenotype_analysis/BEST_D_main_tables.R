#############################################
# BEST-D basline and cytokine tables

# Author: Antonio J Berlanga-Taylor
# Date: 03 March 2017

# Purpose
# 

#############################################


#############################################

#Methods
#=======


#Usage
#=====


#To use type::
#    xxx.R [options] [arguments]
#    xxx.R --help


#Options
#=======

#-I    input file name.
#-S    output file name.
#-L    log file name.


# Input
# Output

# Use docopt, see
#https://github.com/AntonioJBT/various.dir/blob/master/Notes-common-cmds/docopt_argument_parser.txt
#https://github.com/docopt/docopt.R

#############################################


#############################################
# Logging
# TO DO: move to a separate script

##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/tables_and_plots_for_draft/main_tables/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_main_tables",Sys.Date(),".txt", sep=""))
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
# load('../data.dir/R_session_saved_image_cytokines_lm.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_main_tables','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################


#############################
# Import libraries:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
library(xtable)
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)
#############################


#############################################
# Set-up arguments:
cyto_file <- as.character(args[1])
cyto_file <- '../../data.dir/bestd_cytokines.csv'

pheno_file <- as.character(args[4])
pheno_file <- '../../data.dir/BEST-D_phenotype_file_final.tsv'

# TO DO: Read in subset of cytokine transcripts created with 'bestd_cytokines_and_trancripts.R'
# expr_file <- as.character(args[3])
# expr_file <- '../../data.dir/normalised_filtered_annotated.tab'
#############################################


#############################################
# Read files:
cyto_file <- fread(cyto_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
cyto_file
head(cyto_file)
dim(cyto_file)
tail(cyto_file)
summary(cyto_file)

pheno_file <- fread(pheno_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE,
                    na.strings = '-99999')
pheno_file
head(pheno_file)
dim(pheno_file)
tail(pheno_file)
summary(pheno_file)

# TO DO: Add expr file:
#############################################

#############################################
# Join files, these are data.table

# Join pheno and cyto:
pheno_file[, 'pt_id', with = F]
cyto_file[, 'pt_id', with = F]
setkey(pheno_file, 'pt_id')
setkey(cyto_file, 'pt_id')

all_data <- merge(pheno_file, cyto_file, all=TRUE)
all_data
dim(cyto_file)
dim(pheno_file)
dim(all_data)
head(all_data)[, 1:10]
head(all_data)[, (ncol(all_data) - 10):ncol(all_data)]
summary(all_data[, (ncol(all_data) - 10):ncol(all_data)])

# Sanity check:
str(all_data)
all_data$vitd12 <- as.numeric(all_data$vitd12)
summary(all_data$vitd12)

# Correct variable types:
str(all_data$pt_id)
all_data$pt_id <- as.character(all_data$pt_id)
plyr::count(all_data$pt_id)
# 12:41 are ptr_si0 to phosphate12:
str(all_data[, 12:41])
colnames(all_data[, 12:41])
vars_convert <- list('vitd12',
                     'BMI_DEXA12',
                     'grip_strength12',
                     'egfr0',
                     'joint_severity12',
                     'muscle_severity12',
                     'physical_activity12',
                     'corrected_calcium12',
                     'corrected_calcium0',
                     'creatinine0',
                     'creatinine12',
                     "ptr_si0",
                     "ptr_ri0",
                     "ptr_si12",
                     "ptr_ri12",
                     "art_pwv0",
                     "art_AI_aortic0",
                     "art_sbp0",
                     "art_dbp0",
                     "art_hr0",
                     "art_pwv12",
                     "art_AI_aortic12",
                     "art_sbp12",
                     "art_dbp12",
                     "art_hr12",
                     "albumin0",
                     "alk_phosphatase0",
                     "phosphate0",  
                     "vitd0",      
                     "apo_a10",
                     "apo_b0",
                     "tchol0",
                     "ldlc0",
                     "ipth0",
                     "trig0",
                     "hdlc0",
                     "crp0",
                     "vitd6",
                     "albumin12",
                     "alk_phosphatase12",
                     "phosphate12")
vars_convert
all_data <- as.data.frame(all_data)
class(all_data)

# Assign correct types to variables:
for (i in vars_convert){
  print(i)
  # all_data[, i, with = F] <- as.numeric(as.character(all_data[, i, with = F]))
  all_data[, i] <- as.numeric(as.character(all_data[, i]))
}

vars_categorical <- list("male",
                         "incident_fracture",
                         "incident_resp_infection",
                         "diabetes",
                         "heart_disease",
                         "copd_asthma",
                         "basemed_vitamind",
                         "currsmoker",
                         "season_randomisation_2",
                         "arm")
for (i in vars_categorical){
  print(i)
  # all_data[, i, with = F] <- as.numeric(as.character(all_data[, i, with = F]))
  all_data[, i] <- as.factor(as.character(all_data[, i]))
}
str(all_data[, as.character(vars_categorical)])
str(all_data[, as.character(vars_convert)])
dim(all_data)
#############################################

# TO DO:
# Define variables of interest
# Present main summary for cytokines and transcripts (main):
# Mean and SE for cytokine and transcript per group with VD0, 12 and delta values

# Present main summary for baseline characteristics (supplementary):
# By arm: age, gender, VD0, VD12, ?
# Refer to main paper.

#############################################
# Get variables of interest:
colnames(all_data)
colnames(all_data[, c(2, 29, 42, 51, 76, 88:97)])
vars_interest <- c(2, 29, 42, 51, 76, 88:97)
main_table_cyto <- all_data[, vars_interest]
head(main_table_cyto)

# Rename groups for better plotting:
main_table_cyto$arm <- factor(main_table_cyto$arm, levels = c(2, 1, 0),
                              labels = c('Placebo', '2000 IU', '4000 IU'))
plyr::count(main_table_cyto$arm)
main_table_cyto$male <- factor(main_table_cyto$male, levels = c(0, 1),
                               labels = c('Female', 'Male'))
plyr::count(main_table_cyto$male)

# TO DO: Continue here:
# Generate summaries:
options(digits = 2)
# browseVignettes(package = "dplyr")
dplyr::summarise(group_by(main_table_cyto, arm),
                 Mean = mean(vitd0, na.rm = TRUE),
                 SD = sd(vitd0, na.rm = TRUE))

dplyr::summarise(group_by(main_table_cyto, arm),
                 Mean = mean(Ln_TNFalpha12, na.rm = TRUE),
                 SD = sd(Ln_TNFalpha12, na.rm = TRUE))

# Do baseline and 12 months for each, transpose:
main_table_cyto %>%
  group_by(arm) %>%
  select(vitd0, vitd12) %>%
  summarise(
    VD0 = mean(vitd0, na.rm = TRUE),
    VD12 = mean(vitd12, na.rm = TRUE)
    )

###########
# all_data_melt <- melt(all_data, id.vars = vars_interest)
# all_data_melt
# head(all_data_melt)[1:5, 1:5]
# dim(all_data_melt)
# colnames(all_data_melt)
# plyr::count(all_data_melt$variable)
# all_data_melt[1:10, c('value', 'variable')]
# all_data_melt$variable <- factor(all_data_melt$variable, 
#                                  levels = c('arm',
#                                             'vitd0',
#                                             'vitd12',
#                                             'male',
#                                             "calendar_age_ra",
#                                             'Ln_IFNgamma0',
#                                             'Ln_IFNgamma12',
#                                             'Ln_IL10_0',
#                                             'Ln_IL10_12',
#                                             'Ln_IL6_0',
#                                             'Ln_IL6_12',
#                                             'Ln_IL8_0',
#                                             'Ln_IL8_12',
#                                             'Ln_TNFalpha0',
#                                             'Ln_TNFalpha12'),
#                                  labels = c('Randomisation arm',
#                                             '25OHD baseline',
#                                             '25OHD 12 months',
#                                             'Sex',
#                                             'Age at randomisation',
#                                             'IFNg baseline',
#                                             'IFNg 12 months',
#                                             'IL10 baseline',
#                                             'IL10 12 months',
#                                             'IL6 baseline',
#                                             'IL6 12 months',
#                                             'IL8 baseline',
#                                             'IL8 12 months',
#                                             'TNFa baseline',
#                                             'TNFa 12 months')
# )
# plyr::count(all_data_melt$variable)
# plyr::count(all_data_melt$arm)
#############################################


#############################################
# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups
# Extract variables of interest, assign groups, run lm correcting for baseline 
# and other groups
class(all_data_melt)
colnames(all_data_melt)
class(all_data)
head(all_data)
colnames(all_data)

# Sanity check treatment groups and vitamin D levels:
str(all_data$arm)
all_data[1:20, c('arm', 'vitd12')]
summary(all_data[which(all_data$arm == '4000 IU'), 'vitd12']) # 4000IU
summary(all_data[which(all_data$arm == '2000 IU'), 'vitd12']) # 2000IU
summary(all_data[which(all_data$arm == 'Placebo'), 'vitd12']) # Placebo

#############################################


#############################################

#############################################


#############################################
# Save some text:
# cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
# Print/save pretty tables with xtable
# http://blog.revolutionanalytics.com/2010/02/making-publicationready-tables-with-xtable.html
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf

print(xtable(summary(lm(value ~ time + arm, all_data_time_melt))), type = "html", file = 'test_xtable.html')
# Arithmetic mean (SE) shown. Means and SEs are adjusted for baseline values, with missing data imputed using multiple imputation.

#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file = R_session_saved_image, compress='gzip')

sessionInfo()

q()
# Next: run the script for xxx
#############################