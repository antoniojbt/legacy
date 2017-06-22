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
# cyto_file <- as.character(args[1])
# cyto_file <- '../../data.dir/bestd_cytokines.csv'
# 
# pheno_file <- as.character(args[2])
# pheno_file <- '../../data.dir/BEST-D_phenotype_file_final.tsv'

full_pheno_file <- '../data.dir/BEST-D_phenotype_file_final_cytokines_and_transcripts.csv'
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


#############################################
###########
# Get variables of interest:
colnames(all_data)
colnames(all_data[, c(2, 29, 42, 51, 52, 76, 88:97)])
vars_interest <- c(2, 29, 42, 51, 52, 76, 88:97)
all_data_reduced <- all_data[, vars_interest]
head(all_data_reduced)
summary(all_data_reduced)

# Rename groups for better plotting:
all_data_reduced$arm <- factor(all_data_reduced$arm, levels = c(2, 1, 0),
                               labels = c('Placebo', '2000 IU', '4000 IU'))
plyr::count(all_data_reduced$arm)
all_data_reduced$male <- factor(all_data_reduced$male, levels = c(0, 1),
                                labels = c('Female', 'Male'))
plyr::count(all_data_reduced$male)
###########

###########
# Generate summaries
# TO DO: add transcripts:
# Main variables
# Functions for mean and SD for main table:
get_mean <- function(i) round(mean(i, na.rm = TRUE), 2)
get_sd <- function(i) round(sd(i, na.rm = TRUE), 2)

a_mean <- get_mean(all_data_reduced$vitd0)
an_sd <- get_sd(all_data_reduced$vitd0)

nice_print <- sprintf('%s (%s)', a_mean, an_sd)
nice_print

get_mean_sd <- function(i) {
  a_mean <- get_mean(i)
  an_sd <- get_sd(i)
  nice_print <- sprintf('%s (%s)', a_mean, an_sd)
  return(nice_print)
}
get_mean_sd(all_data_reduced$vitd0)

# Get all means and SDs:
get_all_means <- function(df) {
  col_name <- c(
    'Age' = get_mean_sd(df$calendar_age_ra),
    'BMI' = get_mean_sd(df$bmi0),
    '25OHD baseline' = get_mean_sd(df$vitd0),
    '25OHD 12 months' = get_mean_sd(df$vitd12),
    'IFNg baseline' = get_mean_sd(df$Ln_IFNgamma0),
    'IFNg 12 months' = get_mean_sd(df$Ln_IFNgamma12),
    'IL10 baseline' = get_mean_sd(df$Ln_IL10_0),
    'IL10 12 months' = get_mean_sd(df$Ln_IL10_12),
    'IL6 baseline' = get_mean_sd(df$Ln_IL6_0),
    'IL6 12 months' = get_mean_sd(df$Ln_IL6_12),
    'IL8 baseline' = get_mean_sd(df$Ln_IL8_0),
    'IL8 12 months' = get_mean_sd(df$Ln_IL8_12),
    'TNFa baseline' = get_mean_sd(df$Ln_TNFalpha0),
    'TNFa 12 months' = get_mean_sd(df$Ln_TNFalpha12)
  )
  as_df <- as.data.frame(col_name)
  return(as_df)
}

main_table_2 <- data.frame(get_all_means(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                           get_all_means(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                           get_all_means(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
)
colnames(main_table_2) <- c('Placebo', '2000 IU', '4000 IU')
# View(main_table_2)
###########

###########
# TO DO: this should be by complete pairs, then avegared and added to table
# Will error as passing text after calculating values above
# Add deltas:
# main_table_2$`Delta 2000` <- round(as.numeric(as.character(main_table_2$`2000 IU`)) - 
#                                        as.numeric(as.character(main_table_2$Placebo)),
#                                      2)
# main_table_2$`Delta 4000` <- round(as.numeric(as.character(main_table_2$`4000 IU`)) - 
#                                        as.numeric(as.character(main_table_2$Placebo)),
#                                      2)
# main_table_2$`Delta regimens` <- round(as.numeric(as.character(main_table_2$`4000 IU`)) - 
#                                            as.numeric(as.character(main_table_2$`2000 IU`)),
#                                          2)
# # Convert to numeric and round to 2 digits:
# df <- main_table_2
# for (i in colnames(df)) {
#   print(i)
#   df[, i] <- as.numeric(as.character(df[, i]))
#   df[, i] <- round(df[, i], 2)
# }
# main_table_2 <- df
# head(main_table_2)
# # View(main_table_2)
###########

###########
# Get total counts:
head(all_data)
dim(all_data)
counts_by_arm <- plyr::count(all_data_reduced$arm)
rownames(counts_by_arm) <- counts_by_arm$x
counts_by_arm <- t(counts_by_arm)
counts_by_arm <- as.data.frame(counts_by_arm)
counts_by_arm <- counts_by_arm[-1, ]
rownames(counts_by_arm)[1] <- 'Total n'
df <- counts_by_arm
for (i in colnames(df)) {
  print(i)
  df[, i] <- as.numeric(as.character(df[, i]))
}
counts_by_arm <- df
str(counts_by_arm)
counts_by_arm
###########

###########
# Get gender counts, drop males and gender column:
gender_table <- plyr::count(all_data_reduced, vars = c('male', 'arm'))
gender_table$proportion <- (gender_table$freq / sum(gender_table$freq)) * 100
gender_table <- gender_table[-c(4:6), -1]
rownames(gender_table) <- gender_table$arm
gender_table <- t(gender_table)[-1, ]
gender_table <- as.data.frame(gender_table)
rownames(gender_table)[1] <- 'Female n'
rownames(gender_table)[2] <- 'Female %'
df <- gender_table
for (i in colnames(df)) {
  print(i)
  df[, i] <- as.numeric(as.character(df[, i]))
}
gender_table <- df
str(gender_table)
# Gender proportions are when counting each arm per gender, this doesn't make sense when cutting out 'male' rows
# Leaving for now though.
gender_table
###########

###########
# Merge tables
head(main_table_2)
gender_table
# With gender counts:
main_table_2 <- rbind(gender_table, main_table_2)
# With total counts:
main_table_2 <- rbind(counts_by_arm, main_table_2)
# Remove Female %, not useful:
main_table_2 <- main_table_2[-3, ]
# View(main_table_2)

# Save table to file:
print(xtable(main_table_2), type = "html", file = 'BESTD_table_mean_SD.html')
############


############
# Add SEM in separate table
# Calculate standard error of the mean, with sample size to infinity sem approaches 0
# Assumes data have normal distribution
# Functions:
# Standard error of the mean:
get_sem <- function(i) round(sqrt(var(i, na.rm = TRUE) / length(na.omit(i))), 2)
# 95% CIs of the mean:
get_ci95 <- function(i) c(round(mean(i, na.rm = TRUE) - 2 * sem, 2), round(mean(i, na.rm = TRUE) + 2 * sem, 2))
get_ci95up <- function(i) round((mean(i, na.rm = TRUE) + 2 * sem), 2)
get_ci95low <- function(i) round((mean(i, na.rm = TRUE) - 2 * sem), 2)

sem <- get_sem(all_data_reduced$vitd12)
ci95 <- get_ci95(all_data_reduced$vitd12)
ci95up <- get_ci95up(all_data_reduced$vitd12)
ci95low <- get_ci95low(all_data_reduced$vitd12)
sem
ci95
ci95low
ci95up

nice_print <- sprintf('%s (%s, %s)', sem, ci95low, ci95up)
nice_print

get_sem_ci95 <- function(i) {
  sem <- get_sem(i)
  ci95up <- get_ci95up(i)
  ci95low <- get_ci95low(i)
  nice_print <- sprintf('%s (%s, %s)', sem, ci95low, ci95up)
  return(nice_print)
}
get_sem_ci95(all_data_reduced$vitd0)
############

############
get_all_sems <- function(df) {
  col_name <- c(
    '25OHD baseline' = get_sem_ci95(df$vitd0),
    '25OHD 12 months' = get_sem_ci95(df$vitd12),
    'IFNg baseline' = get_sem_ci95(df$Ln_IFNgamma0),
    'IFNg 12 months' = get_sem_ci95(df$Ln_IFNgamma12),
    'IL10 baseline' = get_sem_ci95(df$Ln_IL10_0),
    'IL10 12 months' = get_sem_ci95(df$Ln_IL10_12),
    'IL6 baseline' = get_sem_ci95(df$Ln_IL6_0),
    'IL6 12 months' = get_sem_ci95(df$Ln_IL6_12),
    'IL8 baseline' = get_sem_ci95(df$Ln_IL8_0),
    'IL8 12 months' = get_sem_ci95(df$Ln_IL8_12),
    'TNFa baseline' = get_sem_ci95(df$Ln_TNFalpha0),
    'TNFa 12 months' = get_sem_ci95(df$Ln_TNFalpha12)
  )
  as_df <- as.data.frame(col_name)
  return(as_df)
}

main_table_2_sem <- data.frame(get_all_sems(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                               get_all_sems(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                               get_all_sems(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
)
colnames(main_table_2_sem) <- c('Placebo', '2000 IU', '4000 IU')
# View(main_table_2_sem)

# Save table to file:
print(xtable(main_table_2_sem), type = "html", file = 'BESTD_table_sem_CI95.html')
############

############
# Use bootstrapping (re-sampling technique) for non-normal distributions:
# https://datascienceplus.com/introduction-to-bootstrap-with-applications-to-mixed-effect-models/
# Resample observed data with replacement and compute statistic of interest many times to get a distribution.
# e.g.:
# set.seed(20151101)
# height <- rnorm(100, 175, 6)
# t0 <- median(height)
# t <- sapply(1:1000, function(x) median(sample(x = height, size = 100, replace = TRUE)))
# hist(t)
# abline(v = t0, col = "orange", lwd = 3)

# Use library(boot) for more methods
############
#############################################



#############################################
# Save some text:
# cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
# Arithmetic mean (SE) shown. Means and SEs are adjusted for baseline values, with missing data imputed using multiple imputation.

# Print/save pretty tables with xtable
# http://blog.revolutionanalytics.com/2010/02/making-publicationready-tables-with-xtable.html
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf

# print(xtable(summary(lm(vitd12 ~ vitd0 + . , all_data_reduced))), type = "html", file = 'test_xtable.html')
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