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

pheno_file <- as.character(args[2])
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


#############################################
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


# Generate summaries:
# browseVignettes(package = "dplyr")
dplyr::summarise(group_by(all_data_reduced, arm),
                 Mean = mean(vitd0, na.rm = TRUE),
                 SD = sd(vitd0, na.rm = TRUE))

dplyr::summarise(group_by(all_data_reduced, arm),
                 Mean = mean(Ln_TNFalpha12, na.rm = TRUE),
                 SD = sd(Ln_TNFalpha12, na.rm = TRUE))

main_table <- all_data_reduced %>%
  group_by(arm) %>%
  summarise_if(is.numeric, funs(mean = mean, SD = sd), na.rm = TRUE)
main_table
# Transpose table for better viewing:
colnames(main_table)
rownames(main_table) <- main_table$arm
main_table_t <- t(main_table)[-1, ]
head(main_table_t)
class(main_table_t)
main_table_t <- as.data.frame(main_table_t)
rownames(main_table_t)
colnames(main_table_t)
# Convert to numeric:
for (i in colnames(main_table_t)) {
  print(i)
  main_table_t[, i] <- as.numeric(as.character(main_table_t[, i]))
}
# View(main_table_t)

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
gender_table

# Merge tables:
head(main_table_t)
gender_table
# With gender counts:
main_table_t <- rbind(gender_table, main_table_t)
# With total counts:
main_table_t <- rbind(counts_by_arm, main_table_t)
# View(main_table_t)

# Convert to numeric and round to 2 digits:
for (i in colnames(main_table_t)) {
  print(i)
  main_table_t[, i] <- as.numeric(as.character(main_table_t[, i]))
  main_table_t[, i] <- round(main_table_t[, i], 2)
}

str(main_table_t)
# View(main_table_t)

# Save table to file:
print(xtable(main_table_t), type = "html", file = 'BESTD_main_xtable.html')

# TO DO: add transcripts:
# Second version of the same:
main_table_2 <-
  all_data_reduced %>%
  group_by(arm) %>%
  summarise(
    'Age' = mean(calendar_age_ra, na.rm = TRUE),
    'Age (SD)' = sd(calendar_age_ra, na.rm = TRUE),
    'BMI' = mean(bmi0, na.rm = TRUE),
    'BMI (SD)' = sd(bmi0, na.rm = TRUE),
    '25OHD baseline' = mean(vitd0, na.rm = TRUE),
    '25OHD baseline (SD)' = sd(vitd0, na.rm = TRUE),
    '25OHD 12 months' = mean(vitd12, na.rm = TRUE),
    '25OHD 12 months (SD)' = sd(vitd12, na.rm = TRUE),
    'IFNg baseline' = mean(Ln_IFNgamma0, na.rm = TRUE),
    'IFNg baseline (SD)' = sd(Ln_IFNgamma0, na.rm = TRUE),
    'IFNg 12 months' = mean(Ln_IFNgamma12, na.rm = TRUE),
    'IFNg 12 months (SD)' = sd(Ln_IFNgamma12, na.rm = TRUE),
    'IL10 baseline' = mean(Ln_IL10_0, na.rm = TRUE),
    'IL10 baseline (SD)' = sd(Ln_IL10_0, na.rm = TRUE),
    'IL10 12 months' = mean(Ln_IL10_12, na.rm = TRUE),
    'IL10 12 months (SD)' = sd(Ln_IL10_12, na.rm = TRUE),
    'IL6 baseline' = mean(Ln_IL6_0, na.rm = TRUE),
    'IL6 baseline (SD)' = sd(Ln_IL6_0, na.rm = TRUE),
    'IL6 12 months' = mean(Ln_IL6_12, na.rm = TRUE),
    'IL6 12 months (SD)' = sd(Ln_IL6_12, na.rm = TRUE),
    'IL8 baseline' = mean(Ln_IL8_0, na.rm = TRUE),
    'IL8 baseline (SD)' = sd(Ln_IL8_0, na.rm = TRUE),
    'IL8 12 months' = mean(Ln_IL8_12, na.rm = TRUE),
    'IL8 12 months (SD)' = sd(Ln_IL8_12, na.rm = TRUE),
    'TNFa baseline' = mean(Ln_TNFalpha0, na.rm = TRUE),
    'TNFa baseline (SD)' = sd(Ln_TNFalpha0, na.rm = TRUE),
    'TNFa 12 months' = mean(Ln_TNFalpha12, na.rm = TRUE),
    'TNFa 12 months (SD)' = sd(Ln_TNFalpha12, na.rm = TRUE)
    )
# Transpose table for better viewing:
colnames(main_table_2)
rownames(main_table_2) <- main_table_2$arm
main_table_2_t <- t(main_table_2)[-1, ]

# Merge tables:
head(main_table_2_t)
gender_table
# With gender counts:
main_table_2_t <- rbind(gender_table, main_table_2_t)
# With total counts:
main_table_2_t <- rbind(counts_by_arm, main_table_2_t)
str(main_table_2_t)
# View(main_table_2_t)

# TO DO: this should be by complete pairs, then avegared and added to table
# Add deltas:
main_table_2_t$`Delta 2000` <- round(as.numeric(as.character(main_table_2_t$`2000 IU`)) - 
                                     as.numeric(as.character(main_table_2_t$Placebo)),
                                   2)
main_table_2_t$`Delta 4000` <- round(as.numeric(as.character(main_table_2_t$`4000 IU`)) - 
                                     as.numeric(as.character(main_table_2_t$Placebo)),
                                   2)
main_table_2_t$`Delta regimens` <- round(as.numeric(as.character(main_table_2_t$`4000 IU`)) - 
                                       as.numeric(as.character(main_table_2_t$`2000 IU`)),
                                     2)
str(main_table_2_t)
# Convert to numeric and round to 2 digits:
df <- main_table_2_t
for (i in colnames(df)) {
  print(i)
  df[, i] <- as.numeric(as.character(df[, i]))
  df[, i] <- round(df[, i], 2)
  }
main_table_2_t <- df

# View(main_table_2_t)
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

print(xtable(summary(lm(vitd12 ~ vitd0 + . , all_data_reduced))), type = "html", file = 'test_xtable.html')
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