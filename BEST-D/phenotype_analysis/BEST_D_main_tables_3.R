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
# Print/save pretty tables with xtable
# http://blog.revolutionanalytics.com/2010/02/making-publicationready-tables-with-xtable.html
# https://cran.r-project.org/web/packages/xtable/vignettes/xtableGallery.pdf

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

full_pheno_file <- '../../data.dir/BEST-D_phenotype_file_final_cytokines_and_transcripts.csv'
#############################################


#############################################
# Read files:
full_pheno <- fread(full_pheno_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
full_pheno
head(full_pheno)
dim(full_pheno)
tail(full_pheno)
summary(full_pheno)
str(full_pheno)
full_pheno <- as.data.frame(full_pheno)
all_data <- full_pheno[, -1]
dim(all_data)
#############################################

#############################################
###########
# Get variables of interest:
colnames(all_data)[c(3,4,31,44,53,54,78,88:107)]
vars_interest <- c("pt_id",
                   "arm",
                   "vitd0",
                   "vitd12",
                   "male",
                   "bmi0",
                   "calendar_age_ra",
                   "Ln_IFNgamma0",
                   "Ln_IL10_0",
                   "Ln_IL6_0",
                   "Ln_IL8_0",
                   "Ln_TNFalpha0",
                   "Ln_IFNgamma12",
                   "Ln_IL10_12",
                   "Ln_IL6_12",
                   "Ln_IL8_12",
                   "Ln_TNFalpha12",
                   "transcript_IL10_baseline",
                   "transcript_IL6_baseline",
                   "transcript_TNF_baseline",
                   "transcript_IL8_baseline",
                   "transcript_IFNG_baseline",
                   "transcript_IL10_12months",
                   "transcript_IL6_12months",
                   "transcript_TNF_12months",
                   "transcript_IL8_12months",
                   "transcript_IFNG_12months")
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
get_mean_sd(all_data_reduced$delta_vitd12)

# Get all means and SDs:
get_all_means <- function(df) {
  col_name <- c(
    'Age' = get_mean_sd(df$calendar_age_ra),
    'BMI' = get_mean_sd(df$bmi0),
    '25(OH)D baseline' = get_mean_sd(df$vitd0),
    '25(OH)D 12 months' = get_mean_sd(df$vitd12),
    'IFNg baseline' = get_mean_sd(df$Ln_IFNgamma0),
    'IFNg 12 months' = get_mean_sd(df$Ln_IFNgamma12),
    'IL10 baseline' = get_mean_sd(df$Ln_IL10_0),
    'IL10 12 months' = get_mean_sd(df$Ln_IL10_12),
    'IL6 baseline' = get_mean_sd(df$Ln_IL6_0),
    'IL6 12 months' = get_mean_sd(df$Ln_IL6_12),
    'IL8 baseline' = get_mean_sd(df$Ln_IL8_0),
    'IL8 12 months' = get_mean_sd(df$Ln_IL8_12),
    'TNFa baseline' = get_mean_sd(df$Ln_TNFalpha0),
    'TNFa 12 months' = get_mean_sd(df$Ln_TNFalpha12),
    'IFNg baseline (mRNA)' = get_mean_sd(df$transcript_IFNG_baseline),
    'IFNg 12 months (mRNA)' = get_mean_sd(df$transcript_IFNG_12months),
    'IL10 baseline (mRNA)' = get_mean_sd(df$transcript_IL10_baseline),
    'IL10 12 months (mRNA)' = get_mean_sd(df$transcript_IL10_12months),
    'IL6 baseline (mRNA)' = get_mean_sd(df$transcript_IL6_baseline),
    'IL6 12 months (mRNA)' = get_mean_sd(df$transcript_IL6_12months),
    'IL8 baseline (mRNA)' = get_mean_sd(df$transcript_IL8_baseline),
    'IL8 12 months (mRNA)' = get_mean_sd(df$transcript_IL8_12months),
    'TNFa baseline (mRNA)' = get_mean_sd(df$transcript_TNF_baseline),
    'TNFa 12 months (mRNA)' = get_mean_sd(df$transcript_TNF_12months)
    )
  as_df <- as.data.frame(col_name)
  return(as_df)
}

main_table_2 <- data.frame(get_all_means(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                           get_all_means(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                           get_all_means(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
)
colnames(main_table_2) <- c('Placebo (mean, SD)', '2000 IU (mean, SD)', '4000 IU (mean, SD)')
# colnames(main_table_2) <- c('Placebo', '2000 IU', '4000 IU')
# View(main_table_2)
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
# Change col names for binding later:
colnames(counts_by_arm) <- c('Placebo (mean, SD)', '2000 IU (mean, SD)', '4000 IU (mean, SD)')
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
# Change col names for binding later:
colnames(gender_table) <- c('Placebo (mean, SD)', '2000 IU (mean, SD)', '4000 IU (mean, SD)')
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

###########
# Add deltas in separate table
# By complete pairs, then avegared and added to table
vars_delta_0 <- c("vitd0",
                  "Ln_IFNgamma0",
                  "Ln_IL10_0",
                  "Ln_IL6_0",
                  "Ln_IL8_0",
                  "Ln_TNFalpha0",
                  "transcript_IFNG_baseline",
                  "transcript_IL10_baseline",
                  "transcript_IL6_baseline",
                  "transcript_IL8_baseline",
                  "transcript_TNF_baseline"
)
vars_delta_12 <- c("vitd12",
                   "Ln_IFNgamma12",
                   "Ln_IL10_12",
                   "Ln_IL6_12",
                   "Ln_IL8_12",
                   "Ln_TNFalpha12",
                   "transcript_IFNG_12months",
                   "transcript_IL10_12months",
                   "transcript_IL6_12months",
                   "transcript_IL8_12months",
                   "transcript_TNF_12months"
)

for (i in vars_delta_12) {
  print(i)
  index <- match(i, vars_delta_12)
  print(vars_delta_0[index])
  basal <- vars_delta_0[index]
  delta_i <- sprintf('delta_%s', i)
  print(delta_i)
  all_data_reduced[, delta_i] <- as.numeric(as.character(all_data_reduced[, i])) -
    as.numeric(as.character(all_data_reduced[, basal]))
}

colnames(all_data_reduced)
summary(all_data_reduced[, c(28:38)])
str(all_data_reduced[, c(28:38)])
# View(all_data_reduced[, c(2, 28:38)])

# Use functions from above for means and SDs:
get_mean_sd(all_data_reduced$delta_Ln_IFNgamma12)

get_all_deltas <- function(df) {
  col_name <- c(
    '25(OH)D 12 months' = get_mean_sd(df$delta_vitd12),
    'IFNg 12 months' = get_mean_sd(df$delta_Ln_IFNgamma12),
    'IL10 12 months' = get_mean_sd(df$delta_Ln_IL10_12),
    'IL6 12 months' = get_mean_sd(df$delta_Ln_IL6_12),
    'IL8 12 months' = get_mean_sd(df$delta_Ln_IL8_12),
    'TNFa 12 months' = get_mean_sd(df$delta_Ln_TNFalpha12),
    'IFNg 12 months (mRNA)' = get_mean_sd(df$delta_transcript_IFNG_12months),
    'IL10 12 months (mRNA)' = get_mean_sd(df$delta_transcript_IL10_12months),
    'IL6 12 months (mRNA)' = get_mean_sd(df$delta_transcript_IL6_12months),
    'IL8 12 months (mRNA)' = get_mean_sd(df$delta_transcript_IL8_12months),
    'TNFa 12 months (mRNA)' = get_mean_sd(df$delta_transcript_TNF_12months)
  )
  as_df <- as.data.frame(col_name)
  return(as_df)
  }

main_table_deltas <- data.frame(get_all_deltas(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                           get_all_deltas(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                           get_all_deltas(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
)
colnames(main_table_deltas) <- c('Placebo (delta)', '2000 IU (delta)', '4000 IU (delta)')
# View(main_table_deltas)

# Save table to file:
print(xtable(main_table_deltas), type = "html", file = 'BESTD_table_deltas.html')
###########

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
    'Age' = get_sem_ci95(df$calendar_age_ra),
    'BMI' = get_sem_ci95(df$bmi0),
    '25(OH)D baseline' = get_sem_ci95(df$vitd0),
    '25(OH)D 12 months' = get_sem_ci95(df$vitd12),
    'IFNg baseline' = get_sem_ci95(df$Ln_IFNgamma0),
    'IFNg 12 months' = get_sem_ci95(df$Ln_IFNgamma12),
    'IL10 baseline' = get_sem_ci95(df$Ln_IL10_0),
    'IL10 12 months' = get_sem_ci95(df$Ln_IL10_12),
    'IL6 baseline' = get_sem_ci95(df$Ln_IL6_0),
    'IL6 12 months' = get_sem_ci95(df$Ln_IL6_12),
    'IL8 baseline' = get_sem_ci95(df$Ln_IL8_0),
    'IL8 12 months' = get_sem_ci95(df$Ln_IL8_12),
    'TNFa baseline' = get_sem_ci95(df$Ln_TNFalpha0),
    'TNFa 12 months' = get_sem_ci95(df$Ln_TNFalpha12),
    'IFNg baseline (mRNA)' = get_sem_ci95(df$transcript_IFNG_baseline),
    'IFNg 12 months (mRNA)' = get_sem_ci95(df$transcript_IFNG_12months),
    'IL10 baseline (mRNA)' = get_sem_ci95(df$transcript_IL10_baseline),
    'IL10 12 months (mRNA)' = get_sem_ci95(df$transcript_IL10_12months),
    'IL6 baseline (mRNA)' = get_sem_ci95(df$transcript_IL6_baseline),
    'IL6 12 months (mRNA)' = get_sem_ci95(df$transcript_IL6_12months),
    'IL8 baseline (mRNA)' = get_sem_ci95(df$transcript_IL8_baseline),
    'IL8 12 months (mRNA)' = get_sem_ci95(df$transcript_IL8_12months),
    'TNFa baseline (mRNA)' = get_sem_ci95(df$transcript_TNF_baseline),
    'TNFa 12 months (mRNA)' = get_sem_ci95(df$transcript_TNF_12months)
  )
  as_df <- as.data.frame(col_name)
  return(as_df)
}

main_table_2_sem <- data.frame(get_all_sems(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                               get_all_sems(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                               get_all_sems(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
)
colnames(main_table_2_sem) <- c('Placebo (SEM, CI95)', '2000 IU (SEM, CI95)', '4000 IU (SEM, CI95)')
# View(main_table_2_sem)

# Save table to file:
print(xtable(main_table_2_sem), type = "html", file = 'BESTD_table_sem_CI95.html')
############

############
# Add univariate t-tests
get_t_test <- function(i_basal, i_final) {
  i <- t.test(x = i_basal, y = i_final, paired = T)
  return(i$p.value)
}
t.test(all_data_reduced$vitd0, all_data_reduced$vitd12, paired = T)
get_t_test(all_data_reduced$vitd0, all_data_reduced$vitd12)

get_all_pvalues <- function(df) {
  col_name <- c(
    '25(OH)D 12 months' = get_t_test(df$vitd0, df$vitd12),
    'IFNg 12 months' = get_t_test(df$Ln_IFNgamma0, df$Ln_IFNgamma12),
    'IL10 12 months' = get_t_test(df$Ln_IL10_0, df$Ln_IL10_12),
    'IL6 12 months' = get_t_test(df$Ln_IL6_0, df$Ln_IL6_12),
    'IL8 12 months' = get_t_test(df$Ln_IL8_0, df$Ln_IL8_12),
    'TNFa 12 months' = get_t_test(df$Ln_TNFalpha0, df$Ln_TNFalpha12),
    'IFNg 12 months (mRNA)' = get_t_test(df$transcript_IFNG_baseline, df$transcript_IFNG_12months),
    'IL10 12 months (mRNA)' = get_t_test(df$transcript_IL10_baseline, df$transcript_IL10_12months),
    'IL6 12 months (mRNA)' = get_t_test(df$transcript_IL6_baseline, df$transcript_IL6_12months),
    'IL8 12 months (mRNA)' = get_t_test(df$transcript_IL8_baseline, df$transcript_IL8_12months),
    'TNFa 12 months (mRNA)' = get_t_test(df$transcript_TNF_baseline, df$transcript_TNF_12months)
  )
  as_df <- as.data.frame(col_name)
  return(as_df)
}

main_table_2_pvalues <- data.frame(get_all_pvalues(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), ]),
                                   get_all_pvalues(all_data_reduced[which(all_data_reduced$arm == '2000 IU'), ]),
                                   get_all_pvalues(all_data_reduced[which(all_data_reduced$arm == '4000 IU'), ])
                                   )
colnames(main_table_2_pvalues) <- c('Placebo (p-value)', '2000 IU (p-value)', '4000 IU (p-value)')
# View(main_table_2_pvalues)

t.test(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), 'vitd0'],
       all_data_reduced[which(all_data_reduced$arm == 'Placebo'), 'vitd12'],
       paired = T)

mean(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), 'vitd0'], na.rm = T)
mean(all_data_reduced[which(all_data_reduced$arm == 'Placebo'), 'vitd12'], na.rm = T)

# Save table to file:
# print(xtable(main_table_2_pvalues), type = "html", file = 'BESTD_table_pvalues.html')
write.csv(main_table_2_pvalues, 'BESTD_table_pvalues.csv', quote = FALSE, na = 'NA')
############

############
# Merge tables
# Numeric ID for re-ordering (even with merge(sort = FALSE)):
main_table_2$ID  <- 1:nrow(main_table_2)
# Set common column name:
main_table_2$variable <- rownames(main_table_2)
main_table_2_sem$variable <- rownames(main_table_2_sem)
main_table_deltas$variable <- rownames(main_table_deltas)
main_table_2_pvalues$variable <- rownames(main_table_2_pvalues)
main_table_merged <- merge(main_table_2, main_table_2_sem, by = 'variable', all = TRUE, sort = FALSE)
main_table_merged <- merge(main_table_merged, main_table_deltas, by = 'variable', all = TRUE, sort = FALSE)
main_table_merged <- merge(main_table_merged, main_table_2_pvalues, by = 'variable', all = TRUE, sort = FALSE)

# Original order:
main_table_merged <- main_table_merged[order(main_table_merged$ID), ]
main_table_merged$ID <- NULL
# View(main_table_merged)

# Save table to file:
# print(xtable(main_table_merged),
#       type = "html",
#       file = 'BESTD_table_merged.html',
#       include.rownames = FALSE,
#       NA.string = '-')
write.table(main_table_merged,
          'BESTD_table_merged.tsv',
          sep = '\t',
          quote = F,
          row.names = F,
          na = '-')

#############################################
# Save some text:
title_main_table_merged <- paste(
  'Table 1. Baseline and 12 month values following vitamin D supplementation.',
                                 sep = ''
  )
cat(file = 'title_BESTD_table_merged.tsv', title_main_table_merged, 
    # "\t", xxx_var, '\n', 
    append = FALSE)

legend_main_table_merged <- paste('Arithmetic mean, standard deviation (SD), ',
                                  'standard error of the mean (sem), 95% confidence intervals (CI95%) and two-sided, ',
                                  'univariate, paired t-test p-values shown (baseline versus 12 months within each arm). ',
                                  '\n',
                                  'Values are for observed data only. ',
                                  'Values presented are age in years, ',
                                  'body mass index (BMI) in kg/m^2, 25(OH)D in nmol/L. ',
                                  '\n',
                                  'Circulating cytokine values are natural ',
                                  'logarithm transformed. mRNA values are VSN normalised.',
                                  '\n',
                                  'P-values are shown for completeness and are not adjusted for confounding, ',
                                  'baseline values or multiple testing. ',
                                  '\n',
                                  'Adjusted regression models are shown in supplementary tables with multiple, ',
                                  'imputation for missing data.',
                                  sep = '')
legend_main_table_merged
cat(file = 'legend_BESTD_table_merged.tsv', legend_main_table_merged, 
    # "\t", xxx_var, '\n', 
    append = FALSE)
#############################################
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
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file = R_session_saved_image, compress = 'gzip')

sessionInfo()

q()
# Next: run the script for xxx
#############################