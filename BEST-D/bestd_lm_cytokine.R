#############################################
# BEST-D cytokine data analysis

# Author: Antonio J Berlanga-Taylor
# Date: 03 March 2017

#Purpose
# Regression analysis for BEST-D study for circulating cytokine levels

#############################################
# Cytokine data from BEST-D study, J. Emberson email 23 Jan 2017
# BEST-D cytokine data (IFN gamma, IL-6, IL-8, IL-10, and TNF-alpha). 
# All are measured at 0 and 12 months (variable suffix 0 or 12). Units for all are Ln pg/ml.

# All null from CTSU analysis
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
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/BEST-D_cytokines/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_cytokines_lm",Sys.Date(),".txt", sep=""))
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
# load('../data.dir/R_session_saved_image_cytokines.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokines_lm','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################

# TO DO:
# Check 'normal' levels and ln interpretation, missing lab method.
# Join info with geno, expr, VD, pheno:
# Create one file with pt_id, kit_id, candidate geno, candidate expr, pheno and VD
# This file has pt_id with 0 and 12 months for each variable
# Pull in above data and join
# Plot gene expression with cytokines as main figure

# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups
# Correlation to VD and other measures?
# Association with genotype, SNPs tested only

#############################
# Import libraries:
library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)
library(plyr)
library(dplyr)
#############################


#############################################
# Set-up arguments:

cyto_file <- as.character(args[1])
cyto_file <- 'bestd_cytokines.csv'

pheno_file <- as.character(args[4])
pheno_file <- '../data.dir/BEST-D_phenotype_file_final.tsv'
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
#############################################

#############################################
# Create one file with pt_id, kit_id, candidate geno, candidate expr, pheno and VD
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

#############################################


#############################################
# Explore cytokine data, descriptive analysis

# Rename groups for better plotting:
all_data$arm2[all_data$arm == 0] <- '4000_IU'
all_data$arm2[all_data$arm == 1] <- '2000_IU'
all_data$arm2[all_data$arm == 2] <- 'Placebo'
count(all_data$arm2)

# Get variables of interest:
colnames(all_data)
all_data_melt <- melt(all_data, measure.vars = c(29, 42, 88:97))
all_data_melt
head(all_data_melt)[1:5, 1:5]
dim(all_data_melt)
colnames(all_data_melt)
count(all_data_melt$variable)
all_data_melt[1:10, c('value', 'variable')]
all_data_melt$variable <- factor(all_data_melt$variable, 
                                 levels = c('vitd0',
                                            'vitd12',
                                            'Ln_IFNgamma0',
                                            'Ln_IFNgamma12',
                                            'Ln_IL10_0',
                                            'Ln_IL10_12',
                                            'Ln_IL6_0',
                                            'Ln_IL6_12',
                                            'Ln_IL8_0',
                                            'Ln_IL8_12',
                                            'Ln_TNFalpha0',
                                            'Ln_TNFalpha12'),
                                 labels = c('25OHD baseline',
                                            '25OHD 12 months',
                                            'IFNg baseline',
                                            'IFNg 12 months',
                                            'IL10 baseline',
                                            'IL10 12 months',
                                            'IL6 baseline',
                                            'IL6 12 months',
                                            'IL8 baseline',
                                            'IL8 12 months',
                                            'TNFa baseline',
                                            'TNFa 12 months')
                                 )
count(all_data_melt$variable)
count(all_data_melt$arm2)
group <- factor(all_data_melt$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
count(group)
#############################################


#############################################
# Basic significance tests
# Subgroup:
colnames(all_data)
all_data[, 88:ncol(all_data)]
summary(all_data[, 88:ncol(all_data)])
all_data$arm2 <- as.factor(all_data$arm2)
all_data_placebo <- all_data[which(all_data$arm2 == 'Placebo'), ]
all_data_2000 <- all_data[which(all_data$arm2 == '2000_IU'), ]
all_data_4000 <- all_data[which(all_data$arm2 == '4000_IU'), ]
dim(all_data_placebo)
dim(all_data_2000)
dim(all_data_4000)
# t tests:
t.test(all_data_4000$transcript_IFNG_baseline, all_data_4000$transcript_IFNG_12months)
t.test(all_data_4000$Ln_IFNgamma0, all_data_4000$Ln_IFNgamma12)
t.test(all_data_2000$Ln_IFNgamma0, all_data_2000$Ln_IFNgamma12)
#############################################


#############################################
# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups
# TO DO: extract variables of interest, assign groups, run lm correcting for baseline 
# and other groups
# TO DO: from here
class(all_data_melt)
colnames(all_data_melt)
lm_fit_cyto <- lm.fit(formula = ' ~ .',
                      data = as.matrix(all_data_melt_matrix))

fit <- lm(weight ~ height, data = women)
summary(fit)
fitted(fit)
residuals(fit)
plot(women$height,women$weight, 
     xlab="Height (in inches)",
     ylab="Weight (in pounds)")
abline(fit)
# Plot diagnostics
# Normality, Independence, Linearity, Homoscedasticity, Residual versus Leverage graph (outliers, high leverage values and
# influential observation's (Cook's D))
par(mfrow=c(2,2)) 
plot(fit) 
#############################################


#############################################
# Correlation to VD and other measures?

#############################################


#############################################
# Association with genotype, SNPs tested only

#############################################


#############################################
## Save some text:
# cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
#############################################

#############################################
# 

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
save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx
#############################
