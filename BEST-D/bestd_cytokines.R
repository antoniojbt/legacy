#############################################
# BEST-D cytokine data analysis

# Author: Antonio J Berlanga-Taylor
# Date: 20 January 2017

#Purpose
#=======

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

output_file <- file(paste("R_session_output_cytokines",Sys.Date(),".txt", sep=""))
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
#load('R_session_saved_image_order_and_match.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokines','.RData', sep='')
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
#############################


#############################################
# Set-up arguments:

cyto_file <- as.character(args[1])
# cyto_file <- 'bestd_cytokines.csv'

# See script 
# /cgat_projects/BEST-D/genotype_analysis/plink_output_processing.R
geno_file <- as.character(args[2])
# geno_file <- ''

# TO DO: genes need extracting and processing (see geno file as example)
# See /Users/antoniob/Documents/github.dir/AntonioJBT/Airwave/eQTL_plotting_airwave.R
expr_file <- as.character(args[3])
# expr_file <- ''

pheno_file <- as.character(args[4])
# pheno_file <- 'final_phenotype_file.txt'

#############################################


#############################################
# Read files:
cyto_file <- fread(cyto_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
cyto_file
head(cyto_file)
dim(cyto_file)
tail(cyto_file)
summary(cyto_file)

geno_file <- fread(geno_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
geno_file
head(geno_file)
dim(geno_file)
tail(geno_file)
summary(geno_file)

expr_file <- fread(expr_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
expr_file
head(expr_file)
dim(expr_file)
tail(expr_file)
summary(expr_file)

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
# TO DO: check expr ID as this is kit_id, not pt_id

# Join pheno and cyto:
all_data <- pheno_file[cyto_file]
all_data

# Join expr:
all_data <- all_data[expr_file]
all_data

# Join geno:
all_data <- all_data[geno_file]
all_data

# Sanity check:

#############################################


#############################################
# Explore cytokine data, descriptive analysis

# Plot:
png(paste('qqplot_', SNP_file, '.png', sep = ''))
plot(cytokines_file)
dev.off()
#############################################


#############################################
# Plot gene expression with cytokines as main figure

# Create table with values of interest:


# Save table


# Plot:
png(paste(SNP, '_boxplot_by_genotype_arm_12mo.png', sep=''), width = 6, height = 6, units = 'in', res = 300)
group <- factor(merged_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), labels = c("Placebo", "2000 IU", "4000 IU"))
genotype_order <- factor(merged_data$genotype, levels = c("CC", "AC", "AA"))
ggplot(aes(y = vitd12, x = genotype_order, fill = group), data = merged_data) + 
  labs(x = paste(SNP, '(GC)'), y = '25(OH)D levels (nmol/L) at 12 months') + 
  geom_boxplot(position = position_dodge(1)) + 
  # geom_point(position = position_jitter(width = 0.2)) +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank())
dev.off()

#############################################


#############################################
# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups

#############################################


#############################################
# Correlation to VD and other measures?

#############################################


#############################################
# Association with genotype, SNPs tested only

#############################################


#############################################
## Save some text:
cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
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
