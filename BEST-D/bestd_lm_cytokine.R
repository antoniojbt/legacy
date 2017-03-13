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
str(all_data)
all_data$vitd12 <- as.numeric(all_data$vitd12)
summary(all_data$vitd12)
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
str(all_data)
dim(all_data)
#############################################


#############################################
# Explore cytokine data, descriptive analysis

# Rename groups for better plotting:
all_data$arm2[all_data$arm == 0] <- '4000_IU'
all_data$arm2[all_data$arm == 1] <- '2000_IU'
all_data$arm2[all_data$arm == 2] <- 'Placebo'
plyr::count(all_data$arm2)

# Get variables of interest:
colnames(all_data)
all_data_melt <- melt(all_data, measure.vars = c(29, 42, 88:97))
all_data_melt
head(all_data_melt)[1:5, 1:5]
dim(all_data_melt)
colnames(all_data_melt)
plyr::count(all_data_melt$variable)
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
plyr::count(all_data_melt$variable)
plyr::count(all_data_melt$arm2)
group <- factor(all_data_melt$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)
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
# Some t tests on cytokine levels:
t.test(all_data_placebo$Ln_IFNgamma0, all_data_placebo$Ln_IFNgamma12, paired =  T)
t.test(all_data_2000$Ln_IFNgamma0, all_data_2000$Ln_IFNgamma12, paired =  T)
wilcox.test(all_data_2000$Ln_IFNgamma0, all_data_2000$Ln_IFNgamma12, paired =  T)
t.test(all_data_4000$Ln_IFNgamma0, all_data_4000$Ln_IFNgamma12, paired =  T)
wilcox.test(all_data_4000$Ln_IFNgamma0, all_data_4000$Ln_IFNgamma12, paired =  T)

t.test(all_data_placebo$Ln_IL10_0, all_data_placebo$Ln_IL10_12, paired = T)
t.test(all_data_2000$Ln_IL10_0, all_data_2000$Ln_IL10_12, paired = T)
t.test(all_data_4000$Ln_IL10_0, all_data_4000$Ln_IL10_12, paired = T)

t.test(all_data_placebo$Ln_IL6_0, all_data_placebo$Ln_IL6_12, paired = T)
t.test(all_data_2000$Ln_IL6_0, all_data_2000$Ln_IL6_12, paired = T)
t.test(all_data_4000$Ln_IL6_0, all_data_4000$Ln_IL6_12, paired = T)

t.test(all_data_placebo$Ln_IL8_0, all_data_placebo$Ln_IL8_12, paired = T)
t.test(all_data_2000$Ln_IL8_0, all_data_2000$Ln_IL8_12, paired = T)
t.test(all_data_4000$Ln_IL8_0, all_data_4000$Ln_IL8_12, paired = T)

t.test(all_data_placebo$Ln_TNFalpha0, all_data_placebo$Ln_TNFalpha12, paired = T)
t.test(all_data_2000$Ln_TNFalpha0, all_data_2000$Ln_TNFalpha12, paired = T)
t.test(all_data_4000$Ln_TNFalpha0, all_data_4000$Ln_TNFalpha12, paired = T)

# IL10 4000 IU is borderline significant (without multiple testing)
#############################################


#############################################
# Basic scatterplots of changes in each group:
plot(all_data$Ln_IFNgamma12, all_data$vitd12)
plot(all_data_placebo$Ln_IL10_12, all_data_placebo$Ln_IL10_0)
plot(all_data_2000$Ln_IL10_12, all_data_2000$Ln_IL10_0)
plot(all_data_4000$Ln_IL10_12, all_data_4000$Ln_IL10_0)
boxplot(all_data_placebo$Ln_IL10_0,
        all_data_placebo$Ln_IL10_12,
        all_data_2000$Ln_IL10_0,
        all_data_2000$Ln_IL10_12,
        all_data_4000$Ln_IL10_0,
        all_data_4000$Ln_IL10_12)


##########
group <- factor(all_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)

a1 <- ggplot(data = all_data, aes(x = vitd0, Ln_IFNgamma0, colour = group)) +
      geom_point(shape = 1) + # Hollow circles
      scale_colour_hue(l = 50) + # darker palette
      geom_smooth(method = lm,   # regression line
                  se = FALSE    # exclude confidence region
                  ) + # fullrange = TRUE # Extend regression line
      labs(title = '', x = '') +
      theme_classic() +
      theme(text = element_text(size = 14), 
            legend.title=element_blank()
        )
a1
a2 <- ggplot(data = all_data, aes(x = vitd12, Ln_IFNgamma12, colour = group)) +
      geom_point(shape = 1) + # Hollow circles
      scale_colour_hue(l = 50) + # darker palette
      geom_smooth(method = lm,   # regression line
                  se = FALSE    # exclude confidence region
                  ) + # fullrange = TRUE # Extend regression line
      labs(title = '') +
      theme_classic() +
      theme(text = element_text(size = 14), 
            legend.title=element_blank()
            )
a2
grid.arrange(a1, a2, nrow = 2)
# plots <- arrangeGrob(a1, a2, nrow = 2)
# plots
# ggsave('scatterplots_cytokines_and_VD.png', 
#        plots, height = 10, width = 12)
##########

##########
group <- factor(all_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)

a1 <- ggplot(data = all_data, aes(x = vitd0, Ln_IL10_0, colour = group)) +
  geom_point(shape = 1) + # Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm,   # regression line
              se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '') +
  theme_classic() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank()
  )
a1
a2 <- ggplot(data = all_data, aes(x = vitd12, Ln_IL10_12, colour = group)) +
  geom_point(shape = 1) + # Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm,   # regression line
              se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '25OHD levels') +
  theme_classic() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank()
  )
a2
grid.arrange(a1, a2, nrow = 2)
# plots <- arrangeGrob(a1, a2, nrow = 2)
# plots
# ggsave('scatterplots_cytokines_and_VD.png', 
#        plots, height = 10, width = 12)
##########
#############################################


#############################################
##########
# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups
# Extract variables of interest, assign groups, run lm correcting for baseline 
# and other groups
class(all_data_melt)
colnames(all_data_melt)
class(all_data)
head(all_data)
colnames(all_data)

# Covars used in genotype regression analysis:
# male +
#   vitd0 +
#   incident_fracture +
#   incident_resp_infection +
#   diabetes +
#   heart_disease +
#   copd_asthma +
#   basemed_vitamind +
#   currsmoker +
#   bmi0 +
#   calendar_age_ra +
#   season_randomisation_2 +
#   arm'

# Sanity check treatment groups and vitamin D levels:
str(all_data$arm)
all_data[1:20, c('arm', 'vitd12')]
summary(all_data[which(all_data$arm == 0), 'vitd12']) # 4000IU
summary(all_data[which(all_data$arm == 1), 'vitd12']) # 2000IU
summary(all_data[which(all_data$arm == 2), 'vitd12']) # Placebo

summary(all_data[which(all_data$arm == 0), 'vitd0'])
summary(all_data[which(all_data$arm == 1), 'vitd0'])
summary(all_data[which(all_data$arm == 2), 'vitd0'])

all_data$delta <- (all_data$vitd12 - all_data$vitd0)
summary(all_data$delta)
summary(all_data[which(all_data$arm == 0), 'delta'])
summary(all_data[which(all_data$arm == 1), 'delta'])
summary(all_data[which(all_data$arm == 2), 'delta'])
t.test(all_data_4000$vitd12, all_data_4000$vitd0, paired = T)
##########

##########
# Linear model with all covariates
# Factor arm so that placebo is taken as reference group:
all_data$arm <- factor(all_data$arm, levels = c('2', '1', '0'),
                       labels = c('A_placebo',
                                  'B_2000IU',
                                  'C_4000IU'))
plyr::count(all_data$arm)
str(all_data$arm)

# Code variables for easier referencing downstream:
y <- 'vitd12'
covars <- c('male +
            vitd0 +
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
            arm')

pass_formula <- sprintf('%s ~ %s', y, covars)
pass_formula

lm_cyto <- lm(formula = pass_formula, data = all_data)
summary(lm_cyto)

# Sanity check basemed_vitamind, diabetes, bmi0 and age
# as came out significant:
# TO DO:
summary(all_data[which(all_data$diabetes == 0), 'arm'])
summary(all_data[which(all_data$diabetes == 1), 'arm'])
kruskal.test(all_data$diabetes, all_data$arm)
anova(lm(all_data$diabetes ~ all_data$arm))

plyr::count(all_data$arm)
summary(all_data[which(all_data$arm == 'A_placebo'), 'bmi0'])
summary(all_data[which(all_data$arm == 'B_2000IU'), 'bmi0'])
summary(all_data[which(all_data$arm == 'C_4000IU'), 'bmi0'])
anova(lm(all_data$bmi0 ~ all_data$arm))

summary(all_data[which(all_data$arm == 'A_placebo'), 'calendar_age_ra'])
summary(all_data[which(all_data$arm == 'B_2000IU'), 'calendar_age_ra'])
summary(all_data[which(all_data$arm == 'C_4000IU'), 'calendar_age_ra'])
anova(lm(all_data$calendar_age_ra ~ all_data$arm))

summary(all_data[which(all_data$basemed_vitamind == 0), 'arm'])
summary(all_data[which(all_data$basemed_vitamind == 1), 'arm'])
# tapply(all_data$basemed_vitamind, all_data$arm, plyr::count)
kruskal.test(all_data$basemed_vitamind, all_data$arm)

# TO DO: run with contrasts?
# contrasts(all_data$arm) <- contr.treatment(n = 3, base = 3) # Specifies placebo (arm == 2) is the baseline
# contrasts(all_data$arm)
##########

##########
# Linear model for change in vitd (final minus baseline):
y <- 'delta'
pass_formula
lm_cyto <- lm(formula = pass_formula, data = all_data)
summary(lm_cyto)
# Results as expected, double check but seem fine.
##########


##########
# Linear models for cytokine levels before and after vitD:
# Simple scenario, does x interleukin correlate with vitamin D levels:
lm_cyto <- lm(formula = 'Ln_IL10_0 ~ vitd0', data = all_data)
summary(lm_cyto)
lm_cyto <- lm(formula = 'Ln_IL10_12 ~ vitd12', data = all_data)
summary(lm_cyto)

# Run all cytokines in all groups, use all covariates used in genotype analysis:
cytokines_0 <- c('Ln_IFNgamma0',
               'Ln_IL10_0',
               'Ln_IL6_0',
               'Ln_IL8_0',
               'Ln_TNFalpha0'
               )
cytokines_12 <- c('Ln_IFNgamma12',
               'Ln_IL10_12',
               'Ln_IL6_12',
               'Ln_IL8_12',
               'Ln_TNFalpha12'
               )
covars
plyr::count(all_data$arm)

# i <- 'Ln_TNFalpha12'
for (i in cytokines_12) {
  print(i)
  pass_formula <- sprintf('%s ~ vitd12 + %s', i, covars)
  print(pass_formula)
  df <- all_data
  print(summary(lm(formula = pass_formula, data = df)))
  # # Check results and diagnostic plots:
  # fitted(lm_cyto)
  # residuals(lm_cyto)
  # plot(all_data$vitd12, all_data$Ln_IL10_12)
  # abline(lm_cyto)
  # # Plot diagnostics
  # # Normality, Independence, Linearity, Homoscedasticity, Residual versus Leverage graph (outliers, high leverage values and
  # # influential observation's (Cook's D))
  # par(mfrow=c(2,2)) 
  # plot(lm_cyto)
  }

# Interpretation:
# Running a basic lm adjusted for basic covariates in all data 
# gives null results for all five cytokines at 12 months for vitd12.
# Ln_IL6_12 is assoc. with vitd0, calendar_age_ra and others but not vitd12 and is
# borderline significant for armB_2000IU

# Other 12 month cytokines are associated with age, disease, etc but not arm or
# vitd variables.
##########


##########
# TO DO: run random effects?
library(lme4)
library(nlme)
lm_cyto <- lme('formula = vitd12 ~ vitd0 +
               male +
               bmi0 +
               calendar_age_ra +
               (1 | arm)',
               data = all_data)


lm_cyto_null <- lmer('formula = vitd12 ~ (1 | arm)',
                     data = all_data)
anova(lm_cyto, lm_cyto_null)
summary(lm_cyto)
##########
#############################################


#############################################
# Association with genotype, SNPs tested only

#############################################


#############################################
## Save some text:
# cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
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
