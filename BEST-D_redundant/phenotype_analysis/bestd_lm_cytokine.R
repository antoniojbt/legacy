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
# Basic ANCOVA review:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1121605/pdf/1123.pdf
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
# load('../data.dir/R_session_saved_image_cytokines_lm.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokines_lm','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################

# TO DO:
# Check 'normal' levels and ln interpretation, missing lab method.
# This file has pt_id with 0 and 12 months for each variable

#############################
# Import libraries:
library(ggplot2)
library(data.table)
library(gridExtra)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(svglite)
library(lattice)
library(mice)
library(VIM)
library(miceadds)
# library(car)
# library(gvlma)
# library(biglm)
# library(glmulti)
#############################


#############################################
# Set-up arguments:
cyto_file <- as.character(args[1])
cyto_file <- '../data.dir/bestd_cytokines.csv'

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
# IL6 4000 IU is borderline non-significant (without multiple testing)
#############################################


#############################################
##########
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

##########
# Plot one cytokine:
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
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = 'none'
  )
a2 <- ggplot(data = all_data, aes(x = vitd12, Ln_IL10_12, colour = group)) +
  geom_point(shape = 1) + # Hollow circles
  scale_colour_hue(l = 50) + # darker palette
  geom_smooth(method = lm,   # regression line
              se = FALSE    # exclude confidence region
  ) + # fullrange = TRUE # Extend regression line
  labs(title = '', x = '25OHD levels') +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 24),
        legend.title = element_blank(),
        legend.position = 'bottom'
        # legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
        # plot.margin = margin(b = 0, unit = "pt")
  )
grid.arrange(a1, a2, nrow = 2)
plots <- arrangeGrob(a1, a2, nrow = 2)
ggsave('scatterplot_IL10_and_VD.png', plots, height = 10, width = 12)
##########

##########
# Function to plot all cytkones:
scatter_plot_cyto <- function(df, y_base, y_final, 
                              x_base, x_final, 
                              group, 
                              y_label_base,
                              y_label_final,
                              x_label_base,
                              x_label_final,
                              a1_title) {
  a1 <- ggplot(data = df, aes(x = x_base, y = y_base, colour = group)) +
    geom_point(shape = 1) + # Hollow circles
    scale_colour_hue(l = 50) + # darker palette
    geom_smooth(method = lm,   # regression line
                se = FALSE    # exclude confidence region
    ) + # fullrange = TRUE # Extend regression line
    labs(title = a1_title, x = x_label_base, y = y_label_base) +
    theme_classic() +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 22),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5), # centre title
          legend.position="none"
    )
  # + xlim(0, 300) # fixes the x axis to match the second plot
  a2 <- ggplot(data = df, aes(x = x_final, y = y_final, colour = group)) +
    geom_point(shape = 1) + # Hollow circles
    scale_colour_hue(l = 50) + # darker palette
    geom_smooth(method = lm,   # regression line
                se = FALSE    # exclude confidence region
    ) + # fullrange = TRUE # Extend regression line
    labs(title = '', x = x_label_final, y = y_label_final) +
    theme_classic() +
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 24),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5), # centre title
          legend.position="none",
          plot.margin = margin(b = 0.5, unit = "pt")
    ) 
  # + xlim(0, 300) # fixes x axis, if matching first plot
  return(grid.arrange(a1, a2, nrow = 2))
}
# Test:
scatter_plot_cyto(all_data,
                  all_data$Ln_IFNgamma0, all_data$Ln_IFNgamma12, 
                  all_data$vitd0, all_data$vitd12,
                  group,
                  '',
                  '',
                  '',
                  '',
                  'IFNg')
##########

##########
# All together now
# Create lists for variables of interest at baseline and final:
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
# Create titles for each plot:
labels <- c('IFN gamma',
            'Interleukin 10',
            'Interleukin 6',
            'Interleukin 8',
            'TNF alpha')
# Run the loop to get all plots:
scatter_plot_cyto_loop <- function(base, final, label) {
  scatter_plot_cyto(all_data,
                    all_data[, base], all_data[, final],
                    all_data$vitd0, all_data$vitd12,
                    group,
                    '',
                    '',
                    '',
                    '',
                    label)
}
# scatter_plot_cyto_loop('Ln_IFNgamma0', 'Ln_IFNgamma12', 'IFNG')
# mapply runs one for one for each element, save these as a list:
all_plots <- mapply(scatter_plot_cyto_loop, cytokines_0, cytokines_12, labels)
class(all_plots)
all_plots
# Create legend separately to insert later:
# See: http://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
my_plot <- a2 # Created above for IL10, only used for legend extraction here:
# Extract legend:
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]]
  return(legend)
  }
legend <- g_legend(my_plot)
# Plot in one figure:
grid.arrange(grobs = all_plots, ncol = length(cytokines_0)) # Pass explicitely as a 'grob'
# https://stackoverflow.com/questions/38867430/add-labels-to-a-plot-made-by-grid-arrange-from-multiple-plots
# Plot legend, y-axis:
plots_y_axis <- arrangeGrob(grobs = all_plots,
                     ncol = length(cytokines_0),
                     left = textGrob("Circulating protein levels (natural logarithm)",
                                     rot = 90, vjust = 1, gp = gpar(fontsize = 30)))
plots_legend <- grid.arrange(plots_y_axis, legend,
                             nrow = 2,
                             heights = c(5, 0.2))

# grid.newpage()
# Save to file:
ggsave('scatterplots_cytokines_and_VD.png', plots_legend, height = 15, width = 24)
##########
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
#############################################


#############################################
##########
# Impute missing values
# TO DO: move to separate script
# See:
# https://www.r-bloggers.com/imputing-missing-data-with-r-mice-package/
# https://www.jstatsoft.org/article/view/v045i03
# How to formally check data is missing at random or not?

# Check proportion of missing data, usually <5%:
# cytokines_0
# cytokines_12
# covars_list
covars_list <- c('vitd12', 
                 'male',
                 'vitd0',
                 'incident_fracture',
                 'incident_resp_infection',
                 'diabetes',
                 'heart_disease',
                 'copd_asthma',
                 'basemed_vitamind',
                 'currsmoker',
                 'bmi0',
                 'calendar_age_ra',
                 'season_randomisation_2',
                 'arm')

prop_NA <- function(x) {sum(is.na(x)) / length(x) * 100}
# Individuals with more than 5% of missing variables:
apply(all_data, 1, prop_NA) # by rows
length(all_data[which(apply(all_data, 1, prop_NA) > 10), 'pt_id'])
# By columns:
apply(all_data[, covars_list], 2, prop_NA)
apply(all_data[, cytokines_0], 2, prop_NA)
apply(all_data[, cytokines_12], 2, prop_NA)

# See pattern using VIM and mice libraries
# Plot missing values:
# View(md.pattern(all_data))
png('missing_data_plots_vars_interest.png', width = 7.3, height = 5, units = 'in', res = 600)
aggr_plot <- aggr(all_data[, c(covars_list, cytokines_0, cytokines_12)],
                  only.miss = T, # Plot only missing variables
                  col = c('lightgrey', 'red'), # 1 colour for missing data, 2 observed, 3 imputed
                  numbers = T, sortVars = T,
                  labels = names(all_data[, c(covars_list, cytokines_0, cytokines_12)]),
                  cex.axis = 0.4,
                  gap = 2,
                  ylab = c('Proportion of missing data', 'Pattern'))
dev.off()
##########

###########
# Impute missing data
# Easy tutorial using library mice:
# http://web.maths.unsw.edu.au/~dwarton/missingDataLab.html

names(all_data)
all_data[1:5, 1:5]
all_data_interest <- all_data[, c(covars_list, cytokines_0, cytokines_12)]
# Keep samples IDs but exclude them from imputation (non missing and would be used to estimate
# imputation if left as column):
rownames(all_data_interest) <- all_data[, 'pt_id']
identical(all_data[, 'pt_id'], rownames(all_data_interest))
head(all_data_interest)

# Run imputation
# Roughly one imputation per percent of incomplete data (White et al.,2011),
# but the more the better, 100 can easily be run on small datasets on a laptop
# Roughly 20-30 iterations should be enough, use plot() to check convergence:
# http://stats.stackexchange.com/questions/219013/how-do-the-number-of-imputations-the-maximum-iterations-affect-accuracy-in-mul
imp_all_data <- mice(all_data_interest,
                     m = 50, # Number of imputed datasets, 5 is default
                     maxit = 50, 
                     # meth = 'pmm', # predictive mean matching, leave empty for 
                                    # auto selection depending on variable type
                     diagnostics = T,
                     seed = 500)
summary(imp_all_data)
# Plot convergence of imputed data, only plots the last 3 variables:
png('convergence_plot_imputations.png', width = 7.3, height = 5, units = 'in', res = 600)
plot(imp_all_data)
dev.off()

# Check the imputed data, e.g.:
imp_all_data$method
imp_all_data$imp$vitd12 # Each column is the number of datasets ran
png('missing_data_scatterplots_VD12_cyto12.png', width = 7.3, height = 5, units = 'in', res = 600)
xyplot(imp_all_data, 
       vitd12 ~ Ln_IFNgamma12 + 
         Ln_IL10_12 +
         Ln_IL6_12 +
         Ln_IL8_12 +    
         Ln_TNFalpha12,
       # pch = 1, cex = 1, strip = T, 
       ylab = '25OHD levels at 12 months',
       xlab = '',
       strip = strip.custom(factor.levels = labels),
       type = c('p')) # Magenta are imputed, blue observed
dev.off()
##########

###########
# Further exploratory plots:
png('missing_data_densityplots.png', width = 7.3, height = 5, units = 'in', res = 600)
densityplot(imp_all_data) # Plots all numerical variables with 2 or more missing values 
dev.off()

png('missing_data_bwplots.png', width = 7.3, height = 5, units = 'in', res = 600)
bwplot(imp_all_data)
dev.off()

# Stripplots look better:
png('missing_data_stripplots.png', width = 7.3, height = 5, units = 'in', res = 600)
stripplot(imp_all_data,
          subset = (.imp == 1 | .imp == 2 | .imp == 3 | .imp == 4 | .imp == 5 |
                    .imp == 6 | .imp == 7 | .imp == 8 | .imp == 9 | .imp == 10),
          # col = mdc(1:2), #col = mdc(1:2), pch=20, cex=1.5,
          pch = 1, cex = 0.7,
          strip = strip.custom(par.strip.text = list(cex = 0.7)))
# Magenta are imputed
dev.off()
##########

###########
# Sanity check observed and imputed data
# Variables used as predictors for imputation of each incomplete variable:
imp_all_data$pred
# Implausible results for specific variables:
which(imp_all_data$imp$vitd12 <= 1 | imp_all_data$imp$vitd12 >= 250)
which(imp_all_data$imp$calendar_age_ra <= 60 | imp_all_data$imp$calendar_age_ra >= 100)

# Create dataset with both observed and imputed data:
# imp_all_data_completed <- complete(imp_all_data, action = 'repeat', include = T)
# 'repeated' includes original data (colnames "xxx.0") and all imputations ("xxx.1, etc"). 
imp_all_data_completed <- complete(imp_all_data, include = T)
head(imp_all_data_completed)
dim(imp_all_data_completed)
length(complete.cases(imp_all_data_completed) == TRUE)
# View(imp_all_data_completed)

complete.cases(all_data_interest)
complete.cases(imp_all_data$data)
complete.cases(imp_all_data_completed)
# If mice::complete(action = 'repeat') then this should be TRUE, else FALSE like here:
# action returns the first completed imputed dataset
identical(complete.cases(all_data_interest), complete.cases(imp_all_data_completed))

sapply(all_data_interest, function(x) sum(is.na(x)))
sapply(imp_all_data$data, function(x) sum(is.na(x)))
sapply(imp_all_data_completed, function(x) sum(is.na(x)))

# These should all be TRUE:
identical(all_data$vitd12, all_data_interest$vitd12)
identical(all_data_interest$vitd12, imp_all_data$data$vitd12)
# identical(imp_all_data$data$vitd12, imp_all_data_completed$vitd12.0) # 'xxx.0' is the original data

# Next step: use mice and miceadds libraries to run regression models (bottom of script)
# with imputed dataset and pooled values from multiple imputed datasets and observed data.
##########

##########
# Basic scatterplots of changes in each group with imputed data (see above for more plots)
# TO DO: plot for all cytokines? Make transparent dots, change labels
# To DO: Needs to be separated by arm, or plot only Tx groups as placebo dilutes effect.
# Boxplot and stripchart for pre/post and change in value
# https://www.r-bloggers.com/visualizing-small-scale-paired-data-combining-boxplots-stripcharts-and-confidence-intervals-in-r/
# Similar to line plots in BESTD results paper: fig2, fig3, table2 or 3.

# Set data:
pre <- imp_all_data_completed$Ln_IL6_0
post <- imp_all_data_completed$Ln_IL6_12

# Plot:
png('boxplot_scatterplot_pre_post.png', width = 7.3, height = 5, units = 'in', res = 600)
par(mfrow = c(1, 2))
# First graph
s <- seq(length(pre))
par(bty = "l")
boxplot(pre, post, main = "Raw data", 
        xlab = "Time", ylab = "Measure", 
        names=c("pre", "post"), col = c("lightblue", "lightgreen"))
stripchart(list(pre, post), vertical = T, 
           pch = 16, method = "jitter", 
           cex = 0.5, add = T)
segments(rep(0.95, length(pre))[s], pre[s], rep(2, length(pre))[s],
         post[s], col = 1, lwd = 0.5)
# Second graph
# Confidence intervals:
res <- t.test(post, pre, paired = T, conf.int = T)
# res <- wilcox.test(post, pre, paired = T, conf.int = T)
stripchart(post-pre, vertical = T, pch = 16, 
           method = "jitter", main = "Difference", 
           ylab = "Difference: Post – Pre", xlab = "Median +/- 95% CI")
points(1, res$estimate, col = "red", pch = 16, cex = 2)
arrows(1, res$conf.int[1], 1, res$conf.int[2], col = "red",
       code = 3, lwd = 3, angle = 90)
abline(h = 0, lty = 2) # Zero-effectline
dev.off()

# TO DO: ggplot equivalent: cleanup, correct
# # arm
# # vitd12
# # cytokines_0
# # cytokines_12
# # y is cytokine levels
# # x is time
# all_data_melt_i <- melt(data = all_data, measure.vars = c('Ln_IFNgamma0', 'Ln_IFNgamma12'))
# all_data_melt_i$variable <- factor(all_data_melt_i$variable,
#                                    levels = c('Ln_IFNgamma0',
#                                               'Ln_IFNgamma12'
#                                               ),
#                                    labels = c('IFNg baseline',
#                                               'IFNg 12 months'
#                                               ))
# plyr::count(all_data_melt_i$variable)
# plyr::count(all_data_melt_i$arm2)
# group <- factor(all_data_melt_i$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
#                 labels = c("Placebo", "2000 IU", "4000 IU"))
# plyr::count(group)
# ggplot(all_data_melt_i, aes(y = value, x = variable, fill = group)) +
#   geom_boxplot()
# 
# group <- factor(all_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
#                 labels = c("Placebo", "2000 IU", "4000 IU"))
# 
# all_data$Ln_IFNgamma_delta <- all_data$Ln_IFNgamma12 - all_data$Ln_IFNgamma0
# summary(all_data$Ln_IFNgamma_delta)
# 
# ggplot(data = all_data, aes(x = vitd12, Ln_IFNgamma_delta, colour = group)) +
#   geom_point(shape = 1) + # Hollow circles
#   scale_colour_hue(l = 50) + # darker palette
#   geom_smooth(method = lm,   # regression line
#               se = FALSE    # exclude confidence region
#   ) + # fullrange = TRUE # Extend regression line
#   labs(title = '', x = 'Change in 25OHD levels', y = 'Change in cytokine levels') +
#   theme_classic() +
#   theme(text = element_text(size = 14), 
#         legend.title = element_blank(),
#         legend.position = 'bottom'
#         # legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")
#         # plot.margin = margin(b = 0, unit = "pt")
#   )
##########

###########
# TO DO:
# Consider library(Hmisc) with aregImpute() using additive regression, 
# bootstrapping, and predictive mean matching (pmm)
# It also adapts the method based on variable type automatically
# PMM (with mice or others) for numerical variables
# For categorical in mice use
# polyreg(Bayesian polytomous regression) for factor variables with >= 2 levels
# proportional odds model for ordered variables with >= 2 levels
##########
#############################################


#############################################
##########
# Setup linear model with all covariates
# Factor arm so that placebo is taken as reference group:
all_data$arm <- factor(all_data$arm, levels = c('2', '1', '0'),
                       labels = c('A_placebo',
                                  'B_2000IU',
                                  'C_4000IU'))
imp_all_data_completed$arm <- factor(imp_all_data_completed$arm, levels = c('2', '1', '0'),
                                 labels = c('A_placebo',
                                            'B_2000IU',
                                            'C_4000IU'))
plyr::count(all_data$arm)
plyr::count(imp_all_data_completed$arm)
str(all_data$arm)

# Code variables for easier referencing downstream:
y <- 'vitd12'
# Covars used in genotype regression analysis:
covars <- c('basemed_vitamind +
            male +
            calendar_age_ra +
            season_randomisation_2 +
            currsmoker +
            copd_asthma +
            heart_disease +
            incident_fracture +
            incident_resp_infection +
            diabetes +
            bmi0 +
            vitd0') # Add arm at the end for ANOVA

pass_formula <- sprintf('%s ~ %s', y, covars)
pass_formula
##########

##########
# ANOVA tests
# NOTE:
# In R, the ANOVA tests are sequential. To control for a confounder it must come first in the
# order of the model formula. For the regression output it won't matter but for the ANOVA output
# it does.
# See: 
# http://stats.stackexchange.com/questions/13241/the-order-of-variables-in-anova-matters-doesnt-it
# http://stats.stackexchange.com/questions/212496/why-do-p-values-change-in-significance-when-changing-the-order-of-covariates-in?noredirect=1&lq=1

# Order doesn't matter with type 2 or type 3 sum of squares (e.g. with car's Anova())
# but it does when wrapping as anova(lm())
# Order won't affect summary(lm()) p-values or coefficients as linear models test
# all predictors 'at once' while correcting for all.
# ANOVA requires residuals to be normally distributed
# Use glm() for count variables for instance to take of other residual error distributions
# With balanced designs order doesn't matter.
# The proportion of explained variation attributed to a given variable 
# is eta-squared (or r^2 [effect size] in lm)
# anova is by default type II, for type III car's Anova() is needed
# http://stats.stackexchange.com/questions/20452/how-to-interpret-type-i-type-ii-and-type-iii-anova-and-manova
# Type III anova deals with unbalanced designs, e.g. experimental with randomisation
# See e.g. tutorial:
# https://ww2.coastal.edu/kingw/statistics/R-tutorials/unbalanced.html
# https://ww2.coastal.edu/kingw/statistics/R-tutorials/ancova.html

# Make sure arm and other variables are coded as factors for anova
# otherwise they are run as linear regressions with 0 and 1's
class(imp_all_data)
class(imp_all_data_completed)
plyr::count(imp_all_data_completed$arm)
summary(imp_all_data_completed$arm)
summary(imp_all_data_completed$vitd12)
pass_formula <- sprintf('%s ~ %s + arm', y, covars)
pass_formula

# imputed dataset contains the original observations (imp_all_data$data) and is 
# exactly the same as all_data
class(imp_all_data$data$arm)
sapply(imp_all_data$data[, covars_list], class)

class(all_data$arm)
sapply(all_data[, covars_list], class)
##########

##########
# TO DO: check and complete assumptions tests
# See './BEST-D/eQTL_analysis/run_lm.R' script
# Assumptions linear regression
# Bivariate independent variables
# Continuous dependent variable
# Dependent variable has a normal distribution, 
# with the same variance, σ2, in each group 

# No heteroscedasticity of residuals
# i.e. the variance of residuals should not increase with 
# fitted values of the response variable. 
# Simple tutorial:
# https://www.r-bloggers.com/how-to-detect-heteroscedasticity-and-rectify-it/
# https://datascienceplus.com/linear-regression-predict-energy-output-power-plant/

# Residual analysis, diagnostic plots:
par(mfrow = c(2, 2))
plot(lm(data = imp_all_data_completed, formula = 'vitd0 ~ arm'))
dev.off()
# Without heteroscedastity a random, equal distribution of points 
# throughout the range of X axis and a flat red line should be observed
# (left side plots).

# Breusch-Pagan test for heteroscedasticity:
lmtest::bptest(lm(data = imp_all_data_completed, formula = 'vitd0 ~ arm'))
# ncv test:
car::ncvTest(lm(data = imp_all_data_completed, formula = 'vitd0 ~ arm'))
# If p-values are <0.05, reject null that heteroscedasticity is constant
# among groups for variable of interest
##########

##########
# Run anova on arm with other covariates of interest:
pass_formula <- sprintf('vitd12 ~ %s + arm', covars)
pass_formula
anova(lm(data = imp_all_data_completed, formula = 'vitd12 ~ arm'))
lm_imp_vitd12 <- lm(data = imp_all_data_completed, formula = pass_formula)
anova(lm_imp_vitd12)
# F-tests are significant for:
# arm, vitd0, incident_resp_infection, basemed_vitamind, bmi0
# Marginal for:
# diabetes, calendar_age_ra
# i.e. means between groups are different, if more than two groups, which groups actually
# differ? Run pairwise comparisons.

# Diagnostic plots of normality:
par(mfrow = c(2, 2))
plot(lm_imp_vitd12)
qqnorm(lm_imp_vitd12$residuals)
plot(lm_imp_vitd12$fitted, lm_imp_vitd12$residuals, 
     xlab = "Fitted", ylab = "Residuals")

# Pairwise comparisons and multiple testing
# Default is treatment contrasts with first group treated as baseline (alphanumerical order)
# Run as multiple regression analysis with dummy variable ('0' for baseline group, 
# '1' for others)
# If running with different contrasts do e.g.:
# contrasts(all_data$arm) <- contr.treatment(n = 3, base = 3) # Specifies placebo (arm == 2) is the baseline
# contrasts(all_data$arm)

# t-tests test whether group 0 vs other (pairwise) have the same true mean
summary(lm(data = imp_all_data_completed, formula = 'vitd12 ~ arm'))
pairwise_imp_vitd12 <- summary(lm(data = imp_all_data_completed, formula = pass_formula))
pairwise_imp_vitd12
# So here the probability of the means of arm_placebo vs armB_2000IU
# being the same is p<2e-16
# And arm_placebo vs armB_4000IU is also p<2e-16
# Use relevel() to re-order factor levels to compare e.g. armB vs armC

# All pairwise comparisons with p-value adjustment:
pairwise.t.test(x = imp_all_data_completed$vitd12, 
                g = imp_all_data_completed$arm,
                p.adjust.method = 'bonferroni')

# Plot, e.g.:
names(all_data)
ggplot(imp_all_data_completed, aes(x = arm, y = vitd12, fill = arm)) +
  geom_boxplot() +
  geom_jitter(position = position_jitter(0.2)) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = '')
  # theme(axis.text.x = element_blank())
dev.off()
##########
#############################################


#############################################
##########
# Sanity check significant cofactors:
# Significant: arm, vitd0, incident_resp_infection, basemed_vitamind, bmi0
plyr::count(all_data$arm)
summary(all_data[which(all_data$arm == 'A_placebo'), 'vitd0'])
summary(all_data[which(all_data$arm == 'B_2000IU'), 'vitd0'])
summary(all_data[which(all_data$arm == 'C_4000IU'), 'vitd0'])
kruskal.test(all_data$vitd0, all_data$arm)
anova(lm(all_data$vitd0 ~ all_data$arm))
summary(lm(all_data$vitd0 ~ all_data$arm))
pairwise.t.test(x = imp_all_data_completed$vitd0,
                g = imp_all_data_completed$arm,
                p.adjust.method = 'bonferroni')

summary(all_data[which(all_data$incident_resp_infection == 0), 'arm'])
summary(all_data[which(all_data$incident_resp_infection == 1), 'arm'])
table(all_data$incident_resp_infection, all_data$arm)
chisq.test(x = all_data$incident_resp_infection, y = all_data$arm)

summary(all_data[which(all_data$basemed_vitamind == 0), 'arm'])
summary(all_data[which(all_data$basemed_vitamind == 1), 'arm'])
# tapply(all_data$basemed_vitamind, all_data$arm, plyr::count)
table(all_data$basemed_vitamind, all_data$arm)
chisq.test(x = all_data$basemed_vitamind, y = all_data$arm)

summary(all_data[which(all_data$arm == 'A_placebo'), 'bmi0'])
summary(all_data[which(all_data$arm == 'B_2000IU'), 'bmi0'])
summary(all_data[which(all_data$arm == 'C_4000IU'), 'bmi0'])
kruskal.test(all_data$bmi0, all_data$arm)
anova(lm(all_data$bmi0 ~ all_data$arm))
summary(lm(all_data$bmi0 ~ all_data$arm))
pairwise.t.test(x = imp_all_data_completed$bmi0,
                g = imp_all_data_completed$arm,
                p.adjust.method = 'bonferroni')

# Marginal: diabetes, calendar_age_ra
summary(all_data[which(all_data$diabetes == 0), 'arm'])
summary(all_data[which(all_data$diabetes == 1), 'arm'])
table(all_data$diabetes, all_data$arm)
chisq.test(all_data$diabetes, all_data$arm)
kruskal.test(all_data$vitd12, all_data$diabetes)
anova(lm(vitd12 ~ diabetes + arm, all_data))
summary(lm(vitd12 ~ diabetes + arm, all_data))

summary(all_data[which(all_data$arm == 'A_placebo'), 'calendar_age_ra'])
summary(all_data[which(all_data$arm == 'B_2000IU'), 'calendar_age_ra'])
summary(all_data[which(all_data$arm == 'C_4000IU'), 'calendar_age_ra'])
anova(lm(all_data$calendar_age_ra ~ all_data$arm))
summary(lm(all_data$calendar_age_ra ~ all_data$arm))
##########
#############################################


#############################################
##########
# Linear models for cytokine levels before and after vitD:
# Run all cytokines in all groups, use all covariates used in genotype analysis:
cytokines_0
cytokines_12
covars
plyr::count(all_data$arm)
summary(lm(formula ='Ln_IFNgamma12 ~ vitd12', data = all_data))
summary(lm(formula ='Ln_IFNgamma12 ~ vitd12', data = imp_all_data_completed))
# TO DO: check: Polynomial regression example:
poly_lm_cyto <- lm(formula ='Ln_IFNgamma12 ~ vitd12 + I(vitd12 ^ 2)', data = imp_all_data_completed)
summary(poly_lm_cyto)
# Imputed and observed only give almost the same results
cor.test(y = all_data$Ln_IFNgamma12, x = all_data$vitd12,
         method = 'spearman', #'spearman', #'pearson',
         use = 'pairwise.complete.obs')
cor.test(y = imp_all_data_completed$Ln_IFNgamma12, x = imp_all_data_completed$vitd12,
         method = 'spearman', #'spearman', #'pearson',
         use = 'pairwise.complete.obs')
plot(y = all_data$Ln_IFNgamma12, x = all_data$vitd12)

# Simple scenario, does x interleukin correlate with vitamin D levels:
for (i in cytokines_0) {
  print(i)
  print(summary(lm(formula = sprintf('%s ~ vitd0', i), data = all_data)))
  }
# Only IL6_0 is associated to vitd0
for (i in cytokines_12) {
  print(i)
  print(summary(lm(formula = sprintf('%s ~ vitd12', i), data = all_data)))
}
# No cytokine at 12 months is associated with vitd12, regardless of arm
# TO DO: plot diagnostics for each, e.g.:
# par(mfrow = c(2, 2))
# plot(lm_imp_vitd12)
##########

##########
# Run ANOVA and pairwise for each cytokine at 12 months:
covars_list
covars
# Pairwise t-tests with p-value adjustment:
for (i in cytokines_12) {
  print(i)
  print(
    pairwise.t.test(x = imp_all_data_completed[, i],
                    g = imp_all_data_completed$arm,
                    p.adjust.method = 'bonferroni')
  )
}
# None are significant.

# ANOVAs:
for (i in cytokines_12) {
  print(i)
  pass_formula <- sprintf('%s ~ %s + arm', i, covars)
  print(pass_formula)
  print(
    anova(lm(data = imp_all_data_completed, formula = pass_formula))
    )
}
# None are significant for arm, some with age and disease variables

# Adjusting for baseline cytokine level:
for (i in cytokines_12) {
  print(i)
  index <- match(i, cytokines_12)
  print(cytokines_0[index])
  basal <- cytokines_0[index]
  pass_formula <- sprintf('%s ~ %s + %s + arm', i, covars, basal)
  print(pass_formula)
  print(
    anova(lm(data = imp_all_data_completed, formula = pass_formula))
  )
}
# None are significant for arm, some with age and disease variables
# IL6 and IFNg are borderline or significant for vitD0
##########

##########
# Final results presented are these:
# ANCOVAs
# Some discussion of tests in trials:
# https://stats.stackexchange.com/questions/3466/best-practice-when-analysing-pre-post-treatment-control-designs?noredirect=1&lq=1
# http://www.sciencedirect.com/science/article/pii/S0895435606000813
# BEST-D paper:
# "Comparisons of mean follow-up values between treatment arms involved analysis of 
# covariance (ANCOVA) adjusted, where possible, for the baseline value 
# (with multiple imputation used to impute the few missing data). 
# ANCOVA provides a more powerful test of the null hypothesis than either a comparison 
# of mean follow-up values in isolation or a comparison of mean changes 
# from baseline [38].
# [38] = 
# Borm GF, Fransen J, Lemmens WA (2007) 
# A simple sample size formula for analysis of covariance in randomized clinical trials. 
# J Clin Epidemiol 60(12):1234–1238
# PMID  17998077

# Run for 12 months with arm and adjusting for covariates:
# Include cytokine at baseline in covariates:
for (i in cytokines_12) {
  print(i)
  index <- match(i, cytokines_12)
  print(cytokines_0[index])
  basal <- cytokines_0[index]
  pass_formula <- sprintf('%s ~ %s + %s + arm', i, covars, basal)
  print(pass_formula)
  df <- imp_all_data_completed
  lm_cyto <- lm(formula = pass_formula, data = df)
  print(summary(lm_cyto))
  assign(sprintf('lm_%s', i), lm_cyto)
  # Plot diagnostics:
  # Normality, Independence, Linearity, Homoscedasticity, 
  # Residual versus Leverage graph (outliers, high leverage values and
  # influential observation's (Cook's D))
  svg(sprintf('%s_lm.svg', i))#, width = 7.3, height = 5)
  par(mfrow=c(2,2))
  plot(lm_cyto)
  dev.off()
}

# Diagnostic plots:
# See:
# http://condor.depaul.edu/sjost/it223/documents/resid-plots.gif
# http://rcompanion.org/rcompanion/e_04.html

# Interpretation:
# Running a basic lm adjusted for basic covariates in all data 
# gives null results for all five cytokines at 12 months for arm.
# Ln_IL6_12 is assoc. with vitd0, calendar_age_ra and others but not arm
# borderline significant for armB_2000IU though

# Other 12 month cytokines are associated with age, disease, etc but not arm or
# vitd variables.
# IFNg and IL6 are borderline (0.08) non-significant for armB_2000IU
# Diagnostic plots look like the distribution of residuals is normal 
# and fitted vs residuals show residuals look homoscedastic and unbiased.
##########
#############################################


#############################################
##########
# # TO DO: Plot with corrected values after lm?
# e.g.: (from gex BESTD PC correction)

# # Create dataframe with data (expression value plus PCs)
# gex_1_matched_lm <- data.frame(cbind(gex_1_matched[, 1, with = F], pc_gex_1_matched_to_adjust))
# # gex_2_matched_lm <- data.frame(cbind(gex_2_matched[, 1, with = F], pc_gex_2_matched_to_adjust))
# gex_1_matched_lm[1:5, 1:5]
# colnames(gex_1_matched_lm)[1] <- 'probe'
# 
# # Linear regression corrected for PCs:
# pass_formula <- as.formula(sprintf('%s ~ .', 'probe'))
# pass_formula
# lm_gex_1_matched <- lm(formula = pass_formula, data = gex_1_matched_lm)
# # lm_gex_2_matched <- lm(formula = pass_formula, data = gex_2_matched_lm)
# summary.lm(lm_gex_1_matched)
# 
# # Test assumptions:
# gvmodel <- gvlma(lm_gex_1_matched)
# summary(gvmodel)
# 
# # Model comparison:
# # AIC(lm_1, lm_2)
# # anova(lm_1, lm_2)
# 
# # Get results:
# coefficients(lm_gex_1_matched) # lm_fit_PCs$coefficients
# str(lm_gex_1_matched)
# # Get coefficients (estimate, t-stat, p-value) from lm:
# lm_summary <- summary.lm(lm_gex_1_matched)
# lm_summary$coefficients
# lm_summary_coef <- as.data.frame(lm_summary$coefficients)
# str(lm_summary_coef)
# 
# # Correct expression value:
# gex_1_matched_lm[1:5, 1:5]
# gex_corrected <- gex_1_matched_lm[, 'probe'] - rowSums(coef(lm_gex_1_matched)[c(-1)] * gex_1_matched_lm[, 2:ncol(gex_1_matched_lm)])
# # gex_corrected_2 <- gex_2_matched_lm[, 'probe'] - rowSums(coef(lm_gex_2_matched)[c(-1)] * gex_2_matched_lm[, 2:ncol(gex_2_matched_lm)])
# gex_corrected
##########
#############################################


#############################################
##########
# Pool results from imputed missing data and fit a linear model
# Use imputed set, not complete(), pool() then uses the multiple imputations
# For vitd12:
names(imp_all_data$data)
imp_fit <- with(imp_all_data,
                lm(formula = vitd12 ~ 
                     male +
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
                     arm))
summary(imp_fit)
pool(imp_fit)
summary(pool(imp_fit))
# Manually looked at coefficients and p-values, these are largely the same between pool()
# and lm_vd12:
pass_formula <- sprintf('vitd12 ~ %s + arm', covars)
pass_formula
lm_vd12 <- lm(formula = pass_formula, data = all_data)
summary(lm_vd12)

# TO DO Diagnostic plots of normality with imputed data:
# par(mfrow = c(2, 2))
# mice::plot.mids(pool(imp_fit)) # function doesn't exist
##########

##########
# Run anova on pooled imputations for cytokines:
# Using miceadds to be able to pool the results from multiple imputations
# Type II anova (sum of squares method) is run by default, 
# type III uses library(car) method 'Anova'
# Type III can be used for mildly unbalanced designs, such as experimental and randomised
cytokines_12
pass_formula <- sprintf('Ln_IL6_12 ~ %s + arm', covars)
pass_formula
mi.anova(mi.res = imp_all_data, formula = pass_formula, type = 2)
mi.anova(mi.res = imp_all_data, formula = pass_formula, type = 3)

# pass_formula <- sprintf('%s ~ %s + arm', i, covars)
# pass_formula <- sprintf('%s ~ %s + vitd12', i, covars)

for (i in cytokines_12) {
  print(i)
  index <- match(i, cytokines_12)
  print(cytokines_0[index])
  basal <- cytokines_0[index]
  pass_formula <- sprintf('%s ~ %s + %s + arm', i, covars, basal)
  print(pass_formula)
  # print(
    run_lm <- mi.anova(mi.res = imp_all_data, formula = pass_formula, type = 3)
  # )
  assign(sprintf('anova_%s', i), run_lm)
}

# type II and III seem to give very similar results
# No significant results for arm for type II or III
# No significant results for vitd12 for type II or III
# IFN and IL6 borderline non-significant for arm, p-values 0.0759 and 0.0704, but not vitd12.
# TO DO Diagnostic plots of normality with imputed data
##########
#############################################


#############################################
##########
# Test with interaction term: y ~ + x + z * arm, e.g.
pass_formula <- sprintf('Ln_IFNgamma12 ~ %s + arm + vitd12 * bmi0', covars)
pass_formula
lm_IFN12 <- lm(formula = pass_formula, imp_all_data_completed)
anova(lm_IFN12)
summary(lm_IFN12)
# Plots:
par(mfrow=c(1,2))
plot(Ln_IFNgamma12 ~ vitd12 + arm, data = imp_all_data_completed)
dev.off()
# Convert data frame into long format for time plot:
imp_all_data_completed$pt_id <- rownames(imp_all_data_completed)
imp_all_data_time <- imp_all_data_completed[, c('pt_id', 'Ln_IFNgamma0', 'Ln_IFNgamma12', 'arm')]
head(imp_all_data_time)
imp_all_data_time_melt <- melt(imp_all_data_time, id.vars = c('pt_id', 'arm'))
dim(imp_all_data_time_melt)
head(imp_all_data_time_melt)
# Create column with 0 month and 12 month values:
imp_all_data_time_melt$time <- ifelse(grepl(imp_all_data_time_melt$variable, 
                                            pattern = '12') == TRUE,
                                      '12',
                                      ifelse(grepl(imp_all_data_time_melt$variable, 
                                                   pattern = '0') == TRUE,
                                             '0', NA)
                                      )
imp_all_data_time_melt$time <- factor(imp_all_data_time_melt$time, 
                                      levels = c('0', '12'),
                                      labels = c('Baseline', '12 months'))
# Change response variable name:
response <- 'Ln_IFNgamma'
imp_all_data_time_melt$response <- ifelse(grepl(imp_all_data_time_melt$variable, 
                                        pattern = response) == TRUE,
                                  response, NA)
names(imp_all_data_time_melt)
plyr::count(imp_all_data_time_melt$time)
plyr::count(imp_all_data_time_melt$variable)
plyr::count(imp_all_data_time_melt$response)
plyr::count(imp_all_data_time_melt$arm)
str(imp_all_data_time_melt)
head(imp_all_data_time_melt)
anova(lm(value ~ time + arm, imp_all_data_time_melt))
# Basic summaries:
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$arm == 'A_placebo'), 'value'])
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$arm == 'B_2000IU'), 'value'])
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$arm == 'C_4000IU'), 'value'])

summary(imp_all_data_time_melt[which(imp_all_data_time_melt$time == 'Baseline'), 'value'])
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$time == '12 months'), 'value'])

summary(imp_all_data_time_melt[which(imp_all_data_time_melt$time == '12 months' &
                                   imp_all_data_time_melt$arm == 'A_placebo'), 'value'])
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$time == '12 months' &
                                   imp_all_data_time_melt$arm == 'B_2000IU'), 'value'])
summary(imp_all_data_time_melt[which(imp_all_data_time_melt$time == '12 months' &
                                   imp_all_data_time_melt$arm == 'C_4000IU'), 'value'])
# Interaction plot:
ggplot(imp_all_data_time_melt, aes(x = time, y = value)) +
  geom_boxplot() +
  facet_grid(. ~ arm)
# Actual interaction plot, set na.rm = T, else no error and no plot:
# TO DO: plot all cytokines in facets?
df <- with(imp_all_data_time_melt, aggregate(value, 
                                         list(arm = arm, 
                                              time = time), 
                                         mean, na.rm = T))
df$se <- with(imp_all_data_time_melt, aggregate(value, list(arm = arm, time = time), 
                                      function(x) sd(x)/sqrt(10)))[, 3]
df
gp <- ggplot(df, aes(x = time, y = x, colour = arm, group = arm))
gp + geom_line(aes(linetype = arm), size = 0.6, position = position_dodge(width = 0.2)) + 
  geom_point(aes(shape = arm), size = 2, position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymax = x + se, ymin = x - se), width = 0.1,
                position = position_dodge(width = 0.2)) +
  expand_limits(y = c((df$x - 1), (df$x + 1))) +
  labs(y = 'Cytokine levels (natural logarithm)', x = '')
ggsave('interaction_plot.png', width = 7.3, height = 5, units = 'in')

# Scatterplot:
xyplot(Ln_IFNgamma12 ~ vitd12 | arm, imp_all_data_completed,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.lmline(x, y, ...)
       }
)
# Diagnostic plots:
xyplot(resid(lm_IFN12) ~ fitted(lm_IFN12) | imp_all_data_completed$arm,
       xlab = "Fitted Values",
       ylab = "Residuals",
       main = "Residual Diagnostic Plot",
       panel = function(x, y, ...)
       {
         panel.grid(h = -1, v = -1)
         panel.abline(h = 0)
         panel.xyplot(x, y, ...)
       }
)
# TO DO: Diagnostic plots of normality
##########

##########
# Single cytokine test, then add interaction model:
cytokines_12
imp_fit <- with(imp_all_data,
                lm(formula = Ln_IL6_12 ~
                     # vitd12 +
                     male +
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
                     arm))
imp_fit$call
summary(imp_fit)
pool(imp_fit)
summary(pool(imp_fit))
# IL6 and IFN are significant for 2000 IU vs placebo? Not significant without covariates:
# Observed data only:
anova(with(imp_all_data_completed, lm(Ln_IL6_12 ~ vitd0 + arm)))
anova(with(imp_all_data_completed, lm(Ln_IL6_12 ~ arm)))
# With imputed data:
mi.anova(mi.res = imp_all_data, formula = 'Ln_IL6_12 ~ vitd0 + arm', type = 3)
mi.anova(mi.res = imp_all_data, formula = 'Ln_IL6_12 ~ arm', type = 3)
# With imputed data:
summary(pool(with(imp_all_data, lm(Ln_IL6_12 ~ vitd0 + arm))))
summary(pool(with(imp_all_data, lm(Ln_IL6_12 ~ arm))))
# anova is non-significant with or without correcting for baseline vitd in
# both observed and pooled imputed data anova
# TO DO check:
# pooled imputed lm is significant for 2000 IU (type 2 and 3 anovas)
# without correcting for vitd0 it is non-significant
##########

##########
# With pooled imputations:
pass_formula <- sprintf('Ln_IFNgamma12 ~ %s + arm + vitd12 * bmi0', covars)
pass_formula
mi.anova(mi.res = imp_all_data, formula = pass_formula)
# Quotes cause errors for lm with imputed data in mice (?):
imp_fit_interaction <- with(imp_all_data, 
                            lm(formula = Ln_IFNgamma12 ~ 
                                 # vitd12 + 
                                 male +
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
                                 + arm +
                                 vitd12 * bmi0))                                             
imp_fit_interaction$call
summary(imp_fit_interaction)
pool(imp_fit_interaction)
summary(pool(imp_fit_interaction))
# Not significant for any cytokine for interaction between arms and vitd12.

# Formally test if interaction significantly contribute to the model, e.g.:
# anova(lm.without_interaction, lm.interaction)
# Compare two models with mice::pooled.compare(), pool.r.square(), etc. e.g.:
# Check variance explained:
imp_fit$call
imp_fit_interaction$call
pool.r.squared(imp_fit, adjusted = T)
pool.r.squared(imp_fit_interaction, adjusted = T)
# Wald test and p-value for comparison:
compare_imp_models <- pool.compare(imp_fit_interaction, # larger model first
                                   imp_fit,
                                   # data = imp_all_data,
                                   method = 'Wald') #'likelihood')
compare_imp_models$pvalue
# Results are non-significant, interaction term does not add to the model

# ANOVA for all cytokines at 12 months with pooled imputations and interaction term:
for (i in cytokines_12) {
  pass_formula <- sprintf('%s ~ %s + arm + vitd12 * bmi0', i, covars)
  print(pass_formula)
  print(
    mi.anova(mi.res = imp_all_data, formula = pass_formula)
  )
}
# No significant results
# TO DO Diagnostic plots of normality with imputed data
#############################################


#############################################
# Only for reference, not used:
# Linear model for change in vitd (final minus baseline):
y <- 'delta'
pass_formula <- sprintf('%s ~ %s + arm', y, covars)
pass_formula
lm_delta <- lm(formula = pass_formula, data = all_data)
summary(lm_delta)
# Results as above/expected
#############################################



#############################################
# Run mixed effects?
# # See: http://stackoverflow.com/questions/1169539/linear-regression-and-group-by-in-r
library(lme4)
# library(nlme)
# 
# # Another scatterlot example:
# # xyplot(Ln_IL10_12 ~ vitd12, groups = arm, data = all_data, type = c('p', 'r', 'g'))
# # Run linear regression by group:
# summary(lmList(formula = vitd12 ~
#                  male +
#                  vitd0 +
#                  incident_fracture +
#                  incident_resp_infection +
#                  diabetes +
#                  heart_disease +
#                  copd_asthma +
#                  basemed_vitamind +
#                  currsmoker +
#                  bmi0 +
#                  calendar_age_ra +
#                  season_randomisation_2 | arm,
#                data = all_data))

# Run for all cytokines with covariates, exclude arm:
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

# TO DO: Run with imp_all_data_completed (single imputed, complete dataset)
# Mixed effects model:
lmm_vitd12 <- lmer(formula = vitd12 ~ vitd0 +
               male +
               bmi0 +
               calendar_age_ra +
               (1 | arm),
               data = all_data)

lmm_null <- lmer(formula = vitd12 ~ (1 | arm),
                     data = all_data)
anova(lmm_vitd12, lmm_null)

summary(lmm_vitd12)
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
citation()
citation("ggplot2")
q()
# Next: run the script for xxx
#############################