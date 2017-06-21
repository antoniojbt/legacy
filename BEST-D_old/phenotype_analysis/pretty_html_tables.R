#############################################
# BEST-D cytokine data analysis

# Author: Antonio J Berlanga-Taylor
# Date: 20 January 2017

#Purpose
#=======
# Produce pretty tables
# See:
# http://www.strengejacke.de/sjPlot/sjt.lm/
# http://www.strengejacke.de/sjPlot/sjtbasics/
# https://github.com/strengejacke/sjPlot

# Stargazer:
# https://cran.r-project.org/web/packages/stargazer/vignettes/stargazer.pdf

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
load('../data.dir/R_session_saved_image_cytokines_lm.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokine_lm_table','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################


#############################################
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('sjPlot')
# biocLite('sjmisc')
# biocLite("stargazer")

library(sjPlot)
library(sjmisc)
library(stargazer)
#############################################

#############################################
###################
# With sj libraries, HTML output
# Fit models, generate the table and specify labels:
summary(lm_Ln_IFNgamma12)
sjt.lm(lm_Ln_IFNgamma12,
       lm_Ln_IL10_12,
       lm_Ln_IL6_12,
       lm_Ln_IL8_12,
       lm_Ln_TNFalpha12,
       depvar.labels = c(  'IFNg 12 months',
                           'IL10 12 months',
                           'IL6 12 months',
                           'IL8 12 months',
                           'TNFa 12 months'),
       file = 'lm_table.html',
       group.pred = FALSE,
       pred.labels = c() # For predictor variables
)

# sjt.lm(use.viewer = FALSE) # Suppress output
# separate.ci.col = FALSE) # ci in same cell as estimates
# pred.labels = c() # For predictor variables

legend_lm_table <- paste('Cytokine values are natural logarithm',
                         'COPD = chronic obstructive pulmonary disease',
                         'BMI = body mass index',
                         'Season of randomisation were compared to spring. There were no individuals recruited during summer.',
                         '+ = reference group versus no supplementation use, female, not a smoker, no disease history or no disease incidence accordingly',
                         '* versus placebo as reference',
                         sep = '\n')
legend_lm_table
cat(file = 'legend_lm_table.html', legend_lm_table, 
    # "\t", xxx_var, '\n', 
    append = FALSE)
###################
#############################################

#############################################
###################
# NOT used:
# ANOVA tables with stargazer
# With stargazer, LATEX, HTML outputs
# Needs more work to specify output, more control though
# Example data:
stargazer(attitude)

# Run models, pass to stargazer:
stargazer(lm_Ln_IFNgamma12,
          lm_Ln_IL10_12,
          lm_Ln_IL6_12,
          lm_Ln_IL8_12,
          lm_Ln_TNFalpha12,
          dep.var.labels = c('IFNg 12 months',
                           'IL10 12 months',
                           'IL6 12 months',
                           'IL8 12 months',
                           'TNFa 12 months'),
          summary = TRUE,
          align = TRUE,
          out = 'lm_table_stargazer.html')
###################

###################
# Run models, pass to stargazer:
stargazer(anova_Ln_IFNgamma12$anova.table,
          anova_Ln_IL10_12$anova.table,
          anova_Ln_IL6_12$anova.table,
          anova_Ln_IL8_12$anova.table,
          anova_Ln_TNFalpha12$anova.table,
          dep.var.labels = c('IFNg 12 months',
                             'IL10 12 months',
                             'IL6 12 months',
                             'IL8 12 months',
                             'TNFa 12 months'),
          # summary = TRUE,
          align = TRUE,
          out = 'anova_table.html')

stargazer(anova_Ln_IFNgamma12,
          anova_Ln_IL10_12,
          anova_Ln_IL6_12,
          anova_Ln_IL8_12,
          anova_Ln_TNFalpha12,
          dep.var.labels = c('IFNg 12 months',
                             'IL10 12 months',
                             'IL6 12 months',
                             'IL8 12 months',
                             'TNFa 12 months'),
          # summary = TRUE,
          align = TRUE,
          out = 'anova_table2.html')
###################
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
# Too heavy as loadin normalisation_full RData file
# save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx
#############################