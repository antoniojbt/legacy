#############################
# Post-hoc power analysis of gene expression microarray data
# Antonio J Berlanga-Taylor
# 02 Aug 2016
#############################


#############################
# This script generates simulated datasets using vsn normalisation 
# Used to then pass to differential gene expression scripts for testing
# See:
# http://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/convergence2.pdf

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_3.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/power_analysis.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_vsn_simulation_gene_expression",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific. Check ways of making count comparisons.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
#load('R_session_saved_image_probe_filtering.RData', verbose=T)
# load('R_session_saved_image_pheno_file_check.RData', verbose=T)
# load('R_session_saved_image_diff_expression_full.RData', verbose = T)
#load('R_session_saved_image_diff_expression_3.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_vsn_simulation_gene_expression', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('sizepower')

# library(limma)
# library(ggplot2)
# library(SSPA)
# library(pwr)
# library(sizepower)
library(vsn)

# Get functions from other scripts (eg to add annotations to topTable results):

#############################

#############################
# Simulate data sets:
# Run spike for differentially expressed and null, then pass through limma pipeline.
# Do exploratory power calculation plots for two normalisations, then proper a priori and post-hoc power calc

sim_vsn <- sagmbSimulateData(n = 16700, d = 562, de = 0.20, up = 0.10, nrstrata = 6, miss = 0, log2scale = T)
class(sim_vsn)
head(sim_vsn)
str(sim_vsn)

# Transform/normalise with vsn2:
ny  <- vsn2(sim_vsn$y, strata=sim_vsn$strata)
class(ny)
head(ny)

# 
res <- sagmbAssess(exprs(ny), sim_vsn)
res

#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for xxx.
#############################