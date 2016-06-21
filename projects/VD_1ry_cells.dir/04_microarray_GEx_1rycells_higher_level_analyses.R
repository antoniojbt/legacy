#############################
# To be run after 03 differential gene expression analysis of array data

#############################


#############################
# Illumina gene expression microarray analysis steps (Ritchie et al. 2014 PLoS Comp Bio; 
# Ritchie et al. 2015 Nucleic Acids Res.):

# 1) Data acquisition
# 2) Preprocessing and quality assessment

# 3) Background correction, normalisation, transformation and filtering

# 4) Experimental design matrix specification
# 5) Descriptive statistics
# 6) Differential gene expression analysis

# 7) Higher level analysis: Gene ontology, co-expression, gene set analyses
# 8) Integration with other data

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
setwd("/ifs/projects/proj043/analysis.dir/gene_expression.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_read_and_QC.RData', verbose=T)
load('R_session_saved_image_diff_expression', verbose=T)

# To load multiple .RData files:
#rdata_filenames <- c('.RData')
#lapply(rdata_filenames, load, .GlobalEnv)


# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_higher_level', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_higher_level_full', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("limma")
#install.packages('dendextend')

library(limma)
library(beadarray)
library(lumi)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(genefilter)
library(Biobase)
library(quantro)
library(flashClust)
library(lumiHumanAll.db)
library(vsn)
library(reshape2)
library(gridExtra)
library(plyr)
library(dendextend)
library(gplots)
library(doParallel)

#############################


#############################

# 7) Integration with other data
# a) eQTL analysis
# b) Disease variant overlap
# c) Pathway analysis


#############################



#############################
#The end:
# TO DO: remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes)))

#rm()

# To save R workspace with all objects to use at a later time:
# Also consider dump(), dput() and dget(), see http://thomasleeper.com/Rcourse/Tutorials/savingdata.html
save.image(file=R_session_saved_image_full, compress='gzip')

# Or save specific objects:
# objects_to_save <- (c('normalised_expressed', 'membership_file_cleaned', 'FAILED_QC'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: collect data, plots and tables for a report
#############################