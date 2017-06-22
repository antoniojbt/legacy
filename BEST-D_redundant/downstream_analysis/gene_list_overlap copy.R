####################
# Gene list overlap analysis
# 6 May 2016
####################

####################
# See:
# https://www.bioconductor.org/packages/3.3/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf
# If more sensitive is needed, use rank based methods:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943622/
# http://www.bioconductor.org/packages/release/bioc/vignettes/RRHO/inst/doc/RRHO.pdf
####################

####################
options(echo = TRUE)
##Set working directory and file locations and names of required inputs:
# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('Desktop/scripts_to_upload/GAT_backgrounds_2/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_", Sys.Date(),".txt", sep=""))
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
# load('R_session_saved_image_gene_list_overlap.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_gene_list_overlap','.RData', sep='')
R_session_saved_image
####################

####################
# Libraries:
# source("https://bioconductor.org/biocLite.R")
# biocLite("GeneOverlap")
# install.packages("GeneOverlap")
library(GeneOverlap)
library(data.table)
library(ggplot2)

# Get additional functions needed:
source('gene_expression_functions.R')
# source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene or SNP files:
hits_file <- as.character(args[1])
# hits_file <- 'Diff_exp_BESTD_FDR10.txt'
# hits_file <- 'GAinS_antonio.txt'
# hits_file <- 'full_topTable_pairing_all_treated_probeID_FDR10.txt'
# hits_file <- 'all_treated_baseline_FDR0.05_all_columns.txt'
# hits_file <- 'VD_diff_exp_Hossein-nezhad_Holick_2013_TS4.txt'

background_file <- as.character(args[2])
# background_file <- 'background_genes_BESTD_expressed.txt'
# background_file <- 'background_probeID_BESTD_expressed.txt'

# To test enrichment against a particular set
annotation_file <- as.character(args[3])
# annotation_file <- 'GAinS_Emma_table_B.txt'
# annotation_file <- 'VD_diff_exp_Heikkinen11.txt'
# annotation_file <- 'VD_diff_exp_Hossein-nezhad_Holick_2013_TS4.txt'
# annotation_file <- 'VD_diff_exp_Gerke2014.txt'
# annotation_file <- 'eGenes_BESTD_cis_tx_fdr5_reQTLs_annot.txt'
# annotation_file <- 'all_treated_baseline_FDR0.05_all_columns.txt'
# annotation_file <- 'all_treated_final_FDR0.05_all_columns.txt'
# annotation_file <- 'for_overlap_GO-0071305_cellular_response_to_vitamin_D.txt'

# Ontology term to test against:
# ontology_term <- as.character(args[3])
####################

####################
# Read data:
# hits_data <- fread(hits_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
hits_data <- read.csv(hits_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
class(hits_data)
str(hits_data)
head(hits_data)
dim(hits_data)

background_data <- read.csv(background_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
str(background_data)
head(background_data)
dim(background_data)

annotation_data <- read.csv(annotation_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE,
                            na.string = c('-', '', ' ', 'NA', 'na'))
str(annotation_data)
head(annotation_data)
dim(annotation_data)
####################

####################
# Visual check that whitespaces, names, symbols, etc, are similar in the three files:
head(hits_data)
head(background_data)
head(annotation_data)

length(unique(annotation_data$V1))
dim(annotation_data)

hits_data[1:10, 1]
background_data[1:10, 1]
annotation_data[1:10, 1]

hits_data[(nrow(hits_data)-10):nrow(hits_data), 1]
background_data[(nrow(background_data)-10):nrow(background_data), 1]
annotation_data[(nrow(annotation_data)-10):nrow(annotation_data), 1]
####################

####################
# Run overlap test:
# go.obj <- newGeneOverlap(hits_data$V1, annotation_data$V1, spec = "hg19.gene")
go.obj <- newGeneOverlap(hits_data$V1, annotation_data$V1, genome.size = length(background_data$V1))
# See results of basic overlaps:
go.obj
# Perform statistical test (Fisher's):
go.obj <- testGeneOverlap(go.obj)
go.obj
# Detailed results:
print(go.obj)
# Get specific slots of object (ie overlapping gene set):
head(getIntersection(go.obj))
# Visualise, errors:
# pdf(sprintf('heatmap_overlap_%s_%s_%s.pdf', hits_file, annotation_file, background_file))
# gom.obj <- newGOM(as.list(hits_data$V1), as.list(annotation_data$V1), spec = "hg19.gene")
# drawHeatmap(gom.obj)
# dev.off()
####################

####################
# TO DO: check this package, switch to it:
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2943622/
# http://www.bioconductor.org/packages/release/bioc/html/RRHO.html
# http://www.bioconductor.org/packages/release/bioc/vignettes/RRHO/inst/doc/RRHO.pdf
####################
  
####################
# Interpretation:
# 
# 
####################


####################
# Write results to disk:
file_name <- sprintf('overlap_%s_VS_%s_ON_%s.txt', hits_file, annotation_file, background_file)
sink(file_name, append = FALSE, split = TRUE, type = c("output"))
print(hits_file)
print(annotation_file)
print(background_file)
print(go.obj)
getIntersection(go.obj)
sink(file = NULL)
####################

####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run the script for xxx.
####################
