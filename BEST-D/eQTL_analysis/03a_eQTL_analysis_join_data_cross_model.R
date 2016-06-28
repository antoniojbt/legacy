#############################################
# eQTL file processing script - join files for MODEL_cross in MatrixEQTL
# Antonio J Berlanga-Taylor
# 30 Nov 2015

#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""), open = 'a')
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
# load('R_session_saved_image_processed_files.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
# R_session_saved_image <- paste('R_session_saved_image_order_and_match_2','.RData', sep='')
R_session_saved_image <- paste('R_session_saved_image_order_and_match','.RData', sep='')
R_session_saved_image

#############################################


#############################################
## Update packages if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite('data.table')

#and load them:
packages <-c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse', 'Rsge', 'mePipe', 
             'readr', 'reshape2', 'data.table')
lapply(packages, require, character.only = TRUE)
sessionInfo()
#############################################

#############################################
## Add group membership to PC covariates file:
# 4000 data only for now: 
covar_4000_baseline <- read.csv('covar_PCAs_t_4000_baseline.tsv', sep = '\t', check.names = FALSE)
covar_4000_final <- read.csv('covar_PCAs_t_4000_final.tsv', sep = '\t', check.names = FALSE)

class(covar_4000_baseline)
head(covar_4000_baseline)
tail(covar_4000_baseline)
dim(covar_4000_baseline)
dim(covar_4000_final)

# Add membership as last covariate as MatrixEQTL reads this for interaction analysis:
insertRow2 <- function(existingDF, row_number, newrow){
  existingDF <- rbind(existingDF[1:row_number, ], newrow, existingDF[-(1:row_number), ])
  return(existingDF)
}

covar_4000_baseline <- insertRow2(covar_4000_baseline, 20, '0')
dim(covar_4000_baseline)
head(covar_4000_baseline)
tail(covar_4000_baseline)

covar_4000_final <- insertRow2(covar_4000_final, 20, '1')
dim(covar_4000_final)
head(covar_4000_final)
tail(covar_4000_final)

covar_both <- cbind(covar_4000_baseline, covar_4000_final)
dim(covar_both)
head(covar_both)
tail(covar_both)

# Read genotype files:
genotype_baseline <- 'genotype_data_4000_baseline_dummyline.tsv'
genotype_4000_baseline <- fread(genotype_baseline, sep = 'auto', header = TRUE, stringsAsFactors = FALSE)
dim(genotype_4000_baseline)
genotype_4000_baseline
setkey(genotype_4000_baseline, samples)

genotype_final <- 'genotype_data_4000_final_dummyline.tsv'
genotype_4000_final <- fread(genotype_final, sep = 'auto', header = TRUE, stringsAsFactors = FALSE)
dim(genotype_4000_final)
genotype_4000_final
setkey(genotype_4000_final, samples)
tables()

# Join genotype files:
genotype_4000_both <- merge(genotype_4000_baseline, genotype_4000_final)
dim(genotype_4000_both)
dim(genotype_4000_baseline)
dim(genotype_4000_final)

# Read expression files:
expression_baseline <- 'subset_baseline_4000_dummyline.tsv'
expression_4000_baseline <- fread(expression_baseline, sep = 'auto', header = TRUE, stringsAsFactors = FALSE)
dim(expression_4000_baseline)
expression_4000_baseline
setkey(expression_4000_baseline, samples)

expression_final <- 'subset_finalVisit_4000_dummyline.tsv'
expression_4000_final <- fread(expression_final, sep = 'auto', header = TRUE, stringsAsFactors = FALSE)
dim(expression_4000_final)
expression_4000_final
setkey(expression_4000_final, samples)
tables()

# Join genotype files:
expression_4000_both <- merge(expression_4000_baseline, expression_4000_final)
dim(expression_4000_both)
dim(expression_4000_baseline)
dim(expression_4000_final)

# Check order of samples is the same across covar, expression and genotype:
identical(colnames(genotype_4000_both)[-1], colnames(expression_4000_both)[-1])
identical(colnames(genotype_4000_both)[-1], colnames(covar_both))

# Write to file:
write.table(genotype_4000_both, 'genotype_4000_both_for_cross.tsv', sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(expression_4000_both, 'expression_4000_both_for_cross.tsv', sep='\t', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(covar_both, 'covar_4000_both_for_cross.tsv', sep='\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
#############################################

#############################################
# The end:
# Remove objects that are not necessary to save:
ls()
object_sizes <- sapply(ls(), function(x) object.size(get(x)))
as.matrix(rev(sort(object_sizes))[1:10])

# rm(object_sizes)

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for xxx.
#############################