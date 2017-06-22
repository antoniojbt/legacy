########################
# Script for processing counts from all files after running eQTL_counting.R
# First run eQTL_counting.R and bash_eQTL_counting.sh which will generate one (or many files)
# Output from above will be counts_* files
# Antonio J Berlanga-Taylor
# 1 June 2016
########################


#############################################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
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
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_eQTL_counts_post_processing','.RData', sep='')
R_session_saved_image
####################


########################
# Load packages:
library(data.table)
########################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# eg
# eQTL_counts <- as.character(args[1])
# eQTL_file <- '2000+4000-baseline-1.eQTL_cis'
####################

########################
# Read data:
# Concatenate files:
filename <- paste('eQTL_counts_concatenated_', Sys.Date(), '.txt', sep = '')
cmd_to_run <- sprintf('paste counts_* | column -x > %s', filename)
cmd_to_run
system(cmd_to_run)
# system('head eQTL_counts_concatenated.txt')
# system('head eQTL_counts_concatenated.txt | cut -f1,2')
# system('head eQTL_counts_concatenated.txt | cut -f2')

# File name to read
eQTL_counts_file <- filename

# Read into R:
eQTL_data <- fread(eQTL_counts_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
dim(eQTL_data)
eQTL_data
setkey(eQTL_data)
colnames(eQTL_data)
# View(eQTL_data)

# Every other column contains the headers (odd numbered columns:
# View(eQTL_data)
# system('cat eQTL_counts_paste_column.txt  | cut -f1,3')
# system('cat eQTL_counts_paste_column.txt  | cut -f2,4')

# Delete 3,5,7, etc.:
# Count files to concatenate:
files_to_concat <- as.numeric(system('ls -1 counts_* | wc -l', intern = TRUE))
files_to_concat

# Get odd numbers and delete columns (to delete)except first):
odd_cols <- seq(3, by = 2, length.out = files_to_concat-1)
odd_cols
eQTL_data <- eQTL_data[, !c(odd_cols), with = F]

# ########################
# # Write to file:
write.table(eQTL_data, eQTL_counts_file,  sep = '\t', na = 'NA', quote = FALSE,
            col.names = TRUE, row.names = FALSE)
########################

########################
# Transpose and process if further processing is needed:
# # Transpose file and create headers for easier handling:
# # View(t(eQTL_data))
# eQTL_data <- t(eQTL_data)
# # View(eQTL_data)
# str(eQTL_data[, 1:5])
# eQTL_data[1:5, 1:5]
# colnames(eQTL_data)
# colnames(eQTL_data) <- eQTL_data[1, ]
# eQTL_data <- eQTL_data[-1, ]
# # Convert back to data.frame:
# eQTL_data <- as.data.frame(eQTL_data)
# # # Use row names as first column:
# # eQTL_data$file <- rownames(eQTL_data)
# # rownames(eQTL_data)
# # eQTL_data$file
# # colnames(eQTL_data)
# # eQTL_data[, 1:5]
# # View(eQTL_data)
# ########################
# 
# ########################
# # Process results:
# colnames(eQTL_data)
# rownames(eQTL_data)
# 
# ########################
# 
# ########################
# # Write to file:
# write.table(eQTL_data, eQTL_counts_file,  sep = '\t', na = 'NA', quote = FALSE, 
#             col.names = NA, row.names = TRUE)
########################

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

# Next: run XGR, IPA, etc, cross with GWAS, ENCODE, etc.
####################