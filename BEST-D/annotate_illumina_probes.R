##############
# Antonio Berlanga
# Annotate Illumina probes
# 09 May 2016
##############

####################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D_03_MAR.DIR/eqtl_files.dir/eqtl_data.dir/')

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
# load('R_session_saved_image_gene_list_overlap.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_annotate_illumina_probes','.RData', sep='')
R_session_saved_image
####################

##############
# Load libraries:
library(illuminaHumanv4.db)
library(data.table)
# Get additional functions needed:

##############

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Illumina probes with chr locations:
illumina_file <- as.character(args[1])
# illumina_file <- 'illumina_probes_genomic_locations.txt'

# Read data:
# illumina_data <- fread(illumina_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
illumina_data <- read.csv(illumina_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(illumina_data)
str(illumina_data)
head(illumina_data)
dim(illumina_data)

# Get probe IDs:
# illumina_ids  <-  illumina_data[, 'probe_IDs', with =F]
rownames(illumina_data) <- illumina_data$probe_IDs
head(illumina_data)
illumina_ids  <-  as.character(rownames(illumina_data))
class(illumina_ids)
illumina_ids
# Get Entrez IDs:
probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
probes_by_ENTREZID_RE  <- as.data.frame(probes_by_ENTREZID_RE)
# Get gene symbols:
probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
probes_by_symbol <- as.data.frame(probes_by_symbol)
# Set IDs as columns for merge:
probes_by_symbol$probe_IDs <- row.names(probes_by_symbol)
probes_by_ENTREZID_RE$probe_IDs <- row.names(probes_by_ENTREZID_RE)
# Merge results and annotations:
illumina_annot <- merge(illumina_data, probes_by_symbol)
illumina_annot <- merge(illumina_annot, probes_by_ENTREZID_RE)
head(illumina_annot)
tail(illumina_annot)

# Rename columns:
colnames(illumina_annot)[2] <- 'chr'
colnames(illumina_annot)[3] <- 'chr_start'
colnames(illumina_annot)[4] <- 'chr_end'
colnames(illumina_annot)
####################


####################
# Add plus and minus 100 kb:
illumina_annot$chr_start_minus_100kb <- illumina_annot$chr_start - 100000
illumina_annot$chr_end_plus_100kb <- illumina_annot$chr_start + 100000

# Set negatives to 0:
length(which(illumina_annot$chr_start_minus_100kb < 0))
illumina_annot$chr_start_minus_100kb <- ifelse(illumina_annot$chr_start_minus_100kb < 0, 
                                               0, illumina_annot$chr_start_minus_100kb)
length(which(illumina_annot$chr_start_minus_100kb < 0))
head(illumina_annot)

# TO DO: test if chr ends are longer than actual chr.
####################


####################
# Write results to disk:
write.table(illumina_annot, 'illumina_probes_genomic_locations_annot.txt', row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)
####################

####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run the script for xxx.
####################
