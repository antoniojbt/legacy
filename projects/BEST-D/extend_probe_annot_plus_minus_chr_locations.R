##############
# Antonio Berlanga
# Merge genes/probes to get chr locations and plus/minus extensions
# 09 May 2016
# Add annotations using illuminaHumanv4.db package
# TO DO: clean up as got hacked because of missing locations in original illumina download
##############

####################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D_03_MAR.DIR/GAT_backgrounds/')

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
R_session_saved_image <- paste('R_session_saved_image_gene_plus_minus_chr_locations','.RData', sep='')
R_session_saved_image
####################

##############
# Load libraries:
library(illuminaHumanv4.db)
##############

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Illumina probes with chr locations:
illumina_file <- as.character(args[1])
# illumina_file <- 'illumina_probes_genomic_locations_annot.txt'
# illumina_file <- 'illumina_probes_genomic_locations.txt'

hits_file <- as.character(args[2])
# hits_file <- 'full_topTable_pairing_all_treated_probeID_FDR10.txt'
# hits_file <- 'all_treated_baseline_FDR0.05_all_columns.txt'
# hits_file <- 'all_treated_final_FDR0.05_all_columns.txt'
# hits_file <- 'cis_tx_fdr5_reQTLs_annot.txt'
hits_file <- 'full_topTable_pairing_all_treated.txt'

# TO DO: specify column where illumina probe ids are.
####################

####################
# Read data:
illumina_data <- read.csv(illumina_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(illumina_data)
str(illumina_data)
head(illumina_data)
dim(illumina_data)

hits_data <- read.csv(hits_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
class(hits_data)
str(hits_data)
head(hits_data)
dim(hits_data)

illumina_data[1:10, 1]
hits_data[1:10, 1]

# Rename header to allow merge:
colnames(hits_data)[1] <- 'probe_IDs'
colnames(hits_data)

# Basic checks:
length(unique(hits_data$probe_IDs))
which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs))
length(which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs)))
####################

####################
# Get probe IDs:
rownames(hits_data) <- hits_data$probe_IDs
head(hits_data)
illumina_ids  <-  as.character(rownames(hits_data))
class(illumina_ids)
illumina_ids

# Get Entrez IDs:
probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
probes_by_ENTREZID_RE  <- as.data.frame(probes_by_ENTREZID_RE)
# Get gene symbols:
probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
probes_by_symbol <- as.data.frame(probes_by_symbol)
# Get chr:
probes_by_chr <- unlist(mget(illumina_ids, illuminaHumanv4CHR, ifnotfound = NA))
probes_by_chr <- as.data.frame(probes_by_chr)
# Get chr locations:
probes_by_chrloc <- unlist(mget(illumina_ids, illuminaHumanv4CHRLOC, ifnotfound = NA))
probes_by_chrloc <- as.data.frame(probes_by_chrloc)

# Set IDs as columns for merge:
probes_by_symbol$probe_IDs <- row.names(probes_by_symbol)
probes_by_ENTREZID_RE$probe_IDs <- row.names(probes_by_ENTREZID_RE)
probes_by_chr$probe_IDs <- row.names(probes_by_chr)
probes_by_chrloc$probe_IDs <- row.names(probes_by_chrloc)

# Merge results and annotations:
hits_annot <- merge(hits_data, probes_by_symbol)
hits_annot <- merge(hits_annot, probes_by_ENTREZID_RE)
hits_annot <- merge(hits_annot, probes_by_chr)
hits_annot <- merge(hits_annot, probes_by_chrloc)
head(hits_annot)
tail(hits_annot)
dim(hits_annot)


# Rename columns:
colnames(hits_annot)[2] <- 'chr'
colnames(hits_annot)[3] <- 'chr_start'
colnames(hits_annot)[4] <- 'chr_end'
colnames(hits_annot)
####################


####################
# Add plus and minus 100 kb:
hits_annot$chr_start_minus_100kb <- hits_annot$chr_start - 100000
hits_annot$chr_end_plus_100kb <- hits_annot$chr_start + 100000

# Set negatives to 0:
length(which(hits_annot$chr_start_minus_100kb < 0))
hits_annot$chr_start_minus_100kb <- ifelse(hits_annot$chr_start_minus_100kb < 0, 
                                           0, hits_annot$chr_start_minus_100kb)
length(which(hits_annot$chr_start_minus_100kb < 0))
head(hits_annot)

# TO DO: test if chr ends are longer than actual chr.
####################


####################
# Merge results and annotations:
hits_data_annot <- merge(hits_data, illumina_data)
head(hits_data_annot)
tail(hits_data_annot)
dim(hits_data_annot)

# Rename columns:
colnames(illumina_annot)[2] <- 'chr'
colnames(illumina_annot)[3] <- 'chr_start'
colnames(illumina_annot)[4] <- 'chr_end'
colnames(illumina_annot)
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
