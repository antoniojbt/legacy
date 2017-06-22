##############
# Antonio Berlanga
# Merge genes/probes to get chr locations and plus/minus extensions
# 10 May 2016
# Illumina probe locations are incomplete, run more tests and get locations from other source
##############

##############
# illumina location file is incomplete, ~12800 QC'd BESTD probes match out of ~16700
# so probably many untested eQTLs
# Use biomaRt See: https://support.bioconductor.org/p/75923/
# Illumina download page is:
# http://support.illumina.com/array/array_kits/humanht-12_v4_expression_beadchip_kit/downloads.html
# Missing about one third of genomic locations from differentially expressed genes in BESTD
# Also see another tool:
# http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0139516
# Reasons? Best to leave excluded if low quality but not filtered by illuminaHumanv4.db? Probe annotations are good or
# perfect according to illuminaHumanv4.db though.
# Can also get from illuminaHumanv4.db but erroring
##############

####################
##Set working directory and file locations and names of required inputs:
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
# load('R_session_saved_image_illumina_probes_chr_mismatch_check.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_illumina_probes_chr_mismatch_check','.RData', sep='')
R_session_saved_image
##############

##############
# Load libraries:
library(biomaRt)
library(illuminaHumanv4.db)
library(plyr)
##############

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Illumina probes with chr locations:
illumina_file <- as.character(args[1])
# illumina_file <- 'illumina_probes_genomic_locations_annot.txt'

hits_file <- as.character(args[2])
# hits_file <- 'full_topTable_pairing_all_treated.txt'

illumina_original_file <- 'HumanHT-12_V4_0_R2_15002873_B.txt'
####################


####################
# Read data:
illumina_data <- read.csv(illumina_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(illumina_data)
str(illumina_data)
head(illumina_data)
dim(illumina_data)

hits_data <- read.csv(hits_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
class(hits_data)
str(hits_data)
head(hits_data)
dim(hits_data)

illumina_original <- read.csv(illumina_original_file, sep = '\t', skip = 8, header = TRUE, 
                              stringsAsFactors = FALSE, strip.white = TRUE, na.strings = c('', ' ', 'NA'))
class(illumina_original)
str(illumina_original)
head(illumina_original)
tail(illumina_original)
dim(illumina_original)

illumina_data[1:10, 1]
hits_data[1:10, 2]

# Rename headers:
colnames(hits_data)[2] <- 'probe_IDs'
colnames(hits_data)
####################

####################
# Explore Illumina original file:
colnames(illumina_original)
# View(illumina_original)
length(which(is.na(illumina_original$Chromosome)))
length(which(is.na(illumina_original$ILMN_Gene)))
length(which(is.na(illumina_original$Probe_Coordinates)))
# TO DO: check what processing I did to the original file. For now look at other annotation databases.
####################


####################
# Basic checks:
length(unique(hits_data$probe_IDs))
which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs))
length(which(as.character(hits_data$probe_IDs) %in% as.character(illumina_data$probe_IDs)))
dim(illumina_data)
dim(hits_data)
####################

####################
# Sanity check on QC filtering script, test whether QC'd and normalised probes are of good quality:
# Get probe IDs:
rownames(hits_data) <- hits_data$probe_IDs
head(hits_data)
illumina_ids  <-  as.character(rownames(hits_data))
class(illumina_ids)
illumina_ids

# Get quality annotations:
probes_by_quality <- unlist(mget(illumina_ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
probes_by_quality <- as.data.frame(probes_by_quality)
head(probes_by_quality)
dim(probes_by_quality)
summary(probes_by_quality)
# All probes are of good or perfect quality.
# TO DO: contact Illumina, try to get re-annotations/locations from other pipeline.
####################

####################
# Run biomaRt to get missing locations:
# See: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.pdf

# Check what's available:
listMarts()
# Set the source (database) and dataset to use:
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# Check what information is available
filters <- listFilters(ensembl)
head(filters)
# View(listAttributes(ensembl))
# Search for keyword:
grep('illumina', listAttributes(ensembl), ignore.case = TRUE)

# Obtain values from the database and dataset selected:
# illumina_ids (obtained from QC'd probes in hits_file) and mart were already set:
biomart_attributes <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id', 'entrezgene',
                                         'chromosome_name', 'start_position', 'end_position'), 
                            filters = 'illumina_humanht_12_v4', values = illumina_ids, mart = ensembl)

head(biomart_attributes)
tail(biomart_attributes)
dim(biomart_attributes)
length(illumina_ids)
dim(hits_data)
dim(illumina_data)
####################

####################
# Sanity check: compare to illuminaHumanv4.db package and Illumina's own file:
head(illumina_data)
head(biomart_attributes)

# Get the subset of one file present in the illumina file for direct comparison:
length(which(biomart_attributes$illumina_humanht_12_v4 %in% illumina_data$probe_IDs))
which(biomart_attributes$illumina_humanht_12_v4 %in% illumina_data$probe_IDs)

colnames(biomart_attributes)[1] <- 'probe_IDs'
merged_sets <- merge(biomart_attributes, illumina_data)
dim(merged_sets)
head(merged_sets)

# Check overlap of annotations:
# Percent of biomart gene symbols that match illuminaHumanv4.db:
length(which(as.character(merged_sets$hgnc_symbol) %in% as.character(merged_sets$probes_by_symbol))) / length(merged_sets$probe_IDs)
length(which(is.na(merged_sets$hgnc_symbol)))
length(which(is.na(merged_sets$probes_by_symbol)))
length(which(merged_sets$probes_by_symbol == '<NA>'))
# 86% match, may be due to NAs or ...?

# Percent of biomart Entrez IDs that match illuminaHumanv4.db:
length(which(as.character(merged_sets$entrezgene) %in% as.character(merged_sets$probes_by_ENTREZID_RE))) / length(merged_sets$probe_IDs)

# 95% match

# Percent of biomart chr names that match illuminaHumanv4.db:
# Split illuminaHumanv4.db as reports as 'chr#' while biomaRt as '#' only:
merged_sets$chr_illuminaHumanv4 <- sapply(strsplit(as.character(merged_sets$chr), split = 'chr'), "[[", 2)
head(merged_sets)
count(merged_sets$chr_illuminaHumanv4)
count(merged_sets$chr)
count(merged_sets$chromosome_name)

length(which(as.character(merged_sets$chromosome_name) %in% as.character(merged_sets$chr_illuminaHumanv4))) / length(merged_sets$probe_IDs)

# All match
# Conclusion: use Biomart's output for annotations, run cis eQTLs based on this.
####################


####################
# Return to use Biomart's output:
head(biomart_attributes)

# Add plus and minus 100 kb:
biomart_attributes$chr_start_minus_100kb_bio <- biomart_attributes$start_position - 100000
biomart_attributes$chr_end_plus_100kb_bio <- biomart_attributes$end_position + 100000

# Set negatives to 0:
length(which(biomart_attributes$chr_start_minus_100kb_bio < 0))
biomart_attributes$chr_start_minus_100kb_bio <- ifelse(biomart_attributes$chr_start_minus_100kb_bio < 0, 
                                           0, biomart_attributes$chr_start_minus_100kb_bio)
length(which(biomart_attributes$chr_start_minus_100kb_bio < 0))
head(biomart_attributes)

# TO DO: test if chr ends are longer than actual chr.
####################


####################
# Write results to disk:
write.table(biomart_attributes, 'biomart_QCd_probes_genomic_locations_annot.txt', row.names = FALSE, quote = FALSE, sep = '\t', col.names = TRUE)

# Generate file for MatrixeQTL:
system('cat biomart_QCd_probes_genomic_locations_annot.txt | cut -f1,5,6,7 | awk -v OFS='\t' '$2="chr"$2' - | sed '1d' - > biomart_QCd_probes_genomic_locations_annot_MatrixeQTL.txt')
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
