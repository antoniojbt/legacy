#############################
# Small script to explore illuminaHumanv4.db package for annotation and probe filtering
# Antonio J Berlanga-Taylor
# 29 June 2016
# Requires a (normalised) expression set as input.
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/probe_filtering_check/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_probe_filtering_sanity_check",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Read results from 01_microarrayxxx file (saved as RData object):
# Load a previous R session, data and objects:
#load('R_session_saved_image_read_and_QC.RData', verbose=T)
#load('R_session_saved_image_normalisation_full.RData', verbose=T)
load('R_session_saved_image_normalisation.RData', verbose=T)
# load('R_session_saved_image_normalisation_full_1ry_cells.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_probe_filtering_sanity_check', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

library(limma)
library(illuminaHumanv4.db)
library(plyr)
library(data.table)
#############################


#############################
# Read in Illumina annotation table (data.table errored):
illumina_annotation <- read.csv('HumanHT-12_V4_0_R2_15002873_B.txt', skip = 8, sep = '\t', header = TRUE,
                             stringsAsFactors = TRUE)
head(illumina_annotation)
dim(illumina_annotation)
#############################

#############################
# Once background correction, normalisation and transformation have been performed run probe filtering (ie to exclude multi-mapping probes, probes over
# SNPs, etc.) carry out the following:

## d) Filtering and creating an eset object

# Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays 
# according to a detection p-value of 5% (Limma vignette case study, p. 108):

normalised
dim(normalised)
class(normalised)
summary(normalised)
range(normalised$E)


normalised$other$'Detection Pval'[1:5,1:5]
#normalised$other$Detection[1:5,1:5] # I think these are the inverse of p-values? ie for GAinS data from ArrayExpress
length(normalised$other$'Detection Pval' < 0.05)
length(which(normalised$other$'Detection Pval' < 0.05) == TRUE)
#length(which(normalised$other$Detection < 0.05) == TRUE)
range(normalised$other$'Detection Pval')
dim(normalised$other$'Detection Pval')

#Set filters:
expressed <- which((rowSums(normalised$other$'Detection Pval' < 0.05) >= 3) == TRUE)
#expressed <- which((rowSums(normalised$other$Detection < 0.05) >= 10) == TRUE)
head(expressed, n=10)
length(expressed)

expressed_2 <- rowSums(normalised$other$'Detection Pval' < 0.05) >=3
#expressed_2 <- rowSums(normalised$other$Detection < 0.05) >=10
head(expressed_2, n=10)
head(which(expressed_2 == TRUE), n=10)
length(which(expressed_2 == TRUE))


#Extract data and create an ExpressionSet (eset, or EList object) at the same time (necessary for linear modelling steps):
normalised_expressed <- normalised[expressed,]
#normalised_expressed <- normalised[-expressed,] # I think Array Express data 'Detection' column is the inverse of the p-value

#Explore contents of file:
head(normalised_expressed$E)
dim(normalised_expressed)
class(normalised_expressed)
summary(normalised_expressed)
head(normalised_expressed$source)
head(normalised_expressed$E)[1:5,1:5]
range(normalised_expressed$E)
head(normalised_expressed$genes)
summary(normalised_expressed$other)
head(normalised_expressed$other$'Detection Pval')[1:5,1:5]
range(normalised_expressed$other$'Detection Pval')
normalised_expressed$targets


###################
## Filter based on annotation
# Check http://www.bioconductor.org/packages/release/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
# Could also use Bioconductor's Annotation workflows:
# See http://www.bioconductor.org/help/workflows/annotation/annotation/
# Also download Illumina's annotation file:
# https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanht-12/humanht-12_v4_0_r2_15002873_b.txt.zip
# Merge with file above and exclude all unplaced and non-autosomal probes for example.
# TO DO: The R package illuminaHumanv4.db is known to have issues though:


# Remove probes that aren't annotated with a gene symbol:
annotated_probes <- !is.na(normalised_expressed$genes$SYMBOL)
count(annotated_probes)
head(annotated_probes)
tail(annotated_probes)

un_annotated_probes <- which(annotated_probes == FALSE)
length(un_annotated_probes)
head(un_annotated_probes)
tail(un_annotated_probes)

# Remove un-annotated probes:
normalised_expressed_annotated <- normalised_expressed[annotated_probes,]
#normalised_expressed_annotated <- normalised_expressed # for array data without annotation at this point (eg GAinS)

#Explore contents of file and sanity checks:
head(normalised_expressed_annotated$genes$SYMBOL)
dim(normalised_expressed_annotated)
dim(normalised_expressed)
length(annotated_probes)
count(annotated_probes)

#View(normalised_expressed$genes)
###################


###################
# Remove probes based on quality annotations:
# See: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
# for an example from the BeadArray use cases, p. 27.
illumina_ids  <-  as.character(rownames(normalised_expressed_annotated))
head(illumina_ids)
length(illumina_ids)

probes_by_qual  <- unlist(mget(illumina_ids, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
head(probes_by_qual)
summary(probes_by_qual)
count(probes_by_qual)

# Plot expression signal by probe quality:
aveSignal_probes_by_qual <- rowMeans(normalised_expressed_annotated$E)
png('Average_signal_probes_by_quality.png', width = 4, height = 4, units = 'in', res = 300)
boxplot(aveSignal_probes_by_qual ~ probes_by_qual)
dev.off()

###########
# Sanity check of annotation package
# TO DO: clean up and move to a separate script
# Confirm low quality probes match to multiple locations and should be removed:
all_probes  <- unlist(mget(illumina_ids, illuminaHumanv4REPEATMASK, ifnotfound = NA))
#View(all_probes)
summary(all_probes)
head(all_probes)
length(all_probes)
all_probes_2 <- !is.na(all_probes) # Those marked TRUE are repeats, FALSE are not
#View(all_probes_2)
summary(all_probes_2)

# See which probes have high expression and are of 'Bad' quality:
query_bad_IDs <- names(which(probes_by_qual == "Bad" & aveSignal_probes_by_qual  > 12))
head(query_bad_IDs)
length(query_bad_IDs)

# Identify which of the highly expressed 'Bad' probes are repeats:
repeat_mask_bad <- unlist(mget(query_bad_IDs, illuminaHumanv4REPEATMASK, ifnotfound = NA))
#View(repeat_mask_bad)
head(repeat_mask_bad)
tail(repeat_mask_bad)
length(repeat_mask_bad)

# Counts for repeats:
repeat_mask_bad_2 <- !is.na(repeat_mask_bad) # Those marked TRUE are repeats, FALSE are not
#View(repeat_mask_bad_2)
summary(repeat_mask_bad_2)
head(repeat_mask_bad_2)
length(repeat_mask_bad_2)
# Many classed as Bad are repeats.

# Check which probes multi-map:
multiple_matches_bad <- unlist(mget(query_bad_IDs, illuminaHumanv4SECONDMATCHES, ifnotfound = NA))
head(multiple_matches_bad)
length(multiple_matches_bad)
class(multiple_matches_bad)

# Convert to data frame to check single vs secondary matches:
multiple_matches_bad_df <- as.data.frame(multiple_matches_bad)
multiple_matches_bad_df['probe'] <- row.names(multiple_matches_bad_df)
head(multiple_matches_bad_df)
# View(multiple_matches_bad_df)
dim(multiple_matches_bad_df)
multiple_matches_bad_df_IDs <- row.names(multiple_matches_bad_df)
multiple_matches_bad_df_IDs

# Many classed as Bad multi-map.

# multiple_matches_bad_df <- strsplit(as.character(multiple_matches_bad_df[, 1]), split = 'chr')
# head(multiple_matches_bad_df)
# length(multiple_matches_bad_df)
# multiple_matches_bad_df <- do.call(cbind, multiple_matches_bad_df)
# colnames(multiple_matches_bad_df) <- multiple_matches_bad_df_IDs
# multiple_matches_bad_df <- as.data.frame(multiple_matches_bad_df)
# multiple_matches_bad_df[ multiple_matches_bad_df == "" ] <- NA
# dim(multiple_matches_bad_df)
# head(multiple_matches_bad_df)
# # Count number of mapping locations, each TRUE is one location:
# multiple_matches_bad_df_t <- t(!is.na(multiple_matches_bad_df))
# head(multiple_matches_bad_df_t)
# summary(multiple_matches_bad_df_t)
# count(multiple_matches_bad_df_t)
# count(multiple_matches_bad_df[, 2])
# 
# multiple_matches_bad_2 <- !is.na(multiple_matches_bad) # Those marked TRUE have at least one match
# #View(repeat_mask_bad_2)
# summary(multiple_matches_bad_2)
# head(multiple_matches_bad_2)
# length(multiple_matches_bad_2)

# mget("ILMN_1692145", illuminaHumanv4PROBESEQUENCE)

# Sanity check for 'Perfect' IDs with high expression:
query_perfect_IDs <- names(which(probes_by_qual == "Perfect" & aveSignal_probes_by_qual > 12))
head(query_perfect_IDs)
length(query_perfect_IDs)

repeat_mask_perfect <- unlist(mget(query_perfect_IDs, illuminaHumanv4REPEATMASK, ifnotfound = NA))
head(repeat_mask_perfect)
length(repeat_mask_perfect)
count(repeat_mask_perfect)
# None of the Perfect probes are repeats

# Check which Perfect probes multi-map:
multiple_matches_perfect <- unlist(mget(query_perfect_IDs, illuminaHumanv4SECONDMATCHES, ifnotfound = NA))
head(multiple_matches_perfect)
multiple_matches_perfect_df <- as.data.frame(multiple_matches_perfect)
# View(multiple_matches_perfect_df)
length(multiple_matches_perfect)
# Many have secondary matches.
###########

###########
# Check if VD gene probes are filtered:
# quality, chromosome, multi-mapping, overlapping SNP, symbol, entrezID

# Genes to check:
gene_list <- c('VDR', 
               'GC',
               'CYP2R1',
               'VDBP',
               'CYP27B1',
               'CYP27A1',
               'DHCR7',
               'NADSYN1',
               'CYP24A1')
which(illumina_annotation$Symbol == as.character('VDBP'))
which(illumina_annotation$Symbol == as.character('VDR'))

# Function to get illumina probe ID for each:
get_Probe_id <- function(gene) {
  location <- which(illumina_annotation$Symbol == as.character(gene))
  if (length(location) > 0) {
    info <- illumina_annotation[c(which(illumina_annotation$Symbol == as.character(gene))), 'Probe_Id']
    info <- as.data.frame(as.character(info))
    colnames(info)[1] <- 'Probe_ID'
    info[, 1] <- as.character(info[, 1])
    info['gene'] <- as.character(gene)
    return(info)
    } else {
      info <- data.frame('NA', as.character(gene))
      colnames(info)[1] <- 'Probe_ID'
      colnames(info)[2] <- 'gene'
      return(info)
    }
  }

# Get each and create a dataframe:
vdr_id <- get_Probe_id('VDR')
vdr_id
str(vdr_id)
vdbp <- get_Probe_id('VDBP')
vdbp

# Run loop for gene list:
# Empty data frame:
gene_list_probe_IDs <- data.frame()
gene_list_probe_IDs

# Loop through genes of interest:
for (i in gene_list) {
  probeID <- get_Probe_id(i)
  gene_list_probe_IDs <- rbind(probeID, gene_list_probe_IDs)
}

# Symbols with probe IDs:
head(gene_list_probe_IDs)

# Function to convert to data.frame and merge annotation:
# df_and_merge <- function(df1, df2) {
#   df1 <- as.data.frame(df1)
#   df1['Probe_ID'] <- row.names(df1)
#   # Merge:
#   df_merged <- merge(df2, df1)
#   return(df_merged)
#   }

# Annotate and merge:
probes_by_qual <- as.data.frame(probes_by_qual)
probes_by_qual['Probe_ID'] <- row.names(probes_by_qual)
# Merge:
gene_list_probe_IDs <- merge(gene_list_probe_IDs, probes_by_qual)
gene_list_probe_IDs


# Remove probes that overlap SNPs:
# ?illuminaHumanv4OVERLAPPINGSNP
probes_by_SNPs  <- !is.na(unlist(mget(as.character(gene_list_probe_IDs[, 'Probe_ID']),
                                      illuminaHumanv4OVERLAPPINGSNP, ifnotfound = NA)))

# Annotate and merge:
probes_by_SNPs <- as.data.frame(probes_by_SNPs)
probes_by_SNPs['Probe_ID'] <- row.names(probes_by_SNPs)
# Merge:
gene_list_probe_IDs <- merge(gene_list_probe_IDs, probes_by_SNPs)
gene_list_probe_IDs

# Keep only probes that have Entrez IDs:
probes_by_ENTREZID  <- unlist(mget(as.character(gene_list_probe_IDs[, 'Probe_ID']),
                                   illuminaHumanv4ENTREZID, ifnotfound = NA))

# Annotate and merge:
probes_by_ENTREZID <- as.data.frame(probes_by_ENTREZID)
probes_by_ENTREZID['Probe_ID'] <- row.names(probes_by_ENTREZID)
# Merge:
gene_list_probe_IDs <- merge(gene_list_probe_IDs, probes_by_ENTREZID)
gene_list_probe_IDs


# Remove probes from X and Y chromosomes and unmapped or with missing values:
probes_by_CHR  <- na.omit(unlist(mget(as.character(gene_list_probe_IDs[, 'Probe_ID']),
                                      illuminaHumanv4CHR, ifnotfound = NA)))

# Annotate and merge:
probes_by_CHR <- as.data.frame(probes_by_CHR)
probes_by_CHR['Probe_ID'] <- row.names(probes_by_CHR)
# Merge:
gene_list_probe_IDs <- merge(gene_list_probe_IDs, probes_by_CHR)
gene_list_probe_IDs


# TO DO: remove multi-mapping (these should have been marked as Bad already)

###########

###########
write.table(gene_list_probe_IDs, 'gene_list_probe_IDs_illuminaHuman_annotation.txt', sep='\t', 
            quote=FALSE, na='NA', col.names=NA, row.names=TRUE)
###########

#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image_full, compress='gzip')

# objects_to_save <- (c('normalised_filtered_annotated', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Print session information:
sessionInfo()

q()
#############################