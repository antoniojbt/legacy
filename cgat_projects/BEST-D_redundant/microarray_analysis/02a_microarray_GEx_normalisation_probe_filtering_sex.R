#############################
# To be run after 02 normalisation and filtering script
# Antonio J Berlanga-Taylor

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_probe_filtering",".txt", sep=""), open='a')
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
load('R_session_saved_image_normalisation_full.RData', verbose=T)
# load('R_session_saved_image_normalisation_full_1ry_cells.RData', verbose=T)
# load('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/R_session_saved_image_normalisation_full.RData')
# load('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/R_session_saved_image_probe_filtering.RData')
# load('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/R_session_saved_image_read_and_QC.RData')

# Filename to save current R session, data and objects at the end:
# R_session_saved_image <- paste('R_session_saved_image_probe_filtering', '.RData', sep='')
R_session_saved_image <- paste('R_session_saved_image_probe_filtering_full', '.RData', sep='')
# R_session_saved_image <- paste('R_session_saved_image_probe_filtering_full_noSNP_filter', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:

library(limma)
library(illuminaHumanv4.db)
library(plyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(data.table)
#############################


#############################

# Once background correction, normalisation and transformation have been performed run probe filtering (ie to exclude multi-mapping probes, probes over
# SNPs, etc.) carry out the following:

## d) Filtering and creating an eset object

# Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays 
# according to a detection p-value of 5% (Limma vignette case study, p. 108):

normalised
dim(normalised)
head(normalised$E)
dim(normalised$E)
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
normalised_expressed <- normalised[expressed, ]
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
# TO DO: The R package illuminaHumanv4.db is known to have issues though, complete 'microarray_probe_filtering_sanity_check.R' script.

# TO DO: exlcude this as not necessary at this stage:
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
normalised_expressed_annotated <- normalised_expressed[annotated_probes, ]
#normalised_expressed_annotated <- normalised_expressed # for array data without annotation at this point (eg GAinS)

#Explore contents of file and sanity checks:
head(normalised_expressed_annotated$genes$SYMBOL)
dim(normalised_expressed_annotated)
dim(normalised_expressed)
length(annotated_probes)
count(annotated_probes)
length(which(is.na(normalised_expressed_annotated$genes$SYMBOL)))
length(which(is.na(normalised_expressed$genes$SYMBOL)))

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

# Remove probes of low quality:
remove_bad_probes  <- probes_by_qual == "No match" | probes_by_qual == "Bad"
head(probes_by_qual)
head(remove_bad_probes)
count(remove_bad_probes)
normalised_expressed_annotated_qual <- normalised_expressed_annotated[!remove_bad_probes, ]
dim(normalised_expressed_annotated)
dim(normalised_expressed_annotated_qual)
dim(normalised_expressed_annotated)[1] - dim(normalised_expressed_annotated_qual)[1]
head(normalised_expressed_annotated_qual$E)


# Check counts match:
count(!is.na(normalised_expressed_annotated_qual$genes$SYMBOL))
length(rownames(normalised_expressed_annotated_qual))
###################

###################
# Remove probes that overlap SNPs:
# ?illuminaHumanv4OVERLAPPINGSNP
# Skip if running without the SNP filter (noSNP):
probes_by_SNPs  <- !is.na(unlist(mget(as.character(rownames(normalised_expressed_annotated_qual)), 
                                      illuminaHumanv4OVERLAPPINGSNP, ifnotfound = NA)))
head(probes_by_SNPs)
summary(probes_by_SNPs)
#View(probes_by_SNPs)

normalised_expressed_annotated_qual_noSNPs <- normalised_expressed_annotated_qual[!probes_by_SNPs, ]
dim(normalised_expressed_annotated_qual)
dim(normalised_expressed_annotated_qual_noSNPs)
dim(normalised_expressed_annotated_qual)[1] - dim(normalised_expressed_annotated_qual_noSNPs)[1]
###################


###################
# TO DO: exlcude this as not necessary at this stage:
# Keep only probes that have Entrez IDs:
## Get ENTREZ ID mappings:
dim(normalised_expressed_annotated_qual_noSNPs)
probes_by_ENTREZID  <- unlist(mget(as.character(rownames(normalised_expressed_annotated_qual_noSNPs)),
                                   illuminaHumanv4ENTREZID, ifnotfound = NA))
# If not filtering SNPs for testing, use:
# probes_by_ENTREZID  <- unlist(mget(as.character(rownames(normalised_expressed_annotated_qual)),
#                                    illuminaHumanv4ENTREZID, ifnotfound = NA))

head(probes_by_ENTREZID)
summary(probes_by_ENTREZID)
str(probes_by_ENTREZID)
probes_without_ENTREZID  <- is.na(probes_by_ENTREZID)
count(probes_without_ENTREZID)
length(which(probes_without_ENTREZID))
#View(probes_by_ENTREZID)

# Remove probes without Entrez IDs:
normalised_expressed_annotated_qual_noSNPs_noID <- normalised_expressed_annotated_qual_noSNPs[!probes_without_ENTREZID, ]
# If not filtering SNPs for test use (kept names the same downstream for ease though):
# normalised_expressed_annotated_qual_noSNPs_noID <- normalised_expressed_annotated_qual[!probes_without_ENTREZID, ]
dim(normalised_expressed_annotated_qual_noSNPs)
dim(normalised_expressed_annotated_qual_noSNPs_noID)

# TO DO: Check p. 12 illuminaHumanv4listNewMappings:
# https://www.bioconductor.org/packages/3.3/data/annotation/manuals/illuminaHumanv4.db/man/illuminaHumanv4.db.pdf
# Consider using illuminaHumanv4ENTREZREANNOTATED instead of entrezid re-annotated instead of entrez id.
# Other new mappings are in use but double check.
# The below is currently not used:
probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
head(probes_by_ENTREZID_RE)
summary(probes_by_ENTREZID_RE)
str(probes_by_ENTREZID_RE)
count(!is.na(probes_by_ENTREZID_RE))
#View(probes_by_ENTREZID_RE)

###################


###################
# Remove probes from X and Y chromosomes and unmapped or with missing values:
## Get chromosome mappings:
# illuminaHumanv4CHR can show duplicate mappings, see: 
# http://rgm.ogalab.net/RGM/R_rdfile?f=illuminaHumanv4.db/man/illuminaHumanv4CHR.Rd&d=R_BC
dim(normalised_expressed_annotated_qual_noSNPs_noID)
probes_by_CHR  <- na.omit(unlist(mget(as.character(rownames(normalised_expressed_annotated_qual_noSNPs_noID)),
                                   illuminaHumanv4CHR, ifnotfound = NA)))

#View(probes_by_CHR)
head(probes_by_CHR)
tail(probes_by_CHR)
length(probes_by_CHR)
summary(probes_by_CHR)
str(probes_by_CHR)
count(probes_by_CHR)

# Lengths differ between chromosome annotation and normalised filtered file (multi-mapping probes when illuminaHumanv4 was written).
# Probe IDs aren't duplicated though and it seems like each probe has one chromosome (for each row).

# Convert files to data frames and merge:
probes_by_CHR <- as.data.frame(probes_by_CHR)
probes_by_CHR['Probe_Id'] <- NA
probes_by_CHR['Probe_Id'] <- row.names(probes_by_CHR)
head(probes_by_CHR)
length(which(duplicated(probes_by_CHR$Probe_Id)))
length(which(!duplicated(probes_by_CHR$Probe_Id)))

# Convert elist object to dataframe (this loses some information but retains what is
# needed (probe expression levels, probe ID, illumina ID, annotations, etc.):
# Information lost is source, detection P-values (already used, kept in $other), and $targets (
# which stores the file names read originally for read.illm function):
str(normalised$source)
str(normalised$genes)
str(normalised$other)
str(normalised$targets)

df_normalised_expressed_annotated_qual_noSNPs_noID <- as.data.frame(normalised_expressed_annotated_qual_noSNPs_noID)
#View(df_normalised_expressed_annotated_qual_noSNPs_noID)
head(df_normalised_expressed_annotated_qual_noSNPs_noID)
df_normalised_expressed_annotated_qual_noSNPs_noID['Probe_Id'] <- row.names(df_normalised_expressed_annotated_qual_noSNPs_noID)
head(df_normalised_expressed_annotated_qual_noSNPs_noID$Probe_Id)
length(df_normalised_expressed_annotated_qual_noSNPs_noID$Probe_Id)

# Merge both data frames:
df_normalised_expressed_annotated_qual_noSNPs_noID_chr <- merge(df_normalised_expressed_annotated_qual_noSNPs_noID, probes_by_CHR, by = 'Probe_Id')
df_normalised_expressed_annotated_qual_noSNPs_noID_chr[1:5, 1:5]
head(df_normalised_expressed_annotated_qual_noSNPs_noID_chr)
dim(df_normalised_expressed_annotated_qual_noSNPs_noID_chr)
summary(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$probes_by_CHR)


#########
# TO DO: move this to sanity check script 'microarray_probe_filtering_sanity_check.R'
# Check if annotations between Illumina and illuminaHumanv4.db match:
# Illumina annotations here aret hose from the array itself (ie not the downloaded file from Illumina's webpage, that's another
# option to check).
# This is just sanity as issues have been raised in the past regarding different types 
# of annotations (probe quality, gene symbols, etc.).
illumina_chr <- as.list(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$CHROMOSOME)
head(illumina_chr)
length(illumina_chr)

illuminaHumanv4_chr <- as.list(df_normalised_expressed_annotated_qual_noSNPs_noID_chr$probes_by_CHR)
head(illuminaHumanv4_chr)
length(illuminaHumanv4_chr)

identical(illumina_chr, illuminaHumanv4_chr)
# This test gives false (only takes one mismatch though).

# Check what's happening:
chr_comparison <- df_normalised_expressed_annotated_qual_noSNPs_noID_chr[, c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
head(chr_comparison)
tail(chr_comparison)

# Illumina's annotation has many missing values:
length(which(is.na(chr_comparison$CHROMOSOME)))

# Convert blanks to NAs:
chr_comparison$CHROMOSOME[chr_comparison$CHROMOSOME == ''] <- NA
chr_comparison$probes_by_CHR[chr_comparison$probes_by_CHR == ''] <- NA
#View(chr_comparison)
length(which(is.na(chr_comparison$CHROMOSOME)))
length(which(is.na(chr_comparison$probes_by_CHR)))
count(chr_comparison$CHROMOSOME)
count(chr_comparison$probes_by_CHR)

# TO DO: Check these differences (chromosome counts seem close but don't match...).
# Keeping illuminaHumanv4 package annotations for now (no clear answers on what to use and 
# haven't tested):

# Remove probes mapping to sex chromosomes and unmapped:
# Rename gene expression file to make it manageable:
df_normalised_filtered <- df_normalised_expressed_annotated_qual_noSNPs_noID_chr
remove_sex_probes  <- which(df_normalised_filtered$probes_by_CHR == 'X' | df_normalised_filtered$probes_by_CHR == 'Y' | 
                              df_normalised_filtered$probes_by_CHR == 'Un')
head(remove_sex_probes)

head(which(df_normalised_filtered$probes_by_CHR == 'X'))
head(which(df_normalised_filtered$probes_by_CHR == 'Y'))
head(which(df_normalised_filtered$probes_by_CHR == 'Un'))
# These indexes will cause errors when using other files. Convert or comment out:
#df_normalised_filtered[c(39, 930, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
#df_normalised_filtered[c(17484, 16809, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
#########

#########
#View(remove_sex_probes)
head(remove_sex_probes)
length(remove_sex_probes)
dim(df_normalised_filtered) # This object has the expression values after excluding low quality, probes 
                            # overlapping SNPs, probes without symbol or EntrezID
colnames(df_normalised_filtered)

df_normalised_filtered <- df_normalised_filtered[-remove_sex_probes, ]
dim(normalised)
dim(normalised_expressed)
dim(normalised_expressed_annotated)
dim(normalised_expressed_annotated_qual)
dim(normalised_expressed_annotated_qual_noSNPs)
dim(normalised_expressed_annotated_qual_noSNPs_noID)
dim(df_normalised_filtered) # After removing probes from sex and unannotated chromosomes
colnames(df_normalised_filtered)

which(df_normalised_filtered$probes_by_CHR == 'X')
which(df_normalised_filtered$probes_by_CHR == 'Y')
which(df_normalised_filtered$probes_by_CHR == 'Un')
#df_normalised_filtered[c(39, 930, 16825), c('Probe_Id', 'CHROMOSOME', 'probes_by_CHR')]
#########
###################


##########
# TO DO / CHECK / CLEAN:
#Read Illumina's annotation file downloaded from webpage (alternative to illuminaHumanv4.db package).
# The arrays themselves have columns with annotations. I haven't checked if these are the same and/or which is better
# (illuminaHumanv4 package). All seem to have reported issues.
# Another option is the following:
# http://biorxiv.org/content/biorxiv/early/2015/05/21/019596.full.pdf

# Read Illumina's annotation file:
#illumina_annotation_file <- as.data.frame(read.csv('/ifs/projects/proj043/analysis.dir/HumanHT-12_V4_0_R2_15002873_B.txt', header = TRUE,
#                                     sep = '\t', skip = 8, na.strings=c(""," ","NA")))

#View(illumina_annotation_file)
#class(illumina_annotation_file)
#head(illumina_annotation_file)
#tail(illumina_annotation_file)
#dim(illumina_annotation_file)

#summary(illumina_annotation_file$Chromosome)
#length(which(is.na(illumina_annotation_file$Chromosome)))

#illumina_annotation_file_chr <- illumina_annotation_file[which(!is.na(illumina_annotation_file$Chromosome)), ]
#dim(illumina_annotation_file)
#dim(illumina_annotation_file_chr)
#View(illumina_annotation_file_chr)

#chrs_to_keep <- as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
#                               11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
#                               21, 22))
#chrs_to_keep
#illumina_annotation_file_chr_autosome <- subset(illumina_annotation_file_chr, 
#                                                illumina_annotation_file_chr$Chromosome %in% chrs_to_keep)

#count(illumina_annotation_file_chr$Chromosome)
#length(which(illumina_annotation_file_chr$Chromosome == 'X'))
#length(which(illumina_annotation_file_chr$Chromosome == 'Y'))
#dim(illumina_annotation_file_chr_autosome)

#row.names(illumina_annotation_file_chr_autosome) <- illumina_annotation_file_chr_autosome$Probe_Id
#View(illumina_annotation_file_chr_autosome)
#head(illumina_annotation_file_chr_autosome$Probe_Id)

#normalised_expressed_annotated_qual_noSNPs_noID_nosex <- merge(normalised_expressed_annotated_qual_noSNPs_noID, illumina_annotation_file_chr_autosome, all.x = TRUE)
#Probe_Id <- which(colnames(illumina_annotation_file_chr_autosome) == 'Probe_Id')
#Chromosome <- which(colnames(illumina_annotation_file_chr_autosome) == 'Chromosome')

#probes_to_keep_chr <- illumina_annotation_file_chr_autosome[, c(Probe_Id, Chromosome)]
#head(probes_to_keep_chr)
#dim(probes_to_keep_chr)
############



##############
## Rename object for downstream analysis:
# Rename row IDs:
row.names(df_normalised_filtered) <- df_normalised_filtered$Probe_Id
head(df_normalised_filtered)
#df_normalised_filtered[1:5, 1:5]

# Annotations and expression data:
normalised_filtered_annotated <- df_normalised_filtered

# Keep only expression data for plotting. First 8 columns are annotation, last column is 'probes_by_CHR:
colnames(df_normalised_filtered)
df_normalised_filtered[1:5, 7:12]
normalised_filtered <- df_normalised_filtered[, c(8:(ncol(df_normalised_filtered)-1))]
head(normalised_filtered)

dim(normalised_filtered_annotated)
class(normalised_filtered_annotated)
head(normalised_filtered_annotated$SYMBOL)
colnames(normalised_filtered_annotated)
#View(normalised_filtered_annotated)

dim(normalised_filtered)
class(normalised_filtered)
str(normalised_filtered)
head(normalised_filtered)
#View(normalised_filtered)

# Sanity checks:
count(!is.na(normalised_filtered_annotated$SYMBOL))
count(!is.na(normalised_filtered_annotated$probes_by_CHR))
count(normalised_filtered_annotated$probes_by_CHR)
count(!is.na(normalised_filtered_annotated$Probe_Id))
count(normalised_filtered_annotated$Status)

#############################

###############################
# Save files required for other packages/programmes:
head(normalised_filtered) # expression values only
head(normalised_filtered_annotated) # expression values plus annotations
dim(normalised_filtered)
dim(normalised_filtered_annotated)
#head(membership_file_cleaned)
# Skip if testing without SNP filter:
write.table(normalised_filtered_annotated, 'normalised_filtered_annotated.tab', sep='\t', 
            quote=FALSE, na='NA', col.names=NA, row.names=TRUE)
write.table(normalised_filtered, 'normalised_filtered_expression_values.tab', sep='\t', 
            quote=FALSE, na='NA', col.names=NA, row.names=TRUE)
#col.names=NA, row.names=TRUE is required otherwise it skips the first position and 
#generates a tab file with the first column on rownames (screws up the header).
#write.table(membership_file_cleaned, 'membership_file_cleaned_all.tab', sep='\t', quote=FALSE, na='NA', 
#            col.names=NA, row.names=TRUE)

# fwrite(data.table(normalised$E, keep.rownames = TRUE), 'normalised.csv', na='NA')
write.csv(normalised$E, 'normalised.csv', quote=FALSE, na='NA')
# write.csv(normalised, 'normalised_test.csv', quote=FALSE, na='NA')
# write.csv(normalised_expressed_annotated_qual_noSNPs, 'normalised_expressed_annotated_qual_noSNPs.csv', quote=FALSE, na='NA')
write.csv(normalised_filtered, 'normalised_filtered_annotated.csv', quote=FALSE, na='NA')
write.csv(normalised_filtered, 'normalised_filtered_expression_values.csv', quote=FALSE, na='NA')
###############################

#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
#save.image(file=R_session_saved_image_full, compress='gzip')
objects_to_save <- (c('normalised', 'normalised_filtered_annotated', 'normalised_filtered', 
                      'membership_file_cleaned', 'FAILED_QC_unique'))
save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# Print session information:
sessionInfo()
     
q()

# Next: run the script for PCA, differential gene expression, etc.
#############################
