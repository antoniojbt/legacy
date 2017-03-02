#############################################
# BEST-D cytokine data analysis

# Author: Antonio J Berlanga-Taylor
# Date: 20 January 2017

#Purpose
#=======

#############################################
# Cytokine data from BEST-D study, J. Emberson email 23 Jan 2017
# BEST-D cytokine data (IFN gamma, IL-6, IL-8, IL-10, and TNF-alpha). 
# All are measured at 0 and 12 months (variable suffix 0 or 12). Units for all are Ln pg/ml.

# All null from CTSU analysis
#############################################


#Methods
#=======


#Usage
#=====


#To use type::
#    xxx.R [options] [arguments]
#    xxx.R --help


#Options
#=======

#-I    input file name.
#-S    output file name.
#-L    log file name.


# Input
# Output

# Use docopt, see
#https://github.com/AntonioJBT/various.dir/blob/master/Notes-common-cmds/docopt_argument_parser.txt
#https://github.com/docopt/docopt.R

#############################################


#############################################
# Logging
# TO DO: move to a separate script

##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('~/Desktop/BEST_D.DIR/mac_runs_to_upload/BEST-D_cytokines/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_cytokines",Sys.Date(),".txt", sep=""))
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
# load('R_session_saved_image_cytokines.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokines','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################

# TO DO:
# Check 'normal' levels and ln interpretation, missing lab method.
# Join info with geno, expr, VD, pheno:
# Create one file with pt_id, kit_id, candidate geno, candidate expr, pheno and VD
# This file has pt_id with 0 and 12 months for each variable
# Pull in above data and join
# Plot gene expression with cytokines as main figure

# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups
# Correlation to VD and other measures?
# Association with genotype, SNPs tested only

#############################
# Import libraries:
library(ggplot2)
library(data.table)
library(gridExtra)
# library(GGally)
library(reshape)
library(plyr)
#############################


#############################################
# Set-up arguments:

cyto_file <- as.character(args[1])
cyto_file <- 'bestd_cytokines.csv'

# See script 
# /cgat_projects/BEST-D/genotype_analysis/plink_output_processing.R
geno_file <- as.character(args[2])
# e.g. 'rs1993116_plink_results.raw'
# geno_file <- ''

# See /Users/antoniob/Documents/github.dir/AntonioJBT/Airwave/eQTL_plotting_airwave.R
expr_file <- as.character(args[3])
expr_file <- '../data.dir/normalised_filtered_annotated.tab'

pheno_file <- as.character(args[4])
pheno_file <- '../data.dir/BEST-D_phenotype_file_final.tsv'

# Probe annotation file to get corresponding info for genes:
illumina_file <- as.character(args[5])
illumina_file <- '../annotation_external/HumanHT-12_V4_0_R2_15002873_B.txt'

# Genes of interest:
genes <- list('IFNG', 'IL6', 'IL8', 'IL10', 'TNF')
#############################################


#############################################
# Read files:
cyto_file <- fread(cyto_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
cyto_file
head(cyto_file)
dim(cyto_file)
tail(cyto_file)
summary(cyto_file)

# geno_file <- fread(geno_file, sep = ',', header = TRUE, stringsAsFactors = FALSE)
# geno_file
# head(geno_file)
# dim(geno_file)
# tail(geno_file)
# summary(geno_file)

expr_file <- fread(expr_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
expr_file
head(expr_file)
dim(expr_file)
tail(expr_file)
head(expr_file)[1:5, 1:10]
# summary(expr_file)

pheno_file <- fread(pheno_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE,
                    na.strings = '-99999')
pheno_file
head(pheno_file)
dim(pheno_file)
tail(pheno_file)
summary(pheno_file)

illumina_info <- fread(illumina_file, sep = '\t', 
                       header = TRUE, stringsAsFactors = FALSE,
                       skip = 8, nrows = 47323) # The first rows are metadata, the last rows 
                                                # are additional control probe info
head(illumina_info)
tail(illumina_info)[, 1:10]
dim(illumina_info)
colnames(illumina_info)
#############################################

#############################################
# Extract information for cytokines measured and gene expression

########
# Cytokines measured:
# IFN gamma, IL-6, IL-8, IL-10, and TNF-alpha

# Corresponding genes:
genes
colnames(cyto_file)

# Corresponding Illumina transcripts:
head(illumina_info)
colnames(illumina_info)

# Check how one gene looks like:
illumina_info[which(illumina_info[, 'ILMN_Gene'] == 'IFNG'), ]
# Function to get info from Illumina annotation file:
get_transcript <- function(data_table, colID1, colID2, colID3, gene){
  subsetted <- data_table[which(data_table[, colID1, with = F] == gene), ]
  transcript <- subsetted[, c(colID1, colID2, colID3), with = F]
  return(transcript)
}
# Check:
get_transcript(illumina_info, 'ILMN_Gene', 'Probe_Id', 'Transcript', 'IFNG')

# Get all genes of interest and put them into a dataframe:
transcript_info <- data.frame()
for (i in genes)
{
  print(i)
  transcript <- get_transcript(illumina_info, 'ILMN_Gene', 'Probe_Id', 'Transcript', i)
  transcript_info <- rbind(transcript_info, transcript)
}
transcript_info
########

########
# Check if cytokine transcripts made it to the final expression file:
transcript_info
head(expr_file)[1:5, 1:10]
class(transcript_info)
class(expr_file)

which(transcript_info$Probe_Id %in% expr_file$Probe_Id)
which(expr_file$Probe_Id %in% transcript_info$Probe_Id)
as.character(transcript_info$Probe_Id) %in% as.character(expr_file$Probe_Id)

for (i in transcript_info$Probe_Id)
{
  print(i)
  print(as.character(i) %in% as.character(expr_file$Probe_Id))
}
# There is only one cytokine transcript which passed QC
########

########
# # Check raw, non-QC expression file (before probe filtering, these is individuals QC):
# # Load saved RData analysis from QC stage:
# load('../data.dir/R_session_saved_image_read_and_QC_full.RData', verbose=T)
# head(raw_cleaned_as_matrix)[1:5, 1:5]
# rownames(raw_cleaned_as_matrix)
# which(transcript_info$Probe_Id %in% rownames(raw_cleaned_as_matrix))
# # All transcripts were initially present.
########

########
# # Check why these transcripts were excluded:
# # Load RData file from pre probe-filtering stage:
# # load('../data.dir/R_session_saved_image_normalisation_full.RData', verbose=T)
# # If loading RData file normalised$E (limma format), otherwise read saved file in disk as here:
# normalised <- fread('../data.dir/normalised.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)
# dim(normalised)
# head(normalised)[1:5, 1:5]
# tail(normalised)[1:5, (ncol(normalised)-10):ncol(normalised)]
# normalised[, 'V1']
# which(transcript_info$Probe_Id %in% normalised$V1)
# # All transcripts were present pre-filtering.


# # Check probe filtering script steps in:
# # '02a_microarray_GEx_normalisation_probe_filtering_sex.R'
# load('../data.dir/R_session_saved_image_normalisation.RData', verbose=T)
# which(transcript_info$Probe_Id %in% rownames(normalised_expressed$E))
# which(transcript_info$Probe_Id %in% rownames(normalised_expressed_annotated$E))
# which(transcript_info$Probe_Id %in% rownames(normalised_expressed_annotated_qual$E))
# transcript_info[1, ] # Lost due to low quality illuminaHumanv4.db 'Bad' or 'No match'
# which(transcript_info$Probe_Id %in% rownames(normalised_expressed_annotated_qual_noSNPs$E))
# transcript_info[2, ]
# transcript_info[-2, ] # Lost all but IL10 due to probes matching a SNP
########

########
# Subset expression file to leave only transcripts corresponding to cytokines:
transcripts_present <- expr_file[which(expr_file$Probe_Id %in% transcript_info$Probe_Id), ]
dim(transcripts_present)
transcripts_present[, 1:10]

# Patients IDs are kit ID, each as column, transpose to merge:
transcripts_present_t <- data.table()
transcripts_present_t <- transpose(transcripts_present)
# Expression data has kits as IDs:
transcripts_present_t$kit_id <- colnames(transcripts_present)
transcripts_present_t
# TO DO: loop to name each column according to gene symbol when more than one transcript survives:
setnames(transcripts_present_t, 'V1', 
         sprintf('transcript_%s', transcripts_present[, 'ILMN_GENE']))
dim(transcripts_present_t)
colnames(transcripts_present_t)
str(transcripts_present_t)
# TO DO: set up so that loop renames accordingly
transcripts_present_t$transcript_IL10 <- as.numeric(transcripts_present_t$transcript_IL10)
str(transcripts_present_t)
summary(transcripts_present_t)
########
#############################################


#############################################
# Create one file with pt_id, kit_id, candidate geno, candidate expr, pheno and VD
# Join files, these are data.table

# Join pheno and cyto:
pheno_file[, 'pt_id', with = F]
cyto_file[, 'pt_id', with = F]
setkey(pheno_file, 'pt_id')
setkey(cyto_file, 'pt_id')

all_data <- merge(pheno_file, cyto_file, all=TRUE)
all_data
dim(cyto_file)
dim(pheno_file)
dim(all_data)
head(all_data)[, 1:10]
head(all_data)[, (ncol(all_data) - 10):ncol(all_data)]
summary(all_data[, (ncol(all_data) - 10):ncol(all_data)])

# Join expression data, baseline visit column:
setnames(transcripts_present_t, 'kit_id', 'kit_id_randomisation')

# Make sure variable types are the same:
str(transcripts_present_t)
str(all_data$kit_id_randomisation)
all_data$kit_id_randomisation <- as.character(all_data$kit_id_randomisation)
all_data
transcripts_present_t
# Set keys for merging:
setkey(transcripts_present_t, 'kit_id_randomisation')
setkey(all_data, 'kit_id_randomisation')
all_data <- merge(all_data, transcripts_present_t, all.x = TRUE)
setnames(all_data, 'gene_expr', 'gene_expr_baseline')

# Join expression data, final visit column:
colnames(all_data)
transcripts_present_t
setnames(transcripts_present_t, 'kit_id_randomisation', 'kit_id_finalVisit')
setkey(all_data, 'kit_id_finalVisit')
setkey(transcripts_present_t, 'kit_id_finalVisit')

all_data <- merge(all_data, transcripts_present_t, all.x = TRUE)
all_data
dim(all_data)
setnames(all_data, 'gene_expr', 'gene_expr_12months')
summary(all_data[, (ncol(all_data) - 10):ncol(all_data)])

# # Join geno:
# all_data <- all_data[geno_file]
# all_data

# Sanity check:

#############################################


#############################################
# Explore cytokine data, descriptive analysis

# Rename groups for better plotting:
all_data$arm2[all_data$arm == 0] <- '4000_IU'
all_data$arm2[all_data$arm == 1] <- '2000_IU'
all_data$arm2[all_data$arm == 2] <- 'Placebo'
count(all_data$arm2)

# Get variables of interest:
colnames(all_data)
all_data_melt <- melt(all_data, measure.vars = 88:97)
all_data_melt
head(all_data_melt)[1:5, 1:5]
dim(all_data_melt)
colnames(all_data_melt)

all_data_melt[1:10, c('value', 'variable')]
count(all_data_melt$variable)
count(all_data_melt$arm2)
group <- factor(all_data_melt$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
count(group)

# Plot cytokine distributions:
ggplot(data = as.data.frame(all_data_melt),
       aes(x = variable, y = value, fill = group)) + 
  geom_boxplot(position = position_dodge(1)) +
  labs(title = 'Circulating cytokines') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank())
ggsave('cytokine_boxplots.png')


# Plot IL10 and transcript:
colnames(all_data)
all_data_melt_IL10 <- melt(all_data, measure.vars = c(89, 94, 98, 99))
                           # id.vars = c('Ln_IL10_0', 'Ln_IL10_12',
                           #                        'transcript_IL10.x',
                           #                        'transcript_IL10.y',
                           #                         'arm2')) #c(89, 94, 98, 99))
colnames(all_data_melt_IL10)
all_data_melt_IL10
str(all_data_melt_IL10)
all_data_melt_IL10[1:10, c('value', 'variable')]
count(all_data_melt_IL10$variable)
count(all_data_melt_IL10$arm2)
group <- factor(all_data_melt_IL10$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
count(group)
# Plot:
ggplot(data = as.data.frame(all_data_melt_IL10),
       aes(x = variable, y = value, fill = group)) + 
  geom_boxplot(position = position_dodge(1)) +
  labs(title = 'IL-10 protein and transcript levels') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank())
ggsave('IL10_protein_transcript.png')


# Plot correlations between cytokines:
colnames(all_data)
str(as.data.frame(all_data[, 88:97]))
cormat <- round(cor(as.data.frame(all_data[, 88:97]), 
                    use = "pairwise.complete.obs", 
                    method = 'spearman'), 2)
class(cormat)
cormat
cormat_melted <- melt(cormat)
head(cormat_melted)
str(cormat_melted)
summary(cormat_melted)
count(cormat_melted$X1)
count(cormat_melted$X2)
class(cormat_melted)
# Plot:
ggplot(data = cormat_melted, aes(x = X1, y = X2, fill = value)) + 
  geom_tile() +
  labs(title = 'Circulating cytokines')#, legend(legend = 'Spearman rho'))
ggsave('cytokines_heatmap.png')
# group <- factor(all_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"),
#                 labels = c("Placebo", "2000 IU", "4000 IU"))
#############################################


#############################################
# Plot gene expression with cytokines as main figure

# Create table with values of interest:


# Save table


#############################################


#############################################
# Basic sanity check, see thresholds, basic comparisons between 0, 12 and groups

#############################################


#############################################
# Correlation to VD and other measures?

#############################################


#############################################
# Association with genotype, SNPs tested only

#############################################


#############################################
## Save some text:
# cat(file = 'xxx.txt', xxx_var, "\t", xxx_var, '\n', append = TRUE)
#############################################

#############################################
# 

#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
# ls()
# object_sizes <- sapply(ls(), function(x) object.size(get(x)))
# as.matrix(rev(sort(object_sizes))[1:10])
#rm(list=ls(xxx))
#objects_to_save <- (c('xxx_var'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx
#############################
