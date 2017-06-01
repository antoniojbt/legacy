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
# load('../data.dir/R_session_saved_image_cytokines.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_cytokines','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
#############################################

# TO DO:
# Check 'normal' levels and ln interpretation
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
library(reshape2)
library(plyr)
library(dplyr)
#############################


#############################################
# Set-up arguments:

cyto_file <- as.character(args[1])
cyto_file <- '../data.dir/bestd_cytokines.csv'

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
illumina_info[100:105, c(4:5, 14)]
dim(illumina_info)
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

for (i in transcript_info$Probe_Id) {
  print(i)
  print(as.character(i) %in% as.character(expr_file$Probe_Id))
  }
# There is only one cytokine transcript which passed QC: ILMN_1674167 IL10
########

########
# # Check raw, non-QC expression file (before probe filtering, which has individual samples QC'd but not probes):
# # Load saved RData analysis from QC stage:
# load('../data.dir/R_session_saved_image_read_and_QC_full.RData', verbose=T)
# # head(raw_cleaned_as_matrix)[1:5, 1:5]
# # dim(raw_as_matrix)
# # dim(raw_as_matrix_dedup)
# # dim(raw_cleaned_as_matrix)
# dim(raw_cleaned_as_matrix_dedup) # Use this one
# # dim(raw_minimal_eset)
# # dim(raw_cleaned_minimal_eset)
# 
# rownames(raw_cleaned_as_matrix_dedup)
# which(transcript_info$Probe_Id %in% rownames(raw_cleaned_as_matrix_dedup))
# # All transcripts were initially present.
# rm(list = ls())
########

########
# # Check why these transcripts were excluded:
# # Load RData file from pre probe-filtering stage:
# # load('../data.dir/R_session_saved_image_normalisation_full.RData', verbose=T)
# # If loading RData file normalised$E (limma format), otherwise read saved file in disk as here:
# # normalised <- fread('../data.dir/normalised.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)
# dim(normalised)
# dim(normalised_as_matrix)
# dim(normalised_as_matrix_dedup)
# head(normalised_as_matrix_dedup)[1:5, 1:5]
# dim(normalised_minimal_eset)
# 
# head(normalised$E)[1:5, 1:5]
# tail(normalised$E)[1:5, (ncol(normalised$E)-10):ncol(normalised$E)]
# which(transcript_info$Probe_Id %in% rownames(normalised$E))
# which(as.character(transcript_info$Probe_Id) %in% as.character(rownames(normalised_as_matrix_dedup)))
# transcript_info
# # All transcripts were present pre-filtering.
########


########
# Run the following in order to include cytokine transcripts in plots:
# Check probe filtering script steps in:
# '02a_microarray_GEx_normalisation_probe_filtering_sex.R'
# The object 'normalised_filtered_annotated' contains all probes without the SNP filter from this RData file"
load('../data.dir/R_session_saved_image_probe_filtering_full_noSNP_filter.RData', verbose=T)

which(transcript_info$Probe_Id %in% rownames(normalised_filtered_annotated))
transcript_info
transcript_info[c(3, 5), ] # Lost due to low quality illuminaHumanv4.db 'Bad' or 'No match'
# IL8 ILMN_1666733 ILMN_179575
# IL10 ILMN_2073307   ILMN_9173
transcript_info[-c(3, 5), ] # Lost all but IL10 due to probes matching a SNP
########


########
# Subset expression file to leave only transcripts corresponding to cytokines:
transcripts_present <- expr_file[which(expr_file$Probe_Id %in% transcript_info$Probe_Id), ]
dim(transcripts_present)
transcripts_present[, 1:10]
########

########
# Check how pre-filtered transcripts look like
# normalised_filtered_annotated HAS to be from the RData image 
# R_session_saved_image_probe_filtering_full_noSNP_filter.RData
# as loaded above.
pre_SNP_filter_transcripts <- as.data.frame(normalised_filtered_annotated)
dim(pre_SNP_filter_transcripts)
head(pre_SNP_filter_transcripts)[1:5, 1:5]
colnames(pre_SNP_filter_transcripts)
pre_SNP_filter_transcripts$Probe_Id <- rownames(pre_SNP_filter_transcripts)
# Subset all cytokines in gene expression file:
transcripts_all <- pre_SNP_filter_transcripts[which(pre_SNP_filter_transcripts$Probe_Id %in% transcript_info$Probe_Id), ]
transcripts_all
dim(transcripts_all)
# Overwrite object name for downstream processing only:
transcripts_present <- transcripts_all
dim(transcripts_present)
transcripts_present[, 1:10]
########

########
# Patients IDs are kit ID, each as column, transpose to merge:
transcripts_present_t <- data.table()
transcripts_present_t <- transpose(transcripts_present)
# Expression data has kits as IDs:
transcripts_present_t$kit_id <- colnames(transcripts_present)
head(transcripts_present_t)
colnames(transcripts_present_t)
rownames(transcripts_present_t)
# Check symbols, probe IDs, etc., match:
transcripts_present_t[c(1:4, 10:15, 568), ]
# Name each column according to gene symbol when more than one transcript survives:
setnames(transcripts_present_t, old = c("V1", "V2", "V3", "V4", "V5"), 
         new = sprintf('transcript_%s', transcripts_present[, 'ILMN_GENE']))
colnames(transcripts_present_t)
dim(transcripts_present_t)
str(transcripts_present_t)
head(transcripts_present_t)
# Convert variables to numeric instead of character:
transcripts_present_t$transcript_IL10 <- as.numeric(transcripts_present_t$transcript_IL10)
for(i in c(1, 1:(ncol(transcripts_present_t)-1))) {
  transcripts_present_t[, i] <- as.numeric(as.character(transcripts_present_t[, i]))
  }
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
# convert to data.table:
class(transcripts_present_t)
transcripts_present_t <- as.data.table(transcripts_present_t)
# Set keys for merging:
setkey(transcripts_present_t, 'kit_id_randomisation')
setkey(all_data, 'kit_id_randomisation')
all_data <- merge(all_data, transcripts_present_t, all.x = TRUE)
all_data
colnames(all_data)[98:102]
# Rename columns:
new_names <- list()
for (i in colnames(all_data)[98:102]) {
  print(i)
  new_names[[length(new_names)+1]] <- list(sprintf('%s_%s', i, 'baseline'))
  # c(new_names, i = sprintf('%s_%s', i, 'baseline'))
  }
setnames(all_data, old = colnames(all_data)[98:102], new = unlist(new_names))

# Join expression data, final visit column:
colnames(all_data)
transcripts_present_t
setnames(transcripts_present_t, 'kit_id_randomisation', 'kit_id_finalVisit')
setkey(all_data, 'kit_id_finalVisit')
setkey(transcripts_present_t, 'kit_id_finalVisit')

all_data <- merge(all_data, transcripts_present_t, all.x = TRUE)
all_data
dim(all_data)
# Rename columns:
colnames(all_data)[103:107]
new_names <- list()
for (i in colnames(all_data)[103:107]) {
  print(i)
  new_names[[length(new_names)+1]] <- list(sprintf('%s_%s', i, '12months'))
  # c(new_names, i = sprintf('%s_%s', i, 'baseline'))
}
setnames(all_data, old = colnames(all_data)[103:107], new = unlist(new_names))
colnames(all_data)[98:107]
summary(all_data[, (ncol(all_data) - 10):(ncol(all_data)-1)])
# TO DO: Pretty print:
as.data.frame(summary(all_data[, (ncol(all_data) - 10):(ncol(all_data)-1)]))


# # Join geno:
# all_data <- all_data[geno_file]
# all_data

# Sanity check:
colnames(all_data)
summary(all_data[, c(99, 104)]) # IL6

# Check variable types:
str(all_data)
vars_convert <- c(5:10, 14:52, 54:56, 78, 88:107)
# Assign correct types to variables:
class(all_data)
all_data <- as.data.frame(all_data)
for (i in vars_convert){
  # print(i)
  # all_data[, i, with = F] <- as.numeric(as.character(all_data[, i, with = F]))
  all_data[, i] <- as.numeric(as.character(all_data[, i]))
}
str(all_data)
head(all_data)
#############################################


#############################################
# Explore cytokine data, descriptive analysis

# Rename groups for better plotting:
all_data$arm2[all_data$arm == 0] <- '4000_IU'
all_data$arm2[all_data$arm == 1] <- '2000_IU'
all_data$arm2[all_data$arm == 2] <- 'Placebo'
all_data$arm2 <- factor(all_data$arm2, levels = c('Placebo', '2000_IU', '4000_IU'),
                           labels = c('Placebo', '2000_IU', '4000_IU'))
summary(all_data$arm2)

#########
# Get variables of interest for circulating cytokines:
colnames(all_data)
all_data_melt_cyto <- melt(all_data, measure.vars = 88:97)
all_data_melt_cyto
head(all_data_melt_cyto)[1:5, 1:5]
dim(all_data_melt_cyto)
colnames(all_data_melt_cyto)

all_data_melt_cyto[1:10, c('value', 'variable')]
all_data_melt_cyto$variable <- factor(all_data_melt_cyto$variable, 
                                   levels = c('Ln_IFNgamma0',
                                              'Ln_IFNgamma12',
                                              'Ln_IL10_0',
                                              'Ln_IL10_12',
                                              'Ln_IL6_0',
                                              'Ln_IL6_12',
                                              'Ln_IL8_0',
                                              'Ln_IL8_12',
                                              'Ln_TNFalpha0',
                                              'Ln_TNFalpha12'),
                                   labels = c('IFNg baseline',
                                              'IFNg 12m',
                                              'IL10 baseline',
                                              'IL10 12m',
                                              'IL6 baseline',
                                              'IL6 12m',
                                              'IL8 baseline',
                                              'IL8 12m',
                                              'TNFa baseline',
                                              'TNFa 12m')
                                   )
plyr::count(all_data_melt_cyto$variable)
plyr::count(all_data_melt_cyto$arm2)
group_cyto <- factor(all_data_melt_cyto$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                     labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group_cyto)
# Plot cytokine distributions:
ggplot(data = as.data.frame(all_data_melt_cyto),
       aes(x = variable, y = value, fill = group_cyto)) + 
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'Circulating protein levels (natural logarithm)', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5))
ggsave('cytokine_boxplots.png')
#########

#########
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
plyr::count(cormat_melted$Var1)
plyr::count(cormat_melted$Var2)
class(cormat_melted)
# Plot:
# Improve with: 
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
ggplot(data = cormat_melted, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  labs(title = 'Circulating cytokines', y = '', x = '') +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text('Spearman rho'))
ggsave('cytokines_heatmap.png')

# Add vitamin D levels:
# Plot correlations between cytokines:
colnames(all_data)
head(all_data)
str(as.data.frame(all_data[, c(31, 44, 88:97)]))
cormat <- round(cor(as.data.frame(all_data[, c(31, 44, 88:97)]), 
                    use = "pairwise.complete.obs", 
                    method = 'spearman'), 2)
class(cormat)
cormat
cormat_melted <- melt(cormat)
head(cormat_melted)
str(cormat_melted)
summary(cormat_melted)
plyr::count(cormat_melted$Var1)
plyr::count(cormat_melted$Var2)
class(cormat_melted)
# Plot:
# Improve with: 
# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
ggplot(data = cormat_melted, aes(x = Var1, y = Var2, fill = value)) + 
  geom_tile() +
  labs(title = 'Circulating cytokines', y = '', x = '') +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())
ggsave('cytokines_heatmap_vitd.png')
#########

#########
# Plot all cytokine transcript data (pre-SNP filtering):
colnames(all_data)
all_data_melt_gex <- melt(all_data, measure.vars = 98:107)
all_data_melt_gex
head(all_data_melt_gex)[1:5, 1:5]
dim(all_data_melt_gex)
colnames(all_data_melt_gex)
plyr::count(all_data_melt_gex$variable)

all_data_melt_gex[1:10, c('value', 'variable')]
all_data_melt_gex$variable <- factor(all_data_melt_gex$variable, 
                                     levels = c('transcript_IFNG_baseline', 
                                                'transcript_IFNG_12months',
                                                'transcript_IL10_baseline',
                                                'transcript_IL10_12months',
                                                'transcript_IL6_baseline',
                                                'transcript_IL6_12months',
                                                'transcript_IL8_baseline',
                                                'transcript_IL8_12months',
                                                'transcript_TNF_baseline',
                                                'transcript_TNF_12months'),
                                     labels = c('IFNg baseline',
                                                'IFNg 12m',
                                                'IL10 baseline',
                                                'IL10 12m',
                                                'IL6 baseline',
                                                'IL6 12m',
                                                'IL8 baseline',
                                                'IL8 12m',
                                                'TNFa baseline',
                                                'TNFa 12m')
                                     )
plyr::count(all_data_melt_gex$variable)
plyr::count(all_data_melt_gex$arm2)
group_gex <- factor(all_data_melt_gex$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                    labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group_gex)
# Plot cytokine transcript distributions:
ggplot(data = as.data.frame(all_data_melt_gex),
       aes(x = variable, y = value, fill = group_gex)) + 
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'Transcript levels (VSN normalised, pre-SNP filtering)', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        axis.text.x = element_text(angle=90, vjust = 0.5)
        )
ggsave('cytokine_transcripts_boxplots.png', width = 10, height = 10)
#########


#########
# Facets for all data protein and transcript:
transcripts <- 
  ggplot(data = as.data.frame(all_data_melt_gex),
       aes(x = group_gex, y = value, fill = group_gex)) +
  facet_grid(~variable) +
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'Gene expression levels (VSN normalised)', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position="none",
        strip.background = element_blank()
        # legend.position = 'bottom'
        # axis.text.x = element_text(angle=90, vjust = 0.5)
        )

cytokines <-
  ggplot(data = as.data.frame(all_data_melt_cyto),
         aes(x = group_cyto, y = value, fill = group_cyto)) +
  facet_grid(~variable) +
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'Circulating levels (natural logarithm)', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(text = element_text(size = 14), 
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        strip.text.x = element_blank()
        # axis.text.x = element_text(angle=90, vjust = 0.5)
        )
plots <- arrangeGrob(transcripts, cytokines, nrow = 2)
ggsave('boxplots_cytokine_and_transcripts.png', 
       plots, height = 10, width = 12)
#########


#########
# Plot IL10 protein levels only:
colnames(all_data)
all_data_IL10 <- melt(all_data, measure.vars = c(89, 94))
all_data_IL10$variable <- factor(all_data_IL10$variable, levels = c('Ln_IL10_0', 'Ln_IL10_12'),
                                 labels = c('IL10 baseline', 'IL10 12 months'))
plyr::count(all_data_IL10$variable)
group <- factor(all_data_IL10$arm2, levels=c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
plyr::count(group)
# Plot:
boxplot(all_data$Ln_IL10_0, all_data$Ln_IL10_12)
boxplot(all_data$transcript_IL10_baseline, all_data$transcript_IL10_12months)
ggplot(data = as.data.frame(all_data_IL10),
       aes(x = variable, y = value, fill = group)) + 
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'Circulating protein levels (natural logarithm)', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('IL10_boxplots.png')
#########

#########
# Plot IL10 transcript levels only:
colnames(all_data)
summary(all_data[, c(98, 103)])
all_data_gex <- melt(all_data, measure.vars = c(98, 103))
all_data_gex$variable <- factor(all_data_gex$variable, levels = c('transcript_IL10_baseline', 'transcript_IL10_12months'),
                                 labels = c('IL10 baseline', 'IL10 12 months'))
plyr::count(all_data_gex$variable)
all_data_gex[, c(107, 108)]
summary(all_data_gex[, 108])

# Plot:
ggplot(data = as.data.frame(all_data_gex),
       aes(x = variable, y = value, fill = arm2)) + 
  geom_boxplot(position = position_dodge(1), outlier.alpha = 0.7) +
  labs(title = '', y = 'VSN normalised gene expression levels', x = '') +
  scale_color_brewer(palette = "Dark2") +
  theme_gray() +
  theme(legend.title=element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave('IL10_transcripts_boxplots.png')
#########
#############################################

#############################################
# Save csv file with transcripts info to merge later with pheno and cyto data:
write.csv(all_data, '../data.dir/BEST-D_phenotype_file_final_cytokines_and_transcripts.csv', quote = FALSE, na = 'NA')
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
# Too heavy as loadin normalisation_full RData file
# save.image(file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx
#############################
