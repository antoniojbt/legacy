#############################
#---
# Title: "Gene_expression_Illumina_QC"
# Author: "Antonio Berlanga-Taylor"
# Date: "17 December 2014"
# Updated: "26 March 2015"
# Output: html_document
#---

# Current input: Illumina HumanHT-12 v4 Expression BeadChip summary-level data.
# Outputs: various plots and tables from QC, normalisation and differential expression.
# Requires specifying 
    # phenotype file
    # (will require) array type (affymetrix or illumina)
    # arguments: phenotype file, samples to exclude, ?, check TO DOs

#############################


#############################

# Notes and references:

# Check Limma s users guide, p. 20, Section 9 Single-Channel Experimental Designs on p.40, 
# and case study in Section 17.3 on p.108:
# http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

# Also check:
# - Blog with simple instructions:
# http://gettinggeneticsdone.blogspot.co.uk/2014/12/importing-illumina-beadarray-data-into-r.html?m=1
 
# - BeadArray Expression Analysis Using Bioconductor paper in PLoS Comp Bio:
# www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1002276&representation=PDF
 
# - BeadArray vignette (Dunning et al 2014, same author group as ploscompbio paper above):
# http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf
 
# - Lumi package vignette:
# http://www.bioconductor.org/packages/release/bioc/vignettes/lumi/inst/doc/lumi.pdf

# NOTE! Different functions output different classes (e.g. EList [epxression list], eset [expression set], etc.)
# These may not work downstream! e.g. neqc() (from limma) requires an EListRaw class as input (ie what is output from read.ilmn.targets) and
# produces an EList but normaliseIllumina() (from beadarray) requires an eset as input. To avoid this use the package 'convert':
# http://www.bioconductor.org/packages/release/bioc/manuals/convert/man/convert.pdf
# use the matrices (ie xxx$E data containing probes and samples only) or convert data to other structures (see below).

# - For limma, see p. 13 for data objects created (eg EListRaw, EList, MArrayLM and TestResults) 
# and p.14 for functions used (eg summary, dim, length, ncol, nrow, dimnames, rownames and colnames).

#############################

# Get text file data from the array facility including the following:
# - Control probe file for each sample with:
# ProbeID, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# - Sample probe file for each sample with:
# ProbeID, Symbol, AVG_Signal, BEAD_STDERR, Avg_NBEADS, Detection Pval
# - Annotation columns including:
# SEARCH_KEY, ILMN_GENE, CHROMOSOME, DEFINITION, SYNONYMS
 
# Illumina gene expression microarray analysis steps (Ritchie et al. 2014 PLoS Comp Bio; 
# Ritchie et al. 2015 Nucleic Acids Res.):
 
# 1) Data acquisition
# 2) Preprocessing and quality assessment
 
# 3) Background correction, normalisation, transformation and filtering
 
# 4) Experimental design matrix specification
# 5) Descriptive statistics
# 6) Differential gene expression analysis
 
# 7) Higher level analysis: Gene ontology, co-expression, gene set analyses
# 8) Integration with other data

#############################

# Preliminaries
 
# Script is intended to run from src location (eg "/ifs/devel/antoniob/projects/BEST-D")
# with files containing data and output deposited in separate location 
# (eg /ifs/projects/proj043/analysis.dir/gene_expression.dir).
 
# Things to do before running:
# - Set working directories.
# - Provide file names (if running one file (see read.ilmn function), or multiple files (create a summary file 
# first, see read.ilmn.targets function)).
# - Requires several R packages, see below, mainly limma.
# - Requires Illumina processed files returned as summary level data, no QC for bead level data at the moment. 
# - After QC, arrays will likely need to be removed, this needs intervention 
# (after assessing numbers and plots) and then passing array (Sample IDs) as a list for limma to remove. 
# TO DO: divide this into two scripts.
# - TO DO: Normalisation also needs intervention, set for quantile at the moment.

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
working_dir <- ("/ifs/projects/proj043/analysis.dir/gene_expression.dir")
setwd(working_dir)

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_2015-03-30.RData', verbose=T)

# To load multiple .RData files:
# rdata_filenames <- c('.RData')
# lapply(rdata_filenames, load, .GlobalEnv)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_',Sys.Date(),'.RData', sep='')
R_session_saved_image
#############################


#############################
#File with names of files with data
#Create a tab-delimited file with two columns ("files" header= sample probe profiles and 
#"ctrlfiles" header= control profiles):
targets_file_input <- "targets_file.txt"
targets_file_input

#File with meta-data (group membership, treatment status, replicate number, etc.)
membership_input_file <- 'GEx_BEST-D_targets_file.txt'
membership_input_file

#Lis of samples that fail QC = 
failed_QC_input_file <- c('120005280', '120005133', '120000312', '120005211', '120000131', '120005098', 
                          '120000272', '120005145', 'misload')
failed_QC_input_file

#Normalisation method to call =
#xxx

#Design matrix specification =
#xxx

#Contrast matrix specification =
#xxx


#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("limma")
#biocLite("beadarray")
#biocLite("lumi")
#biocLite("illuminaHumanv4.db")
#biocLite("arrayQualityMetrics")
#biocLite("quantro")
#biocLite("flashClust")
#biocLite('convert')
#install.packages('dendextend')

library(limma)
library(beadarray)
library(lumi)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
library(RColorBrewer)
library(genefilter)
library(Biobase)
library(quantro)
library(flashClust)
library(lumiHumanAll.db)
library(convert)
library(vsn)
library(reshape2)
#library(grid)
library(gridExtra)
library(plyr)
library(dendextend)
library(gplots)
library(doParallel)

#############################

# 1) Data acquisition:
#   Raw data comprises one observation per pixel, per array.
# Bead-level data comprises one observation per bead, per array.
# Summary-level data comprises one observation per probe type, per sample
# (Ritchie et al. PLos Comp Bio, 2014)
  
# a) Bead-level data: intensity and location information for each bead on each BeadArray is 
# available (eg txt files with this information plus other optional/recommended files such as .bab 
# [compressed data from the txt files], raw TIFF image data, sdf files, targets file, metrics file, etc., 
# for analysis with beadarray package). See pdf vignette (Dunning et al 2014) above, analysis not included here.
  
# b) Summary level data (typical): files after processing in the BeadStudio software, for analysis with limma, 
# lumi and beadarray packages. Files can be at different levels of processing (intensities with/without 
# background correction and/or normalised). Files include sample probe profile (text file as data-matrix with 
# 48,000 rowswith non-normalised summary values as output by BeadStudio, required); control probe profile (text 
# file with summarised data for each of the controls of each array used for diagnostics and calibration, recommended);
# targets file (user created text file with information on which sample was hybridised to each array, 
# required if reading in multiple probe profiles). Files with normalised intensities typically end in avg and 
# files with intensites per gene may also be available. Avoid these as well as Illumina background correction and 
# normalisation.
  
# Get probe summary profiles (containing the intensity data). Best to obtain intensities from 
# GenomeStudio or BeadStudio without background correction or normalization. Probe summary files are 
# tab-delimited and usually arrays processed at one time are written to a single file.
# Also obtain profiles for the control probes from BeadStudio or GenomeStudio processed data.
  
# After processing with BeadStudio, for each array probe profile files contain: 
#   summarized expression level (AVG_Signal)
# standard error of the bead replicates (BEAD_STDERR)
# number of beads (Avg_NBEADS) 
# detection p-value (DetectionPval), estimation of the probability of a gene being detected 
# above the background level


## Load and read the Illumina probe files (typically one per experiment containing multiple arrays with 
#a separate control probe profile file output). An EListRaw class is created in limma and should hold ~50,127 
#rows and 6 columns. Each row specificies if it is a control or sample probe. 'xxx'$targets data frame specifies 
#which samples were hybridised to each array. 

#Columns will include source, E (with the expression value for each probe), genes, other, Detection Pval, targets.

#If multiple probe files need to be read use the read.ilmn.targets function, 
#(requires a summary file with probe profiles (one, or two columns if control probe files are available).
#Can also just read one file with corresponding control.

targets_file <- readTargets(targets_file_input)
targets_file


#Then read all files:
# TO DO: check this is true: If 'ctrlfiles=NULL' is specified it will not create the $genes$Status column which
# contains the seven types of control probes (important for QC and normalisation). 

read_files <- read.ilmn.targets(targets=targets_file, probeid="PROBE_ID", 
                                annotation=c("SYMBOL", "SEARCH_KEY", "ILMN_GENE", "CHROMOSOME", "DEFINITION", 
                                             "SYNONYMS"), expr="AVG_Signal", other.columns="Detection Pval", 
                                sep="\t", verbose=TRUE)

#?read.ilmn.targets
# Control probes are labelled: negative, biotin, labeling, cy3_hyb, housekeeping, high_stringency_hyb or low_stringency_hyb.
#?read.ilmn
#length(which(read_files_cleaned_QC$genes$Status == 'negative'))
#length(which(read_files_cleaned_QC$genes$Status == 'regular'))
#length(read_files_cleaned_QC$genes$Status)

#############################


#############################

# xx) Read phenotype file and metadata

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# For BEST-D: Cross kit_ids from 03_lab_kits file with column names (ie kit id labels) from expression data, 
# this to extract only samples measured and have same number of samples.


#Save sample names to cross with phenotype IDs:
membership_file <- read.csv(membership_input_file, header=TRUE, row.names=1, sep='\t')

dim(membership_file)
head(membership_file)
tail(membership_file)
summary(membership_file)

# TO DO: change this to eg grouping1
# Add columns needed for PCA and other plots and groupings:
membership_file$treatment <- ifelse(membership_file$arm == 0 & membership_file$visit_type == 'FinalVisit',  'treated_4000', 
                                    ifelse(membership_file$arm == 1 & membership_file$visit_type == 'FinalVisit',  'treated_2000',
                                           'untreated'))
head(membership_file)
tail(membership_file)

#Get sample names (stored as columns in the limma read files):
#sample_IDs_clean <- colnames(read_files_cleaned_QC)

#Cross IDs:
#matched_IDs <- which(rownames(membership_file) %in% sample_IDs_clean)

#Get file with no NAs and only matched samples to those in the Expression Set object:
#membership_file_cleaned <- membership_file[matched_IDs,]

#head(membership_file_cleaned)
#tail(membership_file_cleaned)

#############################


#############################

# 2) Preprocessing and quality assessment
 
# If available/possible, run per array signal-to-noise value plot (95th percentile of signal divided by 
# the 5th percentile); spatial plots of the intensities across array surface to detect array artefacts; 
# between sample comparison with bloxplots of intensities; and MDS plot to determine between sample differences 
# (biological, batches, etc.) (see Ritchie et al. 2014 Fig 2 for examples)

# Positive controls can be used to identify suspect arrays.
# Negative control probes, which measure background signal on each array, can be used to assess the proportion 
# of expressed probes that are present in a given sample (Beadarray vignette).


#a) Descriptives

#Look at summary information of the probe profiles:

#Array dimensions:
#Number of arrays match number of samples expected? ~570
#Number of rows is expected number of probes? eg ~48,000
dim(read_files)
#All files were included and read?
read_files$targets
read_files$E[1:5, ]

#Check objects class and other attributes:
class(read_files)
class(read_files$E)
str(read_files)
str(read_files$E)

#Number of negative and regular probes. Illumina BeadChip arrays contain 750~1600 negative control probes:
table(read_files$genes$Status)

#View expression values for first 5 columns, first 10 samples:
options(digits=5)
head(read_files$E[1:5,1:10])
#All samples:
head(read_files$other$'Detection Pval')
length(read_files$other$Detection)
length(read_files$other$'Detection Pval')

#Explore contents of file:
head(read_files)
dim(read_files)
class(read_files)
summary(read_files)
head(read_files$source)
head(read_files$E)[1:5,1:5]
range(read_files$E)
median(read_files$E)
head(read_files$genes)
summary(read_files$other)
head(read_files$other$'Detection Pval')[1:5,1:5]
range(read_files$other$'Detection Pval')
read_files$targets


#See p-values for detection, these test whether each probe is more intense than the negative control probes.
#Small values indicate that the probe corresponds to true expression:

range(read_files$other$Detection)
mean(read_files$other$Detection)
median(read_files$other$Detection)


#Boxplots of intensities to assess dynamic range from each sample and identify outliers from signal distributions. 
#The intensities vary from about 5 to 14 on the log2 scale:

#Boxplots for x number of samples (run loop for random sets? Plot all separately?):
intensities_plot <- ("boxplots_of_intensities.png")
png(intensities_plot, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[,400:420]),range=0,ylab="log2 intensity", xlab="Array")
dev.off()
#All samples (too many to visualise for BEST-D):
#boxplot(log2(read_files$E),range=0,ylab="log2 intensity of probe intensities", xlab="Array")


#Separate boxplots of regular probes and control probes to highlight unusual samples:
#Regular probes:
boxplot_intensities_regular_probes <- ("boxplot_intensities_regular_probes.png")
png(boxplot_intensities_regular_probes, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[read_files$genes$Status == "regular", ]), 
        range= 0, las = 2, xlab = "", ylab =expression(log[2](intensity)), main = "Regular probes")
dev.off()

#Negative probes:
length(which(read_files$genes$Status == 'negative'))
boxplot_intensities_negative_probes <- ("boxplot_intensities_negative_probes.png")
png(boxplot_intensities_negative_probes, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files$E[read_files$genes$Status == "negative", ]), 
        range= 0, las = 2, xlab = "", ylab = expression(log[2](intensity)), main = "Negative control probes")
dev.off()


#MDS plots:
intensities_MDS <- ("MDS_of_intensities.png")
png(intensities_MDS, width = 4, height = 4, units = 'in', res = 300)
#'Multidimensional scaling plot of probe intensities')
plotMDS(read_files$E, pch=1)
dev.off()

#The 'propexpr' function estimates the proportion of expressed probes in each array by comparing the empirical 
#intensity distribution of the negative control probes with that of the regular probes. A mixture model is fitted 
#to the data from each array to infer the intensity distribution of expressed probes and estimate the expressed 
#proportion.

#Get proportion of expressed probes and descriptives:
proportion <- propexpr(read_files)
head(proportion)
length(proportion)
range(proportion)
mean(proportion)
median(proportion)
quantile(proportion)

propexpr_plot <- ("boxplot_of_propexpr.png")
png(propexpr_plot, width = 4, height = 4, units = 'in', res = 300)
boxplot(proportion, ylab='Proportion of expressed probes', xlab='All samples')
dev.off()

#Arrays with low expression (which to exclude? what value is min?):
which_low <- which(proportion <= 0.25)
length(which_low)
proportion[which_low]

which_med <- which(proportion > 0.25 & proportion <= 0.75)
length(which_med)
proportion[which_med]

which_high <- which(proportion > 0.75)
length(which_high)
proportion[which_high]

#Which samples to compare if any? Random subsets? Null hypothesis is that 
#t.test(proportion[-which_low], proportion[which_low])


#############################
## Call arrayQualityMetrics
#The arrayQualityMetrics package collates quality assessment plots for summarized data created by beadarray to 
#identify outlier arrays:

# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):
# Different methods, check: 
# ?ExpressionSet
# object<-new(Class=, exprs=as.matrix(normalised_expressed))
# ?new, ?coerce, ?as 
# as(as_matrix, 'ExpressionSet')
# read_files_cleaned_QC_eset <- coerce(read_files_cleaned_QC, to='eset', strict=TRUE)
# convert(read_files_cleaned_QC, 'ExpressionSet')


# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
raw_as_matrix <- as.matrix(read_files, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(raw_as_matrix)
dim(raw_as_matrix)
head(raw_as_matrix)[1:5,1:5]
head(colnames(raw_as_matrix))
head(rownames(raw_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(raw_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
raw_as_matrix_dedup <- raw_as_matrix[-duplicated_probes,]
head(raw_as_matrix_dedup)
dim(raw_as_matrix_dedup)

# Convert expression matrix to expression set:
raw_minimal_eset <- ExpressionSet(assayData=raw_as_matrix_dedup)
class(raw_minimal_eset)
dim(raw_minimal_eset)

# Call arrayQualityMetrics before processing data:

arrayQualityMetrics_preprocessed <- arrayQualityMetrics(expressionset=raw_minimal_eset, 
                                                        outdir=paste('arrayQualityMetrics_preprocessed_data.dir'), 
                                                        force=TRUE, do.logtransform=TRUE)

#############################


#############################
## b) Remove failed samples from EListRaw object before continuing with analysis. These will mainly be from the hybridisation QC 
# and samples failing the arrayQualityMetrics package metrics:
# TO DO: Pass sample labels from arrayQualityMetrics

#Pass samples that failed QC:
FAILED_QC <- failed_QC_input_file

#Get samples IDs and index numbers from EListRaw object:
array_sample_IDs <- colnames(read_files)
to_extract <- match(FAILED_QC, array_sample_IDs)

#Check indexes match ID:
read_files[0,to_extract]
read_files$E[0,to_extract]

#Get clean EListRaw object:
read_files_cleaned_QC <- read_files[,-to_extract]

dim(read_files_cleaned_QC)
dim(read_files)

#Explore contents of file:
head(read_files_cleaned_QC)
class(read_files_cleaned_QC)
summary(read_files_cleaned_QC)
head(read_files_cleaned_QC$source)
head(read_files_cleaned_QC$E)[1:5,1:5]
range(read_files_cleaned_QC$E)
median(read_files_cleaned_QC$E)
head(read_files_cleaned_QC$genes)
summary(read_files_cleaned_QC$other)
head(read_files_cleaned_QC$other$'Detection Pval')[1:5,1:5]
range(read_files_cleaned_QC$other$'Detection Pval')
read_files_cleaned_QC$targets

## Clean membership file. Get sample names (stored as columns in the limma read files):
sample_IDs_clean <- colnames(read_files_cleaned_QC)

#Cross IDs:
matched_IDs <- which(rownames(membership_file) %in% sample_IDs_clean)

#Get file with no NAs and only matched samples to those in the Expression Set object:
membership_file_cleaned <- membership_file[matched_IDs,]

head(membership_file_cleaned)
tail(membership_file_cleaned)


#Separate boxplots of regular probes and control probes to highlight unusual samples in the cleaned file:
#Regular probes:
boxplot_intensities_regular_probes_cleaned_QC <- ("boxplot_intensities_regular_probes_cleaned_QC.png")
png(boxplot_intensities_regular_probes_cleaned_QC, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files_cleaned_QC$E[read_files_cleaned_QC$genes$Status == "regular", ]), 
        range= 0, las = 2, xlab = "", ylab =expression(log[2](intensity)), main = "Regular probes")
dev.off()

#Negative probes:
length(which(read_files_cleaned_QC$genes$Status == 'negative'))
boxplot_intensities_negative_probes_cleaned_QC <- ("boxplot_intensities_negative_probes_clean_QC.png")
png(boxplot_intensities_negative_probes_cleaned_QC, width = 4, height = 4, units = 'in', res = 300)
boxplot(log2(read_files_cleaned_QC$E[read_files_cleaned_QC$genes$Status == "negative", ]), 
        range= 0, las = 2, xlab = "", ylab = expression(log[2](intensity)), main = "Negative control probes")
dev.off()



#############################
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
raw_cleaned_as_matrix <- as.matrix(read_files_cleaned_QC, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(raw_cleaned_as_matrix)
dim(raw_cleaned_as_matrix)
head(raw_cleaned_as_matrix)[1:5,1:5]
head(colnames(raw_cleaned_as_matrix))
head(rownames(raw_cleaned_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(raw_cleaned_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
raw_cleaned_as_matrix_dedup <- raw_cleaned_as_matrix[-duplicated_probes,]
head(raw_cleaned_as_matrix_dedup)
dim(raw_cleaned_as_matrix_dedup)

# Convert expression matrix to expression set:
raw_cleaned_minimal_eset <- ExpressionSet(assayData=raw_cleaned_as_matrix_dedup)
class(raw_cleaned_minimal_eset)
dim(raw_cleaned_minimal_eset)


# Call arrayQualityMetrics before processing data and after cleaning raw files:
arrayQualityMetrics_raw_cleaned <- arrayQualityMetrics(expressionset=raw_cleaned_minimal_eset, 
                                                       outdir=paste('arrayQualityMetrics_raw_cleaned.dir'), 
                                                       force=TRUE, do.logtransform=TRUE)

#############################


#############################

# 3) Background correction, normalisation, transformation and filtering

# Reading the control probe profiles is optional but recommended. If the control probe profiles are available, 
# then the Illumina data can be favorably background corrected and normalized using the neqc or nec functions.

## Test whether quantile normalisation is appropriate.
# The package quantro tests whether quantile normalisation is appropriate for the given dataset:
# http://www.bioconductor.org/packages/release/bioc/vignettes/quantro/inst/doc/quantro-vignette.pdf
# browseVignettes("quantro")

# View the distributions of the samples of interest, matdensity() computes the density for each sample (columns):

matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xlab = " ", ylab = "density", 
           main = "Beta Values")
#legend('top', c("NeuN_neg", "NeuN_pos"), col = c(2,3), lty= 1, lwd = 3)

# matboxplot() orders and colors the samples by a group level variable:
matboxplot(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xaxt = "n", 
           main = "Beta Values")

# Run quantro tests and permutations in parallel with doParallel library:
# TO DO extract number of cores as argument
registerDoParallel(cores=4)
quantro_test <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$treatment, B=10)
quantro_test_by_arm <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$arm, B=10)
quantro_test_visit_type <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$visit_type, B=10)

# Assess statistical significance:
summary(quantro_test)
anova(quantro_test)
quantroStat(quantro_test)
#quantroStat(quantro_test_by_arm)
#quantroStat(quantro_test_visit_type)


# Plot the test statistics of the permuted samples. The plot is a histogram of the null test statistics quantroStatPerm 
# from quantro() and the red line is the observed test statistic quantroStat from quantro():

quantroPlot(quantro_test)
#quantroPlot(quantro_test_by_arm)
#quantroPlot(quantro_test_visit_type)


# a) Background correction
# Non-background corrected, non-normalised, sample and control probe profiles and targets file.
# Check backgroundCorrect and plotFB

class(read_files_cleaned_QC)
summary(read_files_cleaned_QC)
summary(read_files_cleaned_QC$genes$Status)

# Explore foreground vs background intensities and decide which method to use: TO DO: Errors due to no x$Eb present (background intensities)
# plotFB(x=read_files_cleaned_QC, array=1, pch='.')

# Background correct. Single-channel arrays methods incude none, nec or normexp. 
# ?nec, use robust=TRUE? offset set at 16 from suggestion in Ritchie et al PLoS Comp Biol 2011 
# http://www.ploscompbiol.org/article/fetchObject.action?uri=info:doi/10.1371/journal.pcbi.1002276&representation=PDF

range(read_files_cleaned_QC$E)

background_corrected_nec <- nec(x=read_files_cleaned_QC)
background_corrected_nec_robust <- nec(x=read_files_cleaned_QC, robust=TRUE)
background_corrected_nec_16 <- nec(x=read_files_cleaned_QC, offset=16)
background_corrected_nec_16_robust <- nec(x=read_files_cleaned_QC, offset=16, robust=TRUE)

range(background_corrected_nec$E)
range(background_corrected_nec_robust$E)
range(background_corrected_nec_16$E)
range(background_corrected_nec_16_robust$E)

background_corrected_auto <- backgroundCorrect(read_files_cleaned_QC, method='auto', verbose=TRUE)
range(background_corrected_auto$E)

background_corrected_normexp <- backgroundCorrect(read_files_cleaned_QC, method='normexp', verbose=TRUE)
range(background_corrected_normexp$E)


# background_corrected_subtract <- backgroundCorrect(read_files_cleaned_QC, method='subtract', verbose=TRUE)
# range(background_corrected_subtract$E)

# background_corrected_movingmin <- backgroundCorrect(read_files_cleaned_QC, method='movingmin', verbose=TRUE)
# range(background_corrected_movingmin$E)

# background_corrected_edwards <- backgroundCorrect(read_files_cleaned_QC, method='edwards', verbose=TRUE)
# range(background_corrected_edwards$E)


# Visualise and explore the fit of the dataset:
meanSdPlot_background_corrected_nec <- ('mean_SD_plot_background_corrected_nec.png')
png(meanSdPlot_background_corrected_nec, width = 6, height = 12, units = 'in', res = 300)
par(mfrow=c(4,1))
meanSdPlot(background_corrected_nec$E, main='Background corrected expression values - nec')
meanSdPlot(background_corrected_nec_robust$E, main='Background corrected expression values - nec, robust')
meanSdPlot(background_corrected_nec_16$E, main='Background corrected expression values - nec, offset 16')
meanSdPlot(background_corrected_nec_16_robust$E, main='Background corrected expression values - nec, robust, offset 16')
par(mfrow=c(1,1))
dev.off()


meanSdPlot_background_corrected_normexp <- ('mean_SD_plot_background_corrected_normexp_and_auto.png')
png(meanSdPlot_background_corrected_normexp, width = 6, height = 10, units = 'in', res = 300)
par(mfrow=c(2,1))
meanSdPlot(background_corrected_normexp$E, main='Background corrected expression values - normexp')
# meanSdPlot_background_corrected_normexp <- ('mean_SD_plot_background_corrected_normexp.png')
meanSdPlot(background_corrected_auto$E, main='Background corrected expression values - auto')
par(mfrow=c(1,1))
dev.off()


# Mean-difference plots, using limma's plotMA:
# Also check mdplot

plotMA_background_corrected_1 <- ('mean_difference_plots_background_corrected_1.png')
png(plotMA_background_corrected_1, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(read_files_cleaned_QC$E, main='Raw expression values')
limma::plotMA(background_corrected_nec$E, main='Background corrected expression values - nec')
limma::plotMA(background_corrected_nec_robust$E, main='Background corrected expression values - nec, robust')
limma::plotMA(background_corrected_nec_16$E, main='Background corrected expression values - nec, offset 16')
par(mfrow=c(1,1))
dev.off()


plotMA_background_corrected_2 <- ('mean_difference_plots_background_corrected_2.png')
png(plotMA_background_corrected_2, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(background_corrected_nec_16_robust$E, main='Background corrected expression values - nec, robust, offset 16')
limma::plotMA(background_corrected_normexp$E, main='Background corrected expression values - normexp')
limma::plotMA(background_corrected_auto$E, main='Background corrected expression values - auto')
par(mfrow=c(1,1))
dev.off()


# Plot density distributions:
# ?plotDensity

plotDensity_background_corrected <- ('density_distribution_background_corrected.png')
png(plotDensity_background_corrected, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,2))
plotDensity(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density plot of raw expression values')
plotDensity(background_corrected_nec$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec')
plotDensity(background_corrected_nec_16$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, 16')
plotDensity(background_corrected_nec_robust$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, robust')
plotDensity(background_corrected_nec_16_robust$E, logMode=F, addLegend=F, 
            main='Density distribution of background corrected expression values - nec, robust, 16')
plotDensity(background_corrected_normexp$E, logMode=F, addLegend=F, main='Background corrected normexp')

par(mfrow=c(1,1))
dev.off()


# TO DO: extract as parameter:
background_corrected = background_corrected_normexp

class(background_corrected)
head(background_corrected)
str(background_corrected)
summary(background_corrected)
range(background_corrected$E)


#############################
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
background_corrected_as_matrix <- as.matrix(background_corrected, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(background_corrected_as_matrix)
dim(background_corrected_as_matrix)
head(background_corrected_as_matrix)[1:5,1:5]
head(colnames(background_corrected_as_matrix))
head(rownames(background_corrected_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(background_corrected_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
background_corrected_as_matrix_dedup <- background_corrected_as_matrix[-duplicated_probes,]
head(background_corrected_as_matrix_dedup)
dim(background_corrected_as_matrix_dedup)

# Convert expression matrix to expression set:
background_corrected_minimal_eset <- ExpressionSet(assayData=background_corrected_as_matrix_dedup)
class(background_corrected_minimal_eset)
dim(background_corrected_minimal_eset)

# Call arrayQualityMetrics after processing data and background correction:
arrayQualityMetrics_background_corrected <- arrayQualityMetrics(expressionset=background_corrected_minimal_eset, 
                                                                outdir=paste('arrayQualityMetrics_background_corrected.dir'), 
                                                                force=TRUE, do.logtransform=TRUE)


#############################

# b) Transformation
# log2 transformation is run by neqc, see other approaches.

# c) Normalisation
# Quantile approach is run already by limma s neqc function.
# Also check VST and other approaches.


## The neqc function performs normexp background correction using negative controls, then quantile normalizes 
# and finally log2 transforms [see ]. It also automatically removes the control probes, leaving only the regular 
# probes:
#TO DO: When applying quantile normalization, it is assumed that the distribution in signal should be 
# the same from each array. 

run_neqc <- neqc(read_files_cleaned_QC)
dim(run_neqc)
summary(run_neqc)
class(run_neqc)
range(run_neqc$E)

## Other functions to normalise and  transform: 
#?normalizeBetweenArrays # This is run for single channel arrays. For two-colour arrays normaliseWithinArrays must be run first.
# For single-channel data check scale, quantile and cyclicloess 

# TO DO: extract method as argument

between_arrays <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='scale')
class(between_arrays)
head(between_arrays)
dim(between_arrays)
range(between_arrays$E)

between_arrays_loess <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='cyclicloess')
class(between_arrays_loess)
head(between_arrays_loess)
dim(between_arrays_loess)
range(between_arrays_loess$E)

between_arrays_loess_affy <- normalizeBetweenArrays(object=read_files_cleaned_QC, method='cyclicloess', cyclic.method='affy')


# Variance-stabilising normalisation from the 'vsn' package. This background corrects then normalises. 
# Limma's normalizeVSN is an interface to vsn's vsnMatrix function (which run's vsn's vsn2 function, 
# see ?normalizeVSN and ?vsn2 or ?vsnMatrix, vignette('Introduction to vsn'), ?normalizeVSN.

normalize_VSN <- normalizeVSN(read_files_cleaned_QC)
class(normalize_VSN)
head(normalize_VSN)
dim(normalize_VSN)
range(normalize_VSN$E)
head(normalize_VSN$E)

# Other methods for normalisation (require object as ExpressionSet so transform data object):
# ?normaliseIllumina from beadarray

#transform_vst <- normaliseIllumina(read_files_cleaned_QC, transform='vst')
#run_median <- normaliseIllumina(minimal_eset, method='median')
#run_qspline <- normaliseIllumina(minimal_eset, method='qspline')


#Visualise normalisation methods:
# run_neqc
# between_arrays
# between_arrays_loess
# normalize_VSN

meanSdPlot_normalised <- ('meanSdPlot_normalised.png')
png(meanSdPlot_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
meanSdPlot(run_neqc$E, main='Normalised expression values - neqc')
meanSdPlot(between_arrays$E, main='Normalised expression values - scale')
meanSdPlot(between_arrays_loess$E, main='Normalised expression values - cyclic loess')
meanSdPlot(normalize_VSN$E, main='Non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# Mean-difference plots, using limma's plotMA:
# Also check mdplot

plotMA_normalised <- ('mean_difference_plots_normalised.png')
png(plotMA_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
limma::plotMA(run_neqc$E, main='Normalised expression values - neqc')
limma::plotMA(between_arrays$E, main='Normalised expression values - scale')
limma::plotMA(between_arrays_loess$E, main='Normalised expression values - cyclic loess')
limma::plotMA(normalize_VSN$E, main='Non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# Plot density distributions:
# ?plotDensity

plotDensity_normalised <- ('density_distribution_normalised.png')
png(plotDensity_normalised, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
plotDensity(run_neqc$E, logMode=F, addLegend=F, main='Density plot of normalised expression values - neqc')
plotDensity(between_arrays$E, logMode=F, addLegend=F, 
            main='Density distribution of normalised expression values - scale')
plotDensity(between_arrays_loess$E, logMode=F, addLegend=F, 
            main='Density distribution of normalised expression values - cyclic loess')
plotDensity(normalize_VSN$E, logMode=F, addLegend=F, 
            main='Density distribution of non-background corrected, normalised expression values - VSN')
par(mfrow=c(1,1))
dev.off()


# TO DO: Turn into parameter set at the beginning of the script .
# Change object names for downstream analysis:
normalised = normalize_VSN

# Explore contents of the ExpressionSet created (as this is the output from normaliseIllumina):
normalised
dim(normalised)
class(normalised)
summary(normalised)
range(normalised$E)


# If looking at an ExpressionSet:
#head(featureNames(object=normalised))
#head(sampleNames(object=normalised))

# Obtain the actual expression values from the ExpressionSet (this returns a matrix):
#expression_values_normalised <- exprs(normalised)

#class(expression_values_normalised)
#head(expression_values_normalised)
#range(expression_values_normalised[,1])


## TO DO: Plot and explore normalised data
# Extract expression values and convert to matrix (needed for some plots):
# normalised_as_matrix <- as.matrix(normalised$E, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
# class(normalised_as_matrix)
# dim(normalised_as_matrix)
# head(normalised_as_matrix)[1:5,1:5]
# head(colnames(normalised_as_matrix))
# head(rownames(normalised_as_matrix))


## Visualise and explore the fit of the dataset:

# Cumulative density plot

# Pairwise sample correlations

# Heatmap



## Summary plots of raw values and chosen background correction and normalisation methods:
# Density distributions:
summary_plots_density <- ('summary_plots_density_distributions.png')
png(summary_plots_density, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plotDensity(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density distribution of raw expression values')
plotDensity(background_corrected$E, logMode=F, addLegend=F, main='Density distribution of background corrected expression values')
plotDensity(normalised$E, logMode=T, addLegend=F, main='Density distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()

# Plot row standard deviations versus row means:
summary_meanSD_plots <- ('summary_plots_mean_SD.png')
png(summary_meanSD_plots, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
meanSdPlot(read_files_cleaned_QC$E, logMode=F, addLegend=F, main='Density distribution of raw expression values')
meanSdPlot(background_corrected$E, logMode=F, addLegend=F, main='Density distribution of background corrected expression values')
meanSdPlot(normalised$E, logMode=T, addLegend=F, main='Density distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()


# Plot cumulative distributions, ?plotCDF:
plotCDF_summary <- ('summary_cumulative_distribution_plots.png')
png(plotCDF_summary, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plotCDF(read_files_cleaned_QC$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of raw expression values')
plotCDF(background_corrected$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of background corrected expression values')
plotCDF(normalised$E, reverse=FALSE, logMode=TRUE, addLegend=FALSE, main='Cumulative distribution of normalised expression values')
par(mfrow=c(1,1))
dev.off()

#############################
## Call arrayQualityMetrics
# Transform EListRaw class to ExpressionSet class (as more general and accepted by more downstream methods):

# This extracts expression values (with probes as rownames) and sample names (column headers) and returns as a matrix:
normalised_as_matrix <- as.matrix(normalised, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
class(normalised_as_matrix)
dim(normalised_as_matrix)
head(normalised_as_matrix)[1:5,1:5]
head(colnames(normalised_as_matrix))
head(rownames(normalised_as_matrix))


# Remove duplicated probes from expression matrix (these seem to be the negative control probes):
duplicated_probes <- which(duplicated(rownames(normalised_as_matrix)))
duplicated_probes

# If probes have duplicated names but for some reason need to be kept: rownames(as_matrix) = make.names(rownames(as_matrix), unique=TRUE)

# Output clean expression matrix:
normalised_as_matrix_dedup <- normalised_as_matrix[-duplicated_probes,]
head(normalised_as_matrix_dedup)
dim(normalised_as_matrix_dedup)

# Convert expression matrix to expression set:
normalised_minimal_eset <- ExpressionSet(assayData=normalised_as_matrix_dedup)
class(normalised_minimal_eset)
dim(normalised_minimal_eset)

# Call arrayQualityMetrics after processing data and normalisation:
arrayQualityMetrics_normalised <- arrayQualityMetrics(expressionset=normalised_minimal_eset, 
                                                      outdir=paste('arrayQualityMetrics_normalised.dir'), 
                                                      force=TRUE)


#############################

## d) Filtering and creating an eset object

# Filter out probes that are not expressed. Keep probes that are expressed in at least three arrays 
# according to a detection p-value of 5% (Limma vignette case study, p. 108):

#Set filters:
expressed <- rowSums(normalised$other$'Detection Pval' < 0.05) >= 3

#Extract data and create an ExpressionSet (eset, or EList object) at the same time (necessary for linear modelling steps):
normalised_expressed <- normalised[expressed,]

#Explore contents of file:
head(normalised_expressed)
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


#TO DO: Filter based on annotation quality:
# Optional: use Illumina's annotation:
# See http://www.bioconductor.org/help/workflows/annotation/annotation/
# Remove probes that aren't annotated with a gene: #This doesn't work because 
# Biobase expects an eset object while limma outputs EList objects.
# See: http://www.bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf


# annotated_probes <- !is.na((fData(normalised_expressed)$SYMBOL))
# table(annotated)
# eset <- eset[annotated,]
# rm(annotated)


#Plot log2 intensity of regular probes normalised:
boxplot_normalised_expressed_intensities <- ("regular_probes_normalised.png")
png(boxplot_normalised_expressed_intensities, width = 4, height = 4, units = 'in', res = 300)
boxplot(normalised_expressed$E,range= 0, ylab =expression(log[2](intensity)), 
        las = 2, xlab = "", main = "Regular probes, normalised")
dev.off()

#Plot expressed probes in a multi-dimensional scaling plot:
plot_normalised_expressed <- ("MDS_normalised_expressed.png")
png(plot_normalised_expressed, width = 4, height = 4, units = 'in', res = 300)
plotMDS(normalised_expressed, pch=1)
dev.off()

plot_normalised_expressed_by_targets <- ("MDS_normalised_expressed_by_targets.png")
png(plot_normalised_expressed_by_targets, width = 4, height = 4, units = 'in', res = 300)
plotMDS(normalised_expressed, pch=1, labels=normalised_expressed$targets)
dev.off()

#TO DO:
# plot_normalised_expressed_by_sample <- ("MDS_normalised_expressed_by_sample.png")
# png(plot_normalised_expressed_by_sample, width = 4, height = 4, units = 'in', res = 300)
# plotMDS(normalised_expressed, pch=1, labels=normalised_expressed$targets, colors=) 
# dev.off()


# Plot PCA of normalised samples:
# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca_normalised_expressed <- prcomp(t(normalised_expressed$E), center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc <- pca_normalised_expressed$x

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(normalised_expressed$E)
str(pca_normalised_expressed)
summary(pca_normalised_expressed)

# Plot PCA results:
plot_PCA_normalised_expressed <- ('plot_PCA_normalised_expressed.png')
png(plot_PCA_normalised_expressed, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_normalised_expressed, main='Normalised expression values')
# Scatterplot of PC1 and PC2:
biplot(pca_normalised_expressed, main='Normalised expression values')

par(mfrow=c(1,1))
dev.off()

# Check how much of the variance in gene expression values is explained by the first x PCs:
# sum(pca_normalised_expressed$sdev[1:10]^2)/length(normalised_expressed$E[1,])


# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
head(membership_file_cleaned)
tail(membership_file_cleaned)
pc_data <- data.frame(pca_normalised_expressed$x[,1:13])
str(pc_data)
head(pc_data)

pca_by_groups <- data.frame(merge(membership_file_cleaned, pc_data, by='row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)
# View(pc_data)
# View(pca_by_groups)
pc_data['120000222',]
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_normalised_expressed_by_groups_1 <- ('plot_PCA_normalised_expressed_by_groups_1.png')
png(plot_PCA_normalised_expressed_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$arm)) + theme(legend.position="bottom")
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$visit_type)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p4 <- qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()

plot_PCA_normalised_expressed_by_groups_2 <- ('plot_PCA_normalised_expressed_by_groups_2.png')
png(plot_PCA_normalised_expressed_by_groups_2, width = 13, height = 13, units = 'in', res = 300)
p5 <- qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p6 <- qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p7 <- qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p8 <- qplot(x=PC6, y=PC7, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p5, p6, p7, p8, ncol=2)
dev.off()

plot_PCA_normalised_expressed_by_groups_3 <- ('plot_PCA_normalised_expressed_by_groups_3.png')
png(plot_PCA_normalised_expressed_by_groups_3, width = 13, height = 13, units = 'in', res = 300)
p9 <- qplot(x=PC7, y=PC8, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p10 <- qplot(x=PC8, y=PC9, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p11 <- qplot(x=PC9, y=PC10, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p12 <- qplot(x=PC10, y=PC11, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
grid.arrange(p9, p10, p11, p12, ncol=2)
dev.off()


# Plot dendrogram
# Prepare hierarchical cluster:
correlation <- cor(normalised_expressed$E, method='pearson')
head(correlation)
# View(correlation)

hc <- flashClust(dist(correlation))
str(hc)
head(hc)

# Convert to a dendogram object:
dendro <- as.dendrogram(hc)

# Add colours by groups
colourCodes <- c(treated_4000="red", treated_2000="green", untreated="blue")

# Assign labels of dendrogram object with new colors:
labels_colors(dendro) <- colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)]
head(colourCodes[membership_file_cleaned$treatment][order.dendrogram(dendro)])
head(order.dendrogram(dendro))
head(membership_file_cleaned$treatment)

# Plot simple dendrogram, labels at the same level:
plot_dendrogram_normalised_expressed <- ('plot_dendrogram_normalised_expressed.png')
png(plot_dendrogram_normalised_expressed, width = 13, height = 13, units = 'in', res = 300)
par(mfrow=c(1,2), cex = 1)
plot(dendro, hang = -1)
# Zoom in to a sub tree:
plot(dendro[[1]], horiz = TRUE)
par(mfrow=c(1,1))
dev.off()

# Plot correlation between samples in a heatmap:
cor_normalised_expressed <- melt(correlation)
head(cor_normalised_expressed)
summary(cor_normalised_expressed)
plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_expressed.png')
png(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
p1 <- ggplot(cor_normalised_expressed, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
grid.arrange(p1, ncol=1)
dev.off()


# TO DO: 
# Plot expression values in a heatmap:
plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_expressed.png')
png(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(1,2))
#heatmap(normalised_expressed$E[,1:10])
heatmap(normalised_expressed$E[1:100,1:100], ColSideColors=colourCodes)
# p1 <- ggplot(cor_normalised_expressed, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
# grid.arrange(p1, ncol=1)
par(mfrow=c(1,1))
dev.off()

# Try again with heatmap.2:
# http://sebastianraschka.com/Articles/heatmaps_in_r.html

# Create colour palette from red to green:
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition:
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0,0.8,length=100),              # for yellow
               seq(0.8,1,length=100))              # for green

# creates a 5 x 5 inch image
png('heatmap_normalised_expressed_heatmap2.png',    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(normalised_expressed$E[1:100, 1:100], 
          cellnote = normalised_expressed$E[1:100, 1:100],  # same data set for cell labels
          main = "Correlation", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          #          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier 
          #          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="column",     # only draw a column dendrogram
          #          Colv="NA")            # turn off column clustering
)
dev.off()               # close the PNG device


# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# t <- ExpressionSet(e, AnnotatedDataFrame(tab))
# rv <- rowVars(exprs(t))
# idx <- order(-rv)[1:40]
# heatmap(exprs(t)[idx, ], col = hmcol)


#############################


#############################

# 4) Experimental design matrix specification

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

#Check dimensions between annotation file with meta-data (must have the same number of rows, otherwise
#errors downstream):
#TO DO: Raise error if not the same.

dim(membership_file_cleaned)
str(membership_file_cleaned)
dim(normalised_expressed)

#############################


#############################

#5) Descriptive statistics
# 
# 
# 
# 


#############################


#############################

# 6) Differential gene expression analysis
# Linear modelling using weights
 
# Analyze using linear models in limma. A model is fit for every gene and an empirical Bayes method
# moderates the standard errors of the estimated log-fold changes
 
# Requires one or two matrices to be specified. A design matrix which has each array/sample per row, 
# and coefficients that describe the RNA sources in each column (eg group, treatment, batch/plate, etc.).
# A contrast matrix can then be used to group the coeffients for comparisons.
# Input data is the ExpressionSet (eset) or the EList class.
# Main functions: lm(), eBayes(), topTable(), etc.

# Help:
#?lm()
#?eBayes

# For BEST-D:
# arm
# visit type
# plate


## a) Paired (before vs after comparison)

## Test: Compare all samples before vs after according to visit type (baseline vs 12 months):
#Define factors to constrast from annotation file:
arm <- factor(membership_file_cleaned$arm, levels=c('0', '1', '2'))
visit_type <- factor(membership_file_cleaned$visit_type, levels=c("Randomisation","FinalVisit"))

#Define design:
design <- model.matrix(~0+visit_type)
head(design)
tail(design)

#Run linear model and set contrasts:
fit <- lmFit(normalised_expressed, design)
cont.matrix <- makeContrasts(Randomisation_vs_FinalVisit=visit_typeFinalVisit-visit_typeRandomisation, levels=design)

#Obtain differentially expressed genes based on contrasted factors:
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH")

#Get results and plot:
results <-decideTests(fit2)
#vennDiagram(results)

head(fit2)[1]
head(fit2$p.value)
head(fit2$lods)
range(fit2$lods)

# Plot volcano, p-values vs effect size:
# log10_p_values <- -log10(fit2$p.value)
# plot(fit2$lods, log10_p_values, pch='.', xlab='log-ratio', ylab=expression(-log[10]~p))
# abline(h=2)
# Not right... using limma's function:

# volcanoplot(fit2, coef=1, highlight=0, names=fit2$genes$ID, xlab="Log Fold Change", ylab="Log Odds", pch=16, cex=0.35, abline(h=2))
# warnings()

#?volcanoplot 
# Not right... data? design? etc!


# b) Two group comparison: 
#Placebo vs high dose
## Test: Compare all samples before vs after according to visit type (baseline vs 12 months):

#Define factors to constrast from annotation file:
#Define design:
#Run linear model and set contrasts:
#Obtain differentially expressed genes based on contrasted factors:
#Get results and plot:

#Placebo vs low dose
#Define factors to constrast from annotation file:
#Define design:
#Run linear model and set contrasts:
#Obtain differentially expressed genes based on contrasted factors:
#Get results and plot:

#Low dose vs high dose
#Define factors to constrast from annotation file:
#Define design:
#Run linear model and set contrasts:
#Obtain differentially expressed genes based on contrasted factors:
#Get results and plot:

# c) Three group comparison
#Define factors to constrast from annotation file:
#Define design:
#Run linear model and set contrasts:
#Obtain differentially expressed genes based on contrasted factors:
#Get results and plot:

# d) Difference-in-difference estimator
#Define factors to constrast from annotation file:
#Define design:
#Run linear model and set contrasts:
#Obtain differentially expressed genes based on contrasted factors:
#Get results and plot:

##Compare paired samples (eg based on 'SibShip', here using pt_id) for before vs after according to visit 
#type (baseline vs 12 months):
#Define factors to constrast from annotation file (as above plus column that pairs arrays):
#Here visit_type amounts to treament vs baseline

#arm and visit_type defined as above
#Pairing:
pairing <- factor(membership_file_cleaned$pt_id)
head(pairing)

#Define design and set contrasts:
design_all_pairs <- model.matrix(~pairing+visit_type)
head(design_all_pairs)[1:5,1:5]
tail(design_all_pairs)[,-1]

dim(design_all_pairs)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_pairs <- lmFit(normalised_expressed, design_all_pairs)
fit_all_pairs <- eBayes(fit_all_pairs)
head(fit_all_pairs)[1]

topTable(fit_all_pairs, coef="visit_typeFinalVisit", adjust='BH')

#Plot:
results_pairing <- decideTests(fit_all_pairs)
head(results_pairing)
dim(results_pairing)

#vennDiagram(results_pairing)

# volcanoplot(fit_all_pairs, coef="visit_typeFinalVisit", highlight=0, names=fit_all_pairs$genes$ID, 
#             xlab="Log Fold Change", ylab="Log Odds", pch=16, cex=0.35, abline(h=2))


#Compare samples based on treatment arm:
#Define design and set contrasts:
design_by_arm <- model.matrix(~0+arm)
colnames(design_by_arm) <- c('FourIU', 'TwoIU', 'Placebo')
head(design_by_arm)
tail(design_by_arm)

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_by_arm <- lmFit(normalised_expressed, design_by_arm)
head(fit_by_arm)
contrast_matrix_by_arm <- makeContrasts(FourIU-TwoIU, FourIU-Placebo, TwoIU-Placebo, levels=design_by_arm)
head(contrast_matrix_by_arm)

fit_by_arm_2 <- contrasts.fit(fit_by_arm, contrast_matrix_by_arm)
fit_by_arm_2 <- eBayes(fit_by_arm_2)
head(fit_by_arm_2)
head(fit_by_arm_2$coefficients)

topTable(fit_by_arm_2, coef='FourIU - Placebo', adjust='BH', number=10)

#Plot results:
results_by_arm <- decideTests(fit_by_arm_2)
# vennDiagram(results_by_arm)

# volcanoplot(fit_by_arm_2, coef='FourIU - Placebo', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2))


#Factorial design:
#Compare samples based on treatment arm and visit type:
#Define design and set contrasts:

#Check experimental dimensions (factors):
head(membership_file_cleaned)

#Modify to leave only arm and visity type:
membership_file_cleaned_no_pt_id <- subset(membership_file_cleaned, , c(arm, visit_type))
head(membership_file_cleaned_no_pt_id)

#Collect different combinations of factors:
group_combinations <- paste(membership_file_cleaned$arm, membership_file_cleaned$visit_type, sep='.')
head(group_combinations)
tail(group_combinations)

# Define the questions of interest:
# Before and after comparisons for three groups:
# Which genes respond to treatment after stimulation at high dose: 4000IU Final Visit vs 4000IU randomisation (baseline)
#arm.visit_type = 0.FinalVisit vs 0.Randomisation

# Which genes respond to treatment after stimulation at low dose: 2000IU Final Visit vs 2000IU randomisation (baseline)
#arm.visit_type = 1.FinalVisit vs 1.Randomisation

# Which genes are likely due to noise:  Placebo Final Visit vs Placebo randomisation (baseline)
#arm.visit_type = 2.FinalVisit vs 2.Randomisation

# Two group comparison see above:
#

# Compare the difference of differences (ie interaction terms): 
# Differences due to high dose, ie 4000IU minus Placebo: (0.FinalVisit vs 0.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to low dose, ie 2000IU minus Placebo: (1.FinalVisit vs 1.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to high dose only, ie 4000IU minus 2000IU: (0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation) 
# Differences truly due to high dose only, ie (4000IU minus 2000IU) minus (Placebo):
#((0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation)) - (2.FinalVisit vs 2.Randomisation)


#Define the combinations of factors:
factorial_design <- factor(group_combinations, levels=c('0.FinalVisit', '0.Randomisation', 
                                                        '1.FinalVisit', '1.Randomisation', '2.FinalVisit', '2.Randomisation'))
head(factorial_design)

#Extract the comparisons of interest as contrasts:
design_factorial_comparisons <- model.matrix(~0+factorial_design)
head(design_factorial_comparisons)
colnames(design_factorial_comparisons) <- c('FourIU_FinalVisit', 'FourIU_Randomisation', 
                                            'TwoIU_FinalVisit', 'TwoIU_Randomisation', 'Placebo_FinalVisit', 'Placebo_Randomisation')
head(design_factorial_comparisons)  
tail(design_factorial_comparisons)
dim(design_factorial_comparisons)

#Fit a model with a coefficient for each of the factor combinations:
fit_factorial <- lmFit(normalised_expressed, design_factorial_comparisons)
head(fit_factorial)

#Extract the fitted contrasts of interest:
contrast_matrix_factorial <- makeContrasts(High_dose_before_and_after=FourIU_FinalVisit - FourIU_Randomisation, 
                                           Low_dose_before_and_after=TwoIU_FinalVisit - TwoIU_Randomisation,
                                           Placebo_before_and_after=Placebo_FinalVisit - Placebo_Randomisation,
                                           High_dose_minus_Placebo=(FourIU_FinalVisit - FourIU_Randomisation) - 
                                             (Placebo_FinalVisit - Placebo_Randomisation),
                                           Low_dose_minus_Placebo=(TwoIU_FinalVisit - TwoIU_Randomisation) - 
                                             (Placebo_FinalVisit - Placebo_Randomisation),
                                           High_dose_minus_low_dose=(FourIU_FinalVisit - FourIU_Randomisation) -
                                             (TwoIU_FinalVisit - TwoIU_Randomisation),
                                           Truly_high_dose_only=(FourIU_FinalVisit - FourIU_Randomisation) -
                                             (TwoIU_FinalVisit - TwoIU_Randomisation) - (Placebo_FinalVisit - Placebo_Randomisation),                                           
                                           levels=design_factorial_comparisons)
contrast_matrix_factorial

fit_factorial_2 <- contrasts.fit(fit_factorial, contrast_matrix_factorial)
fit_factorial_2 <- eBayes(fit_factorial_2)
head(fit_factorial_2)
head(fit_factorial_2$coefficients)
range(fit_factorial_2$lods)
range(fit_factorial_2$p.value)

topTable_fit_factorial_2 <- topTable(fit_factorial_2, coef='High_dose_before_and_after', sort.by='logFC', number=100)
head(topTable_fit_factorial_2)
range(topTable_fit_factorial_2$logFC)
range(topTable_fit_factorial_2$AveExpr)


topTableF_fit_factorial_2 <- topTable(fit_factorial_2, coef='High_dose_before_and_after', number=100)
head(topTableF_fit_factorial_2)
range(topTableF_fit_factorial_2$logFC)
range(topTableF_fit_factorial_2$AveExpr)

all_probes <- topTable(fit_factorial_2, coef='High_dose_before_and_after', number=Inf, sort.by='p')
head(all_probes)
range(all_probes$logFC)
min_logFC <- min(all_probes$logFC)
max_logFC <- max(all_probes$logFC)
min_FC <- 2^min_logFC
min_FC
max_FC <- 2^max_logFC
max_FC

range(all_probes$AveExpr)
range(all_probes$P.Value)
range(all_probes$adj.P.Val)
range(all_probes$B)


results_by_factorial <- decideTests(fit_factorial_2)
#vennDiagram(results_by_factorial)

#?volcanoplot
#?topTable

#volcano_plot_high_dose_b_vs_a <- ('volcano_plot_high_before_vs_after.png')
#png(volcano_plot_high_dose_b_vs_a, width = 12, height = 12, units = 'in', res = 300)

volcano_plots <- ('volcano_plots.png')
png(volcano_plots, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,3))

volcanoplot(fit_factorial_2, coef='High_dose_before_and_after', highlight=10, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_before_and_after')

volcanoplot(fit_factorial_2, coef='Low_dose_before_and_after', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_before_and_after')

volcanoplot(fit_factorial_2, coef='Placebo_before_and_after', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Placebo_before_and_after')

volcanoplot(fit_factorial_2, coef='High_dose_minus_Placebo', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_Placebo')

volcanoplot(fit_factorial_2, coef='Low_dose_minus_Placebo', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_minus_Placebo')

volcanoplot(fit_factorial_2, coef='High_dose_minus_low_dose', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_low_dose')

volcanoplot(fit_factorial_2, coef='Truly_high_dose_only', highlight=0, names=fit_by_arm_2$genes$ID, 
            xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Truly_high_dose_only')

par(mfrow=c(1,1))
dev.off()

# volcanoplot(fit_factorial_2, coef=4, highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='')

# Plot volcano, p-values vs effect size:
# log10_p_values <- -log10(fit_factorial_2$p.value)
# plot(fit_factorial_2$lods, log10_p_values, pch='.', xlab='log-ratio', ylab=expression(-log[10]~p))
# abline(h=2)

#############################


#############################

# 7) Integration with other data
# a) eQTL analysis
# b) Disease variant overlap
# c) Pathway analysis


#############################


#############################
#The end:
# TO DO: remove objects that are not necessary to save:
#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))
#To save R workspace with all objects to use at a later time:
save.image(R_session_saved_image, compress='gzip')

q()
#############################