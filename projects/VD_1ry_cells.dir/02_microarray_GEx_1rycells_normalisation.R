#############################
# To be run after 01 read and QC array data

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd("/ifs/projects/proj043/analysis.dir/gene_expression.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_normalisation",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Load a previous R session, data and objects:
load('R_session_saved_image_read_and_QC.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_normalisation', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_normalisation_full', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("limma")
#install.packages('dendextend')
#detach("package:pryr", unload=TRUE)

library(limma)
library(beadarray)
library(lumi)
library(illuminaHumanv4.db)
library(arrayQualityMetrics)
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


#############################

# 3) Background correction, normalisation, transformation and filtering

# Reading the control probe profiles is optional but recommended. If the control probe profiles are available, 
# then the Illumina data can be favorably background corrected and normalized using the neqc or nec functions.

## Test whether quantile normalisation is appropriate.
# The package quantro tests whether quantile normalisation is appropriate for the given dataset:
# http://www.bioconductor.org/packages/release/bioc/vignettes/quantro/inst/doc/quantro-vignette.pdf
# browseVignettes("quantro")

# View the distributions of the samples of interest, matdensity() computes the density for each sample (columns):

png('density_plots_raw_cleaned_quantro.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xlab = " ", ylab = "density", 
           main = "Beta Values, by treatment")
#legend('top', c("NeuN_neg", "NeuN_pos"), col = c(2,3), lty= 1, lwd = 3)

matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$arm, col = c(1,2,3), xlab = " ", ylab = "density", 
           main = "Beta Values, by arm")

matdensity(read_files_cleaned_QC$E, groupFactor = membership_file_cleaned$visit_type, col = c(1,2), xlab = " ", ylab = "density", 
           main = "Beta Values, by visit type")

par(mfrow=c(1,1))
dev.off()

# matboxplot() orders and colors the samples by a group level variable:
#png('boxplots_raw_cleaned_subset_quantro.png', width = 12, height = 12, units = 'in', res = 300)
#par(mfrow=c(3,1))
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$treatment, col = c(1,2,3), xaxt = "n", 
#           main = "Beta Values, by treatment")
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$arm, col = c(1,2,3), xaxt = "n",
#           main = "Beta Values, by arm")
#matboxplot(read_files_cleaned_QC$E[1:100,1:200], groupFactor = membership_file_cleaned$visit_type, col = c(1,2), xaxt = "n",
#           main = "Beta Values, by visit type")
#par(mfrow=c(1,1))
#dev.off()

# Run quantro tests and permutations in parallel with doParallel library:
# TO DO extract number of cores as argument; run complete files for permutations
registerDoParallel(cores=4)
quantro_test_treatment <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$treatment, B=1000)
quantro_test_timepoint <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$timepoint, B=1000)
quantro_test_replicate <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$replicate, B=1000)
quantro_test_individual <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$individual, B=1000)
quantro_test_cell_type <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$cell_type, B=1000)
quantro_test_sentrix_ID <- quantro(object=read_files_cleaned_QC$E, membership_file_cleaned$sentrix_ID, B=1000)


quantro_list <- c(quantro_test_treatment, quantro_test_timepoint, quantro_test_replicate, quantro_test_individual, 
                  quantro_test_cell_type, quantro_test_sentrix_ID)

# Assess statistical significance:
print(quantro_list)
lapply(quantro_list, summary)
lapply(quantro_list, anova)
lapply(quantro_list, quantroStat)


# Plot the test statistics of the permuted samples. The plot is a histogram of the null test statistics quantroStatPerm 
# from quantro() and the red line is the observed test statistic quantroStat from quantro():

png('quantro_histogram_permutations.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(1,1))
#lapply(quantro_list, quantroPlot)
quantroPlot(quantro_test_treatment)
#quantroPlot(quantro_test_by_arm) 
#quantroPlot(quantro_test_visit_type)
par(mfrow=c(1,1))
dev.off()

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
background_corrected = background_corrected_nec

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
#arrayQualityMetrics_background_corrected <- arrayQualityMetrics(expressionset=background_corrected_minimal_eset, 
#                                                                outdir=paste('arrayQualityMetrics_background_corrected.dir'), 
#                                                                force=TRUE, do.logtransform=TRUE)


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
class(between_arrays_loess_affy)
head(between_arrays_loess_affy)
dim(between_arrays_loess_affy)
range(between_arrays_loess_affy$E)


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
plotDensity(normalised$E, logMode=F, addLegend=F, main='Density distribution of normalised expression values')
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
#arrayQualityMetrics_normalised <- arrayQualityMetrics(expressionset=normalised_minimal_eset, 
#                                                      outdir=paste('arrayQualityMetrics_normalised.dir'), 
#                                                      force=TRUE)


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
head(pc)

# Explore dimensions and plot first 10 or so components:
head(pc)
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

pca_by_groups <- data.frame(merge(membership_file_cleaned, pc_data, by.x='sample', by.y='row.names'))
#View(pca_by_groups)
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)

# View(pc_data)
# View(pca_by_groups)
pc_data['41',]
head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_normalised_expressed_by_groups_1 <- ('plot_PCA_normalised_expressed_by_groups_1.png')
png(plot_PCA_normalised_expressed_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$individual)) + theme(legend.position="bottom")
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$timepoint)) + theme(legend.position="bottom")
p4 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$replicate)) + theme(legend.position="bottom")
p5 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$cell_type)) + theme(legend.position="bottom")
p6 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(membership_file_cleaned$sentrix_ID)) + theme(legend.position="bottom")
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3)
dev.off()

# plot_PCA_normalised_expressed_by_groups_2 <- ('plot_PCA_normalised_expressed_by_groups_2.png')
# png(plot_PCA_normalised_expressed_by_groups_2, width = 13, height = 13, units = 'in', res = 300)
# p5 <- qplot(x=PC3, y=PC4, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p6 <- qplot(x=PC4, y=PC5, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p7 <- qplot(x=PC5, y=PC6, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p8 <- qplot(x=PC6, y=PC7, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# grid.arrange(p5, p6, p7, p8, ncol=2)
# dev.off()
# 
# plot_PCA_normalised_expressed_by_groups_3 <- ('plot_PCA_normalised_expressed_by_groups_3.png')
# png(plot_PCA_normalised_expressed_by_groups_3, width = 13, height = 13, units = 'in', res = 300)
# p9 <- qplot(x=PC7, y=PC8, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p10 <- qplot(x=PC8, y=PC9, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p11 <- qplot(x=PC9, y=PC10, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# p12 <- qplot(x=PC10, y=PC11, data=pca_by_groups, colour=factor(membership_file_cleaned$treatment)) + theme(legend.position="bottom")
# grid.arrange(p9, p10, p11, p12, ncol=2)
# dev.off()


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
colourCodes <- c(CD14="red", CD19="green", CD4="blue", CD8='yellow')

# Assign labels of dendrogram object with new colors:
labels_colors(dendro) <- colourCodes[membership_file_cleaned$cell_type][order.dendrogram(dendro)]
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

# Plot correlation between samples in a heatmap, TO DO:
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
#plot_heatmap_correlations <- ('plot_heatmap_correlations_normalised_expressed.png')
#png(plot_heatmap_correlations, width = 12, height = 12, units = 'in', res = 300)
#par(mfrow=c(1,2))
#heatmap(normalised_expressed$E[,1:10])
#heatmap(normalised_expressed$E[1:100,1:100], ColSideColors=colourCodes)
# p1 <- ggplot(cor_normalised_expressed, aes(Var1, Var2, fill = value)) + geom_tile() + scale_fill_gradient2(low = 'blue', high = 'yellow')
# grid.arrange(p1, ncol=1)
#par(mfrow=c(1,1))
#dev.off()

# # Try again with heatmap.2:
# # http://sebastianraschka.com/Articles/heatmaps_in_r.html
# 
# # Create colour palette from red to green:
# my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)
# 
# # (optional) defines the color breaks manually for a "skewed" color transition:
# col_breaks = c(seq(-1,0,length=100),  # for red
#                seq(0,0.8,length=100),              # for yellow
#                seq(0.8,1,length=100))              # for green
# 
# # creates a 5 x 5 inch image
# png('heatmap_normalised_expressed_heatmap2.png',    # create PNG for the heat map        
#     width = 5*300,        # 5 x 300 pixels
#     height = 5*300,
#     res = 300,            # 300 pixels per inch
#     pointsize = 8)        # smaller font size
# 
# heatmap.2(normalised_expressed$E[1:100, 1:100], 
#           cellnote = normalised_expressed$E[1:100, 1:100],  # same data set for cell labels
#           main = "Correlation", # heat map title
#           notecol="black",      # change font color of cell labels to black
#           density.info="none",  # turns off density plot inside color legend
#           #          trace="none",         # turns off trace lines inside the heat map
#           margins =c(12,9),     # widens margins around plot
#           col=my_palette,       # use on color palette defined earlier 
#           #          breaks=col_breaks,    # enable color transition at specified limits
#           dendrogram="column",     # only draw a column dendrogram
#           #          Colv="NA")            # turn off column clustering
# )
# dev.off()               # close the PNG device


# hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# t <- ExpressionSet(e, AnnotatedDataFrame(tab))
# rv <- rowVars(exprs(t))
# idx <- order(-rv)[1:40]
# heatmap(exprs(t)[idx, ], col = hmcol)


#############################


#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image_full, compress='gzip')

objects_to_save <- (c('normalised_expressed', 'membership_file_cleaned', 'FAILED_QC_unique'))
save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for differential gene expression.
#############################
