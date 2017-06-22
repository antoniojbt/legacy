###############
# 018 metabolomics
# 16 March 2016
# Metabolomic basis of chronic inflammation and NCD risk factors - differential expression analysis
###############

###############
# To DO's:
#

###############


#############################
##Set working directory and file locations and names of required inputs:
options(echo=TRUE)
# Working directory:
# setwd("/Users/antoniob/Documents/CGAT/Projects_and_work/Pseudomonas-and-CF/018_analysis_to_upload.dir/")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_018_metabolomics",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is very project specific. Check ways of making count comparisons.

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_Airwave_metabolomics.RData', verbose=T) 

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_018_metabolomics.RData')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite('pcaMethods')
# install.packages('pcaMethods')

library(data.table)
library(ggplot2)
library(gridExtra)
library(Hmisc)
library(splines)
library(plyr)
library(metabolomics)
library(limma)

# Get script with functions needed:
# source('functions_for_metabolomics.R')
source('moveme.R')
#############################

#############################################
# Run with command line arguments:
args <- commandArgs(trailingOnly = TRUE)

#  Set variables:
# TO DO: pass to configuration file
phenotype_file <- as.character(args[1])
# phenotype_file <- ''

metabolomics_file <- as.chracter(args[2])
# metabolomics_file <- 'NMR_test_set.txt'

print(args)
#############################################


#############################
## Load data sets:


# Read in metabolomics data:
metabolomics_data <- fread(metabolomics_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE) #check if first cell is not empty
head(names(metabolomics_data))
metabolomics_data[1:5, 1:5, with = F]
dim(metabolomics_data)
setkey(metabolomics_data, 'isolate') # data.table doesn't take column numbers
tables()

##################


##################
# Subset:
names(metabolomics_data)
cols_keep <- c(2, 14:53)
count(metabolomics_data[, 'stage', with = F])
pheno_data <- metabolomics_data[, c(1:13), with = F]
pheno_data
dim(pheno_data)

metabolomics_data <- metabolomics_data[, cols_keep, with = F]
metabolomics_data
dim(metabolomics_data)
# row.names(metabolomics_data) <- metabolomics_data[, 1, with = F]


metabolomics_data_t <- transpose(metabolomics_data)
metabolomics_data_t <- as.data.frame(metabolomics_data_t)
# colnames(metabolomics_data_t) <- metabolomics_data_t[1, ]
head(metabolomics_data_t)
dim(metabolomics_data_t)
metabolomics_data_t <- metabolomics_data_t[-1, ]
colnames(metabolomics_data_t)
head(metabolomics_data_t)
dim(metabolomics_data_t)
str(metabolomics_data_t)
metabolomics_data_t <- sapply(metabolomics_data_t, as.numeric)
metabolomics_data_t <- as.data.frame(metabolomics_data_t)
class(metabolomics_data_t)
head(metabolomics_data_t)
str(metabolomics_data_t)
dim(metabolomics_data_t)
dim(metabolomics_data)

# early <- metabolomics_data[which(metabolomics_data[, stage] == 'early'), cols_keep, with = F]
# late <- metabolomics_data[which(metabolomics_data[, stage] == 'late'), cols_keep, with = F]
# head(early)
# dim(early)
# head(late)
# dim(late)

# Order files so that rows and columns match, both files have rows as sampe IDs and columns as variables:
# pheno_with_metabolomics <- pheno_with_metabolomics[order(pheno_with_metabolomics$BARCODE), ]
# metabolomics_data <- metabolomics_data[order(metabolomics_data[, Row]), ]
# 
# # row.names(pheno_with_metabolomics) <- pheno_with_metabolomics$BARCODE
# 
# # Sanity check:
# # TO DO: Raise error and stop if false:
# identical(pheno_with_metabolomics$BARCODE, metabolomics_data[, Row])
# length(which(pheno_with_metabolomics$BARCODE %in% metabolomics_data[, Row]))
# 
# # Get transposed file as needed downstream:
# metabolomics_data_t <- transpose(metabolomics_data)
##################

##################
# TO DO: continue from here, long computations...
# Get PCs for metabolomics data:
#Plot probes in a multi-dimensional scaling plot:
plot_MDS <- ("MDS_metabolomics.png")
png(plot_MDS, width = 4, height = 4, units = 'in', res = 300)
plotMDS(metabolomics_data[, -1, with = F], pch = 1)
dev.off()

plot_MDS_by_targets <- ("MDS_metabolomics_by_targets.png")
png(plot_MDS_by_targets, width = 4, height = 4, units = 'in', res = 300)
plotMDS(metabolomics_data[, -1, with = F], pch = 1, labels = pheno_data$stage)
dev.off()

dim(metabolomics_data[, -1, with = F])
length(pheno_data[, stage])

# Plot PCA of normalised samples:
# Compute the PCs, use transposed expression values:
dim(metabolomics_data_t)
pca_values <- prcomp(metabolomics_data_t, center=TRUE, scale=TRUE)
head(pca_values$x)
tail(pca_values$x)
dim(pca_values$x)

# Obtain values for all PCs output:
pc <- data.frame(round(pca_values$x, 2))
dim(pc)
head(pc)
pc$sample_id <- row.names(pc)
head(pc$sample_id)

pc <- pc[, moveme(names(pc), 'sample_id first')]
names(pc)[1:10]
class(pc)
write.table(pc, 'principal_components_metabolomics.tsv', quote = FALSE, sep = '\t', row.names = FALSE)

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(metabolomics_data_t)
str(pca_values)
head(pc)
# pc[1:5, 1:5]
summary(pca_values)

# Plot PCA results:
plot_PCA <- ('plot_PCA.png')
png(plot_PCA, width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(2,2))
# Histogram of first x PCs:
plot(pca_values, main='Normalised metabolomics values')
# Scatterplot of PC1 and PC2:
biplot(pca_values, main='Normalised metabolomics values')
par(mfrow=c(1,1))
dev.off()

# Check how much of the variance in gene expression values is explained by the first x PCs:
# sum(pca_normalised_filtered$sdev[1:10]^2)/length(normalised_filtered[1,])

# Run PCA analysis by groups of interest: TO DO: Cross files and IDs first
pc_data_top_20 <- data.frame(pc[, 1:20])
str(pc_data_top_20)
head(pc_data_top_20)
row.names(pc_data_top_20)

pca_by_groups <- data.frame(merge(pheno_data, pc_data, by='row.names'))
head(pca_by_groups)
dim(pca_by_groups)
dim(pc_data)

head(arrange(pc_data, PC1), 10)
head(arrange(pca_by_groups, PC1), 10)

plot_PCA_by_groups_1 <- ('plot_PCA_by_groups_1.png')
png(plot_PCA_by_groups_1, width = 13, height = 13, units = 'in', res = 300)
p1 <- qplot(x=PC1, y=PC2, data=pc, colour=factor(pheno_data$stage)) + theme(legend.position="bottom")
p2 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(pheno_with_metabolomics$TAKES_INTENSE_EXERCISE)) + theme(legend.position="bottom")
p3 <- qplot(x=PC1, y=PC2, data=pca_by_groups, colour=factor(pheno_with_metabolomics$IS_SMOKER)) + theme(legend.position="bottom")
# p4 <- qplot(x=PC2, y=PC3, data=pca_by_groups, colour=factor(pheno_with_metabolomics$)) + theme(legend.position="bottom")
# grid.arrange(p1, p2, p3, p4, ncol=2)
grid.arrange(p1, p2, p3, ncol=2)
dev.off()


##################

##################
# Exploratory analyses of metabolomics data:

# Some stats per sample:
class(metabolomics_data) # data.table object so indexing is DT[i, j, by] with cols requiring names not numbers
range(metabolomics_data[3, ])
range(metabolomics_data[100, ])
range(metabolomics_data[1000, ])

# Some stats per metabolite:
summary(metabolomics_data[, list(Plasma_CPMG_NMR_0_500066)])
summary(metabolomics_data[, 100, with = FALSE])
summary(metabolomics_data[, 1000, with = FALSE])

# TO DO: errors, need to get correct input, see above, merging pbs.
# Some plots per metabolite:
# png('Plasma_CPMG_NMR_0_500066_boxplot.png', width = 6, height = 6, units = 'in', res = 300)
# MetBoxPlots(inputdata = as.data.frame(metabolomics_data), metname = 'Plasma_CPMG_NMR_0_500066', main = 'Plasma_CPMG_NMR_0_500066')
# dev.off()

png('densities_metabolomics.png', width = 6, height = 6, units = 'in', res = 300)
limma::plotDensities(metabolomics_data_t, legend = F)
dev.off()

png('MA_plot_metabolomics.png', width = 6, height = 6, units = 'in', res = 300)
limma::plotMA(as.matrix(metabolomics_data[, -1, with = F]))
# limma::plotMA(as.matrix(metabolomics_data[, c(2,5), with = F]))
dev.off()


#############################

#############################

# 6) Differential gene expression analysis
# Linear modelling using weights
# b) Two group comparison of treated vs untreated: 

# Check NAs:

# Compare samples based on groups:
# Limma will use the first group defined as the reference group:
pheno_data$stage <- factor(pheno_data$stage, levels=c('early', 'late'))

group <- pheno_data$stage
count(group)

#Define design:
design_by_group <- model.matrix(~group)
head(design_by_group)
tail(design_by_group)
count(design_by_group)

#Run linear model and obtain differentially expressed metabolites (?) based on all pairs:
head(colnames(metabolomics_data_t))
head(row.names(metabolomics_data_t))

fit_by_group <- lmFit(metabolomics_data_t, design_by_group)
fit_by_group
names(fit_by_group)
colnames(fit_by_group)
head(row.names(fit_by_group))

fit_by_group_2 <- eBayes(fit_by_group)
fit_by_group_2
head(fit_by_group_2$coefficients)

topTable(fit_by_group_2, adjust = 'BH')
topTable(fit_by_group_2, adjust = 'BH', coef = 'grouplate')

topTable_groups_bmi <- topTable(fit_by_group_2, adjust = 'BH')
topTable_groups_v_less_20 <- topTable(fit_by_group_2, adjust = 'BH', coef = 'groupBMI_less_20', number = Inf)

head(x = topTable_groups_v_less_20)
count(topTable_groups_v_less_20$adj.P.Val < 10e-3)
count(topTable_groups_v_25_30$adj.P.Val < 10e-3)
count(topTable_groups_v_over_30$adj.P.Val < 10e-3)

# count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) > 1.49)
# count(topTable_groups$adj.P.Val < 0.05 & 2^(topTable_groups$logFC) < 0.68)
# count(topTable_groups$adj.P.Val < 0.05 & topTable_groups$logFC > abs(0.57))

colnames(fit_by_group_2)
results_by_groups <- decideTests(fit_by_group_2)
vennDiagram(results_by_groups)
vennDiagram(decideTests(fit_by_group_2[, 3]))

png('volcano_plots_bmi_nmr.png', width = 12, height = 12, units = 'in', res = 300)
par(mfrow=c(3,1))
plot(2^(topTable_groups_v_less_20$logFC), -log10(topTable_groups_v_less_20$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI under 20", abline(h=2, v=c(0.8, 1.2)))
plot(2^(topTable_groups_v_25_30$logFC), -log10(topTable_groups_v_25_30$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI 25 - 30", abline(h=2, v=c(0.8, 1.2)))
plot(2^(topTable_groups_v_over_30$logFC), -log10(topTable_groups_v_over_30$adj.P.Val), pch=20, main="Volcano plot: Lean vs BMI over 30", abline(h=2, v=c(0.8, 1.2)))
# volcanoplot(fit_by_group_2, coef = 4)
# volcanoplot(fit_by_group_2, coef = 2)
# volcanoplot(fit_by_group_2, coef = 3)
par(mfrow=c(1,1))
dev.off()

# Interpretation: Many significant differences between groups for BMI. Lots to check first.
# TO DO: Get random subsets of individuals and run diff. abundance analysis, no differences expected.
#############################

#############################
# Write files to disk:
write.table(x=topTable_groups, sep='\t', quote = FALSE, col.names = NA, row.names = TRUE, 
            file='full_topTable_groups.txt')


#############################

#############################
# Question
# Background
# Test
# Result
# Plot
# Interpretation
# Next step
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for xxx.
#############################