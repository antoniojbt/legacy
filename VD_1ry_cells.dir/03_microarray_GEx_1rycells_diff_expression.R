#############################
# To be run after 02 normalisation of array data

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
#working_dir <- ("/ifs/projects/proj043/analysis.dir/gene_expression.dir")


#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression",".txt", sep=""), open='a')
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
#load('R_session_saved_image_read_and_QC.RData', verbose=T)
load('R_session_saved_image_normalisation.RData', verbose=T)

# To load multiple .RData files:
#rdata_filenames <- c('.RData')
#lapply(rdata_filenames, load, .GlobalEnv)


# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_diff_expression_full', '.RData', sep='')


#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#biocLite("limma")
#install.packages('dendextend')

library(limma)
library(reshape2)
library(gridExtra)
library(plyr)
library(ggplot2)
library(Hmisc)

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
head(membership_file_cleaned)
tail(membership_file_cleaned)
membership_file_cleaned$sample

# Check names/indices match between meta-data file and array file:
head(membership_file_cleaned)
membership_file_cleaned$sample
colnames(normalised_expressed)

dim(normalised_expressed)
normalised_expressed[1:5, 1:5]
class(normalised_expressed)
str(normalised_expressed)
names(normalised_expressed)
head(normalised_expressed$E)

# Make sample IDs row names:
row.names(membership_file_cleaned) <- membership_file_cleaned$sample
row.names(membership_file_cleaned)
colnames(membership_file_cleaned)
head(membership_file_cleaned)

# Re-order so that these match:
identical(row.names(membership_file_cleaned), colnames(normalised_expressed))

normalised_expressed_ordered <- normalised_expressed[, order(colnames(normalised_expressed))]
class(normalised_expressed_ordered)
names(normalised_expressed_ordered)
colnames(normalised_expressed_ordered)

membership_file_cleaned_ordered <- membership_file_cleaned[order(row.names(membership_file_cleaned)), ]

identical(row.names(membership_file_cleaned_ordered), colnames(normalised_expressed_ordered))

# Rename:
membership_file_cleaned <- membership_file_cleaned_ordered
normalised_expressed <- normalised_expressed_ordered

range(normalised_expressed$E[, 1])
range(normalised_expressed$E[, 10])
range(normalised_expressed$E[, 5])
range(normalised_expressed$E[, 17])
range(normalised_expressed$E[1, ])
range(normalised_expressed$E[10, ])
range(normalised_expressed$E[100, ])
range(normalised_expressed$E[1000, ])

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

# Help ?lm()?eBayes

# For VD primary cells:
# individual
# treatment
# timepoint
# replicate
# cell_type
# sentrix_ID



## a) Design:

## Test: 
#Define factors to constrast from annotation file:

individual <- factor(membership_file_cleaned$individual, levels=c('1', '2'))
treatment <- factor(membership_file_cleaned$treatment, levels=c('0', '1'))
timepoint <- factor(membership_file_cleaned$timepoint, levels=c('0', '2', '8', '24'))
cell_type <- factor(membership_file_cleaned$cell_type, levels=c('CD14', 'CD19', 'CD4', 'CD8'))
replicate <- factor(membership_file_cleaned$replicate, levels=c('1', '2'))
sentrix_ID <- factor(membership_file_cleaned$sentrix_ID, levels=c('3999927008', '3999927011'))

#Define design:
design <- model.matrix(~0+individual)
head(design)
tail(design)
length(design)

#Run linear model and set contrasts:
fit <- lmFit(normalised_expressed, design)
cont.matrix <- makeContrasts(Ind1_vs_Ind2=individual1-individual2, levels=design)

#Obtain differentially expressed genes based on contrasted factors:
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
topTable_individuals <- topTable(fit2, adjust="BH", n=Inf)
head(topTable_individuals)

#Get results and plot:
results <- decideTests(fit2)
vennDiagram(results)

head(topTable_individuals)
head(topTable_individuals$adj.P.Val)

volcanoplot(fit2, coef=1, highlight=0, names=fit2$genes$ID, xlab="Log Fold Change", ylab="Log Odds", pch=16, cex=0.35, abline(h=2))

# Interpretation: No difference between individuals.

##################
# b) Two group comparison: 

#Compare samples based on treatment:
#Define design and set contrasts:
design_by_treatment <- model.matrix(~treatment)
head(design_by_treatment)
tail(design_by_treatment)
count(design_by_treatment)
count(membership_file_cleaned$treatment)
design_by_treatment
membership_file_cleaned$description

#Run linear model and obtain differentially expressed genes based on all pairs:
fit_by_treatment <- lmFit(normalised_expressed, design_by_treatment)
fit_by_treatment
names(fit_by_treatment)
colnames(fit_by_treatment)
head(row.names(fit_by_treatment))

#contrast_matrix_by_treatment <- makeContrasts(1-0, levels=design_by_treatment)
#head(contrast_matrix_by_treatment)

fit_by_treatment_2 <- eBayes(fit_by_treatment)
fit_by_treatment_2
head(fit_by_treatment_2$coefficients)

topTable(fit_by_treatment_2, adjust = 'BH')

# Interpretation:
# TO DO: CLEAN UP, check limma examples, straight up stimulated vs basal is negative but PCAs separate cell types and treatments.

# topTable_treatment <- topTable(fit_by_treatment_2, coef='1 - 0', adjust='BH', number=Inf)
# head(topTable_treatment)
# range(topTable_treatment$logFC)
# range(topTable_treatment$AveExpr)
# range(topTable_treatment$adj.P.Val)
# range(topTable_treatment$B)
# dim(topTable_treatment)
# count(topTable_treatment$adj.P.Val < 0.00000000000000000001)

#Plot results:
results_by_treatment <- decideTests(fit_by_treatment_2)
vennDiagram(results_by_treatment)

volcanoplot(fit_by_treatment_2)

#, coef='1 - 0', highlight=10, names=fit_by_treatment_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2))

# Interpreation: Broad comparison of baseline vs stimulated gives large, signifcant changes. 

##########
# TO DO: continue from here, define questions and comparisons of interest though.

# Factorial design:
# Compare samples based on treatment:
# Define design and set contrasts:

#Check experimental dimensions (factors):
head(membership_file_cleaned)
membership_file_cleaned$replicate
count(membership_file_cleaned$cell_type)
count(membership_file_cleaned$treatment)
count(membership_file_cleaned$timepoint)

#Modify to leave only 0 and 2 hour timepoints:
timepoint_8h <- which(membership_file_cleaned$timepoint == 8)
timepoint_24h <- which(membership_file_cleaned$timepoint == 24)
row.names(membership_file_cleaned[c(timepoint_8h, timepoint_24h), ])
sample_names_exclude <- membership_file_cleaned$sample[c(timepoint_8h, timepoint_24h)]
sample_names_exclude

membership_file_cleaned_0h_2h <- membership_file_cleaned[-c(timepoint_8h, timepoint_24h), ]
dim(membership_file_cleaned_0h_2h)
dim(membership_file_cleaned)
which(membership_file_cleaned_0h_2h$timepoint == 8)
which(membership_file_cleaned_0h_2h$timepoint == 24)

# Also from array data:
colnames(normalised_expressed)
membership_file_cleaned$sample
row.names(membership_file_cleaned)

#arrays_exclude <- which(colnames(normalised_expressed) %in% sample_names_exclude)

#normalised_expressed_0h_2h <- normalised_expressed[, -arrays_exclude]
normalised_expressed_0h_2h <- normalised_expressed[, -c(timepoint_8h, timepoint_24h)]
dim(normalised_expressed)
dim(normalised_expressed_0h_2h)
sample_names_exclude %in% colnames(normalised_expressed_0h_2h)

#Collect different combinations of factors:
group_combinations <- paste(membership_file_cleaned_0h_2h$cell_type, 
                            membership_file_cleaned_0h_2h$treatment, sep='.')
group_combinations
count(group_combinations)
length(group_combinations)

# Define the questions of interest:
# Before and after comparisons for three groups:
# Which genes respond to treatment after stimulation at high dose: 4000IU Final Visit vs 4000IU randomisation (baseline)
#arm.visit_type = 0.FinalVisit vs 0.Randomisation

# Which genes respond to treatment after stimulation at low dose: 2000IU Final Visit vs 2000IU randomisation (baseline)
#arm.visit_type = 1.FinalVisit vs 1.Randomisation

# Which genes are likely due to noise:  Placebo Final Visit vs Placebo randomisation (baseline)
#arm.visit_type = 2.FinalVisit vs 2.Randomisation

# Two group comparison see above:


# Compare the difference of differences (ie interaction terms): 
# Differences due to high dose, ie 4000IU minus Placebo: (0.FinalVisit vs 0.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to low dose, ie 2000IU minus Placebo: (1.FinalVisit vs 1.Randomisation) - (2.FinalVisit vs 2.Randomisation)
# Differences due to high dose only, ie 4000IU minus 2000IU: (0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation) 
# Differences truly due to high dose only, ie (4000IU minus 2000IU) minus (Placebo):
#((0.FinalVisit vs 0.Randomisation) - (1.FinalVisit vs 1.Randomisation)) - (2.FinalVisit vs 2.Randomisation)


#Define the combinations of factors:
factorial_design <- factor(group_combinations)#, levels=c('CD14.O', 'CD14.1', 
#                                                         'CD8.0', 'CD8.1', 
#                                                         'CD4.0', 'CD4.1',
#                                                         'CD19.0', 'CD19.1'))
head(factorial_design)
head(membership_file_cleaned_0h_2h$description)
head(membership_file_cleaned_0h_2h$sample)
head(colnames(normalised_expressed_0h_2h))

levels(factorial_design)
length(factorial_design)

#Extract the comparisons of interest as contrasts:
design_factorial_comparisons <- model.matrix(~0+factorial_design)
head(design_factorial_comparisons)
colnames(design_factorial_comparisons)
dim(design_factorial_comparisons)
dim(normalised_expressed_0h_2h)

#Fit a model with a coefficient for each of the factor combinations:
fit_factorial <- lmFit(normalised_expressed_0h_2h, design_factorial_comparisons)
fit_factorial

#Extract the fitted contrasts of interest:
contrast_matrix_factorial <- makeContrasts(CD14VD=factorial_designCD14.1-factorial_designCD14.0,
                                           CD8VD=factorial_designCD8.1-factorial_designCD8.0, 
                                           CD4VD=factorial_designCD4.1-factorial_designCD4.0,
                                           CD19VD=factorial_designCD19.1-factorial_designCD19.0, 
                                           levels=design_factorial_comparisons)

contrast_matrix_factorial

fit_factorial_2 <- contrasts.fit(fit_factorial, contrast_matrix_factorial)
fit_factorial_2 <- eBayes(fit_factorial_2)
colnames(fit_factorial_2$coefficients)
range(fit_factorial_2$lods)
range(fit_factorial_2$p.value)

topTable(fit_factorial_2, adjust = 'BH')
topTable(fit_factorial_2, adjust = 'BH', coef='CD14VD')
topTable(fit_factorial_2, adjust = 'BH', coef='CD4VD')
topTable(fit_factorial_2, adjust = 'BH', coef='CD8VD')
topTable(fit_factorial_2, adjust = 'BH', coef='CD19VD')


topTable_fit_factorial_2_CD14 <- topTable(fit_factorial_2, coef='CD14VD', number=Inf)
head(topTable_fit_factorial_2_CD14)
range(topTable_fit_factorial_2_CD14$logFC)
range(topTable_fit_factorial_2_CD14$AveExpr)
vdr_index <- which(topTable_fit_factorial_2_CD14$SYMBOL == 'VDR')
vdr_index
topTable_fit_factorial_2_CD14[vdr_index, ]
genes_interest <- c('VDR', 'CYP27B1', 'CYP24A1', 'GC', 'VDBP', 'DHCR7', 'PTH', 'FGF23', 'CAMP')

genes_interest_index <- which(topTable_fit_factorial_2_CD14$SYMBOL %in% genes_interest)
genes_interest_index
topTable_fit_factorial_2_CD14[genes_interest_index, ]

topTable_treatment[which(topTable_treatment$SYMBOL %in% genes_interest), ]


results_by_factorial <- decideTests(fit_factorial_2)
#vennDiagram(results_by_factorial)

# Interpretation: Cell-specific comparisons between basal vs stimulated are negative.

# volcano_plots <- ('volcano_plots.png')
# png(volcano_plots, width = 12, height = 12, units = 'in', res = 300)
# par(mfrow=c(3,3))
# 
# volcanoplot(fit_factorial_2, coef='High_dose_before_and_after', highlight=10, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_before_and_after')
# 
# volcanoplot(fit_factorial_2, coef='Low_dose_before_and_after', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_before_and_after')
# 
# volcanoplot(fit_factorial_2, coef='Placebo_before_and_after', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Placebo_before_and_after')
# 
# volcanoplot(fit_factorial_2, coef='High_dose_minus_Placebo', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_Placebo')
# 
# volcanoplot(fit_factorial_2, coef='Low_dose_minus_Placebo', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Low_dose_minus_Placebo')
# 
# volcanoplot(fit_factorial_2, coef='High_dose_minus_low_dose', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='High_dose_minus_low_dose')
# 
# volcanoplot(fit_factorial_2, coef='Truly_high_dose_only', highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='Truly_high_dose_only')
# 
# par(mfrow=c(1,1))
# dev.off()

# volcanoplot(fit_factorial_2, coef=4, highlight=0, names=fit_by_arm_2$genes$ID, 
#             xlab='Log Fold Change', ylab='Log Odds', pch=16, cex=0.35, abline(h=2), main='')

# Plot volcano, p-values vs effect size:
# log10_p_values <- -log10(fit_factorial_2$p.value)
# plot(fit_factorial_2$lods, log10_p_values, pch='.', xlab='log-ratio', ylab=expression(-log[10]~p))
# abline(h=2)

#############################

##########
# KEEP ONLY CD14:
keep_CD14 <- which(membership_file_cleaned_0h_2h$cell_type == 'CD14')
keep_CD14
membership_file_cleaned_0h_2h$description
membership_file_cleaned_0h_2h_CD14 <- membership_file_cleaned_0h_2h[keep_CD14, ]
dim(membership_file_cleaned_0h_2h_CD14)
membership_file_cleaned_0h_2h_CD14$description
membership_file_cleaned_0h_2h_CD14$sample
colnames(normalised_expressed_0h_2h[, keep_CD14])

normalised_expressed_0h_2h_CD14 <- normalised_expressed_0h_2h[, keep_CD14]
dim(normalised_expressed_0h_2h_CD14)

#Collect different combinations of factors:
CD14_Tx <- factor(membership_file_cleaned_0h_2h_CD14$treatment)
CD14_Tx
count(CD14_Tx)

#Extract the comparisons of interest as contrasts:
design_CD14 <- model.matrix(~CD14_Tx)
design_CD14
colnames(design_CD14)
membership_file_cleaned_0h_2h_CD14$description
dim(design_CD14)
dim(normalised_expressed_0h_2h_CD14)

#Fit a model with a coefficient for each of the factor combinations:
fit_CD14 <- lmFit(normalised_expressed_0h_2h_CD14, design_CD14)
fit_CD14
fit_CD14_2 <- eBayes(fit_CD14)

topTable(fit_CD14_2, adjust = 'BH')
# Interpretation: Basal vs stimulated CD14 is negative.
################


#############################
#The end:
# TO DO: remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes)))

#rm()

# To save R workspace with all objects to use at a later time:
# Also consider dump(), dput() and dget(), see http://thomasleeper.com/Rcourse/Tutorials/savingdata.html
save.image(file=R_session_saved_image_full, compress='gzip')

# Or save specific objects:
# objects_to_save <- (c('normalised_expressed', 'membership_file_cleaned', 'FAILED_QC'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################