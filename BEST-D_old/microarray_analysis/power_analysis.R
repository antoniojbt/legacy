#############################
# Power analysis of gene expression microarray data
# Antonio J Berlanga-Taylor
# 02 Aug 2016
# BEST-D project differential expression power analysis
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_3.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/power_analysis.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_power_gene_expression",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific. Check ways of making count comparisons.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
#load('R_session_saved_image_probe_filtering.RData', verbose=T)
load('R_session_saved_image_pheno_file_check.RData', verbose=T)
# load('R_session_saved_image_diff_expression_full.RData', verbose = T)
#load('R_session_saved_image_diff_expression_3.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_power_gene_expression', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('sizepower')

library(limma)
library(ggplot2)
library(SSPA)
library(pwr)
library(sizepower)
#############################

#############################
# See:
# http://www.ncbi.nlm.nih.gov/pubmed/15533245
#############################

#############################
# Power analysis for unpaired t-test using pwr package and the core stats power functions
# http://rcompanion.org/rcompanion/b_02.html
# http://www.statmethods.net/stats/power.html
# ?power.t.test
# ?pwr.t.test

# For delta/effect size use log2 fold-change (Wei, 2004, BMC Genomics)
log(0.5, 2) # This equals 
log(2, 2) # This equals a 2 fold-change
2^0.6
2^0.32
2^0.14
2^4

# Estimate sample size needed:
power.t.test(
  n = NULL,                  # Observations in _each_ group
  d = 1,                 # This is delta, not effect size as such, use mean of the log2 fold change (this should correspond to Cohen's d)
  sig.level = 0.05,         # Type I probability
  power = 0.90,             # 1 minus Type II probability
  type = "two.sample",      # Change for one- or two-sample
  alternative = "two.sided",
  sd = 0.7                  # Include standard deviation, this is arbitrary
)

# Estimate for effect size detectable when sample size is fixed: 
power.t.test(
  n = 100,                  # Observations in _each_ group
  d = NULL,                 # This is delta, not effect size as such, use mean of the log2 fold change (this should correspond to Cohen's d)
  sig.level = 0.05,         # Type I probability
  power = 0.90,             # 1 minus Type II probability
  type = "two.sample",      # Change for one- or two-sample
  alternative = "two.sided",
  sd = 0.7                  # Include standard deviation, this is arbitrary
)
#############################

#############################
# Plot sample size curves for detecting correlations of various sizes
# http://www.evolutionarystatistics.org/document.pdf
# Range of effect sizes:
png('power_calc_curves.png')
deltas <- c(0.14, 0.32, 0.585, 1)
nvals <- seq(2, 200, 5)
plot(nvals, seq(0, 1, length.out = length(nvals)), xlab = 'sample size', ylab = 'power', 
     main = 'Power curves for \nt-test with varying effect size', type = 'n')
for (i in 1:length(deltas)) {
  powvals <- sapply(nvals, function (x) power.t.test(n = x, delta = deltas[i], 
                                                     sig.level = 0.05, type = "two.sample", alternative = "two.sided", sd = 0.7)$power)
  lines(nvals, powvals, lwd = 2, col = i)
}
legend('bottomright', lwd = 2, col = 1:length(deltas), legend = paste(round(2^deltas, 2)))
dev.off()

#############################

#############################
# Power calculations with pwr package (as base power() but with many extended functions)

# pwr package uses d = effect size, such as Cohen's d standardised effect size where SD is taken into account:
# Calculate standardised effect size:
# M1  = 10                      # Mean for sample 1
# M2  = 9                      # Mean for sample 2
# S1  =  0.7                      # Std dev for sample 1
# S2  =  0.7                      # Std dev for sample 2

# Cohen.d = (M1 - M2)/sqrt(((S1^2) + (S2^2))/2) 
# Cohen.d

# Estimate sample size needed (using pwr package):
pwr.t.test(
  n = NULL,                  # Observations in _each_ group
  d = 1,                 # Effect size, use mean of the log2 fold change (this should correspond to Cohen's d)
  sig.level = 0.05,         # Type I probability
  power = 0.90,             # 1 minus Type II probability
  type = "two.sample",      # Change for one- or two-sample
  alternative = "two.sided"
)

# Estimate for effect size detectable when sample size is fixed: 
pwr.t.test(
  n = 100,                  # Observations in _each_ group
  d = NULL,                 # Effect size
  sig.level = 0.05,         # Type I probability
  power = 0.90,             # 1 minus Type II probability
  type = "two.sample",      # Change for one- or two-sample
  alternative = "two.sided"
)


# Paired t-test power calculation for effect size when sample size is known:
pwr.t.test(
  n = 100,                  # Observations in _each_ group
  d = NULL,                 # Effect size
  sig.level = 0.05,         # Type I probability
  power = 0.90,             # 1 minus Type II probability
  type = "paired",      # Change for one- or two-sample
  alternative = "two.sided"
)

#############################


#############################
# sizepower library
# https://bioconductor.org/packages/release/bioc/html/sizepower.html
# Consider a randomized treatment-control design with equal samples per group:
power.multi(
  ER0 = 1, # Mean number of false positives to control
  G0 = 13000, # Number of genes not expected to be differentially expressed
  numTrt = 3, # Number of treatment conditions
  absMu1 =  2, # effect size: difference in expression between treatments in log-intensity scale
  sigma = 0.4, # standard error
  n = 100) # sample size per group

# Results are the estimated power and the non-centrality parameter 

# power.matched()
#############################


#############################
# Use the SSPA library for pilot data estimates

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))

head(membership_file_cleaned)
tail(membership_file_cleaned)
str(membership_file_cleaned)

dim(normalised_filtered)
normalised_filtered[1:5, 1:5]

length(which(complete.cases(normalised_filtered)))
dim(normalised_filtered)

#############################

#############################
# Get statistics (moderated t-tests from limma) for every contrast performed:
head(fit2_groups_before_v_after_time$t)
t_tests_df <- as.data.frame(fit2_groups_before_v_after_time$t)
head(t_tests_df)
dim(t_tests_df)
# View(t_tests_df)
# Total sample size:
total_sample_size <- nrow(design) 
#############################

#############################
# SSPA power calculation use a pilot data set of experimental data
# effectsize calculation
# See: 
# http://bioconductor.org/packages/release/bioc/html/SSPA.html
# http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-10-439
# Functions are chunks from their scripts.
# ?pilotData
pilot_data <- apply(t_tests_df, 2, 
            function(x) pilotData(statistics = x[-1], 
                                  samplesize = total_sample_size,
                                  distribution = "norm"))
head(pilot_data)
class(pilot_data)
str(pilot_data)

# For each column (each set of t-stats for each contrast), plot:
# See: http://bioconductor.org/packages/release/bioc/vignettes/SSPA/inst/doc/SSPA.pdf
# Figure: Exploratory  plots  of  the  pilot-data: 
# Upper  left  panel  shows  a  histogram  of  the  test  statistics,  
# upper  right panel  a  histogram  of  the  p-values.
# Lower  left  panel  shows  a  qq-plot  for  the  test  statistics  using  the  user 
# defined  null distributions.
# Lower right panel shows the sorted p-values against their ranks

# TO DO: plot for arbitrary number of columns/contrasts and insert name
names(pilot_data)
png('power_exploratory_plots.png', width = 12, height = 12, units = 'in', res = 300)
plot(pilot_data[[1]], main = sprintf('%s', names(pilot_data)[1]))
plot(pilot_data[[2]])
plot(pilot_data[[3]])
plot(pilot_data[[4]])
plot(pilot_data[[5]])
dev.off()

qqnorm(pilot_data$UI4000minusplacebo@pvalues)
qqnorm(pilot_data$UI4000minusplacebo@statistics)
qqline(pilot_data$UI4000minusplacebo@statistics)

hist(pilot_data$UI4000minusplacebo@pvalues)
hist(pilot_data$UI2000minusplacebo@pvalues)

# Calculate sample size for each contrast using sampleSize from SSPA library:
sample_size <- lapply(pilot_data, sampleSize, 
                      control = list(pi0Method = "Storey", a = 0, 
                                     resolution = 2^10, verbose = TRUE))
head(sample_size)
names(sample_size)
class(sample_size)

# TO DO: check this is correct:
# Estimated density of effect sizes for the given contrast/condition. True density of effect sizes
# in dotted lines and estimated density of effect sizes in blue.
# plot(sample_size$UI4000minusplacebo)
plot(sample_size$UI4000minusplacebo, panel = function(x, y, ...) {
  panel.xyplot(x, y)
  panel.curve(dbitri(x), lwd = 2, lty = 2, n = 500)
  }, ylim = c(0, 0.6))

# TO DO: plot for each condition:
contrasts <- c("UI4000minusplacebo", "UI2000minusplacebo", "UI4000minus2000", 
               "UI4000minus2000minusplacebo", "UI2000minus4000minusplacebo")


# Estimate the average power for sample sizes other than those of the pilot data provided:
# ?predictpower
Jpred <- seq(10, 300, by = 10)
N <- sqrt(Jpred / 2)
pwrD <- predictpower(ss$UI4000minusplacebo, samplesizes = N, alpha = 0.05, plot = F)
matplot(Jpred, pwrD, type = 'b', pch = 16, ylim = c(0,1),
        ylab = 'Predicted power', xlab = 'Sample size (per group)')
grid()
#############################


#############################
# Web based power calculation:
# http://bioinformatics.mdanderson.org/MicroarraySampleSize/
# Underlying Model
# These sample size computations are based on the assumption that the expression of each gene is 
# normally distributed on the log scale. 
# We assume further that the variance is the same in the two experimental groups, 
# and thus the only difference between the groups is the mean. 
# The acceptable number of false positives is divided by the number of genes to compute a per-gene 
# significance level alpha that would be expected to produce that number of false positives. 
# (As a consequence, asking for 0 false positives will require infinitely many samples.) 
# Finally, we pretend that gene expression measurements are independent, 
# and perform the usual computations to determine the sample size needed 
# for a t-test with an underlying normal distribution. 


#############################

#############################
# TO DO:
# Simulate vs data and pass through pipeline
# http://bioconductor.org/packages/release/bioc/vignettes/vsn/inst/doc/convergence2.pdf

# Run power calculations on post-hoc values and simulated diff. exp. values with SSPA library

#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')
sessionInfo()

q()

# Next: run script for xxx.
#############################