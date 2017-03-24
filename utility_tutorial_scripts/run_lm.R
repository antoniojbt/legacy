# #############################################
# Run linear regression models in R (ordinary least squares) for eQTL studies (candidate SNPs and probes)
# Antonio J Berlanga-Taylor
# 11 June 2016
# Script to test assumptions of linear regression models, run lm, diagnostic plots, test their fit and identify unusual
# observations. 
# Written to test SNPs and genes/probes with PCs as covariates for eQTL studies.
#############################################


#############################################
# See:
# R in Action, Chapters 8, 9, 13, 14.
  # Select model parameters (intercept and slopes) that minimize the difference between actual response values 
  # and those predicted by the model.
# Overview:
  # Define question
  # Prepare data with response, predictors, other covariates and confounders
  # Variable selection
  # Define lm formula and model to run
  # Run model(s)
  # Run diagnostic plots and check statistical assumptions:
    # Normality: For fixed values of the independent variables, the dependent variable is normally distributed.
    # Independence: The Yi values are independent of each other.
    # Linearity: The dependent variable is linearly related to the independent variables.
    # Homoscedasticity: The variance of the dependent variable doesn’t vary with the levels of the independent variables (constant variance).
  # Compare models
  # Exclude outliers, transform data, re-fit models.
  # Model interpretation
  # Rank predictors based on relative importance
  # Evaluate model on new data (test predictions of regression model)
#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/lm_testing.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting//')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_regression",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_regression','.RData', sep='')
R_session_saved_image
####################


####################
# Load packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite("h2o")
library(car) # Run correlations and plots for lm
library(gvlma) # Compare models and fit
library(leaps) # Compare all possible models
library(relaimpo) # Determine importance of variables as predictors
library(effects) # To visualise effects from interactions specified in regression model
library(bootstrap) # Run cross-validation
library(gmodels) # Create cross tables of categorical data
library(dummies)# functions to create dummy variables flexibly using model.matrix returning 
                # them as either matrices or data frames for further analysis. 
library(h2o) # computes parallel distributed machine learning algorithms 
              # https://github.com/h2oai/h2o-3/blob/master/h2o-docs/src/booklets/v2_2015/PDFs/online/R_Vignette.pdf
              # http://docs.h2o.ai
library(ggplot2)
library(data.table)
# TO DO: pass functions to separate script
# source('functions_run_lm.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

#  Variables and files:
lm_file <- as.character(args[1])
# lm_file <- '2000+4000-baseline-1.eQTL_trans'

covariates_file_1 <- as.character(args[2])
# covariates_file_1 <- '2000+4000-baseline-1.eQTL_cis'

covariates_file_2 <- as.character(args[3])
# covariates_file_2 <- ''

response_var <- as.character(args[4])

predictors_var <- as.character(args[5])

lm_formula <- as.character(args[6])

print(args)
####################


####################
# Set up file naming for outputs:
# to create: 2000+4000-baseline-1.trans_in_cis
lm_file_base <- strsplit(lm_file, '[.]')
lm_file_base <- lm_file_base[[1]][1]
lm_file_base
lm_file_ext <- strsplit(lm_file, '_')
lm_file_ext <- lm_file_ext[[1]][2]
lm_file_ext

covariates_file_1_base <- strsplit(covariates_file_1, '[.]')
covariates_file_1_base <- covariates_file_1_base[[1]][1]
covariates_file_1_base
covariates_file_1_ext <- strsplit(covariates_file_1, '_')
covariates_file_1_ext <- covariates_file_1_ext[[1]][2]
covariates_file_1_ext

output_file_name <- sprintf('%s_%s_%s.lm', lm_file, ,)
output_file_name

output_dir <- sprintf('%s_%s.dir', response_var, lm_formula)

# Set-up directory for results:
cmd_run <- sprintf('mkdir %s', output_dir)
system(cmd_run)
cmd_run <- sprintf('cd %s', output_dir)
system(cmd_run)
# Make soft links for data:
# TO DO:
cmd_run <- sprintf('ln -s ../%s .', lm_file)
system(cmd_run)
cmd_run <- sprintf('ln -s ../%s .', covariates_file_1)
system(cmd_run)
cmd_run <- sprintf('ln -s ../%s .', covariates_file_2)
system(cmd_run)
####################


####################
# Read data:
eQTL_data1 <- fread(eQTL_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data1, 'gene', 'Probe_ID')
dim(eQTL_data1)
colnames(eQTL_data1)
setkey(eQTL_data1, SNP, Probe_ID)

eQTL_data2 <- fread(eQTL_file2, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data2, 'gene', 'Probe_ID')
dim(eQTL_data2)
colnames(eQTL_data2)
setkey(eQTL_data2, SNP, Probe_ID)

tables()

eQTL_data1
eQTL_data2
####################


####################
# Variable selection
# TO DO
# Stepwise regression
# All subsets regression
####################


####################
# TO DO: add methods, interpretation, legends for plots and tests

####################


####################
# Run lm
# Examples:
# TO DO: combine example into one for flow and code from there. merge with eQTL_plotting.R
data(women)
fit <- lm(weight ~ height, data = women)
summary(fit)
fitted(fit)
residuals(fit)
plot(women$height,women$weight, 
     xlab="Height (in inches)",
     ylab="Weight (in pounds)")
abline(fit)
# Plot diagnostics
# Normality, Independence, Linearity, Homoscedasticity, Residual versus Leverage graph (outliers, high leverage values and
# influential observation's (Cook's D))
par(mfrow=c(2,2)) 
plot(fit) 

# Polynomial regression example:
fit2 <- lm(weight ~ height + I(height^2), data=women)
summary(fit2)
plot(women$height,women$weight,
     xlab="Height (in inches)",
     ylab="Weight (in lbs)")
lines(women$height,fitted(fit2))
# Scatterplot (car library) example:
scatterplot(weight ~ height, data=women,
            spread=FALSE, lty.smooth=2,
            pch=19,
            main="Women Age 30-39",
            xlab="Height (inches)",
            ylab="Weight (lbs.)")
# Run an interaction term:
data("cars")
fit <- lm(mpg ~ hp + wt + hp:wt, data = mtcars)
summary(fit)
# Interpretation: A significant interaction between two predictor variables tells you that the relationship between 
# one predictor and the response variable depends on the level of the other predictor. Here it means that the
# relationship between miles per gallon and horse power varies by car weight.
# Visualise, library(effects), errors:
# plot(effect("hp:wt", fit, list(wt=c(2.2,3.2,4.2))), multiline=TRUE)

# Explore and summarise results:
# Useful functions:
  # summary()
  # coefficients()
  # confint()
  # fitted()
  # residuals()
  # anova()
  # vcov()
  # AIC()
  # plot()
  # predict()
####################


####################
# Diagnostic plots (car library):
  # qqPlot() example with multiple linear regression:
states <- as.data.frame(state.x77[,c("Murder", "Population", "Illiteracy", "Income", "Frost")])
head(states)
# Get all the bivariate correlations and plot:
cor(states) 
scatterplotMatrix(states, spread = FALSE, lty.smooth = 2, main = "Scatter Plot Matrix") 
# Run linear regression:
fit <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
# Plot to see normality assumption and explore possible outliers: 
  # QQ-plots assess the normality assumption by plotting studentized residuals against a t distribution with n-p-1 
  # degrees of freedom (n = sample size, p = number of regression parameters (including the intercept)).
qqPlot(fit, labels = row.names(states), id.method = "identify", simulate = TRUE, main = "Q-Q Plot")
fitted(fit)['Nevada']
residuals(fit)['Nevada']
rstudent(fit)['Nevada']
# Alternative plot: Histogram of the studentized residuals and superimposed normal curve, 
# kernel density curve, and rug plot:
# residplot() (p.195 R in Action):
residplot <- function(fit, nbreaks=10) {
  z <- rstudent(fit)
  hist(z, breaks=nbreaks, freq=FALSE,
       xlab="Studentized Residual",
       main="Distribution of Errors")
  rug(jitter(z), col="brown")
  curve(dnorm(x, mean=mean(z), sd=sd(z)),
        add=TRUE, col="blue", lwd=2)
  lines(density(z)$x, density(z)$y,
        col="red", lwd=2, lty=2)
  legend("topright",
         legend = c( "Normal Curve", "Kernel Density Curve"),
         lty=1:2, col=c("blue","red"), cex=.7)
}
residplot(fit)

# Test independence of errors:
durbinWatsonTest(fit)
  
# Plot linearity:
crPlots(fit)

# Test and plot homoscedasticity (ie identify non-constant error variance):
  # Score test of the hypothesis of constant error variance against the alternative that the error variance changes 
  # with the level of the fitted values. A significant result suggests heteroscedasticity (nonconstant error variance).
ncvTest(fit)

# The spreadLevelPlot() function creates a scatter plot of the absolute standardized residuals versus the fitted values, 
# and superimposes a line of best fit:
spreadLevelPlot(fit)
  
# Test for multi-colinearity with variance inflaction factor:
vif(fit)
sqrt(vif(fit)) > 2 # TO DO: print warning/error if one is TRUE

# Test unusual observations
  # Outliers:
outlierTest(fit)

# High leverage values:
  # "Observations that have high leverage are outliers with regard to the other predictors.
  # They have an unusual combination of predictor values. The response value isn’t involved in determining leverage.
  # An observation with a hat value greater than 2 or 3 times the average hat value should be examined."
  # Use hat statistic to identify high leverage values (R in Action p.201):
hat.plot <- function(fit) {
  p <- length(coefficients(fit))
  n <- length(fitted(fit))
  plot(hatvalues(fit), main="Index Plot of Hat Values")
  abline(h=c(2,3)*p/n, col="red", lty=2)
  identify(1:n, hatvalues(fit), names(hatvalues(fit)))
}
hat.plot(fit)

# Identify influential observations:
  # Use Cook’s distance (or D statistic), where values greater than 4/(n-k-1), 
  # where n = sample size, k = number of predictor variables, indicate influential observations. 
  # An alternative cut-off of 1 may be more useful and less sensitive.
cutoff <- 4/(nrow(states)-length(fit$coefficients)-2)
plot(fit, which=4, cook.levels=cutoff)
abline(h=cutoff, lty=2, col="red")

# Determine how influential observations affect the model:
# Use added variable plots: 
# For one response variable and k predictor variables, you’d create k added-variable plots:
# For each predictor Xk , plot the residuals from regressing the response variable on the other k-1 predictors versus 
# the residuals from regressing Xk on the other k-1 predictors. 
avPlots(fit, ask = FALSE, onepage = TRUE, id.method = "identify")


# Combine the information from outlier, leverage, and influence plots:
influencePlot(fit, id.method="identify", main="Influence Plot",
              sub="Circle size is proportional to Cook’s distance")
# Plot legend:
# Influence plot. Elements above +2 or below –2 on the vertical axis are considered outliers. 
# Elements above 0.2 or 0.3 on the horizontal axis have high leverage. 
# Circle size is proportional to influence. O
# Observations depicted by large circles may have disproportionate influence on the parameters estimates of the model.

# Global validation of regression model (gvlma):
# TO DO: save output to file:
gvmodel <- gvlma(fit)
summary(gvmodel)

####################


####################
# If regression assumptions are violated:
  # Delete observations
  # Transform variables
  # Add or delete variables
  # Different regression approach
####################

####################
# Compare and select regression models
fit1 <- lm(Murder ~ Population + Illiteracy + Income + Frost, data = states)
fit2 <- lm(Murder ~ Population + Illiteracy, data = states)
anova(fit2, fit1) # Models must be nested

# Aikake's information criterion:
AIC(fit1,fit2)
####################

####################
# Rank predictors according to relative importance
# Use relaimpo package
# TO DO
# relweights()

####################

####################
# Test prediction of regression model:
# TO DO
# Use bootstrap package
# crossval()


####################


####################
# TO DO: run over many rows and columns
# Use dta.table? h2o?
# See example:
# http://stackoverflow.com/questions/22588673/linear-regression-using-data-table-in-r
# http://stackoverflow.com/questions/23947245/use-predict-on-data-table-with-linear-regression?rq=1
# http://stackoverflow.com/questions/11266930/data-table-vs-plyr-regression-output?rq=1

# Example:
  # disAccRegFunc <- function(dt)
  # {
  #   #Compute Discreationary Accrual
  #   model <- lm(ACNew ~ DSALENew + PPEGTNew + ROANew, data = dt);
  #   dt$RES <- residuals(model);
  #   dt$StudRES <- studres(model);  #Calculation of studentized residuals
  #   return(dt)
  # }
#   dt[,disAccRegFunc(.SD),by=.by]
####################


####################
# Save to file:
# TO DO: check feather (H. Wickham) or fwrite on next data.table release 
write.table(shared_pairs, sprintf('shared_pairs_%s', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
####################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run XGR, IPA, etc, cross with GWAS, ENCODE, etc.
####################