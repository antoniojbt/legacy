#############################################
# Multi-tissue MatrixEQTL analysis part 3 of 4.
# 26 November 2015
# Antonio J Berlanga-Taylor
# Code from:
# # http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/
# http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code3.html
# http://arxiv.org/pdf/1311.2948v3.pdf
#############################################


#############################################
# Part 2/4 of the code is in:
# /ifs/devel/antoniob/projects/BEST-D/xxx 
# This produces a large matrix with all the values needed for a MT-eQTL. Run these files here to estimate 
# the parameters for the MT analysis: 

# Input for code 3/4: 
# z-score.matrix.Rdata	—	R data file with the big matrix of z-scores
# df.txt	—	text file with tissue names and the number of degrees of freedom for each tissue

# Output: 
#The code below estimates the model parameters and records the whole sequence of intermediate estimates 
# in the variable paralist and save it in the 1 file:
# paralist.Rdata	—	R data file with the sequence of parameter estimates from EM algorithm

# These are used as input for part 4/4 (see separate script).

#############################################

#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Documents/Work_other/DPhil-09-July-2013/BEST-D/BEST-D_analysis/BEST-D_analysis_eQTL/BEST-D_multi_tissue/MatrixEQTL_results_on_ifs/')  

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
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
#load('R_session_saved_image_order_and_match.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_MatrixQTL_analysis','.RData', sep='')
R_session_saved_image

# Run with cmd arguments:
# args <- commandArgs(trailingOnly = TRUE)
# args[1] = xxxx
# args[2] = 'subset_baseline_placebo.tsv'
# args[3] = 'covar_PCAs_t_placebo_baseline.tsv'

#############################################


#############################################
## Update packages if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite(c('MatrixEQTL', 'devtools', 'trio', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'optparse'))

#and load them:
packages <-c('MatrixEQTL', 'dplyr', 'ggplot2', 'tidyr', 'knitr', 'readr', 'reshape2', 'data.table')
lapply(packages, require, character.only = TRUE)
sessionInfo()
#############################################


#############################################
## Run MatrixEQTL multi-tissue part 3/4:

### Set estimation parameters
maxIterations = 100;

### Load big matrix of z-scores
load(file='z-score.matrix.Rdata');
dim(big.matrix);

### Initialize parameters
{
  param = list();
  
  ### K - the number of tissues
  K = ncol(big.matrix);
  
  ### Delta - null covariance matrix across tissues
  param$Delta = matrix(0.05, K, K);
  diag(param$Delta) = 1;
  
  ### Sigma - signal covariance matrix across tissues
  param$Sigma = matrix(3, K, K) + diag(K);
  
  ### P - the vector of probabilities
  param$P = rep(1/2^K, 2^K);
  
  ### Psubs - the vector of active eQTLs
  ###         for each element of P
  Psubs = vector('list',2^K);
  for( i in 1:2^K) {
    a<-2^((K-1):0);
    b<-2*a;
    Psubs[[i]] = as.double(((i-1) %% b)>=a);
  }
  rm(a,b,i);
  param$Psubs = Psubs;
  rm(Psubs);
  
  ### loglik - the initial likelihood
  param$loglik = -Inf;
  rm(K);
}

ls() # added
class(param) # added
head(param) # added
str(param) # added

### The function does a single iteration of the estimation procedure
DoIteration = function(big.matrix, param) {
  ### extract current model parameters
  K = ncol(big.matrix);
  m = nrow(big.matrix);
  Delta = param$Delta;
  Sigma = param$Sigma;
  P = param$P;
  Psubs = param$Psubs;
  
  ### The function for matrix power
  mat.power = function(mat,pow) {
    e = eigen(mat);
    V = e$vectors;
    return( V %*% diag(e$values^pow) %*% t(V) );
  }
  
  ### Start the timer
  tic = proc.time();
  
  ### variables to accumulate
  ### loglik - likelihood
  ### newP - marginal probabilities
  ### newDelta - the new Delta matrix
  ### newSigmaPlusDelta - Delta+Sigma
  cum.loglik = 0;
  cum.newP = 0;
  cum.newDelta = 0;
  cum.newSigmaPlusDelta = 0;
  
  ### Do calculations in slices of 10000 gene-SNP pairs
  step1 = 1e5L;
  for( j in 1:ceiling(m/step1)) {
    fr = step1*(j-1)+1;
    to = min(step1*j,m);
    X = big.matrix[fr:to, , drop=FALSE];
    
    ### likelihood for the slice
    prob = matrix(0, nrow(X), length(P));
    for( i in 1:length(Psubs)) {
      sigma_star = Delta + Sigma * tcrossprod(Psubs[[i]]);
      sigma_hfiv = mat.power(sigma_star,-0.5);
      sigma_dethfiv = (det(sigma_star))^(-0.5);
      w = (1/(2*pi)^(K/2)) * (P[i]*sigma_dethfiv);
      prob[,i] = exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2),X)^2) ) ;
    }
    
    cum.loglik = cum.loglik + sum(log(rowSums(prob)));
    
    ### Normalize probabilities for each gene-SNP pair
    ### to add up to 1
    prob = prob / rowSums(prob);
    ### new vector of P - tissue specificity probabilities
    cum.newP = cum.newP + colSums(prob);
    
    cum.newDelta = cum.newDelta +
      crossprod(X*sqrt(prob[,1]));
    cum.newSigmaPlusDelta = cum.newSigmaPlusDelta +
      crossprod(X*sqrt(prob[,length(P)]));
  }
  
  {
  ### Calculate Delta from the cumulative sum
  Delta = cum.newDelta / cum.newP[1];
  ### normalize to force the diagonal to 1
  Delta = Delta * tcrossprod(sqrt(1 / diag(Delta)));
  
  ### Same with Sigma
  Sigma = cum.newSigmaPlusDelta / tail(cum.newP,1) - Delta;
  e = eigen(Sigma);
  if( any(e$values<0) ) {
    Sigma = e$vectors %*% diag(pmax(e$values,0)) %*% t(e$vectors);
  }  
  }
  P = cum.newP / sum(cum.newP);
  
  toc = proc.time();
  return(list(Delta = Delta, Sigma = Sigma, P = P,
              Psubs = Psubs, loglik = cum.loglik, time = toc-tic));
}

ls() # added

### The 'paralist' list vector
### will store model estimates at each iteration
paralist = vector('list',maxIterations+1);
paralist # added
paralist[[1]] = param;
paralist # added
str(paralist) # added
str(paralist[[1]]) # added
#rm(param);  # commented out


### Perform up to 'maxIterations' iteration
# for( i in 2:length(paralist) ) {
#   paralist[[i]] = DoIteration(big.matrix=big.matrix, param=paralist[[i-1]]);
#   cat(i, '\t', paralist[[i]]$loglik-paralist[[i-1]]$loglik,
#       '\t', paralist[[i]]$time[3],'\n');
#   if(i>10)
#     if(paralist[[i]]$loglik < paralist[[i-1]]$loglik)
#       break;
# }

### Re-wrote above funcion to avoid error being called when using break (so as to run as an Rscript).
DoMaxIter <- function(x){
  for( i in 2:length(x) ) {
    x[[i]] = DoIteration(big.matrix=big.matrix, param=x[[i-1]]);
    cat(i, '\t', x[[i]]$loglik-x[[i-1]]$loglik,
        '\t', x[[i]]$time[3],'\n');
    if(i>10)
      if(x[[i]]$loglik < x[[i-1]]$loglik)
        break;
  }
}
tryCatch(DoMaxIter(paralist), error=function(e) NULL)

### Continue original script:
paralist = paralist[!sapply(paralist, is.null)];
str(paralist)
class(paralist)
head(paralist)

### Save the results
save(list='paralist', file='paralist.Rdata');

#############################################


#############################################
# The end:
# Remove objects that are not necessary to save:
ls()
object_sizes <- sapply(ls(), function(x) object.size(get(x)))
as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for xxx.
#############################################