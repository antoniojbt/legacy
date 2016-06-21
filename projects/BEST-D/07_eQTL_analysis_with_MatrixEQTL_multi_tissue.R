#############################################
# Multi-tissue MatrixEQTL analysis part 4 of 4.
# 26 November 2015
# Antonio J Berlanga-Taylor
# Code from:
# # http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/
# http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code4.html
# http://arxiv.org/pdf/1311.2948v3.pdf
#############################################


#############################################
# Part 3/4 of the code is in:
# /ifs/devel/antoniob/projects/BEST-D/xxx 
# This estimates the parameters for the MT analysis and is used as input for the last step (eQTL calls, this script): 

# Input for code 4/4: 
# z-score.matrix.Rdata	—	R data file with the big matrix of z-scores
# z-score.matrix.genes.txt	—	text file with gene names for the rows of the big matrix
# z-score.matrix.snps.txt	—	text file with SNP names for the rows of the big matrix
# df.txt	—	text file with tissue names and the number of degrees of freedom for each tissue
# paralist.Rdata	—	R data file with the sequence of parameter estimates from EM algorithm

# Output: 
# The code below uses the parameter estimates to assess which gene-SNP pairs are eQTLs in each of the tissues. 
# The output file MT-eQTLs.txt has contains a record for each significant gene-SNP pair with the following information:
 
# SNP - the SNP name
# gene - the gene name
# isEQTL.* - a set of columns (one per tissue) indicating the most likely tissue specificity profile.
# These columns have 1's for the tissues where the gene-SNP pair is estimated to have an eQTL.
# marginalP.* - a set of columns (one per tissue) with marginal probability of not having an eQTL.
# Small values in these columns indicate high likelihood of the gene-SNP pair being an eQTL in the tissue.


# These are the final results and correspond to the tissue/condition specific eQTLs.

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
args <- commandArgs(trailingOnly = TRUE)
#args[1] = "df_cross_2000.txt"
#args[2] = "_loose_p5_1MB_CROSS_2000.cis"

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
## Run MatrixEQTL multi-tissue part 4/4:

### Parameters
local.FDR.threshold = 0.01;
output.file.name = paste(as.character(args[1]), as.character(args[2]), sep = '');

### Load big matrix of z-scores
load(file="z-score.matrix.Rdata");
dim(big.matrix);
big.matrix

### Load gene names and SNP names
### matching the rows of big.matrix
gnames = readLines("z-score.matrix.genes.txt");
snames = readLines("z-score.matrix.snps.txt");

### Load tissue names
df = read.table(as.character(args[1]), stringsAsFactors=FALSE);
names(df) = c("tissue","df");
show(df);

### Load parameter estimates and pick the last one
load("paralist.Rdata");
#param = tail(paralist,1)[[1]];
param = paralist[[1]]

str(paralist) # added
class(paralist) # added
head(paralist) # added
paralist[[1]]$Delta # added
str(paralist[[1]]) # added

# identical(paralist[[1]]$Delta, paralist[[2]]$Delta) # added

### Number of tissues
K = ncol(big.matrix);
m = nrow(big.matrix);

### The function for matrix power
mat.power = function(mat,pow) {
  e = eigen(mat);
  V = e$vectors;
  return( V %*% diag(e$values^pow) %*% t(V) );
}

### Matrix of possible tissue specificity profiles
Pmat = simplify2array(param$Psubs);
Pmat # added

###
### Call eQTLs and save in a file
###

fid = file(description=output.file.name, open="wt");
writeLines(con=fid, paste(
  "SNP\tgene\t",
  paste("isEQTL.", df$tissue, collapse="\t"),
  "\t",
  paste("marginalP.", df$tissue, collapse="\t")));


ls() # added

### Do calculations in slices of 10000 gene-SNP pairs
step1 = 1e4L;
cumdump = 0;
for( j in 1:ceiling(nrow(big.matrix)/step1)) {
  fr = step1*(j-1)+1;
  to = min(step1*j,nrow(big.matrix));
  X = big.matrix[fr:to, , drop=FALSE];
  
  ### likelihood for the slice
  prob = matrix(0, nrow(X), length(param$P));
  for( i in 1:length(param$Psubs)) {
    sigma_star = param$Delta + param$Sigma * tcrossprod(param$Psubs[[i]]);
    sigma_hfiv = mat.power(sigma_star,-0.5);
    sigma_dethfiv = (det(sigma_star))^(-0.5);
    w = (1/(2*pi)^(K/2)) * (param$P[i]*sigma_dethfiv);
    prob[,i] = exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2),X)^2) ) ;
  }
  prob = prob / rowSums(prob);
  
  ### Select tests with eQTLs significant at local.FDR.threshold level
  keep = (prob[,1] <= local.FDR.threshold);
  if(any(keep)) {
    marginalProb = tcrossprod(prob[keep,,drop=FALSE], 1-Pmat);
    tissueSpecificity = t(Pmat)[
      apply(X=prob[keep,,drop=FALSE], MARGIN=1, FUN=which.max), ];
    
    dump = data.frame(
      snames[(fr:to)[keep]],
      gnames[(fr:to)[keep]],
      tissueSpecificity,
      marginalProb,
      row.names = NULL,
      check.rows = FALSE,
      check.names = FALSE,
      stringsAsFactors = FALSE);
    write.table(dump, file = fid, quote = FALSE,
                sep = "\t", row.names = FALSE, col.names = FALSE);
  }
  cumdump = cumdump + sum(keep);
  cat("Slice",j,"of",ceiling(nrow(big.matrix)/step1)," eQTLs recorded:",cumdump,"\n");
}
close(fid);

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
