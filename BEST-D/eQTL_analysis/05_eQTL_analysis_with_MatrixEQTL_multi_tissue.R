#############################################
# Multi-tissue MatrixEQTL analysis part 2 of 4.
# 26 November 2015
# Antonio J Berlanga-Taylor
# Code from:
# # http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/
# http://www.bios.unc.edu/research/genomic_software/Multi-Tissue-eQTL/code2.html
# http://arxiv.org/pdf/1311.2948v3.pdf
#############################################


#############################################
# Part 1/4 of the code is in:
# /ifs/devel/antoniob/projects/BEST-D/04_eQTL_analysis_with_MatrixEQTL.R
# This produces a single tissue/condition eQTL analysis. Run for each tissue and use the results as input for code 2/4.

# Input: 
# df.txt which contains the names of the tissues/conditions and degrees of freedom for each.
# The results of running MatrixEQTL on each tissue (eg several 'xxx.cis' files). Run cis and trans separately in this 
# multi-tissue (MT_eQTL) analysis.

# Output: 
# a big matrix of z-scores with each row containing informaiton for a particular gene-SNP pair 
# and each column for a particular tissue. It produces 3 files:
# z-score.matrix.Rdata	—	R data file with the big matrix of z-scores
# z-score.matrix.genes.txt	—	text file with gene names for the rows of the big matrix
# z-score.matrix.snps.txt	—	text file with SNP names for the rows of the big matrix

# These are used as input for part 3/4 (see separate script).

# Also check out:
# MT-HESS: 
# http://bioinformatics.oxfordjournals.org/content/early/2015/11/03/bioinformatics.btv568.full.pdf+html
# https://github.com/Bayesian-Omics/CHESS
# Panama: 
# http://pmbio.github.io/envGPLVM/

#############################################

#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Documents/Work_other/DPhil-09-July-2013/BEST-D/BEST-D_analysis/BEST-D_analysis_eQTL/BEST-D_multi_tissue/MatrixEQTL_results_on_ifs/MT_eQTL_3.dir/')

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

# Input name of file with tissue names and degrees of freedom for each:
#args[1] = "df_cross_placebo.txt"
# Input names of MatrixEQTL results to paste below (ie tissue+xxx):
#args[2] = "_MatrixEQTL_loose_p5_1MB_modelLINEAR_CROSS.cis"
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
## Run MatrixEQTL multi-tissue part 2/4:
# Read df.txt for the list of tissues and degrees of freedom of linear models:
#getwd()
#system('ls')

df = read.table(as.character(args[1]), stringsAsFactors=FALSE); # modified
names(df) = c("tissue", "df");
show(df)

# List vector for storing Matrix eQTL results
big.list = vector("list", nrow(df));

# Store gene and SNP names from the first tissue
# for matching with other tissues
genes = NULL;
snps = NULL;

# colClasses for faster reading of Matrix eQTL output
cc.file = NA;

# Loop over tissues
for(t1 in 1:nrow(df) ) {
  
  ### Get tissue name
  tissue = df$tissue[t1];
  
  ### Load Matrix eQTL output for the given tissue
  start.time = proc.time()[3];
  ####### I modified here as my files are called:
  #genotype_data_4000_final.tsv_MatrixEQTL_loose_1MB.trans (or cis)
  #genotype_data_2000_final.tsv
  
  #tbl = read.table(paste0("eQTL_results_AL_",tissue,"_cis.txt"),
  #                 header = T, stringsAsFactors=FALSE, colClasses=cc.file);
  
  tbl = read.table(paste(tissue, as.character(args[2]), sep = ''),
                   header = T, stringsAsFactors=FALSE, colClasses=cc.file);
  
  #######
  end.time = proc.time()[3];
  cat(tissue, "loaded in", end.time - start.time, "sec.", nrow(tbl), "gene-SNP pairs.", "\n");
  
  ### set colClasses for faster loading of other results
  if(any(is.na(cc.file))) {
    cc.file = sapply(tbl, class);
  }
  
  ### Set gene and SNP names for matching
  if(is.null(snps))
    snps = unique(tbl$SNP);
  if(is.null(genes))
    genes = unique(tbl$gene);
  
  ### Match gene and SNP names from Matrix eQTL output
  ### to "snps" and "genes"
  gpos = match(tbl$gene, genes, nomatch=0L);
  spos = match(tbl$SNP,  snps,  nomatch=0L);
  
  ### Assign each gene-SNP pair a unique id for later matching with other tissues
  id = gpos + 2*spos*length(genes);
  
  ### Transform t-statistics into correlations
  r = tbl$t.stat / sqrt(df$df[t1] + tbl$t.stat^2);
  
  ### Record id's and correlations
  big.list[[t1]] = list(id = id, r = r);
  
  ### A bit of clean up to reduce memory requirements
  rm(tbl, gpos, spos, r, id, tissue, start.time, end.time);
  gc();
}
rm(t1, cc.file);

### Find the set of gene-SNP pairs
### present in results for all tissues
keep = rep(TRUE, length(big.list[[1]]$id));
for(t1 in 2:nrow(df)) {
  mch = match(big.list[[1]]$id, big.list[[t1]]$id, nomatch=0L);
  keep[ mch == 0] = FALSE;
  cat(df$tissue[t1], ", overlap size", sum(keep), "\n");
}
final.ids = big.list[[1]]$id[keep];
rm(keep, mch, t1);

### Create and fill in the matrix of z-scores
### Z-scores are calculated from correlations
big.matrix = matrix(NA_real_, nrow=length(final.ids), ncol=nrow(df));
fisher.transform = function(r){0.5*log((1+r)/(1-r))};
for(t1 in 1:nrow(df)) {
  mch = match(final.ids,  big.list[[t1]]$id);
  big.matrix[,t1] = fisher.transform(big.list[[t1]]$r[mch]) * sqrt(df$df[t1]-1);
  cat(t1,"\n");
}
stopifnot( !any(is.na(big.matrix)) )
rm(t1, mch);

### Save the big matrix
save(list="big.matrix", file="z-score.matrix.Rdata", compress=FALSE);

### Save gene names and SNP names for rows of big matrix
writeLines(text=genes[final.ids %% (length(genes)*2)], con="z-score.matrix.genes.txt")
writeLines(text=snps[final.ids %/% (length(genes)*2)], con="z-score.matrix.snps.txt")
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