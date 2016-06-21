####################
# Venn diagram fro reQTLs
# Takes two files and draws a pairwise Venn diagram. Intended for MatrixEQTL data output (baseline vs final).
# 29 May 2016
####################

####################
# See:
# http://mfcovington.github.io/r_club/exercises/2013/05/15/venn-euler-demo/
# https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html
####################

####################
options(echo = TRUE)
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/')

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
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_XGR.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_Venn_reQTLs','.RData', sep='')
R_session_saved_image
####################

####################
# Libraries:
library(VennDiagram)
library(data.table)
library(gridExtra)
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene or SNP files:
hits_file1 <- as.character(args[1])
# hits_file1 <- 'cut_genotype_data_all_treated_baseline.tsv_matched.tsv_MxEQTL_p1_1e+06.cis'
# hits_file1 <- 'trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_baseline.tsv_matched.tsv'

column_file1 <- as.character(args[2])
# column_file1 <- 'gene'

hits_file2 <- as.character(args[3])
# hits_file2 <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv_MxEQTL_p1_1e+06.cis'
# hits_file2 <- 'trans_covSelect_1e-05_0.001_1e+06_cut_genotype_data_all_treated_final.tsv_matched.tsv'

column_file2 <- as.character(args[4])
# column_file2 <- 'gene'

adj_pvalue <- as.integer(args[5])
# adj_pvalue <- 0.05

col_adj_pvalue <- as.character(args[6])
# col_adj_pvalue <- 'FDR'

hits_file3 <- as.character(args[7])
# hits_file3 <- 'full_topTable_pairing_all_treated.txt'

column_file3 <- as.character(args[8])
# column_file3 <- 'probe_ID'

col_adj_pvalue2 <- as.character(args[9])
# col_adj_pvalue2 <- 'adj.P.Val'
####################

####################
# Read data:
hits_data1 <- fread(hits_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
# hits_data <- as.data.frame(hits_data)
hits_data1
str(hits_data1)
hits_data1 <- hits_data1[which(hits_data1[, col_adj_pvalue, with = F] < adj_pvalue), ]
dim(hits_data1)
setkey(hits_data1)


hits_data2 <- fread(hits_file2, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
# hits_data <- as.data.frame(hits_data)
hits_data2
hits_data2 <- hits_data2[which(hits_data2[, col_adj_pvalue, with = F] < adj_pvalue), ]
dim(hits_data2)
setkey(hits_data2)

hits_data3 <- fread(hits_file3, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
# hits_data <- as.data.frame(hits_data)
hits_data3
hits_data3 <- hits_data3[which(hits_data3[, col_adj_pvalue2, with = F] < adj_pvalue), ]
dim(hits_data3)
setkey(hits_data3)
####################

####################
# Remove duplicates and calculate intersection:
area1 <- hits_data1[, column_file1, with = F]
area1
setkey(area1)
dim(area1)
length(which(!duplicated(area1[, column_file1, with = F])))
area1 <- area1[which(!duplicated(area1[, column_file1, with = F])), ]
area1

area2 <- hits_data2[, column_file2, with = F]
area2
setkey(area2)
dim(area2)
length(which(!duplicated(area2[, column_file2, with = F])))
area2 <- area2[which(!duplicated(area2[, column_file2, with = F])), ]
area2

area3 <- hits_data3[, column_file3, with = F]
area3
setkey(area3)
dim(area3)
length(which(!duplicated(area3[, column_file3, with = F])))
area3 <- area3[which(!duplicated(area3[, column_file3, with = F])), ]
area3

intersection12 <- area1[area2, nomatch=0] # data.table INNER join
intersection12

intersection13 <- area1[area3, nomatch=0] # data.table INNER join
intersection13

intersection23 <- area2[area3, nomatch=0] # data.table INNER join
intersection23

intersection123 <- intersection12[area3, nomatch=0] # data.table INNER join
intersection123

####################

####################
# Set labels:
# TO DO: Will break with different file names:
label1 <- strsplit(hits_file1, '_')
label1 <- strsplit(label1[[1]], '[.]')
label1 <- paste(label1[[6]][1], label1[[10]][2], sep = ' ')
label1
label2 <- strsplit(hits_file2, '_')
label2 <- strsplit(label2[[1]], '[.]')
label2 <- paste(label2[[6]][1], label2[[10]][2], sep = ' ')
label2

label3 <- strsplit(hits_file3, '_')
label3 <- strsplit(label3[[1]], '[.]')
label3 <- paste(label3[[2]][1], sep = ' ')
label3

title_Venn <- sprintf('Shared and unique %ss at %s < %s', column_file1,
                      col_adj_pvalue, adj_pvalue)
title_Venn
file1_label <- sprintf('%ss from %s', column_file1, label1)
file1_label
file2_label <- sprintf('%ss from %s', column_file2, label2)
file2_label
# file3_label <- sprintf('%ss from %s', column_file3, label3)
file3_label <- sprintf('genes from %s', label3)
file3_label
####################

####################
# Draw pairwise Venn diagram with VennDiagram package:
png(sprintf('3set_venn_%s_%s.png', column_file1, hits_file1))
grid.newpage()
venn.plot <- draw.triple.venn(area1 = nrow(area1), 
                              area2 = nrow(area2), 
                              area3 = nrow(area3), 
                              n12 = nrow(intersection12),
                              n23 = nrow(intersection23),
                              n13 = nrow(intersection13),
                              n123 = nrow(intersection123),
                              category = c(file1_label, file2_label, file3_label),
                              fill         = c("blue", "red", "green"),
                              alpha        = 0.3,
                              lty          = "blank",
                              cex          = 1,
                              # cat.cex      = 1,
                              cat.pos      = c(340, 20, 180),
                              # cat.dist     = c(1, 1),
                              # cat.just     = list(c(-1, -1), c(1, 1)),
                              # ext.pos      = 0,
                              # ext.dist     = -0.05,
                              # cat.default.pos = 'outer',
                              # ext.length   = 0.85,
                              # ext.line.lwd = 2,
                              # ext.line.lty = "dashed", # Adds external line to circles
                              ind = FALSE, # Do not plot automatically
                              # print.mode = c('percent', 'raw'), sigdigs = 1, # Prints percentage and raw
                              scaled = TRUE, euler.d = TRUE, sep.dist = 0.03) #, rotation.degree = 180)
# grid.draw(venn.plot)
grid.arrange(gTree(children = venn.plot), top = title_Venn)
dev.off()
####################

####################
# # TO DO: run Fisher's test on each set
####################

####################
# TO DO/clean up
####################


####################
# Legend:
# 
# 
####################

####################
# Method:
# 
# 
####################

####################
# Interpretation:
# 
# 
####################


####################
# Write results to disk:

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

# Next: run the script for xxx.
####################
