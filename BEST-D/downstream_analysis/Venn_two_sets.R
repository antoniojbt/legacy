####################
# Venn diagram fro reQTLs
# Takes two files and draws a pairwise Venn diagram. Intended for MatrixEQTL data output (eg baseline vs final).
# Performs Fisher's exact test for overlap between sets.
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
library(GeneOverlap)
library(illuminaHumanv4.db)
# source('/ifs/devel/antoniob/projects/BEST-D/functions_for_MatrixeQTL.R')
source('/Users/antoniob/Desktop/Downloads_to_delete/ifs_scripts_backup.dir/projects/BEST-D/functions_for_MatrixeQTL.R')
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# Gene or SNP files:
hits_file1 <- as.character(args[1])
# hits_file1 <- '2000+4000-baseline-1.eQTL_cis'
# hits_file1 <- '2000+4000-12months-1-VS-2000+4000-baseline-1.reQTL_cis'

column_file1 <- as.character(args[2])
# column_file1 <- 'gene'
# column_file1 <- 'Probe_ID'

hits_file2 <- as.character(args[3])
# hits_file2 <- '2000+4000-12months-1.eQTL_cis'

column_file2 <- as.character(args[4])
# column_file2 <- 'gene'

adj_pvalue <- as.numeric(args[5])
# adj_pvalue <- 0.05

col_adj_pvalue <- as.character(args[6])
# col_adj_pvalue <- 'FDR'
# col_adj_pvalue <- 'FDR.final'

col_adj_pvalue2 <- as.character(args[7])
# col_adj_pvalue2 <- 'FDR' # adj.P.Val for limma topTable

# For Fisher's test:
background_file <- as.character(args[8])
# background_file <- 'background_genes_BESTD_expressed.txt'
# background_file <- 'background_probeID_BESTD_expressed.txt'
# background_file <- 'background_SNPs_MatrixEQTL.txt'

adj_pvalue_file2 <- as.numeric(args[9])
# adj_pvalue_file2 <- 0.05

print(args)

# Make adj_pvalue_file2 equal to adj_pvalue if not provided:
if (is.na(adj_pvalue_file2)) {
  adj_pvalue_file2 <- adj_pvalue
}
adj_pvalue
adj_pvalue_file2
####################


####################
# Read data:
hits_data1 <- fread(hits_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
hits_data1
str(hits_data1)
hits_data1 <- hits_data1[which(hits_data1[, col_adj_pvalue, with = F] < adj_pvalue), ]
dim(hits_data1)
length(which(!duplicated(hits_data1[, column_file1, with = F])))
setkey(hits_data1)

hits_data2 <- fread(hits_file2, sep = '\t', header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
hits_data2
hits_data2 <- hits_data2[which(hits_data2[, col_adj_pvalue2, with = F] < adj_pvalue_file2), ]
dim(hits_data2)
length(which(!duplicated(hits_data2[, column_file2, with = F])))
setkey(hits_data2)

# MatrixEQTL outputs gene as column name but using probes so change for some reason I've forgotten
# (labelling for plot titles):
# TO DO: annotated files (eg reQTLs) have header collisions where I've renamed 'gene' for column with gene symbols
# Pass exact header names in bash file in the meantime, otherwise intersection is zero (comparing Probe IDs vs gene symbols)
if (column_file1 == 'gene') {
  column_file1 <- 'probe'
  setnames(hits_data1, "gene", "probe")
} else {
  column_file1 <- column_file1
}
column_file1
hits_data1

if (column_file2 == 'gene') {
  column_file2 <- 'probe'
  setnames(hits_data2, "gene", "probe")
} else {
  column_file2 <- column_file2
}
column_file2
hits_data2

background_data <- fread(background_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
str(background_data)
background_data
setkey(background_data)
dim(background_data)
length(which(!duplicated(background_data[, V1])))
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

intersection <- area1[area2, nomatch=0] # data.table INNER join
intersection
####################

####################
# Set labels:
title_Venn <- sprintf('Shared and unique %ss at %s < %s and < %s', column_file1, col_adj_pvalue, adj_pvalue, adj_pvalue_file2)
title_Venn
file1_label <- sprintf('%ss from %s', column_file1, hits_file1)
file1_label
file2_label <- sprintf('%ss from %s', column_file2, hits_file2)
file2_label
####################

####################
# Draw pairwise Venn diagram with VennDiagram package:
png(sprintf('venn-%s-%s-VS-%s.png', column_file1, hits_file1, hits_file2))
# eg: plotName-description-filename1-VS-filename2.png
# venn-xxx-2000+4000-baseline-1.eQTL_cis-VS-2000+4000-12months-1.eQTL_cis.png
grid.newpage()
venn.plot <- draw.pairwise.venn(area1 = nrow(area1), 
                                area2 = nrow(area2), 
                                cross.area = nrow(intersection), 
                                category = c(hits_file1, hits_file2),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 1,
                                cat.cex      = 1,
                                cat.pos      = c(170, 190),
                                # cat.dist     = c(1, 1),
                                # cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 0,
                                ext.dist     = -0.05,
                                # ext.length   = 0.85,
                                # ext.line.lwd = 2,
                                # ext.line.lty = "dashed", # Adds external line to circles
                                ind = FALSE, # Do not plot automatically
                                # print.mode = c('percent', 'raw'), sigdigs = 1, # Prints percentage and raw
                                scaled = TRUE, euler.d = TRUE, sep.dist = 0.03, rotation.degree = 180)
# grid.draw(venn.plot)
grid.arrange(gTree(children = venn.plot), top = title_Venn)
dev.off()
####################

####################
# Run Fisher's test / overlap analysis:

# go.obj <- newGeneOverlap(hits_data$V1, annotation_data$V1, spec = "hg19.gene")
str(area1)

go.obj <- newGeneOverlap(area1[, column_file1, with = F], area2[, column_file2, with = F], 
                         genome.size = length(background_data[, V1]))

# Gene overlap package doesn't like data.frame so convert:
area1 <- as.data.frame(area1[, column_file1, with = F])
area2 <- as.data.frame(area2[, column_file2, with = F])

go.obj <- newGeneOverlap(area1[, column_file1], area2[, column_file2], 
                         genome.size = length(background_data[, V1]))
# See results of basic overlaps:
go.obj
# Perform statistical test (Fisher's):
go.obj <- testGeneOverlap(go.obj)
go.obj
# Detailed results:
print(go.obj)
# Get specific slots of object (ie overlapping gene set):
head(getIntersection(go.obj))
length(getIntersection(go.obj))
####################

####################
# Get gene symbols for intersection if running gene expression probes:
# TO DO: do.call (?) as it is in get_illumina_annot generates many duplicates
# if (column_file1 == 'gene' | column_file1 == 'probe' | column_file1 == 'Probe_ID') {
#   get_intersect <- as.data.frame(getIntersection(go.obj))
#   colnames(get_intersect)[1] <- 'probes_intersect'
#   # dim(get_intersect)
#   # head(get_intersect)
#   if (dim(get_intersect)[1] != 0) {
#     gene_names <- get_illumina_annot(get_intersect, 1)
#     gene_names <- gene_names[which(!duplicated(gene_names[, 1])), ]
#     gene_names <- as.data.frame(gene_names)
#   }
# }
# dim(gene_names)
# length(gene_names[, 'Probe_ID'])
# length(unique(gene_names[, 'Probe_ID']))
# class(gene_names)
# head(gene_names)
# tail(gene_names)

# Write results of overlap to disk:
# eg table-xxx-2000+4000-baseline-1.reQTL_trans.txt
file_name <- sprintf('table-Fishers%s-%s-VS-%s.txt', column_file1, hits_file1, hits_file2)
sink(file_name, append = FALSE, split = TRUE, type = c("output"))
print(hits_file1)
print(hits_file2)
print(background_file)
print(go.obj)
length(unique(getIntersection(go.obj)))
getIntersection(go.obj)
# gene_names
# TO DO: if gene symbols were annotated, prin them (ie if intersection wasn't zero and it was GEx probes being tested), 
# this doesn't work:
# exists('gene_names')
# 
# if exists('gene_names') {
#   gene_names
# }
# else {
#   getIntersection(go.obj)
# }
sink(file = NULL)
####################

# TO DO: add legend, method, etc.:
####################
# Legend Venn:
# table/plotName.legend
# - boxplot-xxx-2000+4000-baseline-1.reQTL_trans.legend
####################

####################
# Legend Fisher's:
# y-axis, x-axis, scale, labels and colours if used, stats transformations, background if used, etc.
# table/plotName.legend
# - table-xxx-2000+4000-baseline-1.reQTL_trans.legend
####################

####################
# Method Venn:
# Include packages, versions, etc.
# - boxplot-xxx-2000+4000-baseline-1.reQTL_trans.method
####################

####################
# Method Fisher's:
# Include packages, versions, etc.
# - table-xxx-2000+4000-baseline-1.reQTL_trans.method
####################


####################
# Interpretation Venn:
# Create file as placeholder and then add text manually
# table/plotName.interpretation
# - boxplot-xxx-2000+4000-baseline-1.reQTL_trans.interpretation
####################

####################
# Interpretation Fisher's:
# Create file as placeholder and then add text manually
# table/plotName.interpretation
# - table-xxx-2000+4000-baseline-1.reQTL_trans.interpretation
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
