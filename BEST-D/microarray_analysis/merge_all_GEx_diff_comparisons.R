#############################
# To be run after 03 differential expression analyses scripts
# Antonio J Berlanga-Taylor
# 06 Nov 2016
# BEST-D project differential expression
# Requires all files in one directory (eg from 03 no covar and 03c VD def scripts)
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_3.dir")
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/tables_and_plots_for_draft/final_draft_BEST-D/tables/all_GEx_tables/individual_tables_all_GEx_comparisons/')
# setwd('/Users/antoniob/Desktop/Downloads_to_delete/mePipe_runs_2.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_diff_expression_merge_tables",".txt", sep=""), open='a')
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
# load('R_session_saved_image_pheno_file_check.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_merge_tables', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')
#install.packages("statmod")

library(illuminaHumanv4.db)
library(plyr)

# Get functions from other scripts (eg to add annotations to topTable results):
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/microarray_analysis/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/cgat_projects/BEST-D/moveme.R')
#############################

#############################
# Read and merge all comparisons into one table:
getwd()
pattern <- 'topTable'
# pattern <- '(2000)+.*tsv$'
length(dir(pattern = pattern))
dir(pattern = pattern)
dir(pattern = pattern)[1]
# Set-up first df for merging:
file_name <- dir(pattern = pattern)[1]
initialise_df <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
names(initialise_df)
all_GEx_comparisons <- as.data.frame(initialise_df[, 1])
# names(all_GEx_comparisons)[1] <- 'X'
all_GEx_comparisons <- data.frame(all_GEx_comparisons, logFC = "logFC", AveExpr = "AveExpr", 
                                  t = "t", P.Value = "P.Value", adj.P.Val = "adj.P.Val", B = "B")
dim(all_GEx_comparisons)
class(all_GEx_comparisons)
head(all_GEx_comparisons)
head(initialise_df)
all_GEx_comparisons <- merge(all_GEx_comparisons, initialise_df, by = 'X')
# head(all_GEx_comparisons)

for (i in dir(pattern = pattern)) {
  # print(i)
  df_in <- read.csv(i, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
  i <- strsplit(i, '[.]')
  i <- as.character(i[[1]][1])
  i <- substr(i, 10, nchar(i))
  i <- paste0('.', i)
  # print(i)
  all_GEx_comparisons <- merge(all_GEx_comparisons, df_in, by ='X', suffixes = c('', i))
}
names(all_GEx_comparisons)
head(all_GEx_comparisons)

# Delete dummy columns:
all_GEx_comparisons <- all_GEx_comparisons[, -c(2:7)]

# Columns to expect plus 1 (ID column):
length(dir(pattern = pattern)) * 6
dim(all_GEx_comparisons)
head(all_GEx_comparisons)
# View(all_GEx_comparisons)

# Add gene symbols and FC to each column
names(all_GEx_comparisons)[1] <- 'probe_ID'
# Get probe IDs from topTable results:
illumina_ids  <-  as.character(all_GEx_comparisons[, 'probe_ID'])
# Get Entrez IDs:
probes_by_ENTREZID_RE  <- unlist(mget(illumina_ids, illuminaHumanv4ENTREZREANNOTATED, ifnotfound = NA))
probes_by_ENTREZID_RE  <- as.data.frame(probes_by_ENTREZID_RE)
# Get gene symbols:
probes_by_symbol <- unlist(mget(illumina_ids, illuminaHumanv4SYMBOLREANNOTATED, ifnotfound = NA))
probes_by_symbol <- as.data.frame(probes_by_symbol)
# Set IDs as columns for merge:
probes_by_symbol$probe_ID <- row.names(probes_by_symbol)
probes_by_ENTREZID_RE$probe_ID <- row.names(probes_by_ENTREZID_RE)
# Merge results and annotations:
head(probes_by_symbol)
head(probes_by_ENTREZID_RE)
all_GEx_comparisons <- merge(all_GEx_comparisons, probes_by_symbol, by = 'probe_ID')
all_GEx_comparisons <- merge(all_GEx_comparisons, probes_by_ENTREZID_RE, by = 'probe_ID')
# Move columns:
all_GEx_comparisons <- all_GEx_comparisons[, moveme(names(all_GEx_comparisons), 'probes_by_symbol after probe_ID')]
all_GEx_comparisons <- all_GEx_comparisons[, moveme(names(all_GEx_comparisons), 'probes_by_ENTREZID_RE after probes_by_symbol')]
# Order lists by adj. p-value in increasing order:
# all_GEx_comparisons <- all_GEx_comparisons[order(all_GEx_comparisons[, 'adj.P.Val'], decreasing = FALSE), ]
all_GEx_comparisons[1:10, 1:10]
names(all_GEx_comparisons)[1:10]
# View(all_GEx_comparisons)
names(all_GEx_comparisons)
# Sanity check, placebo paired analyses was run twice, these columns should be the same to 2 or 3 decimals:
# View(all_GEx_comparisons[, c('probe_ID', 'P.Value.joint_pairedtreated_placebo', 'P.Value.pairedtreated_placebo')])


# Write to disk:
write.table(all_GEx_comparisons, 'all_GEx_comparisons.tsv', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)
#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')
sessionInfo()

q()

# Next: run script for higher level analyses.
#############################