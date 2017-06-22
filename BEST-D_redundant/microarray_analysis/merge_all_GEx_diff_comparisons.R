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
# load('R_session_saved_image_diff_expression_merge_tables.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_merge_tables', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:
#install.packages('ellipse')
#install.packages("statmod")

library(illuminaHumanv4.db)
library(plyr)
library(limma)

# Get functions from other scripts (eg to add annotations to topTable results):
# source('/Users/antoniob/Documents/github.dir/AntonioJBT/cgat_projects/BEST-D/microarray_analysis/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/gene_expression_functions.R')
source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('/Users/antoniob/Documents/github.dir/AntonioJBT/cgat_projects/utility_tutorial_scripts//moveme.R')
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
head(initialise_df)
all_GEx_comparisons <- as.data.frame(initialise_df[, 1])
names(all_GEx_comparisons)[1] <- 'X'
names(all_GEx_comparisons)
all_GEx_comparisons <- data.frame(all_GEx_comparisons, logFC = "logFC", AveExpr = "AveExpr", 
                                  t = "t", P.Value = "P.Value", adj.P.Val = "adj.P.Val", B = "B")
dim(all_GEx_comparisons)
class(all_GEx_comparisons)
head(all_GEx_comparisons)
head(initialise_df)
all_GEx_comparisons <- merge(all_GEx_comparisons, initialise_df, by = 'X')
head(all_GEx_comparisons)
dim(all_GEx_comparisons)
# This file now has 13 columns (probe ID once, the rest twice). All but probe ID (here 'X') are dummy columns
# to be deleted later.

for (i in dir(pattern = pattern)) {
  print(i)
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
dim(all_GEx_comparisons)
# This object now has the 13 initial columns plus 6 additional for each file read
# Columns to expect plus 1 (ID column):
length(dir(pattern = pattern)) * 6
dim(all_GEx_comparisons)
dim(all_GEx_comparisons) - length(dir(pattern = pattern)) * 6 # Should give 13
head(all_GEx_comparisons)
head(all_GEx_comparisons)[, 1:10]
# View(all_GEx_comparisons)

# Delete dummy columns, these are from the dummy variables first created (2 to 7):
head(all_GEx_comparisons[, c(1:8)])
head(all_GEx_comparisons[, c(9:16)])
head(all_GEx_comparisons[, c(17:24)])
head(all_GEx_comparisons[, c(25:32)])
all_GEx_comparisons <- all_GEx_comparisons[, -c(2:7)]
head(all_GEx_comparisons[, c(1:7)]) # To delete still, these are the extra VDdef25 from first file read
head(all_GEx_comparisons[, c(8:13)]) # These now correspond to VDdef25, see:
grep(pattern = 'VDdef25', colnames(all_GEx_comparisons))
grep(pattern = '-0.102632226534278', all_GEx_comparisons[, 2]) # Should give 1188, ILMN_1664922
grep(pattern = '-0.102632226534278', all_GEx_comparisons[, 8]) # Should also give 1188, ILMN_1664922
head(all_GEx_comparisons[, c(14:19)]) # These now correspond to VDdef50
# Delete second set of dummy variables:
colnames(all_GEx_comparisons)[1:20]
all_GEx_comparisons <- all_GEx_comparisons[, -c(2:7)]
# Rename first column set back to all_baseline_VDdef25
colnames(all_GEx_comparisons)[2] <- "logFC.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[3] <- "AveExpr.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[4] <- "t.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[5] <- "P.Value.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[6] <- "adj.P.Val.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[7] <- "B.all_baseline_VDdef25"
colnames(all_GEx_comparisons)[1:20]

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

# Delete columns which are not needed
# All ave. exp. columns
dim(all_GEx_comparisons)
pattern <- 'AveExpr'
ave_expr_cols <- grep(pattern = pattern, colnames(all_GEx_comparisons))
colnames(all_GEx_comparisons)
colnames(all_GEx_comparisons)[ave_expr_cols]
head(all_GEx_comparisons[, c(ave_expr_cols)])
all_GEx_comparisons <- all_GEx_comparisons[, -c(ave_expr_cols)]
dim(all_GEx_comparisons)
colnames(all_GEx_comparisons)
head(all_GEx_comparisons)
grep(pattern = pattern, colnames(all_GEx_comparisons))
#############################


#############################
# Calculate average expression per gene per arm and place in the first columns:

############
# This is from '03c_microarray_GEx_diff_expression_VD_def_all_no_covar.R'

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
load('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/R_session_saved_image_pheno_file_check.RData', verbose=T)
# load('R_session_saved_image_diff_expression.RData', verbose=T) # Has subsetted objects for array data.

# Read in the file with the sample membership information (group membership, replicate number, 
# treatment status, etc.) to be able to create a design and contrast (done above for PCA plots).

#Check dimensions between annotation file with meta-data (must have the same number of rows, otherwise
#errors downstream):
#TO DO: Raise error if not the same.

dim(membership_file_cleaned)
str(membership_file_cleaned)
head(membership_file_cleaned)
tail(membership_file_cleaned)
head(normalised_filtered)
dim(normalised_filtered)
str(normalised_filtered)
dim(normalised_filtered_annotated)

# Sanity check:
# TO DO: Raise error and stop if false:
identical(row.names(membership_file_cleaned), colnames(normalised_filtered))
length(which(row.names(membership_file_cleaned) %in% colnames(normalised_filtered)))
############


############
# Check GEx and cross match file (membership_file), then merge membership_file with GEx (normalised_filtered)
dim(membership_file_cleaned)
head(membership_file_cleaned)
dim(normalised_filtered)
head(normalised_filtered)
plyr::count(membership_file_cleaned$arm)
plyr::count(membership_file_cleaned$group_membership)

# New table to hold summary stats per probe, per arm, per timepoint:
by_group_GEx <- list()

# Subset kit_id array file for placebo at baseline:
by_group_GEx$placebo_baseline <- normalised_filtered[,
                                                             which(colnames(normalised_filtered) 
                                                                   %in% row.names(
                                                                     membership_file_cleaned[which(
                                                                       membership_file_cleaned$group_membership == 'baseline_placebo'),
                                                                       ])
                                                                   )
                                                             ]
# For placebo at 12 months:
by_group_GEx$placebo_final <- normalised_filtered[,
                                                             which(colnames(normalised_filtered) 
                                                                   %in% row.names(
                                                                     membership_file_cleaned[which(
                                                                       membership_file_cleaned$group_membership == 'final_placebo'),
                                                                       ])
                                                             )
                                                             ]
# Subset kit_id array file for 2000 IU at baseline:
by_group_GEx$baseline_2000 <- normalised_filtered[,
                                                          which(colnames(normalised_filtered) 
                                                                %in% row.names(
                                                                  membership_file_cleaned[which(
                                                                    membership_file_cleaned$group_membership == 'baseline_2000'),
                                                                    ])
                                                          )
                                                          ]
# Subset kit_id array file for 2000 IU at 12 months:
by_group_GEx$final_2000 <- normalised_filtered[,
                                                          which(colnames(normalised_filtered) 
                                                                %in% row.names(
                                                                  membership_file_cleaned[which(
                                                                    membership_file_cleaned$group_membership == 'final_2000'),
                                                                    ])
                                                          )
                                                          ]
# Subset kit_id array file for 4000 IU at baseline:
by_group_GEx$baseline_4000 <- normalised_filtered[,
                                                       which(colnames(normalised_filtered) 
                                                             %in% row.names(
                                                               membership_file_cleaned[which(
                                                                 membership_file_cleaned$group_membership == 'baseline_4000'),
                                                                 ])
                                                       )
                                                       ]
# Subset kit_id array file for 4000 IU at 12 months:
by_group_GEx$final_4000 <- normalised_filtered[,
                                                       which(colnames(normalised_filtered) 
                                                             %in% row.names(
                                                               membership_file_cleaned[which(
                                                                 membership_file_cleaned$group_membership == 'final_4000'),
                                                                 ])
                                                       )
                                                       ]
# Check subsets:
sapply(by_group_GEx, dim)
plyr::count(membership_file_cleaned$group_membership)
dim(normalised_filtered)

# Sanity:
identical(row.names(by_group_GEx$placebo_baseline), row.names(by_group_GEx$placebo_final))
identical(row.names(by_group_GEx$placebo_baseline), row.names(by_group_GEx$baseline_2000))
identical(row.names(by_group_GEx$placebo_baseline), row.names(by_group_GEx$final_2000))
identical(row.names(by_group_GEx$placebo_baseline), row.names(by_group_GEx$baseline_4000))
identical(row.names(by_group_GEx$placebo_baseline), row.names(by_group_GEx$final_4000))


# Get means and SD for each probe:
summary_by_group_GEx <- data.frame()

summary_by_group_GEx <- data.frame('baseline_placebo_mean' = rowMeans(by_group_GEx$placebo_baseline))
summary_by_group_GEx$baseline_placebo_SD <- apply(by_group_GEx$placebo_baseline, 1, sd, na.rm = FALSE)
summary_by_group_GEx$baseline_placebo_var <- apply(by_group_GEx$placebo_baseline, 1, var, na.rm = FALSE)

summary_by_group_GEx$final_placebo_mean <- rowMeans(by_group_GEx$placebo_final)
summary_by_group_GEx$final_placebo_SD <- apply(by_group_GEx$placebo_final, 1, sd, na.rm = FALSE)
summary_by_group_GEx$final_placebo_var <- apply(by_group_GEx$placebo_final, 1, var, na.rm = FALSE)

summary_by_group_GEx$baseline_2000_mean <- rowMeans(by_group_GEx$baseline_2000)
summary_by_group_GEx$baseline_2000_SD <- apply(by_group_GEx$baseline_2000, 1, sd, na.rm = FALSE)
summary_by_group_GEx$baseline_2000_var <- apply(by_group_GEx$baseline_2000, 1, var, na.rm = FALSE)

summary_by_group_GEx$final_2000_mean <- rowMeans(by_group_GEx$final_2000)
summary_by_group_GEx$final_2000_SD <- apply(by_group_GEx$final_2000, 1, sd, na.rm = FALSE)
summary_by_group_GEx$final_2000_var <- apply(by_group_GEx$final_2000, 1, var, na.rm = FALSE)

summary_by_group_GEx$baseline_4000_mean <- rowMeans(by_group_GEx$baseline_4000)
summary_by_group_GEx$baseline_4000_SD <- apply(by_group_GEx$baseline_4000, 1, sd, na.rm = FALSE)
summary_by_group_GEx$baseline_4000_var <- apply(by_group_GEx$baseline_4000, 1, var, na.rm = FALSE)

summary_by_group_GEx$final_4000_mean <- rowMeans(by_group_GEx$final_4000)
summary_by_group_GEx$final_4000_SD <- apply(by_group_GEx$final_4000, 1, sd, na.rm = FALSE)
summary_by_group_GEx$final_4000_var <- apply(by_group_GEx$final_4000, 1, var, na.rm = FALSE)

dim(summary_by_group_GEx)
head(summary_by_group_GEx)
tail(summary_by_group_GEx)
colnames(summary_by_group_GEx)
row.names(summary_by_group_GEx)[1:10]

# Sanity:
identical(row.names(by_group_GEx$placebo_baseline), row.names(summary_by_group_GEx))

# Basic distribution plots for a few
hist(as.numeric(by_group_GEx$placebo_baseline[1, ]))
hist(as.numeric(by_group_GEx$placebo_baseline[2, ]))
hist(summary_by_group_GEx$baseline_placebo_mean)
hist(summary_by_group_GEx$baseline_placebo_SD)
hist(summary_by_group_GEx$baseline_placebo_var)
############

############
# Merge tables
head(all_GEx_comparisons)
head(summary_by_group_GEx)
dim(all_GEx_comparisons)
dim(summary_by_group_GEx)
colnames(all_GEx_comparisons)
colnames(summary_by_group_GEx)
row.names(all_GEx_comparisons)[1:10]
all_GEx_comparisons$probe_ID[1:10]
row.names(summary_by_group_GEx)[1:10]

# Sanity:
identical(as.character(all_GEx_comparisons$probe_ID),
          as.character(row.names(summary_by_group_GEx)))
row.names(all_GEx_comparisons) <- as.character(all_GEx_comparisons$probe_ID)
summary_by_group_GEx$probe_ID <- row.names(summary_by_group_GEx)

identical(row.names(all_GEx_comparisons), row.names(summary_by_group_GEx))
identical(as.character(all_GEx_comparisons$probe_ID), as.character(summary_by_group_GEx$probe_ID))

# Merge to get final table
# Sort out variable types and column to merge by:
head(summary_by_group_GEx$probe_ID)
head(all_GEx_comparisons$probe_ID)
length(unique(summary_by_group_GEx$probe_ID))
dim(summary_by_group_GEx)

all_GEx_comparisons$probe_ID <- as.character(all_GEx_comparisons$probe_ID)
all_GEx_comparisons$probes_by_symbol <- as.character(all_GEx_comparisons$probes_by_symbol)
all_GEx_comparisons$probes_by_ENTREZID_RE <- as.character(all_GEx_comparisons$probes_by_ENTREZID_RE)
str(summary_by_group_GEx)
str(all_GEx_comparisons)

# Merge:
final_table <- merge(summary_by_group_GEx, all_GEx_comparisons, by = 'probe_ID')
dim(final_table)
dim(summary_by_group_GEx)
dim(all_GEx_comparisons)
head(final_table)[, c(1:5, 15:20)]
tail(final_table)[, c(1:5, 15:20)]

# Check order of columns after merging
colnames(final_table)
# Move 
final_table <- final_table[, moveme(names(final_table), "probes_by_symbol after probe_ID")]
final_table <- final_table[, moveme(names(final_table), "probes_by_ENTREZID_RE after probes_by_symbol")]
colnames(final_table)
colnames(final_table)[1:10]
############

############
# Create separate table for ST1B with 'main comparisons' only:
main_comps <- c(
  'joint_pairedtreated_2000+4000',
  'joint_pairedtreated_placebo',
  'treatment_joint_vd_def_50treated_2000+4000',
  'treatment_joint_vd_def_50treated_placebo',
  'treatment_joint_vd_deltatreated_2000+4000',
  'UI2000minusplacebo',
  'UI4000minusplacebo')

get_main_comps <- function(pattern, df, df2) {
  cols <- grep(pattern = pattern, colnames(df), fixed = TRUE)
  colnames(df)[cols]
  df_temp <- df[, c(1, cols)]
  df2 <- merge(df2, df_temp, by = 'probe_ID')
  return(df2)
}

colnames(final_table)[1:5]
ST1B_GEx_main_comparisons <- as.data.frame(final_table[, c('probe_ID',
                                                           'probes_by_symbol',
                                                           'probes_by_ENTREZID_RE')])
head(ST1B_GEx_main_comparisons)
for (i in main_comps) {
  print(i)
  ST1B_GEx_main_comparisons <- get_main_comps(i, final_table, ST1B_GEx_main_comparisons)
}

dim(ST1B_GEx_main_comparisons)
head(ST1B_GEx_main_comparisons)
tail(ST1B_GEx_main_comparisons)
colnames(ST1B_GEx_main_comparisons)
############

############
# Get version without variance columns for final_table (which is ST1C):
cols <- grep(pattern = '_var', fixed = TRUE, colnames(final_table))
colnames(final_table)[cols]
dim(final_table)
final_table_no_var <- final_table[, -c(cols)]
dim(final_table_no_var)
grep(pattern = '_var', fixed = TRUE, colnames(final_table_no_var))
############

#############################

#############################
# Write files to disk:
# write.table(final_table, '../ST1C_GEx_all_comparisons.tsv', sep = '\t', 
#             quote = F, na = 'NA', col.names = NA, row.names = TRUE)
write.table(final_table_no_var, '../ST1C_GEx_all_comparisons.tsv', sep = '\t', 
            quote = F, na = 'NA', col.names = NA, row.names = TRUE)
write.table(ST1B_GEx_main_comparisons, '../ST1B_GEx_main_comparisons.tsv', sep = '\t', 
            quote = F, na = 'NA', col.names = NA, row.names = TRUE)
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
# saved_image name gets overwritten as loading a different image above, renaming:
R_session_saved_image <- paste('R_session_saved_image_diff_expression_merge_tables', '.RData', sep='')
save.image(file = R_session_saved_image, compress = 'gzip')
sessionInfo()

q()

# Next: run script for higher level analyses.
#############################