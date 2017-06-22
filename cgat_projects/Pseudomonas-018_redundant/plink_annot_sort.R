###########
# Script to process plink's output and sort by p-values 
# Antonio B.
###########

###########
# Set up command line arguments:
args <- commandArgs(trailingOnly = TRUE)

infile <- as.character(args[1])
# infile <- 'all_awk_loose_processed_filtered_MAF_PASS_IDs_AT.qfam.total.perm'
adj_p_value <- as.numeric(args[2])
# adj_p_value <- 10e-3
column_to_sort <- as.numeric(args[3])
# column_to_sort <-  9 # For FDR BH in plink2 assoc adjusted output
# column_to_sort <-  6 # For EMP1 in plink2 qfam output
###########

###########
# Clean up plink's output containing multiple spaces, etc.:
system_call <- sprintf("cat %s | tr -s ' ' | cut -d ' ' -f2- | tr ' ' '\t' | sort -nk%d 1<> %s.tr3.p-sorted",
                       infile, column_to_sort, infile)
system_call
system(system_call)
###########

###########
# Read the new file now sorted by p-value:
infile2 <- sprintf("%s.tr3.p-sorted", infile)

assoc_data <- read.csv(infile2, sep = '\t', header = F, na.strings = c('NA', '', ' '))
head(assoc_data)
###########

###########
# Get significant results, NAs (from variants which couldn't be tested) come up on top, so clean up:
data_hits <- assoc_data[which(!is.na(assoc_data[, column_to_sort])), ]
colnames(data_hits) # Colnames get messed up with unix sorting, they'll appear at the end of NAs.
data_hits[, column_to_sort] <- as.numeric(as.character(data_hits[, column_to_sort])) # Colnames get coerced into NA.
head(data_hits)
dim(data_hits)
head(data_hits[, column_to_sort])
class(data_hits[, column_to_sort])

# Keep only variants passing p-value specified at args:
data_hits_adjp <- data_hits[which(data_hits[, column_to_sort] < adj_p_value), ]
head(data_hits_adjp)
###########

###########
# Write to file:
write.table(data_hits_adjp, paste(infile2, '_adj_P_', adj_p_value, '.tsv', sep = ''),
            sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
system(sprintf('head %s_adj_P_%G.tsv', infile2, adj_p_value))
###########

###########
q()
###########
