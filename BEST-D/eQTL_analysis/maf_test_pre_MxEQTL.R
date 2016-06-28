########################
# Script to check mior allele frequencies from 0, 1 and 2 genotype codings
# Antonio J Berlanga-Taylor
# 8 June 2016
# The second allele should correspond to the minor.
# Claculate MAF and check if any is > 0.50
########################

########################
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting')
getwd()

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
########################


########################
# Get file:
snp_file_name <- as.character(args[1])
# snp_file_name <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv'

# Rows should be SNP IDs and columns sample IDs
########################

########################
# Read file:
snp_file <- fread(snp_file_name, sep = '\t', header = TRUE, stringsAsFactors = F)
setkey(snp_file, FID)
snp_file
dim(snp_file)
str(snp_file)
snp_file[1:5, 1:5, with = F]
########################


########################
# Check minor allele frequencies:
snp_file[1:5 ,1:5, with = F]
head(snp_file[, -1, with = F])
# snp_file_df <- as.data.frame(snp_file)
str(snp_file)
length(which(complete.cases(snp_file)))
maf <- rowMeans(snp_file[, -1, with = F], na.rm = TRUE)/2 # First column has the SNPs IDs
head(maf)
class(maf)
length(maf)
length(which(maf < 0.50))
length(which(maf > 0.50))

# Flip to get correct minor allele frequencies printed:
maf <- pmin(maf, 1-maf)
head(maf)
length(which(maf < 0.50))
length(which(maf > 0.50))

# Match back to SNP IDs:
maf_IDs <- cbind(snp_file[, 1, with = F], maf)
maf_IDs
dim(maf_IDs)
str(maf_IDs)
setnames(maf_IDs, 'FID', 'rsID')
setnames(maf_IDs, 'maf', 'MAF')
########################

########################
# Save to file:
options(digits = 2)
write.table(maf_IDs, sprintf('MAF_%s', snp_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
########################

########################
q()

########################