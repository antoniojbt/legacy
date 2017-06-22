#############################
# Script to process Illumina gene expression probe annotations and to subset SNP locations for input into MatrixeQTL
# Antonio J Berlanga-Taylor
# 29 Sept 2015

#############################


#############################
# Requires a probe file as input (from Illumina arrays) and a SNP file (e.g. from dbSNP)
# 

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/")
# setwd("/Users/antoniob/Desktop/BEST_D_03_MAR.DIR/")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_",".txt", sep=""), open='a')
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
#load('R_session_saved_image_read_and_QC.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_old_illumina_probe_locs', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:

# vignette("lumi")
# source("http://bioconductor.org/biocLite.R")
# biocLite('qdap')

library(illuminaHumanv4.db)
# library(lumiHumanAll.db)
# library(qdap)
library(reshape2)
library(data.table)

# update.packages('illuminaHumanv4.db')
# sessionInfo()

#############################


#############################

# Get the probe identifiers:
gen_loc <- illuminaHumanv4GENOMICLOCATION
str(gen_loc)
gen_loc_mapped_probes <- mappedkeys(gen_loc)
probe_IDs <- (unlist(gen_loc_mapped_probes))
head(probe_IDs)
length(probe_IDs)

# Get locations of probes:
genomic_locs <- unlist(mget(probe_IDs, illuminaHumanv4GENOMICLOCATION, ifnotfound = NA))
head(genomic_locs)
str(genomic_locs)
length(genomic_locs)
# View(genomic_locs)


# Identify probes of poor quality:
probes_by_qual  <- unlist(mget(probe_IDs, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
head(probes_by_qual)
bad_probes  <- probes_by_qual == "No match" | probes_by_qual == "Bad"
head(bad_probes)
str(bad_probes)
length(which(bad_probes == TRUE))


# Eliminate probes of poor quality:

# View(bad_probes)
bad_probes_df <- as.data.frame(bad_probes)
bad_probes_df['probe_IDs'] <- row.names(bad_probes_df)
# bad_probes_df <- bad_probes_df[which(bad_probes_df == TRUE), ]
# head(bad_probes_df)
# dim(bad_probes_df)
# length(bad_probes_df[, 2])

good_probes_df <- bad_probes_df[which(bad_probes_df == FALSE), ]
head(good_probes_df)
dim(good_probes_df)


genomic_locs_df <- as.data.frame(genomic_locs)
genomic_locs_df['probe_IDs'] <- row.names(genomic_locs_df)
head(genomic_locs_df)
# View(genomic_locs)

# Multi-mapping probes:
multi_mapping <- unlist(mget(probe_IDs, illuminaHumanv4SECONDMATCHES, ifnotfound = NA))
head(multi_mapping)
length(which(is.na(multi_mapping) == TRUE))
uniquely_mapping_df <- as.data.frame(multi_mapping[is.na(multi_mapping)])
head(uniquely_mapping_df)
uniquely_mapping_df['probe_IDs'] <- row.names(uniquely_mapping_df)
head(uniquely_mapping_df)

multi_mapping_2 <- unlist(mget(probe_IDs, illuminaHumanv4OTHERGENOMICMATCHES, ifnotfound = NA))
head(multi_mapping_2)
length(which(is.na(multi_mapping_2) == TRUE))
uniquely_mapping_2_df <- as.data.frame(multi_mapping_2[is.na(multi_mapping_2)])
head(uniquely_mapping_2_df)
uniquely_mapping_2_df['probe_IDs'] <- row.names(uniquely_mapping_2_df)
head(uniquely_mapping_2_df)
dim(uniquely_mapping_2_df)


repeat_masked <- unlist(mget(probe_IDs, illuminaHumanv4REPEATMASK, ifnotfound = NA))
head(repeat_masked)
length(which(is.na(repeat_masked) == TRUE))
repeat_masked_df <- as.data.frame(repeat_masked[is.na(repeat_masked)])
head(repeat_masked_df)
repeat_masked_df['probe_IDs'] <- row.names(repeat_masked_df)
head(repeat_masked_df)
dim(repeat_masked_df)


# Probes to keep that don't multi-map or are of bad quality:
probes_keep <- merge(good_probes_df, genomic_locs_df, by.x = 'probe_IDs')
head(probes_keep)
dim(probes_keep)

probes_keep <- merge(probes_keep, uniquely_mapping_df, by = 'probe_IDs')
head(probes_keep)
dim(probes_keep)
# View(probes_keep)

probes_keep <- merge(probes_keep, uniquely_mapping_2_df, by = 'probe_IDs')
head(probes_keep)
dim(probes_keep)
# View(probes_keep)

probes_keep <- merge(probes_keep, repeat_masked_df, by = 'probe_IDs')
head(probes_keep)
dim(probes_keep)
# View(probes_keep)

gene_locations_MatrixEQTL <- probes_keep[, c(1,3)]
head(gene_locations_MatrixEQTL)
dim(gene_locations_MatrixEQTL)
# View(gene_locations_MatrixEQTL)

# Clean probes that still have more than one genomic location:
nchar(as.character(gene_locations_MatrixEQTL[150, 2]))
length(which(nchar(as.character(gene_locations_MatrixEQTL[, 2])) > 40))
gene_locations_MatrixEQTL[2102, ]
gene_locations_MatrixEQTL[25300, ]
gene_locations_MatrixEQTL[22367, ]
still_multi_mapping <- which(nchar(as.character(gene_locations_MatrixEQTL[, 2])) > 40)
still_multi_mapping

gene_locations_MatrixEQTL_clean <- gene_locations_MatrixEQTL[-still_multi_mapping, ]
head(gene_locations_MatrixEQTL_clean)
dim(gene_locations_MatrixEQTL_clean)

# Split locations and process to match MatrixEQTL format
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html#own

gene_locations_MatrixEQTL_split <- transform(gene_locations_MatrixEQTL_clean, FOO = colsplit(
  gene_locations_MatrixEQTL_clean[, 2], ':', names = c('chr', 'start', 'end', 'strand')))

head(gene_locations_MatrixEQTL_split)
dim(gene_locations_MatrixEQTL_split)

gene_locations_MatrixEQTL_split_MatrixEQTL <- gene_locations_MatrixEQTL_split[, c(1, 3, 4, 5)]
head(gene_locations_MatrixEQTL_split_MatrixEQTL)
dim(gene_locations_MatrixEQTL_split_MatrixEQTL)

# Write to file:
write.table(gene_locations_MatrixEQTL_split_MatrixEQTL, 'old_illumina_probes_genomic_locations.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
#############################


#############################
## Save all except multi-mapping instead:
# Clean probes that still have more than one genomic location and split locations and process to match MatrixEQTL format (no multi-mapping or repeatmask)
# Identify probes of poor quality:
probes_by_qual_loose  <- unlist(mget(probe_IDs, illuminaHumanv4PROBEQUALITY, ifnotfound = NA))
head(probes_by_qual_loose)
bad_probes_loose  <- probes_by_qual_loose == "No match" #| probes_by_qual == "Bad"
head(bad_probes_loose)
str(bad_probes_loose)
length(which(bad_probes_loose == TRUE))

# Eliminate probes of poor quality:

# View(bad_probes)
bad_probes_loose_df <- as.data.frame(bad_probes_loose)
bad_probes_loose_df['probe_IDs'] <- row.names(bad_probes_loose_df)
# bad_probes_df <- bad_probes_df[which(bad_probes_df == TRUE), ]
# head(bad_probes_df)
# dim(bad_probes_df)
# length(bad_probes_df[, 2])

good_probes_loose_df <- bad_probes_loose_df[which(bad_probes_loose_df == FALSE), ]
head(good_probes_loose_df)
dim(good_probes_loose_df)


genomic_locs_loose_df <- as.data.frame(genomic_locs)
genomic_locs_loose_df['probe_IDs'] <- row.names(genomic_locs_loose_df)
head(genomic_locs_loose_df)
# View(genomic_locs)

# Probes to keep that don't multi-map or are of bad quality:
probes_keep_loose <- merge(good_probes_loose_df, genomic_locs_loose_df, by.x = 'probe_IDs')
head(probes_keep_loose)
dim(probes_keep_loose)

gene_locations_MatrixEQTL_loose <- probes_keep_loose[, c(1,3)]
head(gene_locations_MatrixEQTL_loose)
dim(gene_locations_MatrixEQTL_loose)
# View(gene_locations_MatrixEQTL)

# Clean probes that still have more than one genomic location:
nchar(as.character(gene_locations_MatrixEQTL_loose[150, 2]))
length(which(nchar(as.character(gene_locations_MatrixEQTL_loose[, 2])) > 40))
still_multi_mapping_loose <- which(nchar(as.character(gene_locations_MatrixEQTL_loose[, 2])) > 40)
still_multi_mapping_loose

gene_locations_MatrixEQTL_clean_loose <- gene_locations_MatrixEQTL_loose[-still_multi_mapping_loose, ]
head(gene_locations_MatrixEQTL_clean_loose)
dim(gene_locations_MatrixEQTL_clean_loose)

# Split locations and process to match MatrixEQTL format
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html#own

gene_locations_MatrixEQTL_split_loose <- transform(gene_locations_MatrixEQTL_clean_loose, FOO = colsplit(
  gene_locations_MatrixEQTL_clean_loose[, 2], ':', names = c('chr', 'start', 'end', 'strand')))

head(gene_locations_MatrixEQTL_split_loose)
dim(gene_locations_MatrixEQTL_split_loose)

gene_locations_MatrixEQTL_split_MatrixEQTL_loose <- gene_locations_MatrixEQTL_split_loose[, c(1, 3, 4, 5)]
head(gene_locations_MatrixEQTL_split_MatrixEQTL_loose)
dim(gene_locations_MatrixEQTL_split_MatrixEQTL_loose)

gene_locations_MatrixEQTL_split_MatrixEQTL_loose$length <- gene_locations_MatrixEQTL_split_MatrixEQTL_loose$FOO.end -
  gene_locations_MatrixEQTL_split_MatrixEQTL_loose$FOO.start

summary(gene_locations_MatrixEQTL_split_MatrixEQTL_loose$length)

# Write to file:
write.table(gene_locations_MatrixEQTL_split_MatrixEQTL_loose[, 1:4], 'old_illumina_probes_genomic_locations_loose.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')

#############################


#############################
# Subset SNP file to those in the genotype file:

genotyped_SNPs <- fread('P140343-Results_FinalReport_clean-individuals_and_SNPs.map', sep = '\t', header = FALSE, stringsAsFactors = FALSE)
head(genotyped_SNPs)
dim(genotyped_SNPs)
str(genotyped_SNPs)
genotyped_SNPs_IDs <- genotyped_SNPs[, 2, with = FALSE]
head(genotyped_SNPs_IDs)
dim(genotyped_SNPs_IDs)
colnames(genotyped_SNPs_IDs) <- 'ID'
head(genotyped_SNPs_IDs)

all_SNPs <- fread('snp142_MatrixEQTL_snp_pos.txt', sep = ' ', header = FALSE, stringsAsFactors = FALSE)
colnames(all_SNPs)[1] <- 'ID'
str(all_SNPs)
head(all_SNPs)

# Set keys for indexing and later merging with data.table:
setkey(all_SNPs, ID)
setkey(genotyped_SNPs_IDs, ID)

# Merge to get SNP locations for genotyped data:
SNP_locations <- all_SNPs[genotyped_SNPs_IDs]
# merge(all_SNPs, genotyped_SNPs_IDs, by = 'ID')
SNP_locations

# Eliminate SNPs without information:
SNPs_to_remove <- which(is.na(SNP_locations[, 2, with = FALSE]))
length(SNPs_to_remove)

SNP_locations <- SNP_locations[!SNPs_to_remove, ]
SNP_locations
dim(SNP_locations)

head(SNP_locations, 200)
dim(SNP_locations)

# Write to file:
write.table(SNP_locations, 'SNP_genomic_locations.txt', row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')

#############################


#############################
#The end:
# TO DO: remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes)))

#rm()

# To save R workspace with all objects to use at a later time:
# Also consider dump(), dput() and dget(), see http://thomasleeper.com/Rcourse/Tutorials/savingdata.html
save.image(file=R_session_saved_image_full, compress='gzip')

# Or save specific objects:
# objects_to_save <- (c('normalised_expressed', 'membership_file_cleaned', 'FAILED_QC'))
# save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: collect data, plots and tables for a report
#############################