#############################################
# Bacterial GWAS
# Antonio J Berlanga-Taylor
# 10 Aug 2015

# Current input: genotyping data, annotation files


# Outputs: various plots and tables for 
# See libraries and packages required below

# Inputs are:
# Genotype file
# Covariates file

# e.g.
# ln -s 

#############################################


#############################################
# Check 
# http://adegenet.r-forge.r-project.org/files/simGWAS/MSc-GWAS-1.1.pdf
# http://www.r-phylo.org/wiki/Main_Page
# https://github.com/thibautjombart/adegenet/wiki/Tutorials

# Other sources to check:
# 

#############################################


#############################################
## :

# 1) 
# 2) 
# 3) 

#############################################


#############################################
##Set working directory and file locations and names of required inputs:

# Set arguments to run from the command line:
# e.g.
# args <- commandArgs(trailingOnly = TRUE)
# rnorm(n=as.numeric(args[1]), mean=as.numeric(args[2]))

args <- commandArgs(trailingOnly = TRUE)
# args[1] is for genotype input file
# args[2] ...

# Working directory:
# working_dir <- ('/ifs/projects/proj018/runs_second_round.dir/GWAS.dir/adegenet.dir')
# setwd(working_dir)


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
# load('R_session_saved_image.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_load_data','.RData', sep='')
R_session_saved_image

#############################################


#############################################
## Update packages if necessary:
# source("http://bioconductor.org/biocLite.R")
# biocLite()
# install.packages("ade4", dep=TRUE)
# install.packages("adegenet", dep=TRUE)

#and load them:
# library()
packages <-c('dplyr', 'ggplot2', 'reshape', 'tidyr', 'knitr', 'readr', 'ade4', 'adegenet', 'gtools')
lapply(packages, require, character.only = TRUE)

#############################################


#############################################


#############################################
# Read data in:
# TO DO: pass to configuration file
# help.search("Hardy-Weinberg")

# Genetic data:

read_geno <- read.table(as.character(args[1]))
# TO DO:This doesn't have the locations (ie the .map file is needed, important as adegenet result will then require to map back 
# column number to chromosome location).

#head(read_table)
str(read_geno)
read_geno[1:5, 1:5]

# Split sample names for easier viewing (labels are dragged from fastq, trimming, mapping, variant calling, etc.). 
# These should match the phenotype file (i.e. the 'unique_ID_fastq' column name):

# TO DO: Rename fastqs/bams/etc before to avoid all this...
new_names <- strsplit(as.character(read_geno[, 1]), fixed = TRUE, split='_')
new_names <- as.data.frame(rapply(new_names, function(x) x[[2]]))
head(new_names)
str(new_names)

# Replace sample names:
read_geno[, 1] <- new_names
read_geno[1:5, 1:5]
read_geno[, ncol(read_geno)]
str(read_geno)

# TO DO: Pass as args:
# Read in phenotype file:
file_to_read <- 'metadata_strain_ids_pseudomonas_final_13Aug2015.csv.tr'
read_pheno <- read.csv(file_to_read, header = TRUE)
head(read_pheno)

# Match phenotype and genotype file order:
read_pheno_matched <- read_pheno[match(read_geno[, 1], read_pheno[, 1]), ]
head(read_pheno_matched)
tail(read_pheno_matched)
as.data.frame(read_pheno_matched[,1])
str(read_pheno_matched)

# Groupings for time since sampling:
# six_months <- which(read_pheno_matched[, column_of_interest] <= 0.5)
# one_year <- which(read_pheno_matched[, column_of_interest] > 0.5 & 
#                     read_pheno_matched[, column_of_interest] <= 1)
# six_years <- which(read_pheno_matched[, column_of_interest] > 1 &
#                      read_pheno_matched[, column_of_interest] <= 6)
# more_than_6_years <-which(read_pheno_matched[, column_of_interest] > 6)
# first_sample <- which(read_pheno_matched[, 20] == 1)
# last_sample <- which(read_pheno_matched[, 20] == 2)
# 
# length(six_months)
# date_bins <- list(six_months, one_year, six_years, more_than_6_years, 
#                   first_sample, last_sample)
# date_bins
# sapply(date_bins, length)

# Recode and add as a variable (first and last are coded in separate variable (column 20, otherwise levels are overwritten:
read_pheno_matched$time_bins[read_pheno_matched[, 10] <= 0.5] <- 'six_months'
head(read_pheno_matched$time_bins)
length(which(read_pheno_matched$time_bins == 'six_months'))

read_pheno_matched$time_bins[read_pheno_matched[, 10] > 0.5 & read_pheno_matched[, 10] <= 1] <- 'one_year'
head(read_pheno_matched$time_bins)
length(which(read_pheno_matched$time_bins == 'one_year'))

read_pheno_matched$time_bins[read_pheno_matched[, 10] > 1 & read_pheno_matched[, 10] <= 6] <- 'six_years'
read_pheno_matched$time_bins[read_pheno_matched[, 10] > 6] <- 'more_than_6_years'

summary(as.factor(read_pheno_matched$time_bins))

# Recode morphology column for case/control type analysis:
summary(as.factor(read_pheno_matched$morphology))
head(read_pheno_matched)
colnames(read_pheno_matched)

read_pheno_matched$morphology_binary[read_pheno_matched[, 3] == 'colony'] <- 1
read_pheno_matched$morphology_binary[read_pheno_matched[, 3] == 'mucoid'] <- 2
read_pheno_matched$morphology_binary[is.na(read_pheno_matched[, 25])] <- 0
read_pheno_matched$morphology_binary
colnames(read_pheno_matched)

# Dump to file for plink:
order_columns <- c(5,1,3,9,20,21,22,23,24,25)
unlist(colnames(read_pheno_matched))[order_columns]
read_pheno_matched_plink <- read_pheno_matched[unlist(colnames(read_pheno_matched))[order_columns]]
head(read_pheno_matched_plink)
str(read_pheno_matched_plink)

# Rename first two columns:
names(read_pheno_matched_plink)[1] <- 'FID'
names(read_pheno_matched_plink)[2] <- 'IID'

colnames(read_pheno_matched_plink)
head(read_pheno_matched_plink)

# Save as table in file:
write.table(read_pheno_matched_plink, 'read_pheno_matched_plink.tsv', sep='\t', 
            quote = FALSE, col.names = TRUE, row.names = FALSE)

# TO DO: Pass as args:
# Identify column with phenotype of interest to test:
colnames(read_pheno_matched)
morphology <- 3
time_of_sampling <- 10
patient_ID <- 5
time_bins <- 24
# column_of_interest <- time_of_sampling
# column_of_interest <- patient_ID
column_of_interest <- time_bins
column_of_interest
phenotype_name <- colnames(read_pheno_matched)[column_of_interest]
phenotype_name

# For time since sampling, round numbers:
# read_pheno_matched[, column_of_interest] <- round(read_pheno_matched[, column_of_interest], 2)
# head(read_pheno)
# str(read_pheno)
# head(read_pheno_matched[, column_of_interest])


# Replace sample names in the phenotype file:
# Unnecesary as the 'unique_ID_fasq' is already like this:
# read_pheno_matched[, 1] <- new_names
# read_pheno_matched[1:5, ]
# 
# read_geno[1:5, 1:10]

# Check sample names match across files:
which(read_geno[, 1] %in% read_pheno_matched[, 1])
cbind(as.character(read_geno[, 1]), as.character(read_pheno_matched[, 1]))

# identical() nees variables to be of the same type (?):
identical(as.character(read_geno[, 1]), as.character(read_pheno_matched[, 1]))


# Eliminate duplicates, technical replicates, etc. according to analysis being done:
# TO DO: Check time of sampling and what it actually means (may be sample sample from petri dish but different morphotype, e.g.
# S04 and S05 are the same time and patient but different morphotypes):
# TO DO: pass this as argument to run from the command line:
samples_remove <- c('S15b')
samples_remove
rows_to_remove <- which(read_pheno_matched[, 1] %in% samples_remove)
rows_to_remove
# Check in genotype file:
which(read_geno[, 1] %in% samples_remove)

# Remove from pheno and geno files:
read_pheno_matched <- read_pheno_matched[-rows_to_remove, ]
str(read_pheno_matched)

read_geno <- read_geno[-rows_to_remove, ]
str(read_geno)

# Sanity check for sample names and order:
identical(as.character(read_geno[, 1]), as.character(read_pheno_matched[, 1]))


# Convert data to adegent format:
# ?df2genind
data_to_adegent <- df2genind(read_geno, ploidy=1, sep='', NA.char = c('0'))
data_to_adegent

validObject(data_to_adegent)
data_to_adegent$tab[1:5, 1:5]
head(data_to_adegent$loc.n.all)
head(data_to_adegent$loc.fac)
head(data_to_adegent$all.names)
data_to_adegent$ploidy
data_to_adegent$call
summary(data_to_adegent)

#############################


#############################
# The end:
# Remove objects that are not necessary to save:
ls()
object_sizes <- sapply(ls(), function(x) object.size(get(x)))
as.matrix(rev(sort(object_sizes))[1:20])

# rm()

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for xxx.
#############################

