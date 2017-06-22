################
# Simple counts of number of individuals with genotype data pre and post QC
# For CONSORT diagram figure
# Antonio Berlanga
# BEST-D project
# 01 Nov 2016
################


################
# setwd('Downloads_to_delete/')
################

################
library(plyr)
################

################
# Phenotype data after merging all datasets sent by CTSU
# /ifs/projects/proj043/analysis.dir/association_analysis_6.dir/final_phenotype_data.tab
# /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport.fam
# /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.fam


pheno_data <- read.csv('final_phenotype_data.tab', sep = '\t', stringsAsFactors = FALSE,
                       header = TRUE, na.strings = c('-99999'))

geno_preQC <- read.csv('P140343-Results_FinalReport.fam', sep = ' ', stringsAsFactors = FALSE,
                       header = FALSE, na.strings = c('-99999'))

geno_postQC <- read.csv('P140343-Results_FinalReport_clean_SNPs_autosome_individuals.fam', sep = ' ', 
                        stringsAsFactors = FALSE, header = FALSE, na.strings = c('-99999'))
################

################
head(pheno_data)
dim(pheno_data)
# View(pheno_data)
count(pheno_data$placebo_v_Tx)
count(pheno_data$placebo_v_4000)
count(pheno_data$placebo_v_2000)
count(pheno_data$arm)
summary(pheno_data$vitd12_placebo_only)

# PreQC genotype data from plink fam files:
head(geno_preQC)
dim(geno_preQC)
names(geno_preQC)[1] <- 'FID'

# PostQC genotype data from plink fam files:
head(geno_postQC)
dim(geno_postQC)
names(geno_postQC)[1] <- 'FID'
################

################
# Genotype data corresponds to randomisation kit IDs, so all should be present in 
# phenotype file:
length(which(geno_postQC$FID %in% geno_preQC$FID))
length(which(geno_preQC$FID %in% pheno_data$FID))
length(which(geno_postQC$FID %in% pheno_data$FID))
################


################
# Merge sets so as to have arm and count for CONSORT diagram:
geno_preQC <- merge(pheno_data[, c('FID', 'arm', 'visit_type')], geno_preQC)
head(geno_preQC)
dim(geno_preQC)

geno_postQC <- merge(pheno_data[, c('FID', 'arm', 'visit_type')], geno_postQC)
head(geno_postQC)
dim(geno_postQC)

# All samples here are from randomisation kit IDs, so only count arm:
# 0 = 4000 ; 1 = 2000 ; 2 = placebo
count(pheno_data$arm)
count(geno_preQC$arm)
count(geno_postQC$arm)
################

################
sessionInfo()
q()
################