########################
# Script for counting SNPs and probes/genes following MatrixEQTL analysis
# Antonio J Berlanga-Taylor
# 30 Sept 2015

########################


########################
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')

library(ggplot2)
# library(data.table)

args <- commandArgs(trailingOnly = TRUE)
# args[1] = 'genotype_data_2000_baseline.tsv_MatrixEQTL_loose_1MB.cis'

########################


########################
# Basic counting:

file_eqtl_name <- as.character(args[1])
file_eqtl <- read.table(file_eqtl_name, sep = '\t', header = TRUE)
head(file_eqtl)
dim(file_eqtl)

tissue <- file_eqtl_name
tissue

order <- c(order(file_eqtl[, 6], decreasing = FALSE))
file_eqtl_FDR <- file_eqtl[order, ]
head(file_eqtl_FDR)

total_snps <- length(file_eqtl$SNP)
total_snps

total_genes <- length(file_eqtl$gene)
total_genes

unique_snps <- unique(file_eqtl$SNP)
head(unique_snps)
length(unique_snps)

unique_genes <- unique(file_eqtl$gene)
head(unique_genes)
length(unique_genes)

# Save results to file:
cat(file = 'unique_eQTL_counts_all.txt', tissue, 
    '\t', 'unique_snps', '\t', length(unique_snps), 
    '\t', 'unique_genes', '\t', length(unique_genes), 
    '\t', 'total_snps_(snp/gene pairs)', '\t', total_snps, 
    '\t', 'total_genes/probes_(snp/gene pairs)',  '\t', total_genes, 
    '\n', append = TRUE)

########################

########################
q()

########################
