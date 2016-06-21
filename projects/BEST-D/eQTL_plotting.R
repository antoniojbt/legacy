########################
# Script for plotting eQTLs following MatrixEQTL analysis
# Antonio J Berlanga-Taylor
# 30 Sept 2015

########################

########################
# setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis.dir/')
getwd()

# source("http://bioconductor.org/biocLite.R")
# biocLite('venneuler')
library(ggplot2)
library(data.table)
library(VennDiagram)

args <- commandArgs(trailingOnly = TRUE)
# args[1] = 'genotype_data_4000_final.tsv_MatrixEQTL_loose_1MB.cis'
# args[2] = 'genotype_data_4000_final.tsv'
# args[3] = 'subset_finalVisit_4000.tsv'
# args[4] = 'rs11150882'
# args[5] = 'ILMN_1707137'
# e.g.
# /Library/Frameworks/R.framework/Resources/bin/Rscript \
# eQTL_plotting.R \
# cut_genotype_data_all_treated_final.tsv_matched.tsv_MatrixEQTL_loose_p5_1MB.cis \
# cut_genotype_data_all_treated_final.tsv_matched.tsv \
# cut_GEx_treated_4000_and_2000.tsv_matched.tsv \
# rs4236435 \
# ILMN_1651819
########################


########################
## Plot expression levels against genotypes
# Set-up ggplot for plotting:

# Get eQTL file:
# TO DO: automate to plot top ten SNPs/probes:
eQTL_file_name <- as.character(args[1])
# eQTL_file_name <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv_MatrixEQTL_loose_p5_1MB.cis'

# Get SNP file and read:
snp_file_name <- as.character(args[2])
# snp_file_name <- 'cut_genotype_data_all_treated_final.tsv_matched.tsv'
snp_file <- fread(snp_file_name, sep = '\t', header = TRUE, stringsAsFactors = F)
setkey(snp_file, FID)
snp_file
dim(snp_file)
str(snp_file)
snp_file[1:5, 1:5, with = F]

# Get probe file and read:
probe_file_name <- as.character(args[3])
# probe_file_name <- 'cut_GEx_treated_4000_and_2000.tsv_matched.tsv'
probe_file <- fread(probe_file_name, sep = '\t', header = T, stringsAsFactors = F)
setkey(probe_file, FID)
dim(probe_file)
probe_file
probe_file[1:5, 1:5, with = F]

# Specify SNP to plot:
snp <- as.character(args[4])
# snp <- 'rs4236435'
snp_data <- as.data.frame(snp_file[snp, ])
head(snp_data)
dim(snp_data)

# Specify probe to plot:
probe <- as.character(args[5])
# probe <- 'ILMN_1651819'
probe_data <- as.data.frame(probe_file[probe, ])
head(probe_data)
dim(probe_data)

# Get label from eQTL file name, this breaks with changes to name. TO DO:
group_label <- strsplit(eQTL_file_name, split = '_')
# group_label <- paste(group_label[[1]][3], 
#                      strsplit(group_label[[1]][4], '\\.')[[1]][1], 
#                      group_label[[1]][7],
#                      sep = ' ')
group_label <- eQTL_file_name
group_label

identical(colnames(snp_data), colnames(probe_data))

snp_probe_data <- rbind(snp_data, probe_data)
head(snp_probe_data)
snp_probe_data <- t(snp_probe_data) #, check.names = FALSE)
snp_probe_data <- as.data.frame(snp_probe_data)
head(snp_probe_data)
class(snp_probe_data)
str(snp_probe_data)
colnames(snp_probe_data)[1] <- as.character(snp)
colnames(snp_probe_data)[2] <- as.character(probe)
colnames(snp_probe_data)
snp_probe_data <- snp_probe_data[-1, ]
head(snp_probe_data)
snp_probe_data[, 2] <- as.numeric(as.character(snp_probe_data[, 2]))
snp_probe_data[, 1] <- as.integer(as.character(snp_probe_data[, 1]))
snp_probe_data[, 1] <- as.factor(as.character(snp_probe_data[, 1]))
# snp_probe_data[, 1] <- as.factor(snp_probe_data[, 1])
str(snp_probe_data)
head(snp_probe_data)
colnames(snp_probe_data)

ggplot(snp_probe_data, aes(snp_probe_data[, 1], snp_probe_data[, 2])) + 
  geom_jitter(colour='darkgrey', position=position_jitter(width=0.25)) + 
  geom_boxplot(outlier.size=0, alpha=0.6, fill='grey') + 
  xlab(names(snp_probe_data)[1]) + 
  ylab(names(snp_probe_data)[2]) + ggtitle(paste(group_label)) + theme_bw() +
  theme(text = element_text(size = 6))

ggsave(plot = last_plot(), filename = paste(snp, '_', probe, '_', eQTL_file_name, '.png', 
                                            sep = ''))
dev.off()

########################


########################
# TO DO: move to another script:
## Venn diagrams:
# placebo_baseline = read.table('genotype_data_placebo_baseline.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# placebo_final =	read.table('genotype_data_placebo_final.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# t2000_baseline =	read.table('genotype_data_2000_baseline.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# t2000_final =	read.table('genotype_data_2000_final.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# t4000_baseline =	read.table('genotype_data_4000_baseline.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# t4000_final =	read.table('genotype_data_4000_final.tsv_MatrixEQTL_loose_1MB.cis', sep = '\t', check.names = FALSE, header = TRUE)
# 
# placebo_baseline[1:5, 1:5]
# 
# placebo_baseline['file'] <- 'placebo_baseline'
# placebo_final['file'] <- 'placebo_final'
# t2000_baseline['file'] <- 't2000_baseline'
# t2000_final['file'] <- 't2000_final'
# t4000_baseline['file'] <- 't4000_baseline'
# t4000_final['file'] <- 't4000_final'
# 
# head(placebo_baseline)
# 
# merged_eqtls <- rbind(placebo_baseline, placebo_final,
#                       t2000_baseline, t2000_final,
#                       t4000_baseline, t4000_final)
# head(merged_eqtls)
# tail(merged_eqtls)
# dim(merged_eqtls)
# colnames(merged_eqtls)
# merged_eqtls_subset <- merged_eqtls[, c(1,4)]
# venn.diagram(x = 'file', filename=merged_eqtls_subset)
# plot(venneuler(merged_eqtls_subset))

########################


########################
q()

########################