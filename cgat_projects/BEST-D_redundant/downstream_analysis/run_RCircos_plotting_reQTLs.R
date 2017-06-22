#################
# Script to run RCircos plot
# Antonio J Berlanga-Taylor
# 18 May 2016
# Requires processed MatrixEQTL output (or files with chr, start, end format)
# Outputs plots to disk
# Many parameters, currently set to plot from MatrixEQTL output to globally visualise links between 
# eSNPs and eGenes
#################


#################
# RCircos facilitates circular plotting and interface with R
# See:
# https://cran.r-project.org/web/packages/RCircos/vignettes/Using_RCircos.pdf
# Check the demo:
# demo("RCircos.Demo.Human")
# Copied to file 'RCircos_demo.R'. The code is the same as the vignette (copied/paraphrased in this file).

# RCircos implements Circos concepts/tools in R, Krzywinski et al 2009:
# http://genome.cshlp.org/content/19/9/1639.full
# http://circos.ca/documentation/

# Other R circular plotting tools are:
# https://bioinformatics.oxfordjournals.org/content/30/19/2811.full
# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3921174/
# http://bioinformatics.oxfordjournals.org/content/31/17/2912.full.pdf+html

# See other scripts with examples from RCircos demo and vignette:
# 'RCircos_plotting.R'
# 'RCircos_demo.R'
#################


#################
# RCircos needs data frame as input with chr, start and end 
# Plotting is carried out by:
# Initialising core components
# Initialising graphics device
# Specifying parameters
# Specifying plot data and locations
# Adding tracks sequentially
# Saving to disk
#################

#################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/scripts_to_upload/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))
options(echo = TRUE)

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_RCircos','.RData', sep='')
R_session_saved_image
####################

#################
# Load libraries
# source("https://bioconductor.org/biocLite.R")
# biocLite('RCircos')
# biocLite('cowplot')

library(RCircos)
library(data.table)
# library(cowplot) # Complements ggplot2 for publication ready figures (multi-panels, labels, etc.)
# See: https://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
#################

####################
# TO DO: change highlighted genes to cis, trans, or pathways
# TO DO: heatmap legend
# TO DO: scatter plot legend and looks ugly
# TO DO: increase margins
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)
# Data input as data frame (processed later for chr, strat, end):

# Cis eQTL file to get eSNPs and eGenes (as link lines for link plot?) 
# and p-values (as outer ring, eg Manhattan type plot):
cis_file <- as.character(args[1])
# cis_file <- 'biomart_annot_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_cis.txt.txt'

# eQTL trans file to place as ribbon lines with link plot:
trans_file <- as.character(args[2])
# trans_file <- 'biomart_annot_cis_tx_fdr5_reQTLs_annot_all_Tx_joint_trans.txt.txt'

# Gene expression levels as heatmap for cis genes, using FC file:
Gex_file <- as.character(args[3])
# Gex_file <- 'biomart_annot_full_topTable_pairing_all_treated.txt.txt'

# Selected genes (from pathways, disease, etc):
highlight_genes_file <- as.character(args[4])
# highlight_genes_file <- 'genes_in_pathway_of_interest.txt'
# highlight_genes_file <- 'RCircos_biomart_annot_GO_VD_proteins.txt.txt'
# highlight_genes_file <- cis_file

# TO DO: track numbering isn't clear: vignette skips tracks 3 and 4? Ideogram is counted?
# Initialise track number
inside_track_number <- as.integer(args[5])
# inside_track_number <- as.integer(5) #Counting ideogram, gene names, connectors, link lines and ribbon lines.

outside_track_number <- as.integer(args[6])
# outside_track_number <- as.integer(1) # For Manhattan type plot

chr_exclude_file <- as.character(args[7])
# chr_exclude <- c("chrX", "chrY", 'chrM')
# chr_exclude_file <- 'chrs_to_exclude_RCircos.txt'
####################


####################
# Read data:
cis_data <- fread(cis_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
cis_data
dim(cis_data)
# Convert to data.frame:
cis_data <- as.data.frame(cis_data)
head(cis_data)
colnames(cis_data)
# Select columns for use:
# links_plot_file <- cis_data[, c(13, 14, 15, 19, 20, 21)]  # For reQTL files
links_plot_file <- cis_data[, c("SNP_chr_name", "SNP_chrom_start", "SNP_chrom_end", 
                     "Probe_ID_chromosome_name", "Probe_ID_start_position", "Probe_ID_end_position")]
head(links_plot_file)
dim(links_plot_file)

# cis_data_pvalues:
colnames(cis_data) # columns differ between reQTLs and MxEQTL files
# cis_data_pvalues <- cis_data[, c(13, 14, 15, 9)] # For reQTL files
cis_data_pvalues <- cis_data[, c("SNP_chr_name", "SNP_chrom_start", "SNP_chrom_end", "p.value.final")] # Different to MxEQTL files
head(cis_data_pvalues)
dim(cis_data_pvalues)
# Convert to -log10
cis_data_pvalues$log10_pvalues <- -log10(cis_data_pvalues[, c("p.value.final")])
summary(cis_data_pvalues$log10_pvalues)
head(cis_data_pvalues)

trans_data <- fread(trans_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
trans_data
dim(trans_data)
# Convert to data.frame:
trans_data <- as.data.frame(trans_data)
head(trans_data)
colnames(trans_data)
# Select columns for use:
links_plot_file_trans <- trans_data[, c(13, 14, 15, 19, 20, 21)] # For reQTLs file
head(links_plot_file_trans)
dim(links_plot_file_trans)


Gex_data <- fread(Gex_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
Gex_data
dim(Gex_data)
# Convert to data.frame:
Gex_data <- as.data.frame(Gex_data)
head(Gex_data)
colnames(Gex_data)
# Select columns for use:
Gex_data <- Gex_data[, c(14, 15, 16, 4, 3)]
head(Gex_data)
dim(Gex_data)


highlight_genes_data <- fread(highlight_genes_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
highlight_genes_data
dim(highlight_genes_data)
# Convert to data.frame:
highlight_genes_data <- as.data.frame(highlight_genes_data)
head(highlight_genes_data)
colnames(highlight_genes_data)
# Sort to get the top hits by pvalue if using eQTL list:
highlight_genes_data <- highlight_genes_data[order(highlight_genes_data$FDR.final, highlight_genes_data$beta.final), ]
head(highlight_genes_data)
# Select columns:
head(highlight_genes_data)
colnames(highlight_genes_data)
highlight_genes_data <- highlight_genes_data[, c("Probe_ID_chromosome_name", "Probe_ID_start_position", 
                                                 "Probe_ID_end_position", "hgnc_symbol_biomart")]
head(highlight_genes_data)
# Remove duplicates:
highlight_genes_data <- highlight_genes_data[which(!duplicated(highlight_genes_data$hgnc_symbol_biomart)), ]
head(highlight_genes_data)
dim(highlight_genes_data)

chr_exclude <- read.csv(chr_exclude_file, header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
chr_exclude
class(chr_exclude)
chr_exclude <- as.list(as.data.frame(chr_exclude))
chr_exclude
class(chr_exclude)
chr_exclude <- chr_exclude$V1
chr_exclude
#################


#################
# Plot with RCircos
# Initialise core components and plot ideogram:

# Load ideograms from RCircos itself or use a text file manually. 
# Download cytoBandIdeo.txt.gz files from e.g.:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBandIdeo.txt.gz
# head(UCSC.HG19.Human.CytoBandIdeogram)

hg38 <- read.csv(file = 'cytoBandIdeo.txt', sep = '\t', header = FALSE, 
                 col.names = c("Chromosome", "ChromStart", "ChromEnd", "Band", "Stain"))
head(hg38)
dim(hg38)

# Check parameters:
RCircos.List.Parameters()
# Modify with:
# rcircos.params <- RCircos.Get.Plot.Parameters();
# rcircos.params$base.per.unit <- 30000;
# rcircos.params$radius.len <- 2.0;
# rcircos.params$track.padding <- 0.10;
# RCircos.Reset.Plot.Parameters(rcircos.params);
# RCircos.List.Parameters()
# TO DO: changing these gives errors

#   Setup RCircos core components:
#   1.  Chromosome ideogram plot information
#   2.  x and y coordinates for a circular line and degrees of the text rotation at each point
#   3.  Plot parameters for image plot control
RCircos.Set.Core.Components(cyto.info = hg38, chr.exclude = chr_exclude, 
                            tracks.inside = inside_track_number, 
                            tracks.outside = outside_track_number)

# Initialize Graphics Device:
out_file <- sprintf('rcircos_%s.pdf', cis_file) # Also supports pdf, tiff, svg, png, etc.
pdf(file = out_file, height = 14, width = 14, compress=TRUE);
RCircos.Set.Plot.Area()
title('Genome wide SNP/gene associations')

# RCircos.Set.Plot.Area() will setup plot area base on total number of tracks 
# inside and outside of chromosome ideogram. Modify this with:
# par(mai=c(0.25, 0.25, 0.25, 0.25));
# plot.new();
# plot.window(c(-2.5,2.5), c(-2.5, 2.5));

# Plot Chromosome Ideogram:
RCircos.Chromosome.Ideogram.Plot()
#################

#################
# Plot connectors and gene names:
# Use selected genes (from pathways, disease, etc):
# TO DO: create file based on XGR results and process to get chr coords.

# Connectors and gene names. Requires file with chr, start, end and gene.
# For best visualization, cex should be >= 0.4 when drawing gene labels
# head(RCircos.Gene.Label.Data)
RCircos.Gene.Connector.Plot(genomic.data = highlight_genes_data, track.num = 1, side = "in")
RCircos.Gene.Name.Plot(gene.data = highlight_genes_data, name.col = 4, track.num = 2, side = "in")
#################

#################
# Cis eQTL data p-values (as outer ring, eg Manhattan type plot):
# Scatterplot:
head(cis_data_pvalues)
RCircos.Scatter.Plot(scatter.data = cis_data_pvalues, data.col = 5, track.num = 1,
                     side = "out", by.fold = 2) # -log10(0.01) equals 2, p-value to show

# Histogram plot
# RCircos.Histogram.Plot(hist.data = cis_data_pvalues, data.col = 5, track.num = 1, side = "out")
#################

#################
# Gene expression levels as heatmap for cis genes, using FC file:
# Heatmap:
head(Gex_data)
RCircos.Heatmap.Plot(heatmap.data = Gex_data, data.col = 5, track.num = 3, side = "in")
#################

#################
# Links plots:
# Link lines are used for presenting relationship of two small genomic regions.
# Ribbons are plotted for bigger genomic regions. Link lines are always drawn inside of chromosome ideogram.
# Examples from RCircos:
# head(RCircos.Link.Data)
# head(RCircos.Ribbon.Data)

# Add link lines to the center of plot area, eQTL cis data:
# RCircos.Link.Plot(link.data = links_plot_file, track.num = 4, by.chromosome = FALSE)

# Add ribbon lines to the center of plot area, eQTL trans data:
head(links_plot_file_trans)
RCircos.Ribbon.Plot(ribbon.data = links_plot_file_trans,  track.num = 4, by.chromosome = FALSE, twist = FALSE)
#################

#################
# Other plots which can be added:
# # Line plot.
# RCircos.Line.Plot(line.data=RCircos.Line.Data, data.col=5, track.num=7, side="in")

# # Tile plot. Note: tile plot data have only chromosome locations and each data file is for one track
# RCircos.Tile.Plot(tile.data=RCircos.Tile.Data, track.num=9,  side="in")
#################

#################
# Close the graphic device:
dev.off()
#################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next run:
####################