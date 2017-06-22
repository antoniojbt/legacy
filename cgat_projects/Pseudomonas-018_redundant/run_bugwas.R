#############################################
# GWAS analysis of bacteria using bugwas
# Antonio J Berlanga-Taylor
# 8 June 2016
#############################################

#############################################
## General steps:
# Outline of bacterial GWAS (also see tutorial for pipeline bacterial GWAS):
  # Set up - assemble data files, check dependencies (see below for links and libraries to install).
  # Run basic GWAS analysis: kmer GWAS using chi-square test. 
  # Generate kmers fo each genome to be tested (from bam or fasta files)
  # Generate file with path to each kmer file and a 2nd column with phenotype (eg 0/1 for resistance/sensitivity)
  # Run kmerAnalysis.R script
  # Visualise significant kmers
  
  # Control for population structure with a linear mixed model (running Gemma)
  
  # Identify lineage-associated as well as locus-specific variants with bugwas
  
  # Annotate significant/top kmers with BLAST.  

# TO DO / add:
  # Correct for multiple testing 
  # Assess low-frequency variants
  # Pan genome GWAS:
    # https://github.com/jessiewu/bacterialGWAS/blob/master/manual/pangenomeGWAS_user_manual.pdf
  # SNP GWAS: 
    # https://github.com/jessiewu/bacterialGWAS/blob/master/manual/GWAS_for_bacteria_user_manual.pdf
    # Input:
      # Requires mapped fasta files (genomic sequences in fasta format and that the genomic sequences have been
      # all mapped/aligned to a "reference" such that all the sequences in the fasta format are the same length 
      # (as the reference) and that the nucleotides at each site are homologous (the nth position in the all your genomes 
      # corresponds to the nth position in the reference). 
      # This requires an alignment including all your genomes and the reference genome.
      # Phylogeny is created on the fly but best to run before, check and pass as argument
  # Kmer GWAS (same manual as SNP):
      # https://github.com/jessiewu/bacterialGWAS/blob/master/manual/GWAS_for_bacteria_user_manual.pdf
      # kmers are created on the fly
#############################################


#############################################
# Requires:
# From https://github.com/jessiewu/bacterialGWAS
  # R scripts: kmerAnalysis.R, kmerAnnotation.R, kmerLMM.R, LMM_kmerAnnotation.R
# From https://github.com/sgearle/bugwas
  # bugwas library in R
# NCBI blast
# And see dependencies below for R packages

# Files neeed:
  # genotypes
  # phenotype
  # relatedness matrix 
  # path to the GEMMA software.
  
  # Kmer files - 1 genome per file. These were generated from raw sequencing reads using the kmer counting software dsk.
    # https://gatb.inria.fr/software/dsk/
    # https://github.com/GATB/dsk
  # List of kmer file paths linked to phenotypes.
  # Relatedness matrix calculated from genomic SNPs (this is created by GEMMA).
  # Maximum-likelihood phylogeny calculated from genomic SNPs (we use PhyML for < 100 genomes and RAxML for > 100 genomes)
  # BLAST databases (nt)
  # Gene look-up table for annotation for genome of interest
#############################################

#############################################
# See paper:
  # http://www.nature.com/articles/nmicrobiol201641
  # http://www.danielwilson.me.uk/virulogenomics.html
  # https://github.com/sgearle/bugwas
  # https://github.com/sgearle/bugwas/blob/master/manual/bugwas.pdf

# Pipeline that uses bugwas, other scripts (needed here) and tutorial:
  # https://github.com/jessiewu/bacterialGWAS
  # https://github.com/janepipistrelle/bacterial_GWAS_tutorial/blob/master/tutorial.rmd
#############################################


#############################################
### Tutorial introduction:
# https://github.com/janepipistrelle/bacterial_GWAS_tutorial/blob/master/tutorial.rmd
# Bacterial GWAS problems:
# Different species vary greatly in how often their genomes recombine and most species have only a single, 
# circular chromosome. This means that even variants that are not in close physical proximity on the chromosome are 
# in linkage disequilibrium. 
# Bacterial populations often exhibit strong signals of structure, 
# due to expansion of ecologically successful clones in free-living species and isolation in 
# different hosts in pathogens. 
# Finally, individuals of the same bacterial species often vary quite dramatically in gene 
# content, so just looking at SNPs called relative to a reference genome risks ignoring much interesting variation. 

# The pipeline uses variation in SNPs, gene presence/absence and 31nt kmers to carry out a GWAS which accounts for the above.
#############################################

#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd('/ifs/projects/proj018/')
# setwd('/Users/antoniob/Desktop/018_25_FEB.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_bugwas",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output"))

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
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_bugwas','.RData', sep='')
R_session_saved_image
####################


####################
# Get packages:
# Get bugwas and bacterialGWAS pipeline:
# system('git clone https://github.com/jessiewu/bacterialGWAS.git')
# system('git clone https://github.com/sgearle/bugwas.git')
# sudo apt-get install -y ncbi-blast+
# For Mac OSX, follow:
# https://www.blaststation.com/intl/members/en/howtoblastmac.html
# Download nt database, these are pre-formatted only inflate with tar:
# ftp://ftp.ncbi.nih.gov/blast/db/
# Create db if using a fasta file

# R dependencies:
# source('http://bioconductor.org/biocLite.R')
# biocLite('genoPlotR')
# biocLite("ape")
# biocLite("phangorn")
# install.packages("/Users/antoniob/Applications/apps.dir/bugwas/build/bugwas_1.0.tar.gz", repos = NULL, type="source")

# Load packages:
library('genoPlotR')
library("ape")
library("phangorn")
library('bugwas')
####################


####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# MatrixeQTL files:
eQTL_file1 <- as.character(args[1])
# eQTL_file1 <- '2000+4000-baseline-1.eQTL_trans'

eQTL_file2 <- as.character(args[2])
# eQTL_file2 <- '2000+4000-baseline-1.eQTL_cis'

total_SNPs_tested <- as.numeric(args[3])
# total_SNPs_tested <- 477422
####################

####################
# Set up file naming for outputs:
# to create: 2000+4000-baseline-1.trans_in_cis
eQTL_file1_base <- strsplit(eQTL_file1, '[.]')
eQTL_file1_base <- eQTL_file1_base[[1]][1]
eQTL_file1_base
eQTL_file1_ext <- strsplit(eQTL_file1, '_')
eQTL_file1_ext <- eQTL_file1_ext[[1]][2]
eQTL_file1_ext

eQTL_file2_base <- strsplit(eQTL_file2, '[.]')
eQTL_file2_base <- eQTL_file2_base[[1]][1]
eQTL_file2_base
eQTL_file2_ext <- strsplit(eQTL_file2, '_')
eQTL_file2_ext <- eQTL_file2_ext[[1]][2]
eQTL_file2_ext

output_file_name <- sprintf('%s_in_%s_%s.eQTL', eQTL_file1_ext, eQTL_file2_ext, eQTL_file1_base)
output_file_name
####################

####################
# Read data:
eQTL_data1 <- fread(eQTL_file1, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data1, 'gene', 'Probe_ID')
dim(eQTL_data1)
colnames(eQTL_data1)
setkey(eQTL_data1, SNP, Probe_ID)

eQTL_data2 <- fread(eQTL_file2, sep = '\t', header = TRUE, stringsAsFactors = FALSE)
setnames(eQTL_data2, 'gene', 'Probe_ID')
dim(eQTL_data2)
colnames(eQTL_data2)
setkey(eQTL_data2, SNP, Probe_ID)

tables()

eQTL_data1
eQTL_data2
####################

####################
## Run basic GWAS analysis: kmer GWAS using chi-square test. 
  # Generate kmers fo each genome to be tested (from bam or fasta files)
  # Generate file with path to each kmer file and a 2nd column with phenotype (eg 0/1 for resistance/sensitivity)


## Run kmerAnalysis.R script:
  # Use kmerAnalysis.R to test for patterns of presence or absence of each kmer in each genome 
  # and test each kmer for a significant association with the phenotype using a Chi-Square test:

#Bash Shell
# Rscript bacterialGWAS/kmerGWAS/kmerAnalysis.R -dataFile \
# ~/data/fus300.kmerfiles.txt -prefix fus300 \
# -removeKmerTxt FALSE -minCov 5 -externalSoftware dependencies.txt

## Visualise significant kmers:
# This generates p-values from the chi-square test for each kmer and two plots:
# Empirical cumulative distribution function (ECDF) plot (p-values for each kmer, ordered by significance)
# Quantile-quantile (QQ) plot comparing the distribution of -log10 p-values for the kmers in our dataset 
# to a theoretical distribution. Red line is the expectation.
####################


####################
## Control for population structure with a linear mixed model (running Gemma).
# Run our GWAS using a control for population structure with a linear mixed model (LMM) 
# implemented in the software GEMMA (Zhou & Stephens 2012). 
# This assigns all variants a background significance level and tests significance of each variant against the background level. 
# This acts to remove variants that are associated with specific lineages.

# This required the output of kmerAnalysis.R and an additional file:
# relatedness_matrix.txt 
# calculated from genomic SNPs from the same dataset.

#Bash Shell
# Rscript bacterialGWAS/kmerGWAS/kmerLMM.R \
#   -chisqStat fus300.gwaskmer-out.chisqStat.txt \
#   -patternKey fus300.gwaskmer-out.patternKey.txt \
#   -patternIndex fus300.gwaskmer-out.patternIndex.txt  \
#   -signif 5000 \
#   -relateMatrix ~/data/fus300_gemma_relmatrixout.cXX.txt \
#   -phenotype ~/data/fus300.pheno.txt \
#   -prefix fus300 \
#   -externalSoftware dependencies.txt

# Outputs a plot comparing the p-values of our kmers before and after controlling for population structure:
# Top kmers in unadjusted and adjust should be the same with drop in significance.
####################


####################
## Identify lineage-associated as well as locus-specific variants with bugwas
# Define variables:

gem.path = "bugwas/gemma/gemma.0.93b"
output.dir = "./"

# Call the function “lin_loc” which tests for both lineage and locus effects and generates several plots.
# This function needs the genotypes, phenotype, the relatedness matrix we used earlier and the path to the GEMMA software.
lin_loc(gen = "data/fus300_bugwas_gemma_gen_format.txt", pheno = "data/fus300_bugwas.pheno.txt",
        phylo = "data/RAxML_bestTree.fus300", prefix = "fus300", gem.path = gem.path,
        var.matrix = "data/fus300_bugwas.var.matrix.txt", relmatrix = "data/fus300_gemma_relmatrixout.cXX.txt",
        output.dir = output.dir)

## Plots output:
# Lineages that are associated with phenotype of interest:
# eg check xxx_tree_branchescolouredbyPC.png
# Lineages are defined by principal components (PCs) which can split groups of resistant and sensitive isolates.

# Variants are associated with the phenotype:
# xxx_genVar1_ManhattanLMMPvalues.png
# Lineage-specific variants associated with phenotype (in green) 
# and variants that are associated with the phenotype (high -log10 p-value) but which are not significantly associated 
# with any lineage (not shaded). 

# TO DO: map kmers corresponding to the locus-specific effects and annotate
####################


####################
# Annotate significant/top kmers with BLAST.
# Align kmers to a database of reference genomes using BLAST and check positions against a look-up table of annotated genes
# Kmers that do not have a good hit to reference genome are BLASTed against the whole nucleotide database from NCBI. 

#Bash Shell
# # Rscript bacterialGWAS/kmerGWAS/kmerAnnotation.R \
#   -chisq_results fus300.gwaskmer-out.chisqStat.txt \
#   -kmer_results fus300.gwaskmer-out.kmer.txt \
#   -present_ctrl fus300.gwaskmer-out.nPresentCtrl.txt  \
#   -present_case  fus300.gwaskmer-out.nPresentCase.txt \
#   -blastdb1 ~/data/staphdb3.blast \
#   -blastdb2 ~/data/nt_blastdb \
#   -ncbi_db ~/data/geneannot.allstaph.txt \
#   -signif 5000 \
#   -nproc 2 \
#   -prefix fus300 \
#   -externalSoftware dependencies.txt

# Annotate the significant kmers as before:

#Bash Shell
# Rscript LMM_kmerAnnotation.R \
#   -chisq_results fus300.gwaskmer-out.chisqStat.txt \
#   -kmer_results fus300.gwaskmer-out.kmer.txt \
#   -present_ctrl fus300.gwaskmer-out.nPresentCtrl.txt \
#   -present_case fus300.gwaskmer-out.nPresentCase.txt ~/data/dbs/staphdb3.blast \
#   -blastdb2 ~/data/dbs/nt_blastdb \
#   -ncbi_db ~/data/dbs/geneannot.allstaph.txt \
#   -signif 5000 \
#   -nproc 2 \
#   -prefix fus300 \
#   -LMM_kmers fus300_lmm_kmers_used.txt \
#   -LMM_output fus300_lmm_LMM_allkmers_out.txt \
#   -externalSoftware dependencies.txt

# TO DO: Annotate in other ways: blast against reference genome, database of interest (onown resistance genes, etc.)
####################


####################
# TO DO / add:
# Correct for multiple testing 
# Assess low-frequency variants
####################


####################
# Save to file:
write.table(xxx, sprintf('xxx_%s', output_file_name), 
            sep = '\t', na = 'NA', quote = FALSE, col.names = TRUE, row.names = FALSE)
####################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run downstream analysis, GAT, etc. Check bacterial specific annotations, etc.
####################
