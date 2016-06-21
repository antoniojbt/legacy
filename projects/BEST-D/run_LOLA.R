#############################################
# Enrichment analysis using LOLA
# Antonio J Berlanga-Taylor
# 4 June 2016
# Requires:
# region database
# Regions of interest
# Universe
#############################################

#############################################
# See:
# http://bioinformatics.oxfordjournals.org/content/32/4/587.long
# http://databio.org/lola/
# https://github.com/sheffien/LOLA#lola-core-database
# http://databio.org/simplecache/
# Vignettes:
# https://github.com/sheffien/LOLA/blob/master/vignettes/gettingStarted.Rmd
# https://github.com/sheffien/LOLA/blob/master/vignettes/usingLOLACore.Rmd
# https://github.com/sheffien/LOLA/blob/master/vignettes/choosingUniverse.Rmd
#############################################

# TO DO: this vs GAT?
# TO DO: Create folder with universe(s) and databases to test against to re-use against multiple regions of interest.
# Pass GAT backgrounds and files to this:
# '/Users/antoniob/Desktop/Enrichment_test_databases'
# It contains raw text files which could be passed to GAT
# Create a frequency matched universe (bin variants by frequency and only pass those regions?)
# Ancient alleles? 
# Also see regioneR:
# http://bioinformatics.oxfordjournals.org/content/32/2/289.full.pdf+html
# Helper functions may be useful but doesn't seem different to GAT (provides sampling strategies though?)
# For isochores, see:
# https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
# Usually, isochores are defined as genomic regions of homogeneous G+C content. 
# In GAT, isochores can mean any genomic regions with a shared property. Isochores are used to eliminate the effect of 
# confounders in an analysis.

#############################################
##Set working directory and file locations and names of required inputs:

# Working directory:
#setwd('/ifs/projects/proj043/analysis.dir/eqtl_analysis_5.dir/')
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/scripts_to_upload/counting_and_plotting/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_LOLA",Sys.Date(),".txt", sep=""))
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
R_session_saved_image <- paste('R_session_saved_image_LOLA','.RData', sep='')
R_session_saved_image
####################


####################
# Load packages:
# source("https://bioconductor.org/biocLite.R")
# biocLite("LOLA")
# Use simpleCache to download databases:
# https://github.com/sheffien/simpleCache
# require(devtools)
# install_github("sheffien/simpleCache") #public

library(LOLA)
library(simpleCache)
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
# Run LOLA
# Load a region database:
dbPath = system.file("extdata", "hg19", package="LOLA")
regionDB = loadRegionDB(dbPath)
names(regionDB)
# dbLocation: A string recording the location of the database folder you passed to loadRegionDB().
# collectionAnno: A data.table annotating the collections, with rows corresponding to the rows in your collection annotation files in the database.
# regionAnno: A data.table annotating each region set, with rows corresponding to bed files in the database (there is also a collection column recording which collection each region set belongs to).
# regionGRL: A GRangesList object holding the actual regions, with one list element per region set, ordered as in regionAnno.

# Load data for regions of interest and universe:
data("sample_input", package="LOLA") # load userSets
data("sample_universe", package="LOLA") # load userUniverse

# Other universes have been created and curated, such as DHS sites (to use vs TF binding data for example):
activeDHS = readBed("lola_vignette_data/activeDHS_universe.bed")
# Created by doing:
# activeDHS = unlist(regionDB$regionGRL[which(regionDB$regionAnno$collection == "sheffield_dnase")])
# activeDHS = disjoin(activeDHS)
# activeDHS

# userSets and userUniverse are GRanges objects, required to run the enrichment calculation: 
# calcLocEnrichment() ?? Function not available
# which will test the overlap between userSets and each region set in the regionDB.
# ??

# Test for pairwise overlap between each user set and each region set in regionDB
# using a Fisher's exact test to assess signficance of the overlap:
locResults = runLOLA(userSet, userUniverse, regionDB, cores=1)

# Results are a data.table with several columns:
colnames(locResults)
head(locResults)


# Columns userSet and dbSet are indexes into the respective GRangeList objects, identifying each pairwise comparison. 
# Columns describing the results: pValueLog, logOdds, and the actual values from the contingency table 
# (support is the overlap, and b, c, and d complete the 2x2 table). 
# Rank columns simply rank the tests by pValueLog, logOdds, or support; 
# following these are a series of columns annotating the database regions, depending on how you populated the index table 
# in the regionDB folder.

# Rank with different orders:
locResults[order(support, decreasing=TRUE), ]

# Order by one of the rank columns:
locResults[order(maxRnk, decreasing=TRUE), ]

# Write to file:
writeCombinedEnrichment(locResults, outFolder= "lolaResults")

# Use includeSplits parameter to also print out additional tables that are subsetted by userSet
# so that each region set you test has its own result table:
writeCombinedEnrichment(locResults, outFolder= "lolaResults", includeSplits=TRUE)

# Get he regions responsible for enrichment:
oneResult = locResults[2,]
extractEnrichmentOverlaps(oneResult, userSets, regionDB)

#Extract specific regions, eg regions from the "vistaEnhancers" region set:
getRegionSet(regionDB, collections="ucsc_example", filenames="vistaEnhancers.bed")

# Or just give the path to the database and LOLA will only load the specific region set(s) of interest:
getRegionSet(dbPath, collections="ucsc_example", filenames="vistaEnhancers.bed")
# Multiple filenames and collections can be passed


####################

####################
# Save to file:
write.table(shared_pairs, sprintf('shared_pairs_%s', output_file_name), 
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

# Next: run XGR, IPA, etc, cross with GWAS, ENCODE, etc.
####################