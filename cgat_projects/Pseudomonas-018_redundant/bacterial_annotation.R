#############################
# SNP, gene and genome annotation scripts for bacterial strains

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
setwd('/ifs/projects/proj018/runs_second_round.dir/annotation.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output",".txt", sep=""), open='a')
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Load a previous R session, data and objects:
#load('R_session_saved_image.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
#R_session_saved_image <- paste('R_session_saved_image_normalisation', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image_annotation', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages('')
#detach("package:pryr", unload=TRUE)
#biocLite(c("GenomicFeatures", "AnnotationDbi"))
#biocLite("MeSH.Pae.PAO1.eg.db")
#biocLite("AnnotationHub")
#biocLite("FunctSNP")
#install.packages("FunctSNP",dependencies=TRUE)


library(AnnotationHub)
library(GenomicFeatures)
library(AnnotationDbi)
library(MeSH.Pae.PAO1.eg.db)
library(ggplot2)


#############################


#############################
# Get annotations for PAO1

ah <- AnnotationHub()
ah
#?AnnotationHub

snapshotVersion(ah)
snapshotDate(ah)
## how many resources?
length(ah)

## list currently active filters
filters(ah)

## list values that can be used to filter on:
columns(ah)
keytypes(ah)
head(ah)

## list possible values for one of these filter types
head(keys(ah, keytype="Species"))

## OR retrieve metadata values about several keys at once
## (This approach may not always scale the way you want it to)
metadata(ah, columns = c("Species","RDataPath"))

## create and apply a new filter to only include people
filters(ah) <- list(Species="Pseudomonas aeruginosa")
filters(ah)
display(ah)

## now how many resources are there?
length(ah)

## what are the names for these resources?
head(names(ah))

## What are the URLs for these resources?
head(snapshotUrls(ah))

## what web service is this AnnotationHub pointing to?
hubUrl()

## and more explicitly
snapshotUrl()

## Where are the files that get downloaded being cached?
## (there is also a setter if you wish to assign this to another location)
hubCache(ah)

## Download a resource (using tab completion) and put it into "res"
head(ah)
res <- ah$inparanoid8.Orthologs.hom.Pseudomonas_aeruginosa.inp8.sqlite
res
help('select')

## query and subset
query(ah, "GRCh37")    # 'GRCh37' anywhere in the metadata
subset(ah, Species == "Homo sapiens" & any(Tags == "FASTA"))

## For brevity, lets define a resource pathname we are interested in 
pathName <- "goldenpath.hg19.encodeDCC.wgEncodeUwTfbs.wgEncodeUwTfbsMcf7CtcfStdPkRep1.narrowPeak_0.0.1.RData"

## We can get metadata information about that resource:
ahinfo(ah, pathName)

## Or we can extract it by name like this (using a list semantic):
res <- ah[[pathName]]

## And we can also use "[" restrict the things that are in the
## AnnotationHub object.  Either by position:
subHub <- ah[1:3]

## Or by name
subHub <- ah[pathName]

unique(ah$dataprovider)

#############################









#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image_full, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

#############################