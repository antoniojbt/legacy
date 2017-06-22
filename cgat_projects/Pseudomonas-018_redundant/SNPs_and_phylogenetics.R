#############################
# Phylogenetics scripts of bacterial strains

#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
setwd('/ifs/projects/proj018/runs_second_round.dir/samtools_bwa.dir')

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
R_session_saved_image_full <- paste('R_session_saved_image', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages('')
#detach("package:pryr", unload=TRUE)
#biocLite('SNPRelate')
#biocLite("gdsfmt")

library(SNPRelate)
library(ggplot2)

# SNPRelate webpages 
# https://github.com/zhengxwen/SNPRelate
# http://corearray.sourceforge.net/tutorials/SNPRelate/

# https://www.biostars.org/p/83232/


#############################

# vcf to GDS
# snpgdsVCF2GDS('trimmomatic_S02_TP_D7_001_TP_D5_001-1-1.bowtie.vcf.gz', 'trimmomatic_S02_TP_D7_001_TP_D5_001-1-1.gds')
# snpgdsSummary('trimmomatic_S02_TP_D7_001_TP_D5_001-1-1.gds')
# genofile <- openfn.gds('trimmomatic_S02_TP_D7_001_TP_D5_001-1-1.gds')

snpgdsVCF2GDS('all-samples-merged.mod.vcf', 'all-samples-merged.gds') #, method='biallelic.only')
snpgdsSummary('all-samples-merged.gds')
genofile <- snpgdsOpen('all-samples-merged.gds')


#Sample information
sample_info <- read.csv('SNPRelate_ids.mod.txt', header=T, sep='\t')
head(sample_info)

sample.id <- read.gdsn(index.gdsn(genofile, 'sample.id'))

#pop_code <- read.gdsn(index.gdsn(genofile, 'pop.id'))
#pop_code <- scan('SNPRelate_ids.txt', what=character(), sep='\t')


# dendogram
dissMatrix  <-  snpgdsDiss(genofile , sample.id=NULL, snp.id=NULL, autosome.only=F,remove.monosnp=TRUE, maf=NaN, missing.rate=NaN, 
                           num.thread=10, verbose=TRUE)
snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)
cutTree <- snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, n.perm = 5000, samp.group=NULL,col.outlier="red", col.list=NULL, pch.outlier=4, 
                         pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)

#pca
pca <- snpgdsPCA(genofile, autosome.only=F)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# make a data.frame when no prior population information is available:
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

png('pca_1_SNPRelate.png', width = 12, height = 12, units = 'in', res = 300)
plot(tab$EV2, tab$EV1, xplab="eigenvector 2", ylab="eigenvector 1")
dev.off()

#With population information:
color_list <- c("gray", "blue", "green", "red", "blue", "yellow",
                 'black', 'orange', 'pink', 'lightblue', 
                'cyan', 'brown', 'gold', 'limegreen', 'sienna', 'navy')

tab2 <- data.frame(sample.id = pca$sample.id, 
                   pop = factor(sample_info$pop.id)[match(pca$sample.id, sample_info$sample.id)], 
                   EV1 = pca$eigenvect[,1], 
                   EV2 = pca$eigenvect[,2],
                   stringsAsFactors = FALSE)

head(tab2)
head(sample_info)
png('pca_2_SNPRelate.png', width = 12, height = 12, units = 'in', res = 300)
plot(tab2$EV2, tab2$EV1, col=color_list[as.integer(tab2$pop)], xlab="eigenvector 2", ylab="eigenvector 1")
legend("topleft", legend=levels(tab2$pop), pch="o", col=1:nlevels(tab2$pop))
dev.off()


# Other way with IBS, ??
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, autosome.only=F, num.thread=2))
rv <- snpgdsCutTree(ibs.hc)

png('dendrogram_SNPRelate_IBS.png', width = 12, height = 12, units = 'in', res = 300)
plot(rv$dendrogram, leaflab='perpendicular', main='P. aeruginosa')
#legend("topright", legend=levels(tab2$pop), col=1:nlevels(tab2$pop))
dev.off()

#Plot the principal component pairs for the first four PCs:
png('First_4_PCAs_SNPRelate.png', width = 12, height = 12, units = 'in', res = 300)
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab2$pop, labels=lbls)
dev.off()

#Given two or more populations, Fst can be estimated by the method of Weir & Cockerham (1984).

# Two populations:
head(sample_info)
flag <- sample_info$pop.id %in% c("P001", "P002")
samp.sel <- sample_info[flag,]
head(samp.sel)
snpgdsFst(genofile, sample.id=samp.sel$sample.id, autosome.only=F, 
          population=as.factor(samp.sel$pop.id), method="W&C84")

snpgdsClose(genofile)
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

# Next: run the script for xxx.
#############################
