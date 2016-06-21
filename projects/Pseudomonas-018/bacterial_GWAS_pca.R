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

# Working directory:
# working_dir <- ('/ifs/projects/proj018/runs_second_round.dir/GWAS.dir/adegenet.dir/test.dir')
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
# load('R_session_saved_image_run_analysis.RData', verbose=T)
load('R_session_saved_image_load_data.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_pca','.RData', sep='')
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
# Check data and re-set phenotype variables if needed:
# TO DO: Pass as args and move from load_data.R file to here:
# Check phenotype file:
head(read_pheno_matched)

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
read_pheno_matched$time_bins[read_pheno_matched[, 10] > 6] <- 'more_than_6_years '

summary(as.factor(read_pheno_matched$time_bins))

# Change phenotypes of interest if necessary:
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

head(read_pheno_matched)

# Adjust colour palette to grouping according to phenotype of interest:
colours <- rownames(as.table(summary(as.factor(read_pheno_matched[, column_of_interest]))))
# View(as.table(summary(as.factor(read_pheno_matched[, column_of_interest]))))
str(colours)
#colours <- as.factor(colours)
# Order for time bins:
colours <- factor(colours, ordered = TRUE, labels = c("six_months", "one_year",  "six_years", "more_than_6_years"))
str(colours)
colours
colour_palette <- rainbow(length(colours), alpha = 0.5)
colour_palette
palette(colour_palette)

head(read_pheno_matched)
tail(read_pheno_matched)
#############################################


#############################################
## High-level assessment of genetic diversity
# Run PCA, this displays a barplot of eigenvalues and asks for number of axis to retain (keep the largest). You can specify how many in the function, otherwise
# it asks interactively once plotted. There is no answer as to how many axes to keep. Typically axes that sum up to >80%, 95%, etc. of variance. 
# Others quote ~40% or whenever the cumulative variance flattens:
# ??dudi.pca

run_pca <- assign(paste('pca_', phenotype_name, sep = ''), 
                  dudi.pca(data_to_adegent, scale = FALSE, scannf = FALSE, nf = 12))
run_pca


# Kept 12 axes here.

# Eigenvalues (amount of information contained in each PC):
run_pca$eig
# Principal components:
run_pca$li
# Principal axes (loadings of the variables):
run_pca$c1

# Visualise data representing PCs:
# TO DO: run PCs kept for all groupings available and plot. 
png(paste('pca_1_slabel_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.label(run_pca$li, sub="PCA - PC 1 and 2")
add.scatter.eig(run_pca$eig,4,1,2, ratio=.3, posi="topleft")
dev.off()
# Genetic relationships between isolates? 
# Distinct lineages of bacteria?  How many?

# Get quantitative differences of the clustering with squared Euclidean distances (dist):
D <- dist(run_pca$li[,1:4])^2
D

# Use hierarchical clustering with complete linkage (hclust) to 
# define tight clusters:
clust <- hclust(D, method="complete")
clust

# Visualise number of clusters:
png(paste('clust_ade4_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
plot(clust, main="Clustering (complete linkage) based on the first 12 PCs", cex=.4)
dev.off()

# TO DO: define where to cut:
# Define clusters based on the dendrogram:
pop <- factor(cutree(clust, k=11))
pop

# Show groups on top of the PCs (clusters as colours and ellipses):
png(paste('populations_on_PCs_1_and_2_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, fac=pop, col=transp(funky(5)), cpoint=2, sub="PCA - axes 1 and 2")
dev.off()

# Repeat for PCs 3 and 4:
png(paste('populations_on_PCs_3_and_4_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=3, yax=4, fac=pop, col=transp(funky(5)),
        cpoint=2, sub="PCA - axes 3 and 4")
dev.off()

# Run through all axes kept from PCA or until all groups can be visualised.
png(paste('populations_on_PCs_5_and_6_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=5, yax=6, fac=pop, col=transp(funky(5)),
        cpoint=2, sub="PCA - axes 5 and 6")
dev.off()

png(paste('populations_on_PCs_7_and_8_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=7, yax=8, fac=pop, col=transp(funky(5)),
        cpoint=2, sub="PCA - axes 7 and 8")
dev.off()

png(paste('populations_on_PCs_9_and_10_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=9, yax=10, fac=pop, col=transp(funky(5)),
        cpoint=2, sub="PCA - axes 9 and 10")
dev.off()

png(paste('populations_on_PCs_11_and_12_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=11, yax=12, fac=pop, col=transp(funky(5)),
        cpoint=2, sub="PCA - axes 11 and 12")
dev.off()

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
