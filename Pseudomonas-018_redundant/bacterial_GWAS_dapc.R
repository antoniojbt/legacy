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
# working_dir <- ('/ifs/projects/proj018/runs_second_round.dir/GWAS.dir')
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
load('R_session_saved_image_load_data.RData', verbose=T)
# load('R_session_saved_image.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_dapc','.RData', sep='')
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

# Change phenotypes of interest if necessary:
#TO DO:
#unique_ID_fastq, antibiotics
#Done:
#morphology,Patient_ID,time_of_inf_years,Phenotype_TFAM_first_and_last,RAPD_type,time by bins,
#Phenotype_TFAM_chronic_1y,Phenotype_TFAM_chronic_6y,Phenotype_TFAM_6months

colnames(read_pheno_matched)
# phenotype_of_interest <- 'time_of_inf_years'
# For time, round to one digit:
# read_pheno_matched[, column_of_interest] <- round(read_pheno_matched[, column_of_interest], 0)
# read_pheno_matched[1:5, keep_columns]

# phenotype_of_interest <- 'Patient_ID'
# phenotype_of_interest <- 'time_bins'
# phenotype_of_interest <- 'Phenotype_TFAM_first_and_last'
# phenotype_of_interest <- 'RAPD_type'
# phenotype_of_interest <- 'Phenotype_TFAM_chronic_1y'
# phenotype_of_interest <- 'Phenotype_TFAM_chronic_6y'
phenotype_of_interest <- 'Phenotype_TFAM_6months'
column_of_interest <- which(colnames(read_pheno_matched) == phenotype_of_interest)
column_of_interest
phenotype_name <- as.character(colnames(read_pheno_matched)[column_of_interest])
phenotype_name

head(read_pheno_matched)

# Adjust colour palette to grouping according to phenotype of interest:
colours <- rownames(as.table(summary(as.factor(read_pheno_matched[, column_of_interest]))))
# View(as.table(summary(as.factor(read_pheno_matched[, column_of_interest]))))
str(colours)
#colours <- as.factor(colours)
# Order for time bins:
colours <- factor(colours)#, ordered = TRUE, labels = c("six_months", "one_year",  "six_years", "more_than_6_years"))
str(colours)
colours
colour_palette <- rainbow(length(colours), alpha = 0.5)
colour_palette
palette(colour_palette)

head(read_pheno_matched)
tail(read_pheno_matched)
#############################################


#############################
## Identify SNPs associated to phenotype of interest

# TO DO: plot all groups/colours of phenotypes by PCs.

# Overlay phenotype labels onto PCA and visualise.
# Set fac option in s.class() as factor for phenotype of interest:
fac <- as.factor(read_pheno_matched[, column_of_interest])

# Plot first axes of PCA analysis:
png(paste('labels_on_PCA_1_and_2_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, fac=fac, col=colour_palette, cpoint=2, sub="PCA - axes 1 and 2")
dev.off()

# Second pair of axes:
png(paste('labels_on_PCA_3_and_4_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=3, yax=4, fac=fac, col=colour_palette, cpoint=2, sub="PCA - axes 3 and 4")
dev.off()

# Run through axes that were kept as above.
png(paste('labels_on_PCA_5_and_6_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
s.class(run_pca$li, xax=5, yax=6, fac=fac, col=colour_palette, cpoint=2, sub="PCA - axes 5 and 6")
dev.off()

# TO DO: run plots for remaining axes. Run a plotting function to loop through any PCs kept...

# TO DO: run apply with chi2 through top SNPs (ie top 0.01%)?
# Run statistical test on this (based on the populations identified using PCA), 
# use chi-square:
table(read_pheno_matched[, 2], pop)
chisq.test(table(read_pheno_matched[, 2], pop), simulate=TRUE)

# Note: PCA optimizes the representation of the overall genetic 
# diversity and does not explicitly look for distinctions between 
# pre-defined groups of isolates.
# If only a few loci are correlated, they will be overlooked, 
# especially  if  stronger  structures  such  as  separate
# lineages or populations are present.  


## Use DAPC to look for combinations of SNPs correlated to a phenotype.
# in a given group of individuals, 
# ?dapc

# TO DO, run through all phenotypes:

# TO DO: Decide PCs and DAs to keep based on something objective, currently keeping maximum possible:
n.pca_keep <- length(read_pheno_matched[, 1]) - 1
percent_pca <- 40
n.da_keep <- length(colours) - 1
n.pca_keep
n.da_keep
percent_pca
colours

png(paste('dapc_plot_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
run_dapc <- assign(paste('dapc_', phenotype_name, sep = ''), 
                   dapc(data_to_adegent, read_pheno_matched[, column_of_interest], perc.pca = percent_pca, 
                        pca.select = 'percVar', n.da = n.da_keep))
dev.off()
run_dapc

# dapc_to_plot <- dapc_time_bins
# dapc_to_plot

# DAPC retains PCs for the dimension reduction step (PCA, retain 30 PCs)
# and for the subsequent discriminant analysis (DA). 
# For the latter, only one axis can be retained (the maximum number of 
# axes in DA is always the number of groups minus 1).

# Visualise with a scatter plot:
png(paste('scatter_run_dapc_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
scatter(run_dapc, bg="white", scree.da=TRUE, scree.pca=TRUE, posi.pca="topright", col=colour_palette, legend=TRUE, posi.leg="topleft")
dev.off()

# TO DO: run through all DAs kept to visualise separation:
# Visualise a single discriminant function as a density plot:
png(paste('scatter_one_da_run_dapc_', phenotype_name, '.png', sep=''), width = 10, height = 10, units = 'in', res = 300)
# par(mfrow=c(1, 2))
scatter(run_dapc, 1, 1, scree.da = F, legend=TRUE, posi.leg="topright")
dev.off()

# The  contribution  of  each  variable  to  the  separation  of  the  
# groups is stored in run_dapc$var.contr, visualise with:
head(run_dapc$var.contr)
tail(run_dapc$var.contr)
dim(run_dapc$var.contr)
range(run_dapc$var.contr)

# Run loadingplot and use the temp object returned to get the alleles with 
# the largest contribution to group separation

png(paste('loading_plot_top0.01pct_SNPs_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
run_loadingplot <- assign(paste('loadingplot_', phenotype_name, sep = ''), 
                          loadingplot(run_dapc$var.contr, main = 'Loading plot - top variants by contribution to group separation',
                                      threshold = quantile(run_dapc$var.contr, 0.999), lab.jitter=1))
dev.off()

run_loadingplot$threshold
head(run_loadingplot$var.names)
head(run_loadingplot$var.idx)
head(run_loadingplot$var.values)
length(run_loadingplot$var.names)
head(run_loadingplot$var.values)
str(run_loadingplot$var.values)
# qplot(y = run_loadingplot$var.values, x = run_loadingplot$var.idx, label = labels)

# Plot with ggplot loadings (contributions of each allele to group separation) by index (position) and label:
df <- as.data.frame(run_loadingplot)
y <- df$var.values
x <- df$var.idx
labels <- df$var.names

png(paste('loading_plot_top0.01pct_SNPs_ggplot_', phenotype_name, '.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
ggplot(data = df, aes(y = y, x = x)) + geom_point() + geom_text(label = labels) + ggtitle('Loading plot from DAPC') + ylab('Loadings') + xlab('Index')
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

