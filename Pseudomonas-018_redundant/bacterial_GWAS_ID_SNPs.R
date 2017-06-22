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
load('R_session_saved_image_dapc.RData', verbose=T)
# load('R_session_saved_image.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_ID_SNPs','.RData', sep='')
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


#############################

# TO DO: extract variants automatically:
# Save the top SNPs based on visual inspection of loadingplot:
# The population can be defined as the groupings, epi-cluster, etc.:
# ?genind2genpop
# data(H3N2)
# H3N2$other$epid # Here the populations are the years of epidemic

# variant_position <- c('V60107') # This was from morphology
head(data_to_adegent$all.names)
# data_to_adegent$loc.n.all[variant_position]
# head(data_to_adegent$tab[variant_position])

# variant_position <- c('V52421') # This was from time_bins
# variant_position <- c('V247549') # This was from time_bins
variant_position
variant_data <- data_to_adegent[loc=variant_position]
variant_data

# Get the allelic profiles for the most strongly correlated variants:
freq_allele1 <- tab(genind2genpop(variant_data, pop = read_pheno_matched[, column_of_interest]), freq=TRUE)
freq_allele1

pop_allele1 <- genind2genpop(variant_data, pop = read_pheno_matched[, column_of_interest])
pop_allele1
summary(pop_allele1)
# TO DO: clean up and extract without intervention:
# sel.profiles <- apply(data_to_adegent[, top_SNPs$var.idx]$tab, 1, paste, collapse="-")
# head(sel.profiles)
# table(sel.profiles)
# 
# sel.profiles_top_2 <- apply(data_to_adegent[, top_2_SNPs$var.idx]$tab, 1, paste, collapse="-")
# head(sel.profiles_top_2)
# table(sel.profiles_top_2)

sel.profiles_test <- apply(variant_data$tab, 1, paste, collapse="-")
head(sel.profiles_test)
table(sel.profiles_test)
df_sel.profiles_test <- as.data.frame(table(sel.profiles_test))
df_sel.profiles_test

# Bind profiles of alleles per sample to phenotype and sample name:
columns_chi2 <- c(1, column_of_interest)
columns_chi2
pheno_geno <- cbind.data.frame(read_pheno_matched[, columns_chi2], sel.profiles_test)
head(pheno_geno, 5)
dim(pheno_geno)

# Create a contingency table between  phenotype  and  SNP  profiles:
allele_by_pop_table <- table(read_pheno_matched[, column_of_interest], sel.profiles_test)
chisq.test(allele_by_pop_table, simulate=TRUE)

# TO DO: plot top alleles:
# Plot top SNPs:  
# png('loadingplot_top_SNPs.png', width = 4, height = 4, units = 'in', res = 300)
# loadingplot(dapc_to_plot$var.contr, threshold=2.0e-05)
# dev.off()
# 
# png('loadingplot_top_2_SNPs.png', width = 4, height = 4, units = 'in', res = 300)
# loadingplot(dapc_to_plot$var.contr, threshold=2.7e-05)
# dev.off()

# TO DO: check assign function to test whether dapc assigned to group of 
# interest correctly and plot with posterior prob. see p. 23.

#############################


#############################
# Plot time course:
# data(H3N2)
# pop(H3N2) <- H3N2$other$epid
# dapc.flu <- dapc(H3N2, n.pca=30,n.da=10)
# set.seed(4)
# contrib <- loadingplot(dapc.flu$var.contr, axis=2, thres=.07, lab.jitter=1)
# freq399 <- tab(genind2genpop(H3N2[loc=c("399")]),freq=TRUE)
# freq399
# matplot(freq399, pch=c("c","t"), type="b", xlab="year",
#         ylab="allele frequency", xaxt="n", cex=1.5,
#         main="SNP # 399")

#############################

# variant_data$all.names
# matplot(allele_by_pop_table, pch=c("a","c"), type="b", xlab="Time since first sample", ylab="Allele frequency", 
#         xaxt="n", cex=1.5, main="SNP # ")
# 
# matplot(freq_allele1, pch=c('C','T', '0'), type="b", xlab="Time since first sample", ylab="Allele frequency", 
#         xaxt="n", cex=1.5, main=paste('SNP # ', variant_position))
# axis(side=1, at=1:4, lab=colours)

# matplot(freq399, pch=c("c","t"), type="b", xlab="year",
#         ylab="allele frequency", xaxt="n", cex=1.5,
#         main="SNP # 399")
# axis(side=1, at=1:6, lab=2001:2006)

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
