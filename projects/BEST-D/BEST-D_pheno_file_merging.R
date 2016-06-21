#############################
# Antonio J Berlanga-Taylor
# 24 Feb 2016
# BEST-D project phenotype/metadata file processing/merging
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file("R_session_output_pheno_file_merging.txt", open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific.

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_pheno_file_merging', '.RData', sep='')

#############################

#############################
## Update packages if necessary and load them:

library(lubridate)
library(plyr)
#############################

#############################
## Read in files:

pheno_first_set<- read.csv('JULIAN_sheet_1.txt', sep = '\t', header = TRUE, na.string = c('-99999', "", " ", "NA", "NaN"))
pheno_new_variables <- read.csv('bestd_antonio.txt', sep = '\t', header = TRUE, na.string = c('-99999', "", " ", "NA", "NaN"))
lab_kit_IDs <- read.csv('03_form_lab_kits.csv.tr', header = TRUE, na.string = c('-99999', "", " ", "NA", "NaN"))

head(pheno_first_set)
head(pheno_new_variables)
head(lab_kit_IDs)
pheno_first_set[1:5, 1:5]
pheno_new_variables[1:5, 1:5]
dim(pheno_first_set)
dim(pheno_new_variables)
dim(lab_kit_IDs)

## Check missing/mis-labeled samples in kit ids file:
# pt_id 105281: Kit ID 120005281 should be 120005145, see email from Jon 23/03/2015:
lab_kit_IDs[which(lab_kit_IDs$pt_id == '105281'), ]
lab_kit_IDs[which(lab_kit_IDs$kit_id == '120005281'), ] # wrong label
lab_kit_IDs[which(lab_kit_IDs$kit_id == '120005145'), ] # correct label
lab_kit_IDs[640, 4] <- 120005145
lab_kit_IDs[640, 4]

# pt_id 103485: Kit ID 120000270 should be 120000272
lab_kit_IDs[which(lab_kit_IDs$pt_id == '103485'), ]
lab_kit_IDs[which(lab_kit_IDs$kit_id == '120000270'), ] # wrong label
lab_kit_IDs[which(lab_kit_IDs$kit_id == '120000272'), ] # correct label
lab_kit_IDs[447, 4] <- 120000272
lab_kit_IDs[447, 4]

# Merge phenotype files:
length(which(row.names(pheno_first_set) %in% row.names(pheno_new_variables)) == TRUE)

pheno_all <- merge(pheno_first_set, pheno_new_variables)
dim(pheno_all)
names(pheno_all)
which(names(pheno_first_set) %in% names(pheno_new_variables))
names(pheno_first_set)[c(1, 31, 45)]

head(pheno_all)
tail(pheno_all)
summary(pheno_all)
class(pheno_all)

# Cross with lab kit IDs for genotype and gene expression array data matching (which is based on kit IDs):
head(lab_kit_IDs)
dim(lab_kit_IDs)
which(names(pheno_all) == 'arm')
which(names(pheno_all) == 'visit_type')

head(t(lab_kit_IDs))
head(t(pheno_all))

# Subset baseline and final samples to make merging easier:
lab_kit_IDs_Randomisation <- subset(lab_kit_IDs, subset = (lab_kit_IDs$visit_type == 'Randomisation'))
head(lab_kit_IDs_Randomisation)
dim(lab_kit_IDs_Randomisation)

lab_kit_IDs_FinalVisit <- subset(lab_kit_IDs, subset = (lab_kit_IDs$visit_type == 'Final Visit'))
head(lab_kit_IDs_FinalVisit)
dim(lab_kit_IDs_FinalVisit)

lab_kit_IDs_first_and_last <- merge(lab_kit_IDs_Randomisation, lab_kit_IDs_FinalVisit, by = 'pt_id', all.x = TRUE)
dim(lab_kit_IDs_first_and_last)
head(lab_kit_IDs_first_and_last)
summary(lab_kit_IDs_first_and_last) 
names(lab_kit_IDs_first_and_last)

# Clean up names and unnecessary columns:
lab_kit_IDs_first_and_last <- lab_kit_IDs_first_and_last[, -5]
colnames(lab_kit_IDs_first_and_last)[2] <- 'arm'
colnames(lab_kit_IDs_first_and_last)
colnames(lab_kit_IDs_first_and_last)[4] <- 'kit_id_randomisation'
colnames(lab_kit_IDs_first_and_last)[6] <- 'kit_id_finalVisit'
lab_kit_IDs_first_and_last <- lab_kit_IDs_first_and_last[, -c(3, 5)]
colnames(lab_kit_IDs_first_and_last)
head(lab_kit_IDs_first_and_last)
dim(lab_kit_IDs_first_and_last)

# Merge with phenotype file:
which(names(pheno_all) %in% names(lab_kit_IDs_first_and_last))
names(pheno_all)[c(1, 8)]

pheno_all_kit_ids <- merge(pheno_all, lab_kit_IDs_first_and_last)
dim(pheno_all)
dim(lab_kit_IDs_first_and_last)
dim(pheno_all_kit_ids)
head(pheno_all_kit_ids)
names(pheno_all_kit_ids)
summary(pheno_all_kit_ids)

# Sanity check for corrected labels (see at start of script):
which(pheno_all_kit_ids$kit_id_finalVisit == 120000272)
which(pheno_all_kit_ids$kit_id_randomisation == 120000272)

# Also subset one file with baseline and final kit id to write to file for array QC steps:
lab_kits_for_array_QC <- subset(lab_kit_IDs, subset = (lab_kit_IDs$visit_type == 'Randomisation' | 
                                                         lab_kit_IDs$visit_type == 'Final Visit'))
dim(lab_kits_for_array_QC)
head(lab_kits_for_array_QC)
tail(lab_kits_for_array_QC)
which(lab_kits_for_array_QC$kit_id == 120000272)

# Date (randomised_date) needs to be recoded to set as season:
pheno_all_kit_ids['randomised_date_2'] <- NA
pheno_all_kit_ids['randomised_date_2'] <- pheno_all_kit_ids['randomised_date']
summary(pheno_all_kit_ids['randomised_date_2'])
summary(pheno_all_kit_ids['randomised_date'])
which(is.na(pheno_all_kit_ids['randomised_date']))

pheno_all_kit_ids$randomised_date_2 <- as.POSIXlt(pheno_all_kit_ids$randomised_date_2, tz = 'GMT', '%d/%m/%y')
summary(pheno_all_kit_ids['randomised_date_2'])
head(pheno_all_kit_ids['randomised_date_2'])
head(pheno_all_kit_ids['randomised_date'])
str(pheno_all_kit_ids['randomised_date_2'])

# Use lubridate to extract month:
pheno_all_kit_ids['month_randomisation'] <- NA
pheno_all_kit_ids['month_randomisation'] <- month(ymd(pheno_all_kit_ids$randomised_date_2))
head(pheno_all_kit_ids$month_randomisation)
head(pheno_all_kit_ids$randomised_date_2)
tail(pheno_all_kit_ids$month_randomisation)
tail(pheno_all_kit_ids$randomised_date_2)
count(pheno_all_kit_ids$month_randomisation)

# Recode to seasons:
pheno_all_kit_ids['season_randomisation'] <- NA
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 1] <- 'winter'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 2] <- 'winter'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 3] <- 'spring'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 4] <- 'spring'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 5] <- 'spring'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 6] <- 'summer'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 7] <- 'summer'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 8] <- 'summer'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 9] <- 'autumn'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 10] <- 'autumn'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 11] <- 'autumn'
pheno_all_kit_ids['season_randomisation'][pheno_all_kit_ids['month_randomisation'] == 12] <- 'winter'
summary(pheno_all_kit_ids['season_randomisation'])
head(pheno_all_kit_ids['season_randomisation'])
count(pheno_all_kit_ids['season_randomisation'])
count(pheno_all_kit_ids$month_randomisation)

# Plink errors with character variables, recode season:
pheno_all_kit_ids['season_randomisation_2'] <- NA
pheno_all_kit_ids['season_randomisation_2'][pheno_all_kit_ids['season_randomisation'] == 'winter'] <- 4
pheno_all_kit_ids['season_randomisation_2'][pheno_all_kit_ids['season_randomisation'] == 'spring'] <- 1
pheno_all_kit_ids['season_randomisation_2'][pheno_all_kit_ids['season_randomisation'] == 'summer'] <- 2
pheno_all_kit_ids['season_randomisation_2'][pheno_all_kit_ids['season_randomisation'] == 'autumn'] <- 3
summary(pheno_all_kit_ids['season_randomisation_2'])
head(pheno_all_kit_ids['season_randomisation'])
count(pheno_all_kit_ids['season_randomisation_2'])
#############################

#############################
# Write files to disk:
write.table(pheno_all_kit_ids, 'BEST-D_phenotype_file_final.tsv', sep='\t', quote=F, 
            na='NA', col.names = TRUE, row.names = FALSE)

write.table(lab_kits_for_array_QC, 'BEST-D_lab_kit_IDs_for_array_QC.tsv', sep='\t', quote=F, 
            na='NA', col.names = TRUE, row.names = FALSE)

#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for differential expression analysis.
#############################

