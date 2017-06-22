##########################
# Antonio Berlanga-Taylor 
# BEST-D file processing
# 08 Feb 2016
# Add new variables Jon sent end of Jan 2016 (age, disease, smoking, etc.) 
##########################


##########################
#source('http://bioconductor.org/biocLite.R')
#biocLite('lubridate')
library(lubridate)
library(plyr)

# Existing phenotype data:
phenotype_data <- read.csv('BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups.txt', sep = '\t', header = TRUE, na.string = '-99999')
head(phenotype_data)
summary(phenotype_data)

# Variables to add:
new_variables <- read.csv('bestd_antonio.txt', sep = '\t', header = TRUE, na.string = c('', ' ', NA))
head(new_variables)
summary(new_variables)

# Date (randomised_date) needs to be recoded to set as season:
new_variables['randomised_date_2'] <- NA
new_variables['randomised_date_2'] <- new_variables['randomised_date']
summary(new_variables['randomised_date_2'])
summary(new_variables['randomised_date'])
which(is.na(new_variables['randomised_date']))

new_variables$randomised_date_2 <- as.POSIXlt(new_variables$randomised_date_2, tz = 'GMT', '%d/%m/%y')
summary(new_variables['randomised_date_2'])
head(new_variables['randomised_date_2'])
head(new_variables['randomised_date'])
str(new_variables['randomised_date_2'])

# Use lubridate to extract month:
new_variables['month_randomisation'] <- NA
new_variables['month_randomisation'] <- month(ymd(new_variables$randomised_date_2))
head(new_variables$month_randomisation)
head(new_variables$randomised_date_2)
tail(new_variables$month_randomisation)
tail(new_variables$randomised_date_2)
count(new_variables$month_randomisation)

# Recode to seasons:
new_variables['season_randomisation'] <- NA
new_variables['season_randomisation'][new_variables['month_randomisation'] == 1] <- 'winter'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 2] <- 'winter'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 3] <- 'spring'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 4] <- 'spring'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 5] <- 'spring'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 6] <- 'summer'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 7] <- 'summer'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 8] <- 'summer'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 9] <- 'autumn'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 10] <- 'autumn'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 11] <- 'autumn'
new_variables['season_randomisation'][new_variables['month_randomisation'] == 12] <- 'winter'
summary(new_variables['season_randomisation'])
head(new_variables['season_randomisation'])
count(new_variables['season_randomisation'])
count(new_variables$month_randomisation)

# Plink errors with character variables, recode season:
new_variables['season_randomisation_2'] <- NA
new_variables['season_randomisation_2'][new_variables['season_randomisation'] == 'winter'] <- 4
new_variables['season_randomisation_2'][new_variables['season_randomisation'] == 'spring'] <- 1
new_variables['season_randomisation_2'][new_variables['season_randomisation'] == 'summer'] <- 2
new_variables['season_randomisation_2'][new_variables['season_randomisation'] == 'autumn'] <- 3
summary(new_variables['season_randomisation_2'])
head(new_variables['season_randomisation'])
count(new_variables['season_randomisation_2'])
count(new_variables$month_randomisation)


# Merge datasets:
# Rename column for merge:
names(new_variables)[names(new_variables) == 'pt_id'] <- 'IID_pt_id'
full_phenotype_data <- merge(phenotype_data, new_variables, by = 'IID_pt_id', all.x = TRUE)
# Re-order first three columns for plink:
full_phenotype_data <- full_phenotype_data[, c(2, 3, 1, 4:ncol(full_phenotype_data))]
head(full_phenotype_data)
tail(full_phenotype_data)
#View(full_phenotype_data)
summary(full_phenotype_data)

dim(phenotype_data)
dim(new_variables)
dim(full_phenotype_data)


# Check variables to pass to plink:
plink_adjust_full_phenotype_data <- full_phenotype_data[, c("IID", "FID", "IID_pt_id", "incident_fall", "incident_fracture", "incident_resp_infection", "bmi0", "diabetes", "hypertension", "heart_disease", 
                                                            "stroke", "fracture", "falls", "copd_asthma", "muscle_pain0", "joint_pain0", "physical_activity12", 
                                                            "muscle_pain12", "joint_pain12", "physical_activity0", "basemed_antihypertensive", 
                                                            "basemed_statin", "basemed_antithrombotic", "basemed_vitamind", 
                                                            "basemed_calcium", "calendar_age_ra", "currsmoker", "exsmoker", "BMI_DEXA12", "season_randomisation_2")]

head(plink_adjust_full_phenotype_data)
summary(plink_adjust_full_phenotype_data)


# Write results to file:
write.table(full_phenotype_data, 'BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt', quote = FALSE, sep = '\t',
            na = '-99999', row.names = FALSE, col.names = TRUE)

q()

