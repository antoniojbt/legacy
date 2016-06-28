##########################
# Antonio Berlanga-Taylor 
# BEST-D file processing
# 20 Jan 2016
#Recode arm into new variables for sub-group analysis (because plink2's coding takes 0 or -9 = missing, 1 = unaffected, 2 = affected). 
# I've changed this below in the same way as here for the fam file
# Here I'll add new columns to create further labels of cases and controls within the phenotype file for sub-group tests 
#(ie placebo vs 2000+4000, 4000 vs placebo, etc.).
# Recode into new variables from original 'arm':
#0=4000 IU, 1=2000 IU, 2=Placebo
#to new
# -99999 = missing ; 1 = untreated or 2000UI ; 2 = 4000 or both treatments groups

phenotype_data <- read.csv('BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt', sep = ' ', header = TRUE)
head(phenotype_data)
summary(phenotype_data)

# Change placebo to 1 and others (2000 and 4000) to 2:
phenotype_data['placebo_v_Tx'] = NA
phenotype_data$placebo_v_Tx[phenotype_data$arm == 2] = 1
phenotype_data$placebo_v_Tx[phenotype_data$arm == 0] = 2
phenotype_data$placebo_v_Tx[phenotype_data$arm == 1] = 2
head(phenotype_data$placebo_v_Tx)
summary(phenotype_data$placebo_v_Tx)


phenotype_data['placebo_v_4000'] = NA
phenotype_data$placebo_v_4000[phenotype_data$arm == 2] = 1 # Change placebo to 1.
phenotype_data$placebo_v_4000[phenotype_data$arm == 0] = 2 # Change 4000 to 2.
phenotype_data$placebo_v_4000[phenotype_data$arm == 1] = -99999 # Change 2000 to missing.
head(phenotype_data$placebo_v_4000)

phenotype_data['placebo_v_2000'] = NA
phenotype_data$placebo_v_2000[phenotype_data$arm == 2] = 1 # Change placebo to 1.
phenotype_data$placebo_v_2000[phenotype_data$arm == 1] = 2 # Change 2000 to 2.
phenotype_data$placebo_v_2000[phenotype_data$arm == 0] = -99999 # Change 4000 to missing.
head(phenotype_data$placebo_v_2000)

phenotype_data['IU4000_v_2000'] = NA
phenotype_data$IU4000_v_2000[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
phenotype_data$IU4000_v_2000[phenotype_data$arm == 1] = 1 # Keep 2000 as 1.
phenotype_data$IU4000_v_2000[phenotype_data$arm == 0] = 2 # Change 4000 to 2.
head(phenotype_data$IU4000_v_2000)

head(phenotype_data[,c('arm', 'placebo_v_Tx', 'placebo_v_4000', 'placebo_v_2000',
                       'IU4000_v_2000')])

tail(phenotype_data[,c('arm', 'placebo_v_Tx', 'placebo_v_4000', 'placebo_v_2000',
                       'IU4000_v_2000')])

#View(phenotype_data)


# Add alternate phenotype codings for vitamin D results so that only placebo is tested, only 2000, only 4000 and only treated in vitD gen. assoc. analysis:
# Recode into new variables from original 'arm':
#0=4000 IU, 1=2000 IU, 2=Placebo
#to new variables where only sub-groups are left for vitD measurements
# -99999 = missing

# For vitd6:
phenotype_data['vitd6_placebo_only'] = phenotype_data$vitd6
phenotype_data$vitd6_placebo_only[phenotype_data$arm == 0] = -99999 # Change 4000 to missing.
phenotype_data$vitd6_placebo_only[phenotype_data$arm == 1] = -99999 # Change 2000 to missing.
head(phenotype_data$vitd6_placebo_only)

phenotype_data['vitd6_2000_only'] = phenotype_data$vitd6
phenotype_data$vitd6_2000_only[phenotype_data$arm == 0] = -99999 # Change 4000 to missing.
phenotype_data$vitd6_2000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
head(phenotype_data$vitd6_2000_only)

phenotype_data['vitd6_4000_only'] = phenotype_data$vitd6
phenotype_data$vitd6_4000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
phenotype_data$vitd6_4000_only[phenotype_data$arm == 1] = -99999 # Change 2000 to missing.
head(phenotype_data$vitd6_4000_only)

phenotype_data['vitd6_2000_4000_only'] = phenotype_data$vitd6
phenotype_data$vitd6_2000_4000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
head(phenotype_data$vitd6_2000_4000_only)

head(phenotype_data[,c('arm', 'vitd6', 'vitd6_placebo_only', 'vitd6_2000_only', 'vitd6_4000_only', 'vitd6_2000_4000_only')])

tail(phenotype_data[,c('arm', 'vitd6', 'vitd6_placebo_only', 'vitd6_2000_only', 'vitd6_4000_only', 'vitd6_2000_4000_only')])


# For vitd12:
phenotype_data['vitd12_placebo_only'] = phenotype_data$vitd12
phenotype_data$vitd12_placebo_only[phenotype_data$arm == 0] = -99999 # Change 4000 to missing.
phenotype_data$vitd12_placebo_only[phenotype_data$arm == 1] = -99999 # Change 2000 to missing.
head(phenotype_data$vitd12_placebo_only)

phenotype_data['vitd12_2000_only'] = phenotype_data$vitd12
phenotype_data$vitd12_2000_only[phenotype_data$arm == 0] = -99999 # Change 4000 to missing.
phenotype_data$vitd12_2000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
head(phenotype_data$vitd12_2000_only)

phenotype_data['vitd12_4000_only'] = phenotype_data$vitd12
phenotype_data$vitd12_4000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
phenotype_data$vitd12_4000_only[phenotype_data$arm == 1] = -99999 # Change 2000 to missing.
head(phenotype_data$vitd12_4000_only)

phenotype_data['vitd12_2000_4000_only'] = phenotype_data$vitd12
phenotype_data$vitd12_2000_4000_only[phenotype_data$arm == 2] = -99999 # Change placebo to missing.
head(phenotype_data$vitd12_2000_4000_only)

head(phenotype_data[,c('arm', 'vitd12', 'vitd12_placebo_only', 'vitd12_2000_only', 'vitd12_4000_only', 'vitd12_2000_4000_only')])

tail(phenotype_data[,c('arm', 'vitd12', 'vitd12_placebo_only', 'vitd12_2000_only', 'vitd12_4000_only', 'vitd12_2000_4000_only')])


# Write results to file:
write.table(phenotype_data, 'BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups.txt', quote = FALSE, sep = '\t',
            na = '-99999', row.names = FALSE, col.names = TRUE)

q()
