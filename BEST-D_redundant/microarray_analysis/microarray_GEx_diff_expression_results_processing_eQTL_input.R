#####################################
# TO DO: Process files containing logFC (lmfit object) so that baseline kit_id is the 
# identifier and can be used to match downstream for eQTL.
# TO DO: Use log odds from eBayes object or coefficient. 

# Use pairing output from limma:
head(topTable_pairing_high_dose)
head(topTable_pairing_low_dose)

dim(fit_all_pairs_2)
class(fit_all_pairs_2)
names(fit_all_pairs_2)
str(fit_all_pairs_2)
summary(fit_all_pairs_2)
head(fit_all_pairs_2)
#View(fit_all_pairs)

# Convert eBayes object with lods to dataframe for lods:
df_fit_all_pairs_2_lods <- as.data.frame(fit_all_pairs_2$lods[, -1 ]) # Remove intercept column
head(df_fit_all_pairs_2_lods)[1:5, 1:5]
tail(df_fit_all_pairs_2_lods)[1:5, 1:5]
dim(df_fit_all_pairs_2_lods)


# Split string (delete 'pairing'):
split_colnames <- as.data.frame(strsplit(colnames(df_fit_all_pairs_2_lods), split = 'g', fixed = TRUE))
head(split_colnames)[, 1:5]
dim(split_colnames)
split_colnames[2, 1:5]
split_colnames <- t(split_colnames[2, ])
#View(t(split_colnames[2, ]))
head(split_colnames)

# Rename column headers with split strings:
names(df_fit_all_pairs_2_lods) <- split_colnames
#View(df_fit_all_pairs_2_lods)
dim(df_fit_all_pairs_2_lods)
df_fit_all_pairs_2_lods <- df_fit_all_pairs_2_lods[, 1:(ncol(df_fit_all_pairs_2_lods)-2)] #Remove last two columns
head(df_fit_all_pairs_2_lods)

# Tranpose so that it can be merged:
df_fit_all_pairs_2_lods_t <- as.data.frame(t(df_fit_all_pairs_2_lods))
head(df_fit_all_pairs_2_lods_t)[1:5, 1:5]
#View(df_fit_all_pairs_2_lods_t)
class(df_fit_all_pairs_2_lods_t)
row.names(df_fit_all_pairs_2_lods_t)
colnames(df_fit_all_pairs_2_lods_t)
dim(df_fit_all_pairs_2_lods_t)
df_fit_all_pairs_2_lods_t['pt_id'] <- NA
df_fit_all_pairs_2_lods_t['pt_id'] <- row.names(df_fit_all_pairs_2_lods_t)
dim(df_fit_all_pairs_2_lods_t)
head(df_fit_all_pairs_2_lods_t)#[1:5, 1:5]
tail(df_fit_all_pairs_2_lods_t)#[1:5, 1:5]

# Subset membership file so that it only contains baseline kit_id:
#120005004 101242 finalVisit
head(membership_file_cleaned)
tail(membership_file_cleaned)
membership_file_cleaned_finalVisit_only <- membership_file_cleaned[which(membership_file_cleaned$visit_type == 'FinalVisit'), ]
membership_file_cleaned_finalVisit_only['kit_id'] <- row.names(membership_file_cleaned_finalVisit_only)
head(membership_file_cleaned_finalVisit_only)
tail(membership_file_cleaned_finalVisit_only)
dim(membership_file_cleaned_finalVisit_only)
dim(membership_file_cleaned)

# Merge files: 
df_fit_all_pairs_2_lods_t_finalVisit <- merge(df_fit_all_pairs_2_lods_t, membership_file_cleaned_finalVisit_only, by = 'pt_id')
head(df_fit_all_pairs_2_lods_t_finalVisit)[1:5, 1:5]
dim(df_fit_all_pairs_2_lods_t_finalVisit)
# Keep kit_id (corresponding to final visit) as identifier:
df_fit_all_pairs_2_lods_t_finalVisit[1:5, c(1, 16728:16731)] # columns to exclude (pt_id, arm, treatment, visit_type)
df_fit_all_pairs_2_lods_t_finalVisit <- df_fit_all_pairs_2_lods_t_finalVisit[, -c(1, 16728:16730)]
dim(df_fit_all_pairs_2_lods_t_finalVisit)
row.names(df_fit_all_pairs_2_lods_t_finalVisit) <- df_fit_all_pairs_2_lods_t_finalVisit$kit_id
df_fit_all_pairs_2_lods_t_finalVisit <- df_fit_all_pairs_2_lods_t_finalVisit[, -ncol(df_fit_all_pairs_2_lods_t_finalVisit)]
dim(df_fit_all_pairs_2_lods_t_finalVisit)
head(df_fit_all_pairs_2_lods_t_finalVisit)#[1:5, 1:5]

# Transpose FC/lods/coefficients file: 
df_fit_all_pairs_2_lods_kit_id_finalVisit <- as.data.frame(t(df_fit_all_pairs_2_lods_t_finalVisit))
head(df_fit_all_pairs_2_lods_kit_id_finalVisit)
dim(df_fit_all_pairs_2_lods_kit_id_finalVisit)
#View(df_fit_all_pairs_2_lods_kit_id_finalVisit)

# Keep only patient pairs from treated groups:
membership_file_cleaned_finalVisit_treated <- membership_file_cleaned_finalVisit_only[-which(membership_file_cleaned_finalVisit_only$arm == 2), ]
head(membership_file_cleaned_finalVisit_treated)
tail(membership_file_cleaned_finalVisit_treated)
dim(membership_file_cleaned_finalVisit_treated)
dim(membership_file_cleaned)
# Merge files for treated only samples: 
df_fit_all_pairs_2_lods_t_finalVisit_treated <- merge(df_fit_all_pairs_2_lods_t, membership_file_cleaned_finalVisit_treated, by = 'pt_id')
head(df_fit_all_pairs_2_lods_t_finalVisit_treated)[1:5, 1:5]
dim(df_fit_all_pairs_2_lods_t_finalVisit_treated)
# Keep kit_id (corresponding to final visit) as identifier for treated only samples: 
df_fit_all_pairs_2_lods_t_finalVisit_treated[1:5, c(1, 16728:16731)] # columns to exclude (pt_id, arm, treatment, visit_type)
df_fit_all_pairs_2_lods_t_finalVisit_treated <- df_fit_all_pairs_2_lods_t_finalVisit_treated[, -c(1, 16728:16730)]
dim(df_fit_all_pairs_2_lods_t_finalVisit_treated)
row.names(df_fit_all_pairs_2_lods_t_finalVisit_treated) <- df_fit_all_pairs_2_lods_t_finalVisit_treated$kit_id
df_fit_all_pairs_2_lods_t_finalVisit_treated <- df_fit_all_pairs_2_lods_t_finalVisit_treated[, -ncol(df_fit_all_pairs_2_lods_t_finalVisit_treated)]
dim(df_fit_all_pairs_2_lods_t_finalVisit_treated)
head(df_fit_all_pairs_2_lods_t_finalVisit_treated)#[1:5, 1:5]
# Transpose FC/lods/coefficients file for treated only samples: 
df_fit_all_pairs_2_lods_kit_id_finalVisit_treated <- as.data.frame(t(df_fit_all_pairs_2_lods_t_finalVisit_treated))
head(df_fit_all_pairs_2_lods_kit_id_finalVisit_treated)
dim(df_fit_all_pairs_2_lods_kit_id_finalVisit_treated)
#View(df_fit_all_pairs_2_lods_kit_id_finalVisit_treated)


# Write to disk:
write.table(df_fit_all_pairs_2_lods_kit_id_finalVisit, 'lods_paired_tests_all.tab', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

write.table(df_fit_all_pairs_2_lods_kit_id_finalVisit_treated, 'lods_paired_tests_treated_only.tab', sep='\t', 
            quote=F, na='NA', col.names=NA, row.names=TRUE)

#####################################
