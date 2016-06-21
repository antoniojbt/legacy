# Do the same for the covariates file (leave only those with data in the genotye and expression files):
length(which(colnames(covar_transposed2) %in% colnames(geno2_matched)))

covar_transposed2_to_match <- which(colnames(covar_transposed2) %in% colnames(geno2_matched))
covar_transposed2_to_match
covar_transposed2[1:5,1:10]
covar_transposed2_matched <- covar_transposed2[ , covar_transposed2_to_match]
covar_transposed2_matched[1:5,1:10]
dim(covar_transposed2_matched)

length(which((colnames(covar_transposed2_matched) %in% colnames(covar_transposed2))))
length(which(colnames(geno2_matched) %in% colnames(covar_transposed2_matched)))

# Check final files match:
dim(geno2_matched)
dim(expr2_matched)
dim(covar_transposed2_matched)

geno2_matched[1:5,1:10]
expr2_matched[1:5,1:10]
covar_transposed2_matched[1:5,1:10]

# Check colnames match exactly (label and order):
identical(colnames(geno2_matched), colnames(expr2_matched))
identical(colnames(expr2_matched), colnames(covar_transposed2_matched))
identical(colnames(geno2_matched), colnames(covar_transposed2_matched))
#############################################

minus_column <- subset_baseline_2000[,-1]
to_match <- which(colnames(minus_column) %in% colnames(geno2_matched))
matched <- minus_column[ , to_match]
ordered1 <- order(names(matched))
ordered2 <- matched[, ordered1]
plus_column <- cbind(subset_baseline_2000[,1], ordered2)
minus_column[1:5, 1:5]
to_match
matched[1:5, 1:5]
ordered2[1:5, 1:5]
plus_column[1:5, 1:5]
subset_baseline_2000  <- plus_column

minus_column <- subset_baseline_4000[,-1]
to_match <- which(colnames(minus_column) %in% colnames(geno2_matched))
matched <- minus_column[ , to_match]
ordered1 <- order(names(matched))
ordered2 <- matched[, ordered1]
plus_column <- cbind(subset_baseline_4000[,1], ordered2)
minus_column[1:5, 1:5]
to_match
matched[1:5, 1:5]
ordered2[1:5, 1:5]
plus_column[1:5, 1:5]
subset_baseline_4000  <- plus_column

# Final visits:
head(kit_ids_matched)
kit_ids_matched$FID

minus_column <- subset_finalVisit_placebo[,-1]
to_match <- which(colnames(minus_column) %in% kit_ids_matched$FID)
matched <- minus_column[ , to_match]
ordered1 <- order(names(matched))
ordered2 <- matched[, ordered1]
plus_column <- cbind(subset_finalVisit_placebo[,1], matched)
minus_column[1:5, 1:5]
to_match
matched[1:5, 1:5]
ordered2[1:5, 1:5]
plus_column[1:5, 1:5]
subset_finalVisit_placebo  <- plus_column

minus_column <- subset_finalVisit_2000[,-1]
to_match <- which(colnames(minus_column) %in% kit_ids_matched$FID)
matched <- minus_column[ , to_match]
ordered1 <- order(names(matched))
ordered2 <- matched[, ordered1]
plus_column <- cbind(subset_finalVisit_2000[,1], ordered2)
minus_column[1:5, 1:5]
to_match
matched[1:5, 1:5]
ordered2[1:5, 1:5]
plus_column[1:5, 1:5]
subset_finalVisit_2000 <- plus_column

minus_column <- subset_finalVisit_4000[,-1]
to_match <- which(colnames(minus_column) %in% kit_ids_matched$FID)
matched <- minus_column[ , to_match]
ordered1 <- order(names(matched))
ordered2 <- matched[, ordered1]
plus_column <- cbind(subset_finalVisit_4000[,1], ordered2)
minus_column[1:5, 1:5]
to_match
matched[1:5, 1:5]
ordered2[1:5, 1:5]
plus_column[1:5, 1:5]
subset_finalVisit_4000  <- plus_column
dim(subset_finalVisit_4000)
subset_finalVisit_4000[1:5, 1:5]

sapply(subset_list, dim)
