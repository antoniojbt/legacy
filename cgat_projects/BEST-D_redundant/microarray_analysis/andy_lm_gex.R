########################
# Compare joint 2000+4000 vs baseline and placebo with pairing.
# Define factors to constrast from annotation file:
#Pairing:
head(membership_file_cleaned)
tail(membership_file_cleaned)
# kit_id appears once with corresponding pt_id twice (baseline + final visit):
which(membership_file_cleaned == 104018)
membership_file_cleaned[c(275, 560), ]

pairing_joint <- factor(membership_file_cleaned$pt_id)
head(pairing_joint)
length(pairing_joint)
pairing_joint

treatment_joint <- factor(membership_file_cleaned$treatment, levels = c('untreated', 'treated_2000', 'treated_4000', 'treated_placebo'))
head(treatment_joint)
str(treatment_joint)
summary(treatment_joint)
treatment_joint

#Define design and set contrasts:
design_all_treated_pairs <- model.matrix(~pairing_joint+treatment_joint)
head(design_all_treated_pairs)[1:5, 1:5]
tail(design_all_treated_pairs)[, -1]
dim(design_all_treated_pairs)

head(design_all_treated_pairs[,299:301])
colSums(design_all_treated_pairs[,299:301])
quantile( rowSums(design_all_treated_pairs[,299:301]) )

newdesign <- cbind(
  design_all_treated_pairs[,-(299:301)],
  design_all_treated_pairs[,299] + 2*design_all_treated_pairs[,300], # 1* for  x, y, y , 2* for y, x+y, 2x+y
  # design_all_treated_pairs[,299],  # indicator of placebo-time; 299:301 to account for time (y)
  rowSums( design_all_treated_pairs[,299:301] ) # indicator of time > 0;  299:301 to account for time (y)
)
newdesign


head(newdesign[,299:300])
sum(newdesign[,300])
quantile( rowSums(newdesign[,299:300]) )
colnames( newdesign ) <- c( colnames(newdesign)[1:298], c( 'VitaminD', 'Time' ) )


colnames(newdesign)
#Run linear model and obtain differentially expressed genes based on all pairs:
fit_all_treated <- lmFit(normalised_filtered, newdesign)
fit_all_treated_2 <- eBayes(fit_all_treated)
dim(fit_all_treated_2)
str(fit_all_treated_2)
summary(fit_all_treated_2)
names(fit_all_treated_2)
head(fit_all_treated_2$coefficients)#[1:5, 1:5]
head(fit_all_treated_2$lods)#[1:5, 1:5]
head(fit_all_treated_2$cov.coefficients)#[1:5, 1:5]
fit_all_treated_2$coefficients
colnames(fit_all_treated_2)

# Get results:
topTable(fit_all_treated_2, adjust='BH')
topTable(fit_all_treated_2, coef="VitaminD", adjust='BH')
topTable(fit_all_treated_2, coef="Time", adjust='BH')
