#######################################
# Exploratory analysis of phenotypes P. aeruginosa in CF project (018)
# Antonio J Berlanga-Taylor
# 21 Sept 2015
#######################################


#######################################
# install.packages("Kendall")
library(ggplot2)
library(Kendall)
library(Hmisc)
library(boot)
library(gridExtra)
library(lubridate)

source('/ifs/devel/antoniob/projects/BEST-D/moveme.R')
# source('moveme.R')
# setwd('/ifs/projects/proj018/runs_second_round.dir/phenotype_analysis.dir/')
getwd()
options(digits = 10)
#######################################


#######################################
# Phenotypes to test:
# Increased antibiotic resistance as time goes by?
# Change in morphotypes according to time?
# Are morphotypes associated to antibiotic resistance?
# 

#######################################


#######################################
# Get antibiotics data and add EUCAST cut-offs:
antibiotics_data <- read.csv('antibiotic_resistance_final_13Aug2015.csv.tr', header = TRUE, stringsAsFactors = F)

# Antibiotic resistance keys:
# TM: Tobramycin
# IP: Imipenem
# AT: Aztreonam
# CI: Ciprofloxacin
# CO: Colistin
# http://www.eucast.org/fileadmin/src/media/PDFs/EUCAST_files/Breakpoint_tables/v_5.0_Breakpoint_Table_01.pdf

# View(antibiotics_data)
head(antibiotics_data)
antibiotics_data$TM_EUCAST <- as.factor(ifelse(antibiotics_data$TM < 4, 1, 2))#'Sensitive', 'Resistant'))
antibiotics_data$IP_EUCAST <- as.factor(ifelse(antibiotics_data$IP < 4, 1,  ifelse(antibiotics_data$IP > 8, 2, 3)))
                                               #'Sensitive', ifelse(antibiotics_data$IP > 8, 'Resistant', 'Intermediate')))
antibiotics_data$AT_EUCAST <- as.factor(ifelse(antibiotics_data$AT < 1, 1, ifelse(antibiotics_data$AT > 16, 2, 3)))
                                               #'Sensitive', ifelse(antibiotics_data$AT > 16, 'Resistant', 'Intermediate')))
antibiotics_data$CI_EUCAST <- as.factor(ifelse(antibiotics_data$CI < 0.5, 1, ifelse(antibiotics_data$CI > 1, 2, 3)))
                                               #'Sensitive', ifelse(antibiotics_data$CI > 1, 'Resistant', 'Intermediate')))
antibiotics_data$CO_EUCAST <- as.factor(ifelse(antibiotics_data$CO < 4, 1, 2))#'Sensitive', 'Resistant'))

class(antibiotics_data)
str(antibiotics_data)
head(antibiotics_data)
tail(antibiotics_data)
# View(antibiotics_data)
summary(antibiotics_data)
#######################################


#######################################
# Explore correlations between antibiotic classes:

cor.test(as.matrix(antibiotics_data[, 2]), as.matrix(antibiotics_data[, 3]), method = 'pearson', use = "pairwise.complete.obs")
cor.test(as.matrix(antibiotics_data[, 2]), as.matrix(antibiotics_data[, 3]), method = 'spearman', use = "pairwise.complete.obs", exact = FALSE)

cor.test(as.matrix(antibiotics_data[, 2]), as.matrix(antibiotics_data[, 4]), method = 'pearson', use = "pairwise.complete.obs")
cor.test(as.matrix(antibiotics_data[, 2]), as.matrix(antibiotics_data[, 5]), method = 'pearson', use = "pairwise.complete.obs")
cor.test(as.matrix(antibiotics_data[, 2]), as.matrix(antibiotics_data[, 6]), method = 'pearson', use = "pairwise.complete.obs")

# Subset data to run correlations all v. all:
antibiotics_data_subset <- antibiotics_data[, 2:6]
str(antibiotics_data_subset)
head(antibiotics_data_subset)
cor(antibiotics_data_subset[, 2], antibiotics_data_subset[, 3], method = 'pearson', use = "pairwise.complete.obs")

# Run all correlations of antimicrobials:
amb_pearson <- rcorr(as.matrix(antibiotics_data_subset), type = 'pearson')
amb_spearman <- rcorr(as.matrix(antibiotics_data_subset), type = 'spearman')
output_file <- file(paste('amb_correlation_results_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))
print('amb_spearman')
amb_spearman
print('amb_pearson')
amb_pearson
sink(file = NULL)

# Plots of individual correlations:
# qplot(antibiotics_data_subset[, 1], antibiotics_data_subset[, 2])
# qplot(antibiotics_data_subset[, 1], antibiotics_data_subset[, 3])
# qplot(antibiotics_data_subset[, 1], antibiotics_data_subset[, 4])
# qplot(antibiotics_data_subset[, 1], antibiotics_data_subset[, 5])

# Plot correlations of all antimicrobials:
png('antimicrobials_scatterplots.png')
plot(antibiotics_data_subset)
dev.off()

# Define multi-reistant phenotype (any sample resistant to two or more antibiotics): 
antibiotics_data$resistance_count <- rowSums(x = antibiotics_data[, 7:11] == 2, na.rm = FALSE)
head(antibiotics_data)
tail(antibiotics_data)
summary(antibiotics_data[, 7:12])
length(which(antibiotics_data[, 12] == 0))
length(which(antibiotics_data[, 12] == 1))
length(which(antibiotics_data[, 12] > 1))

antibiotics_data$multi_resistance <- as.factor(ifelse(antibiotics_data$resistance_count > 1, 2, 1))#'Multi-resistant', 'Sensitive'))
antibiotics_data$resistance <- as.factor(ifelse(antibiotics_data$resistance_count >= 1, 2, 1))#'Resistant', 'Sensitive'))
head(antibiotics_data)
tail(antibiotics_data)
summary(antibiotics_data)

#######################################


#######################################
# Determine if there is a trend towards resistance according to time:
# Merge with metadata file to get all other phenotypes and metadata:
metadata_file <- read.csv('metadata_strain_ids_pseudomonas_21Sept2015.csv.tr', header = TRUE, stringsAsFactors = F)
head(metadata_file)
str(metadata_file)

metadata_antibiotics <- merge(metadata_file, antibiotics_data, by = 'unique_ID_fastq')

head(metadata_antibiotics)
# View(metadata_antibiotics)
summary(metadata_antibiotics)


# Bin samples by quartiles from length of infection:
# boxplot(metadata_antibiotics$time_of_inf_years)
summary(metadata_antibiotics$time_of_inf_years)
length(which(metadata_antibiotics$time_of_inf_years < 2.1))
length(which(metadata_antibiotics$time_of_inf_years >= 2.1  & metadata_antibiotics$time_of_inf_years < 5.47))
length(which(metadata_antibiotics$time_of_inf_years >= 5.47 & metadata_antibiotics$time_of_inf_years < 11.75))
length(which(metadata_antibiotics$time_of_inf_years >= 11.75))

metadata_antibiotics$time_bins <- as.factor(ifelse(metadata_antibiotics$time_of_inf_years < 2.1, '1', 
                                         ifelse(metadata_antibiotics$time_of_inf_years >= 2.1 & metadata_antibiotics$time_of_inf_years < 5.47, '2',
                                         ifelse(metadata_antibiotics$time_of_inf_years >= 5.47 & metadata_antibiotics$time_of_inf_years < 11.75, '3',
                                         ifelse(metadata_antibiotics$time_of_inf_years >= 11.75, '4', 'NA')))))

summary(metadata_antibiotics$time_bins)

# Round to nearest year since infection:
metadata_antibiotics$time_discrete <- round(metadata_antibiotics$time_of_inf_years, 0)
summary(metadata_antibiotics$time_discrete)
#######################################


#######################################
# Summarise data of resistance according to time bins:
# View(metadata_antibiotics)

output_file <- file(paste('tables_018_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))

print('summary(metadata_antibiotics$time_of_inf_years)')
summary(metadata_antibiotics$time_of_inf_years)

print('
table(metadata_antibiotics$time_bins, metadata_antibiotics$multi_resistance)
table(metadata_antibiotics$time_bins, metadata_antibiotics$resistance_count)
table(metadata_antibiotics$time_discrete, metadata_antibiotics$multi_resistance)
table(metadata_antibiotics$time_discrete, metadata_antibiotics$resistance_count)
table(metadata_antibiotics$replicate_number_by_date, metadata_antibiotics$multi_resistance)
')

table(metadata_antibiotics$time_bins, metadata_antibiotics$multi_resistance)
table(metadata_antibiotics$time_bins, metadata_antibiotics$resistance_count)
table(metadata_antibiotics$time_discrete, metadata_antibiotics$multi_resistance)
table(metadata_antibiotics$time_discrete, metadata_antibiotics$resistance_count)
table(metadata_antibiotics$replicate_number_by_date, metadata_antibiotics$multi_resistance)

sink(file = NULL)

which(colnames(metadata_antibiotics) == 'replicate_number_by_date')
which(colnames(metadata_antibiotics) == 'multi_resistance')
which(colnames(metadata_antibiotics) == 'unique_ID_fastq')
which(colnames(metadata_antibiotics) == 'Patient_ID')
which(colnames(metadata_antibiotics) == 'morphology')
which(colnames(metadata_antibiotics) == 'time_discrete')
which(colnames(metadata_antibiotics) == 'original_date')
which(colnames(metadata_antibiotics) == 'time_of_inf_years')
which(colnames(metadata_antibiotics) == 'RAPD_type_short')

# TO DO: set date variable for R:
# dates <- as.Date(metadata_antibiotics$original_date, "%d/%B/%Y")
# dates
# 
# head(metadata_antibiotics[, c(1, 3, 5, 19, 37, 38, 39, 40)])
# View(metadata_antibiotics[order(metadata_antibiotics$Patient_ID, metadata_antibiotics$time_of_inf_years), c(1, 3, 5, 10, 19, 38, 40)])

# View other time bins (binary for first and last, <6m, <1y, <6y):
first_last <- table(metadata_antibiotics$resistance_count, metadata_antibiotics$Phenotype_TFAM_first_and_last) # '1' = less than x time (one year), 0 = NA, 2 = last or more than x time
first_last <- first_last[, 2:3]
first_last
time_6m <- table(metadata_antibiotics$resistance_count, metadata_antibiotics$Phenotype_TFAM_6months)
time_6m
time_1y <- table(metadata_antibiotics$resistance_count, metadata_antibiotics$Phenotype_TFAM_chronic_1y)
time_1y
time_1y_by_row <- table(metadata_antibiotics$Phenotype_TFAM_chronic_1y, metadata_antibiotics$resistance_count)
time_1y_by_row
time_1y_multi_resistance <- table(metadata_antibiotics$Phenotype_TFAM_chronic_1y, metadata_antibiotics$multi_resistance)
time_1y_multi_resistance
time_6y <- table(metadata_antibiotics$resistance_count, metadata_antibiotics$Phenotype_TFAM_chronic_6y)
time_6y

# Proportion test:
#table_to_test <- table(metadata_antibiotics$time_bins, metadata_antibiotics$multi_resistance)
output_file <- file(paste('proportions_test_resistance_time_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))

table_to_test <- table(metadata_antibiotics$time_bins, metadata_antibiotics$resistance)
#table_to_test <- time_1y_by_row # first_last, time_6m, time_1y, time_6y
table_to_test

successes <- table_to_test[, 1]
# successes <- table_to_test[, '2'] # when table is by columns (eg first and last are columns)
# successes <- table_to_test['1', ] # when table is by rows (eg less than 1 y = row 1)
successes

total_events <- margin.table(table_to_test, 1) # 1 for rows, 2 for columns, arrange table accordingly
# total_events <- margin.table(table_to_test, 1)
# total_events <- margin.table(table_to_test, 2)
total_events

# This tests whther the proportions of k-groups differ:
prop.test(successes, total_events)

# Test of trend in the proportions. This is a weighted linear regression of the proportions on the group scores where the test is against a zero slope.
prop.trend.test(successes, total_events)
# The proportion of antimicrobial resistant samples increases or decreases (should increase for late samples, decrease for early samples). Not appropriate like this as outcome should be resistant vs sensitive.
# Better to test whether a simple count of resistant samples (at least to one antibiotic) increases over time (binned by quartiles).

sink(file = NULL)
#######################################


#######################################
# Basic plots of phenotype variables:
# Resistance counts by time bins:
ggplot(metadata_antibiotics, aes(time_bins, resistance_count)) + geom_boxplot()
ggsave('boxplots_resistance_count_by_time_bins.png')

# Proportions of resistance by time bins:
png('resistance_time_bins_proportion.png')
plot(metadata_antibiotics$time_bins, metadata_antibiotics$resistance, xlab = 'time by quartiles', ylab = 'Proportion resistant to at least one antimicrobial')
dev.off()

# Resistance counts by patient ID, time bins and RAPD type:
ggplot(metadata_antibiotics, aes(y = resistance_count, x = Patient_ID, colour = RAPD_type_short, size = time_bins)) + 
  geom_point(position = position_jitter(width = 0.2, height = 0.3), alpha = 0.6) +
  theme_bw() + theme(legend.position = 'right', axis.text.x = element_text(angle = 90, hjust = 0)) +
  scale_size_discrete(range = c(2, 7))
ggsave('resistance_counts_by_time_bins_patientID_and_RAPDtype.png')

# Resistance counts by RAPD type: 
ggplot(metadata_antibiotics, aes(x = RAPD_type_short, fill = resistance)) + geom_dotplot(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Set1")
ggsave('Resistance_counts_by_RAPD_type.png')

# Resistance counts by patient ID: 
ggplot(metadata_antibiotics, aes(x = Patient_ID, fill = resistance)) + geom_dotplot(position = position_jitter(width = 0.2, height = 0), alpha = 0.7) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_fill_brewer(palette="Set1")
ggsave('Resistance_counts_by_patient_ID.png')

# Resistance by patient IDs and time bins:
ggplot(metadata_antibiotics, aes(y = Patient_ID, x = time_bins)) + 
  geom_point(shape = 'o', size = 3, position = position_jitter(height = 0.35, width = 0.35), aes(colour = resistance)) + theme_bw()
ggsave('Resistance_by_patient_IDs_and_time_bins.png')

# Resistance by RAPD types and time bins:
ggplot(metadata_antibiotics, aes(y = RAPD_type_short, x = time_bins)) + 
  geom_point(shape = 'o', size = 3, position = position_jitter(height = 0.35, width = 0.35), aes(colour = resistance)) + theme_bw()
ggsave('Resistance_by_RAPD_types_and_time_bins.png')

# Resistance by patient IDs, RAPD types and time bins:
ggplot(metadata_antibiotics, aes(y = Patient_ID, x = time_bins, colour = RAPD_type_short, shape = factor(resistance))) + 
  geom_point(size = 3, position = position_jitter(height = 0.35, width = 0.35)) + theme_bw()
ggsave('Resistance_by_patient_IDs_RAPD_types_and_time_bins.png')

# Resistance by patient IDs, morphotypes and time bins:
ggplot(metadata_antibiotics, aes(y = Patient_ID, x = time_bins, colour = resistance, shape = factor(morphology))) + 
  geom_point(size = 3, position = position_jitter(height = 0.35, width = 0.35)) + theme_bw()
ggsave('Resistance_by_patient_IDs_morphotypes_and_time_bins.png')

# Resistance by patient IDs, morphotypes and time as discrete:
ggplot(metadata_antibiotics, aes(y = Patient_ID, x = time_discrete, colour = resistance, shape = factor(morphology))) + 
  geom_point(size = 3, position = position_jitter(height = 0.35, width = 0.35)) + theme_bw()
ggsave('Resistance_by_patient_IDs_morphotypes_and_time_discrete.png')

#######################################


#######################################
#Further tests:
# Antibiotic resistance accumulates over time:
table(metadata_antibiotics$multi_resistance, metadata_antibiotics$Phenotype_TFAM_chronic_1y)
table(metadata_antibiotics$multi_resistance, metadata_antibiotics$time_bins)
# Basic tests:
kruskal.test(metadata_antibiotics$multi_resistance ~ metadata_antibiotics$time_bins)
kruskal.test(metadata_antibiotics$multi_resistance ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y)
fisher.test(table(metadata_antibiotics$multi_resistance, metadata_antibiotics$Phenotype_TFAM_chronic_1y))
fisher.test(table(metadata_antibiotics$multi_resistance, metadata_antibiotics$time_bins))


## Individual antibiotics and reistance:
# Ciprofloxacin shows the highest resistance (simple counts):
summary(antibiotics_data)
table(metadata_antibiotics$TM_EUCAST, metadata_antibiotics$time_bins)
table(metadata_antibiotics$AT_EUCAST, metadata_antibiotics$time_bins)
table(metadata_antibiotics$IP_EUCAST, metadata_antibiotics$time_bins)
table(metadata_antibiotics$CO_EUCAST, metadata_antibiotics$time_bins)
table(metadata_antibiotics$CI_EUCAST, metadata_antibiotics$time_bins)

# Basic tests:
fisher.test(table(metadata_antibiotics$TM_EUCAST, metadata_antibiotics$time_bins))
fisher.test(table(metadata_antibiotics$AT_EUCAST, metadata_antibiotics$time_bins))
fisher.test(table(metadata_antibiotics$IP_EUCAST, metadata_antibiotics$time_bins))
fisher.test(table(metadata_antibiotics$CI_EUCAST, metadata_antibiotics$time_bins))
fisher.test(table(metadata_antibiotics$CO_EUCAST, metadata_antibiotics$time_bins))

fisher.test(table(metadata_antibiotics$TM_EUCAST, metadata_antibiotics$Phenotype_TFAM_chronic_1y))
fisher.test(table(metadata_antibiotics$AT_EUCAST, metadata_antibiotics$Phenotype_TFAM_chronic_1y))
fisher.test(table(metadata_antibiotics$IP_EUCAST, metadata_antibiotics$Phenotype_TFAM_chronic_1y))
fisher.test(table(metadata_antibiotics$CI_EUCAST, metadata_antibiotics$Phenotype_TFAM_chronic_1y))
fisher.test(table(metadata_antibiotics$CO_EUCAST, metadata_antibiotics$Phenotype_TFAM_chronic_1y))

# Test variances and homoskedasticity before doing ANOVA
# http://www.r-bloggers.com/analysis-of-variance-anova-for-multiple-comparisons/
bartlett.test(x = metadata_antibiotics$CI, g = metadata_antibiotics$Phenotype_TFAM_chronic_1y)
fligner.test(x = metadata_antibiotics$CI, g = metadata_antibiotics$Phenotype_TFAM_chronic_1y)
# variances are homogeneous if p-value > 0.05

fit <- lm(metadata_antibiotics$CI ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y)
anova(fit)

t.test(metadata_antibiotics$TM ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y, var.equal = FALSE)
t.test(metadata_antibiotics$AT ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y, var.equal = FALSE)
t.test(metadata_antibiotics$IP ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y, var.equal = FALSE)
t.test(metadata_antibiotics$CI ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y, var.equal = FALSE)
t.test(metadata_antibiotics$CO ~ metadata_antibiotics$Phenotype_TFAM_chronic_1y, var.equal = FALSE)

# If p-value > 0.05, accept the null hypothesis H0: the four means are statistically equal.


# Basic plots of morphotypes and resistance for each antibiotic:
png('morphotypes_and_resistance_per_antibiotic.png', height = 15, width = 15, units = 'in', res = 300)
par(mfrow = c(4,2), las = 2)
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$TM_EUCAST))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$AT_EUCAST))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$IP_EUCAST))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$CI_EUCAST))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$CO_EUCAST))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$multi_resistance))
plot(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$resistance_count))
dev.off()

# Basic plots of Patient IDs and resistance for each antibiotic:
png('Patient_IDs_and_resistance_per_antibiotic.png', height = 15, width = 15, units = 'in', res = 300)
par(mfrow = c(2,4), las = 2)
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$TM_EUCAST))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$AT_EUCAST))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$IP_EUCAST))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$CI_EUCAST))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$CO_EUCAST))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$multi_resistance))
plot(as.factor(metadata_antibiotics$Patient_ID), as.factor(metadata_antibiotics$resistance_count))
dev.off()

# Basic plots of RAPD types and resistance for each antibiotic:
png('RAPD_types_and_resistance_per_antibiotic.png', height = 15, width = 15, units = 'in', res = 300)
par(mfrow = c(2,4), las = 2)
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$TM_EUCAST))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$AT_EUCAST))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$IP_EUCAST))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$CI_EUCAST))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$CO_EUCAST))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$multi_resistance))
plot(as.factor(metadata_antibiotics$RAPD_type_short), as.factor(metadata_antibiotics$resistance_count))
dev.off()
#######################################


#######################################
# MannKendall trend tests for morphotypes, antibiotic resistance over time, etc.

# Test trends for antibiotic resistance accumulation over time. See MannKendall package?:
# See https://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&ved=0CCMQFjABahUKEwjBhKbVpYjIAhULtBQKHcXrDOs&url=http%3A%2F%2Fwww.ghement.ca%2FMann-Kendall%2520Trend%2520Test%2520in%2520R.doc&usg=AFQjCNEO0yK7QxSgQFIwzs0WsObjQMCDNg&cad=rja
# data(PrecipGL)
# View(PrecipGL)
# ?PrecipGL

# Test whether resistance count to antimicrobials increases over time:
head(metadata_antibiotics[-13, c(11, 37, 40)])
resistance_count_sorted_noNAs <- metadata_antibiotics[-13, c(11, 37, 40)]
head(resistance_count_sorted_noNAs)
resistance_count_sorted_noNAs <- resistance_count_sorted_noNAs[order(resistance_count_sorted_noNAs$time_of_inf_years), 2]
head(resistance_count_sorted_noNAs)

# Check assumptions of correlations:
png('resistance_count_sorted_noNAs_MannK_assumptions_of_autocorrelations.png')
par(mfrow = c(2,1))
acf(resistance_count_sorted_noNAs)
pacf(resistance_count_sorted_noNAs)
dev.off()
# Seems OK according to example

png('resistance_count_sorted_noNAs_MannK.png')
plot(resistance_count_sorted_noNAs)
lines(lowess(time(resistance_count_sorted_noNAs), resistance_count_sorted_noNAs), col="blue", lwd=2)
dev.off()


# Run test:
output_file <- file(paste('MannK_test_resistance_count_time_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))
print('Mann Kendall test of trends for time series data: antimicrobial resistance counts dorted by time of infections years (original variable coding)')
trend_MannK_test <- MannKendall(resistance_count_sorted_noNAs)
print(trend_MannK_test)
summary(trend_MannK_test)

sink(file = NULL)

# TO DO, bootstrap gives error: 
# Compute bootstrap and CI:
#Use block bootstrap 
#data(PrecipGL)
#MannKendall(PrecipGL)
MKtau <- function(z) MannKendall(z)$tau
#tsboot(PrecipGL, MKtau, R=500, l=5, sim="fixed")
#?tsboot
boot_test <- tsboot(resistance_count_sorted_noNAs, MKtau, R=1000, sim='model')
boot_test

# Get 95% CI from bootstrap:
#?boot.ci
boot.ci(boot_test, type='all')


## Check each antibiotic, raw values sorted by time from first sample:
head(metadata_antibiotics[-13, c(11, 27:31, 37:41)])
amb_all_sorted_noNAs <- metadata_antibiotics[-13, c(11,27:31, 37:41)]
head(amb_all_sorted_noNAs)
amb_sorted_noNAs <- amb_all_sorted_noNAs[order(amb_all_sorted_noNAs$time_of_inf_years), 6]
head(amb_sorted_noNAs)

#qplot(y = amb_EUCAST_sorted_noNAs) + geom_point()
# plot(amb_sorted_noNAs)
# lines(lowess(time(amb_sorted_noNAs), amb_sorted_noNAs), col="blue", lwd=2)
#dev.off()

# Check assumptions of correlations:
# acf(amb_sorted_noNAs)
# pacf(amb_sorted_noNAs)
# Seems OK according to example

# Run test:
trend_MannK_test <- MannKendall(amb_sorted_noNAs)
print(trend_MannK_test)
summary(trend_MannK_test)

# TO DO, bootstrap gives error/warnings: 
# Compute bootstrap and CI:
MKtau <- function(z) MannKendall(z)$tau
boot_test <- tsboot(amb_sorted_noNAs, MKtau, R=1000, sim='model')
boot_test <- tsboot(amb_sorted_noNAs, MKtau, R=1000, sim='fixed', l = 5)
boot_test
# Get 95% CI from bootstrap:
boot.ci(boot_test, type='all')

print('# TM, AT, IP, CI, are significant for tests of upward trend with MannK (increased acquisition of 
      resistance over time for time bins. CO is not. ')
#######################################



#######################################
# Test morphotypes:
# Morphotypes aren't double counted, if there was more than one sampling per time-point this was because more than one morphotype was present in the dish:
# First and last variable keeps the same morphotype where there is more than one (this is the case for P003 and P008):

table(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_first_and_last)
table(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_chronic_6y)
table(metadata_antibiotics$morphology, metadata_antibiotics$time_discrete)

# Colony and mucoid accumulate in later time-points but not dwarf:
# Proportion test:
output_file <- file(paste('proportions_test_morphotypes_time_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))

table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$time_bins)
table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_chronic_1y)
print('table(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_chronic_1y)')
table_to_test

successes <- table_to_test[5, ]
successes

total_events <- margin.table(table_to_test, 2) # 1 for rows, 2 for columns, arrange table accordingly
total_events

# This tests whther the proportions of k-groups differ:
prop.test(successes, total_events)

# Test of trend in the proportions. This is a weighted linear regression of the proportions on the group scores where the test is against a zero slope.
prop.trend.test(successes, total_events)
# Use only chi square for two proportions (ie time cut-offs at 6m, 1y, etc.)
# Only mucoid at 1y (only tested time bins and 1y) has proportions that are significantly different, dwarf is borderline. Mucoid is increased late and dwarf early.

sink(file = NULL)

kruskal.test(metadata_antibiotics$morphology, metadata_antibiotics$time_bins)
kruskal.test(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_chronic_1y)

fisher.test(table(metadata_antibiotics$morphology, metadata_antibiotics$time_bins))
fisher.test(table(metadata_antibiotics$morphology, metadata_antibiotics$Phenotype_TFAM_chronic_1y))

# Dwarf has more resistant types but not mucoid or colony:
# Basic tests:
# Proportions test:
output_file <- file(paste('proportions_test_morphotypes_resistance_', Sys.Date(), '.txt', sep=''))
sink(output_file, append = FALSE, split = TRUE, type = c("output", "message"))

table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$multi_resistance)
table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$resistance_count)
table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$resistance)
print('table_to_test <- table(metadata_antibiotics$morphology, metadata_antibiotics$resistance)')
table_to_test

successes <- table_to_test[5, ]
successes

total_events <- margin.table(table_to_test, 2) # 1 for rows, 2 for columns, arrange table accordingly
total_events

# This tests whther the proportions of k-groups differ:
prop.test(successes, total_events)

# Test of trend in the proportions. This is a weighted linear regression of the proportions on the group scores where the test is against a zero slope.
prop.trend.test(successes, total_events)
# Use only chi square for two proportions (ie time cut-offs at 6m, 1y, etc.)
print('# None of the morphotypes are significantly different for amb resistance to any one amb count(resistance_count). 
# Dwarf is significantly more multi-resistant but not any of the others.
')
sink(file = NULL)

kruskal.test(as.factor(metadata_antibiotics$multi_resistance) ~ as.factor(metadata_antibiotics$morphology))
kruskal.test(as.factor(metadata_antibiotics$resistance_count) ~ as.factor(metadata_antibiotics$morphology))
fisher.test(table(as.factor(metadata_antibiotics$morphology), as.factor(metadata_antibiotics$multi_resistance)))
#######################################


#######################################
# Write antibiotic classes and additional phenotypes to file:
getwd()
#write.table(antibiotics_data, paste('antibiotic_resistance_EUCAST_', Sys.Date(), '.tsv', sep = ''), row.names = FALSE, 
#		col.names = TRUE, quote = FALSE, sep = '\t')
#write.table(metadata_antibiotics, paste('metadata_and_antibiotics_data_final_', Sys.Date(), '.tsv', sep = ''), row.names = FALSE, 
#		col.names = TRUE, quote = FALSE, sep = '\t')
#######################################


#######################################
# Metadata update: 
# Updated IDs for plink manually, read here: 
metadata_plink <- read.csv('metadata_strain_ids_pseudomonas_12Mar2016_plink_IDs.tsv', sep = '\t', header = TRUE, 
                     na.string = c('-99999', "", " ", "NA"), stringsAsFactors = F)

metadata_plink <- metadata_plink[, moveme(names(metadata_plink), 'unique_ID_fastq first')]
# Keep only relevant columns:
metadata_plink <- metadata_plink[, c('unique_ID_fastq', 'IID', 'FID')]
head(metadata_plink)
which(is.na(metadata_plink$IID))
which(is.na(metadata_plink$FID))

# Merge IDs with full phenotype/metadata file:
metadata_final <- merge(metadata_plink, metadata_antibiotics)
colnames(metadata_final)
dim(metadata_final)
dim(metadata_antibiotics)
head(metadata_final)

# Prep for plink:
metadata_final <- metadata_final[, moveme(names(metadata_final), 'Patient_ID first')]
colnames(metadata_final)
colnames(metadata_final)[4] <- 'FID_old' 
colnames(metadata_final)[3] <- 'IID_old'
colnames(metadata_final)[2] <- 'IID'
colnames(metadata_final)[1] <- 'FID'
colnames(metadata_final)
# Remove S18 which has NAs and no sequencing data, TO DO: update labels:
which(is.na(metadata_final[, 3]))
# View(metadata_final)
metadata_final <- metadata_final[-16, ]

summary(metadata_final$TM)
which(is.na(metadata_final$TM))
metadata_final[which(is.na(metadata_final$TM)), 'TM']

summary(as.factor(metadata_final$TM_EUCAST))
summary(as.factor(metadata_final$IP_EUCAST))

# write.table(metadata_final, 'metadata_final_pseudomonas.ssv', row.names = FALSE,
#                 col.names = TRUE, quote = FALSE, sep = ' ', na = '-99999')

# This didn't write out decimals despite options(), format() ensures this but then writes NAs 
# as NA, not as specified string. 

# Plink errors (not scalar variable) if --linear is used but variable doesn't
# have decimals (eg 1 vs 1.0001, doesn't seem to work with 1.000 either). So adding 0.001 to each antibiotic variable:
metadata_final$TM_new <- metadata_final$TM + 0.001
metadata_final$IP_new <- metadata_final$IP + 0.001
metadata_final$AT_new <- metadata_final$AT + 0.001
metadata_final$CI_new <- metadata_final$CI + 0.001
metadata_final$CO_new <- metadata_final$CO + 0.001

names(metadata_final)

# Fix date ('original_date') as prints with spaces and messes delimiters:
class(metadata_final$original_date)
head(metadata_final$original_date)
metadata_final$original_date_2 <- as.POSIXlt(metadata_final$original_date, tz = 'GMT', '%d %B %y')
metadata_final$original_date_2 <- format(metadata_final$original_date_2, '%d/%B/%y')
class(metadata_final$original_date_2)
head(metadata_final$original_date_2)
which(colnames(metadata_final) == 'original_date')
metadata_final <- metadata_final[, -15]
colnames(metadata_final)

# Format data frame and write to disk:  
metadata_final_formatted <- format(metadata_final, digits = 4, na.encode = T)

write.table(metadata_final_formatted, 'metadata_final_pseudomonas.ssv', row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = ' ', na = '-99999')
system("sed -i 's/NA/-99999/g' metadata_final_pseudomonas.ssv")

write.table(metadata_final_formatted, 'metadata_final_pseudomonas.tsv', 
            na = '-99999',
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = '\t')
system("sed -i 's/NA/-99999/g' metadata_final_pseudomonas.tsv")

#######################################


#######################################
# The end
q()

#######################################

