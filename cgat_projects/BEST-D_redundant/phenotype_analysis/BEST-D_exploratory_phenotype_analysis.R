#############################
# Antonio J Berlanga-Taylor
# 22 Feb 2016
# BEST-D project phenotype/metadata 
  # File processing
  # Multi-co-linearity check
  # Determine variables to pass to plink as covariates
  # Exploratory analysis
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/ifs/projects/proj043/analysis.dir/gene_expression_2.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_exploratory_analysis",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is project specific.

# Load results from 02_microarrayxxx file, saved as RData object:
# Re-load a previous R session, data and objects:
#load('R_session_saved_image_probe_filtering.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_exploratory_analysis', '.RData', sep='')

#############################


##########################
#source('http://bioconductor.org/biocLite.R')
#biocLite('lubridate')
library(lubridate)
library(plyr)
library(ggplot2)
library(Hmisc)

#############################



#############################
# TO DO: clean this section up, duplicated with 03bis_2_quant_assoc file, also messy, better to re-run 
# from sqlite3 steps. Move covariates adjustment tests to different file. Keep merging and pheno file 
# processing separate.
# TO DO: NA handling, missing cases, create functions to test all.

# Load full phenotype data for covariates adjustment:
# TO DO: Load pheno file with the 305 individuals, this has 298 (which are the ones with genotype date)

#phenotype_data <- read.csv('BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt', sep = '\t', 
#                           header = TRUE, na.string = c('-99999', "", " "))
phenotype_data <- read.csv('BEST-D_phenotype_file_final.tsv', sep = '\t', 
                           header = TRUE, na.string = c('-99999', "", " ", "NA"))
head(phenotype_data)
tail(phenotype_data)
summary(phenotype_data)
str(phenotype_data)
class(phenotype_data)
names(phenotype_data)
length(which(complete.cases(phenotype_data)))
#View(phenotype_data)

# Determine variable types:
# Keep two data frames as NA handling differs when changing variable data types (NA vs <NA> in numeric vs factor for example):
col_names_character <- c("pt_id")
phenotype_data_vars <- phenotype_data
phenotype_data_vars[, col_names_character] <- as.character(unlist((phenotype_data_vars[, col_names_character])))
str(phenotype_data_vars$pt_id)

col_names_factors <- c("arm", "incident_fall", "incident_fracture", "incident_resp_infection", "male", 
                       "diabetes", "hypertension", "heart_disease", "stroke", "fracture", "falls", 
                       "copd_asthma", "basemed_antihypertensive", "basemed_statin", "basemed_antithrombotic", 
                       "basemed_vitamind", "basemed_calcium", "calendar_age_ra", 
                       "currsmoker", "exsmoker", "muscle_pain0", "joint_pain0","muscle_pain12", "joint_pain12")

#col_names_factors_2 <- c("placebo_v_Tx", "placebo_v_4000", "placebo_v_2000", "IU4000_v_2000")
# These last variables come from the plink genotype pheno file but since then I updated the full pheno file 
# to include all samples (n=305).

head(phenotype_data_vars[, col_names_factors])

# This returns chr instead of factor because of conversion to matrix, not data.frame (?):
#phenotype_data[, col_names_factors] <- as.factor(unlist(phenotype_data[, col_names_factors]))
# See:
# http://stackoverflow.com/questions/2392216/why-does-as-factor-return-a-character-when-used-inside-apply
# Use:
phenotype_data_vars[, col_names_factors] <- lapply(phenotype_data_vars[, col_names_factors], as.factor)

# Convert <NA> back to NA:
str(phenotype_data_vars[, col_names_factors])
phenotype_data_vars$arm
summary(phenotype_data_vars[, col_names_factors])
head(phenotype_data_vars[, col_names_factors])

# Rename file back:
phenotype_data <- phenotype_data_vars

# Some plots:
ggplot(aes(x = male, y = vitd0), data = phenotype_data) + geom_boxplot()
ggplot(aes(x = vitd12, y = vitd0), data = phenotype_data) + geom_point()

summary(phenotype_data$vitd0)
summary(phenotype_data$vitd6)
summary(phenotype_data$vitd12)
boxplot(phenotype_data$vitd0, phenotype_data$vitd6, phenotype_data$vitd12)
qplot(phenotype_data$vitd0)
qplot(log2(phenotype_data$vitd0))
qplot(phenotype_data$vitd6)
qplot(phenotype_data$vitd12)

# Variables to adjust for in association analysis and differential gene expression:
# age, gender, season, and BMI as in GWASs
# Also correct for vitD intake, baseline vitD, smoking, disease, medications used (?) and
# sun exposure (not measured but could use physical activity as proxy)
# cell count (not measured) needs to be accounted for in gene expression (use PCs?)
# incident_x refers to incidents during follow-up

# Check variables to pass to plink:
names(phenotype_data)
plink_adjust_phenotype_data <- phenotype_data[, c("pt_id", "incident_fall", "incident_fracture", "incident_resp_infection", # incidents
                                                       "bmi0", "BMI_DEXA12", "vitd0", "vitd6", "vitd12", #numerics
                                                       "diabetes", "hypertension", "heart_disease", "stroke", "fracture", "falls", "copd_asthma", #medical history
                                                       "muscle_pain0", "joint_pain0", "muscle_pain12", "joint_pain12", # pain
                                                       "physical_activity12", "physical_activity0", # activity
                                                       "basemed_antihypertensive", "basemed_statin", "basemed_antithrombotic", "basemed_vitamind", "basemed_calcium", # medications
                                                       "calendar_age_ra", 'male', "season_randomisation_2", # demographics
                                                       "currsmoker", "exsmoker" # smoking
                                                       )]

incidents <- phenotype_data_vars[, c("pt_id", "incident_fall", "incident_fracture", "incident_resp_infection")]
med_history <-phenotype_data_vars[, c("pt_id", "diabetes", "hypertension", "heart_disease", "stroke", "fracture", "falls", "copd_asthma")] #medical history
bmi_vitd <- phenotype_data_vars[, c("pt_id", "bmi0", "BMI_DEXA12", "vitd0", "vitd6", "vitd12")] #numerics
pain <- phenotype_data_vars[, c("pt_id", "muscle_pain0", "joint_pain0", "muscle_pain12", "joint_pain12")] # pain
activity <- phenotype_data_vars[, c("pt_id", "physical_activity12", "physical_activity0")] # activity
medications <- phenotype_data_vars[, c("pt_id", "basemed_antihypertensive", "basemed_statin", "basemed_antithrombotic", "basemed_vitamind", "basemed_calcium")] # medications
demographics <- phenotype_data_vars[, c("pt_id", "calendar_age_ra", 'male', "season_randomisation_2")] # demographics
smoking <- phenotype_data_vars[, c("pt_id", "currsmoker", "exsmoker")] # smoking
  
head(plink_adjust_phenotype_data)
summary(plink_adjust_phenotype_data)
str(plink_adjust_phenotype_data)
count(plink_adjust_phenotype_data)

# These have many missing values:
# vitd6 has 8 NAs
many_missing <- c('BMI_DEXA12', 'creatinine12', 'grip_strength12', 'joint_severity12', 'physical_activity12', 'vitd12')
summary(phenotype_data_vars[, many_missing])
head(phenotype_data_vars[, many_missing])

# These are scales (0 to 10, 10 worse) = ("physical_activity0",  "physical_activity12", "muscle_severity12", , "joint_severity12")
str(phenotype_data_vars[, c("physical_activity0",  "physical_activity12", "muscle_severity12", "joint_severity12")])

# TO DO: Check correlation between variables:
rcorr(incidents[, -1])
rcorr(as.matrix(incidents[, -1]))
rcorr(as.matrix(med_history[, -1]))
rcorr(as.matrix(bmi_vitd[, -1]))
rcorr(as.matrix(pain[, -1]))
rcorr(as.matrix(activity[, -1]))
rcorr(as.matrix(medications[, -1]))
rcorr(as.matrix(demographics[, -1]))
rcorr(as.matrix(smoking[, -1]))

# Keep only baselines?
# Keep incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, muscle_pain0, joint_pain0, physical_activity0, 
# basemed_antihypertensive, basemed_vitamind, currsmoker, exsmoker, bmi0, calendar_age_ra, season_randomisation_2
# Also male, vitd0
# Check: fracture falls

vars_to_keep <- c('pt_id', 'incident_fracture', 'incident_resp_infection', 'diabetes', 'heart_disease', 'copd_asthma', 'muscle_pain0', 'joint_pain0', 'physical_activity0', 
             'basemed_antihypertensive', 'basemed_vitamind', 'currsmoker', 'exsmoker', 'bmi0', 'calendar_age_ra', 'season_randomisation_2',
             'male', 'vitd0', 'vitd6', 'vitd12',
             'fracture', 'falls')

phenotype_data_vars_to_keep <- phenotype_data_vars[, vars_to_keep]
summary(phenotype_data_vars_to_keep)
rcorr(as.matrix(phenotype_data_vars_to_keep[, -1]))

# Keep incidents
# Exclude: fracture, muscle_pain0, falls, exsmoker, joint_pain0, physical_activity0, basemed_antihypertensive
# bmi0, male, age, correlate with many variables...
vars_to_keep_2 <- c('pt_id', 'incident_fracture', 'incident_resp_infection', 
                    'diabetes', 'heart_disease', 'copd_asthma',
                    'basemed_vitamind', 'currsmoker', 
                    'bmi0', 'calendar_age_ra', 'season_randomisation_2',
                    'male', 'vitd0', 'vitd6', 'vitd12')

phenotype_data_vars_to_keep_2 <- phenotype_data_vars[, vars_to_keep_2]
summary(phenotype_data_vars_to_keep_2)
rcorr(as.matrix(phenotype_data_vars_to_keep_2[, -1]), type = 'spearman')

# Pass to plink: 
#vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2
# male, 'vitd6', 'vitd12'

#heatmap_P <- heatmap(as.matrix(rcorr_plink_vars_df_P))
#######################


#######################
# Explore likely variables that can confound vitd levels:
# TO DO: Adjusted for in assoc analysis in plink: 
# vitd0, bmi0, calendar_age_ra
# incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker
# season_randomisation_2, male, 'vitd6', 'vitd12'

# TO DO: run non-parametric

#rcorr(na.omit(as.matrix(phenotype_data$vitd0)), na.omit(as.matrix(phenotype_data[, 5])), type = 'spearman')

# Normality test:
shapiro.test(phenotype_data$vitd0)
shapiro.test(log2(phenotype_data$vitd0))
# Plot normality:
qqnorm(phenotype_data$vitd0); qqline(phenotype_data$vitd0, col = 2)
qqnorm(log2(phenotype_data$vitd0)); qqline(log2(phenotype_data$vitd0), col = 2)


rcorr(phenotype_data$vitd0, phenotype_data$bmi0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$bmi0, type = 'pearson')

rcorr(phenotype_data$vitd0, phenotype_data$calendar_age_ra, type = 'pearson')
rcorr(phenotype_data$vitd0, phenotype_data$calendar_age_ra, type = 'spearman')

ggplot(phenotype_data, aes(y = vitd0, bmi0)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

ggplot(phenotype_data, aes(y = vitd0, calendar_age_ra)) + 
  geom_point(shape=1) +    # Use hollow circles
  geom_smooth(method=lm)   # Add linear regression line, (by default includes 95% confidence region)

t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$male))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$male))
ggplot(phenotype_data, aes(log2(vitd0), fill = as.factor(male))) + 
  geom_histogram(binwidth=.5, alpha=.5, position="identity")

ggplot(phenotype_data, aes(as.factor(male), vitd0)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

summary(aov(phenotype_data$vitd0 ~ as.factor(phenotype_data$season_randomisation_2)))
kruskal.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$season_randomisation_2))

ggplot(phenotype_data, aes(as.factor(season_randomisation_2), vitd0)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggplot(phenotype_data, aes(y = vitd0 , as.factor(month_randomisation))) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

summary(aov(phenotype_data$vitd0 ~ as.factor(phenotype_data$arm)))
ggplot(phenotype_data, aes(as.factor(arm), vitd0)) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

kruskal.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$arm))

ggplot(phenotype_data, aes(log2(vitd0), fill = as.factor(arm))) + 
  geom_histogram(binwidth=.5, alpha=.5, position="identity")

ggplot(phenotype_data, aes(log2(vitd0), fill = as.factor(arm))) + 
  geom_density(binwidth=.5, alpha=.5, position="identity")

ggplot(phenotype_data, aes(vitd0, fill = as.factor(arm))) + 
  geom_density(binwidth=.5, alpha=.5, position="identity")

# incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker
# TO DO: What to use? Incident or prevalent data? Fracture varies for example:
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$incident_fracture))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$fracture))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$incident_fracture))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$fracture))

ggplot(phenotype_data, aes(y = vitd0, as.factor(incident_fracture))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggplot(phenotype_data, aes(y = vitd0, as.factor(fracture))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$incident_resp_infection))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$incident_resp_infection))
ggplot(phenotype_data, aes(y = vitd0, as.factor(incident_resp_infection))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$diabetes))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$heart_disease))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$copd_asthma))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$currsmoker))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$exsmoker))

wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$diabetes))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$heart_disease))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$copd_asthma))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$currsmoker))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$exsmoker))

t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_vitamind))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_vitamind))
ggplot(phenotype_data, aes(y = vitd0, as.factor(basemed_vitamind))) + 
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# Others not in plink assoc adjustment:
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$hypertension))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$stroke))

wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$hypertension))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$stroke))

summary(aov(phenotype_data$vitd0 ~ as.factor(phenotype_data$physical_activity0)))
kruskal.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$physical_activity0))
ggplot(phenotype_data, aes(y = vitd0, physical_activity0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$falls))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$muscle_pain0))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$joint_pain0))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_statin))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_antihypertensive))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_antithrombotic))
t.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_calcium))

wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$falls))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$muscle_pain0))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$joint_pain0))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_statin))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_antihypertensive))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_antithrombotic))
wilcox.test(phenotype_data$vitd0 ~ as.factor(phenotype_data$basemed_calcium))

rcorr(phenotype_data$vitd0, phenotype_data$egfr0)
rcorr(phenotype_data$vitd0, phenotype_data$egfr0, type = 'spearman')
ggplot(phenotype_data, aes(egfr0, vitd0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(phenotype_data$vitd0, phenotype_data$foot_tscore)
rcorr(phenotype_data$vitd0, phenotype_data$foot_tscore, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$forearm_tscore)
rcorr(phenotype_data$vitd0, phenotype_data$forearm_tscore, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$albumin0)
rcorr(phenotype_data$vitd0, phenotype_data$albumin0, type = 'spearman')
ggplot(phenotype_data, aes(albumin0, vitd0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(phenotype_data$vitd0, phenotype_data$alk_phosphatase0)
rcorr(phenotype_data$vitd0, phenotype_data$albumin0, type = 'spearman')
ggplot(phenotype_data, aes(alk_phosphatase0, vitd0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(phenotype_data$vitd0, phenotype_data$phosphate0)
rcorr(phenotype_data$vitd0, phenotype_data$phosphate0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$apo_a10)
rcorr(phenotype_data$vitd0, phenotype_data$apo_a10, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$apo_b0)
rcorr(phenotype_data$vitd0, phenotype_data$apo_b0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$tchol0)
rcorr(phenotype_data$vitd0, phenotype_data$tchol0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$creatinine0)
rcorr(phenotype_data$vitd0, phenotype_data$creatinine0, type = 'spearman')

rcorr(phenotype_data$vitd0, phenotype_data$ipth0)
rcorr(phenotype_data$vitd0, phenotype_data$ipth0, type = 'spearman') 
ggplot(phenotype_data, aes(ipth0, vitd0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(phenotype_data$vitd0, phenotype_data$trig0)
rcorr(phenotype_data$vitd0, phenotype_data$ldlc0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$ldlc0, type = 'pearson')
rcorr(phenotype_data$vitd0, phenotype_data$hdlc0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$hdlc0, type = 'pearson')

rcorr(phenotype_data$vitd0, phenotype_data$crp0, type = 'pearson')
rcorr(phenotype_data$vitd0, phenotype_data$crp0, type = 'spearman')

rcorr(phenotype_data$vitd0, phenotype_data$corrected_calcium0)
rcorr(phenotype_data$vitd0, phenotype_data$corrected_calcium0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$grip_strength12)
rcorr(phenotype_data$vitd0, phenotype_data$grip_strength12, type = 'spearman')
ggplot(phenotype_data, aes(grip_strength12, vitd0)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(phenotype_data$vitd0, phenotype_data$art_sbp0, type = 'pearson')
rcorr(phenotype_data$vitd0, phenotype_data$art_sbp0, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$art_dbp0, type = 'pearson')
rcorr(phenotype_data$vitd0, phenotype_data$art_dbp0, type = 'spearman')

rcorr(phenotype_data$vitd0, phenotype_data$ptr_si0)#, type = 'spearman')
rcorr(phenotype_data$vitd0, phenotype_data$ptr_ri0)#, type = 'spearman')

## Significant results were:
# bmi0, ipth0, age, season, basemed_vitamind, fracture, physical_activity0, albumin0
# alk_phosphatase0, incident_resp_infection, grip_strength12

#############################


#############################
# TO DO: run comparisons for vitd6, vitd12, subset for paired tests
# vitd12

## Subset groups:
dim(phenotype_data)
UI4000_group <- subset(phenotype_data, subset = (phenotype_data$arm == 0))
head(UI4000_group)
dim(UI4000_group)

UI2000_group <- subset(phenotype_data, subset = (phenotype_data$arm == 1))
dim(UI2000_group)

placebo_group <- subset(phenotype_data, subset = (phenotype_data$arm == 2)) 
dim(placebo_group)

# Check normality for vitd12:
shapiro.test(UI4000_group$vitd12)
shapiro.test(UI2000_group$vitd12)

ggplot(phenotype_data, aes(vitd12)) +
  geom_histogram() #binwidth=.5, alpha=.5, position="identity")

ggplot(phenotype_data, aes(vitd12, fill = as.factor(arm))) +
  geom_density(binwidth=.5, alpha=.5, position="identity") +
  scale_fill_brewer(palette = 3)

ggplot(phenotype_data, aes(vitd12, fill = as.factor(arm))) +
  geom_histogram() +
  scale_fill_brewer(palette = 2)

# Check relationship with other variables:
# TO DO: create function and run systematically for all variables.
rcorr(UI4000_group$vitd12, UI4000_group$crp12)
rcorr(UI4000_group$vitd12, UI4000_group$crp12, type = 'spearman')
rcorr(UI2000_group$vitd12, UI2000_group$crp12, type = 'spearman')
rcorr(UI4000_group$vitd12, UI4000_group$calendar_age_ra)
rcorr(UI4000_group$vitd12, UI4000_group$calendar_age_ra, type = 'spearman')
rcorr(UI4000_group$vitd12, UI4000_group$BMI_DEXA12)
ggplot(UI4000_group, aes(BMI_DEXA12, vitd12)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(UI4000_group$vitd12, UI4000_group$ipth12, type = 'spearman')
ggplot(UI4000_group, aes(ipth12, vitd12)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

rcorr(UI2000_group$vitd12, UI2000_group$ipth12, type = 'spearman')
ggplot(UI2000_group, aes(ipth12, vitd12)) +
  geom_point(shape=1) +
  geom_smooth(method=lm)

# Respiratory infection seems significant:
wilcox.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$incident_resp_infection))
t.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$incident_resp_infection))
ggplot(UI4000_group, aes(as.factor(incident_resp_infection), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

wilcox.test(UI2000_group$vitd12 ~ as.factor(UI2000_group$incident_resp_infection))
t.test(UI2000_group$vitd12 ~ as.factor(UI2000_group$incident_resp_infection))
ggplot(UI2000_group, aes(as.factor(incident_resp_infection), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

wilcox.test(placebo_group$vitd12 ~ as.factor(placebo_group$incident_resp_infection))
t.test(placebo_group$vitd12 ~ as.factor(placebo_group$incident_resp_infection))
ggplot(placebo_group, aes(as.factor(incident_resp_infection), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)


# Continue with other variables:
wilcox.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$incident_fracture))
ggplot(UI4000_group, aes(as.factor(incident_fracture), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

wilcox.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$incident_fall))
ggplot(UI4000_group, aes(as.factor(incident_fall), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

rcorr(UI4000_group$vitd12, UI4000_group$albumin12)
rcorr(UI4000_group$vitd12, UI4000_group$albumin12, type = 'spearman')

rcorr(UI4000_group$vitd12, UI4000_group$physical_activity12)
rcorr(UI4000_group$vitd12, UI4000_group$physical_activity12, type = 'spearman')

rcorr(UI4000_group$vitd12, UI4000_group$alk_phosphatase12)
rcorr(UI4000_group$vitd12, UI4000_group$alk_phosphatase12, type = 'spearman')

wilcox.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$basemed_vitamind))
t.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$basemed_vitamind))
ggplot(UI4000_group, aes(as.factor(basemed_vitamind), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggplot(UI2000_group, aes(as.factor(basemed_vitamind), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

kruskal.test(UI4000_group$vitd12 ~ as.factor(UI4000_group$season_randomisation_2))
summary(aov(UI4000_group$vitd12 ~ as.factor(UI4000_group$season_randomisation_2)))
ggplot(UI4000_group, aes(as.factor(season_randomisation), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggplot(UI2000_group, aes(as.factor(season_randomisation), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)
ggplot(UI2000_group, aes(as.factor(month_randomisation), vitd12)) +
  geom_boxplot() + geom_point(position = position_jitter(width = 0.2), alpha = 0.5)

# Significant results were:
# BMI_DEXA12, ipth12, basemed_vitamind, incident_resp_infection (4000 group only and higher rate in the placebo group)
# TO DO: run multi-variate analysis, test systematically, check trial papers to avoid duplicating (BMI was the only significant there).
#############################

#############################
# TO DO: write report to file with plots, etc.
# Write results to file:
#write.table(x = xx, 'xxx.txt', quote = FALSE, sep = '\t',
#            na = '-99999', row.names = FALSE, col.names = TRUE)
#############################


#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run association analysis with genotypes and also GEx considering covariates.
#############################