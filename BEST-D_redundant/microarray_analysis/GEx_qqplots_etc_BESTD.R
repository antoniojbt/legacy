###########################
# Exploration of variance and distributions

# F tests, qqplotting, etc.
###########################


###########################
# Re-load a previous R session, data and objects:
# setwd('/Users/antoniob/Desktop/BEST_D.DIR/mac_runs_to_upload/tables_and_plots_for_draft/final_draft_BEST-D/tables/all_GEx_tables/individual_tables_all_GEx_comparisons/')
load('R_session_saved_image_diff_expression_merge_tables.RData', verbose=T)
###########################



###########################
#########
# Regression tests distributions:
qqnorm(final_table$t.UI4000minusplacebo)
qqnorm(final_table$t.UI2000minusplacebo)
qqnorm(final_table$t.joint_pairedtreated_placebo)

qqnorm(final_table$P.Value.UI4000minusplacebo)
qqnorm(final_table$P.Value.UI2000minusplacebo)
qqnorm(final_table$P.Value.joint_pairedtreated_placebo)

qqnorm(final_table$adj.P.Val.UI4000minusplacebo)
qqnorm(final_table$adj.P.Val.UI2000minusplacebo)
qqnorm(final_table$adj.P.Val.joint_pairedtreated_placebo)
qqnorm(final_table$`adj.P.Val.joint_pairedtreated_2000+4000`)

qqnorm(final_table$logFC.UI4000minusplacebo)
qqnorm(final_table$logFC.UI2000minusplacebo)
qqnorm(final_table$logFC.joint_pairedtreated_placebo)

qqnorm(final_table$B.UI4000minusplacebo)
qqnorm(final_table$B.UI2000minusplacebo)
qqnorm(final_table$B.joint_pairedtreated_placebo)

qqplot(y = final_table$`t.joint_pairedtreated_2000+4000`,
       final_table$t.joint_pairedtreated_placebo)
abline(coef = c(0, 1))

qqplot(y = final_table$`P.Value.joint_pairedtreated_2000+4000`,
       final_table$P.Value.joint_pairedtreated_placebo)
abline(coef = c(0, 1))

qqplot(y = final_table$`adj.P.Val.joint_pairedtreated_2000+4000`,
       final_table$adj.P.Val.joint_pairedtreated_placebo)
abline(coef = c(0, 1))
#########

#########
# Get log10 Pvalues for all raw p-value columns
cols <- grep(pattern = 'P.Value', fixed = TRUE, colnames(final_table))
log10s <- final_table[, c(1:3, cols)]
head(log10s)
colnames(log10s)
dim(log10s)

colnames(log10s)[4:ncol(log10s)]
for (i in colnames(log10s)[4:ncol(log10s)]){
  new_col <- sprintf('neglog10.%s', i)
  # print(new_col)
  log10s[, new_col] <- -log10(log10s[, i])
}
dim(log10s)
colnames(log10s)
colnames(log10s)[1:33]

# Bonferroni:
bonf <- -log10(0.05 / nrow(log10s))

# Plots:
qqnorm(log10s$neglog10.P.Value.UI4000minusplacebo)
qqnorm(log10s$neglog10.P.Value.UI2000minusplacebo)
qqnorm(log10s$neglog10.P.Value.joint_pairedtreated_placebo)
qqnorm(log10s$`neglog10.P.Value.joint_pairedtreated_2000+4000`)

# No difference between doses:
qqplot(y = log10s$neglog10.P.Value.UI4000minusplacebo,
       log10s$neglog10.P.Value.UI2000minusplacebo)

# Possibly some difference between treated paired and placebo paired:
pdf('../qqplot_pairedjoint_vs_pairedplacebo.pdf')
qqplot(y = log10s$`neglog10.P.Value.joint_pairedtreated_2000+4000`,
       x = log10s$neglog10.P.Value.joint_pairedtreated_placebo,
       ylim = c(0, 10)
       )
abline(coef = c(0, 1))#, h = 6, v =6)
abline(h = bonf, col = c('red'))
dev.off()

# The qqplot for raw p-values shows that only ~20 transcripts are more significant in the paired treated group vs
# placebo paired. This is without account for time (i.e. placebo itself as a control) and 
# that the joint group is roughly twice the size of the placebo group (as the 2000 and 4000 
# arms have been placed together here).
# However, there is a consistent deviation from around 10^-3 to lower p-values in the treated groups.


qqplot(y = final_table$`t.joint_pairedtreated_2000+4000`,
       final_table$t.joint_pairedtreated_placebo)
abline(coef = c(0, 1))
# The distribution of t-statistics seems almost identical between joint treated and placebo paired.
# 
# Together with SF7 plot (histograms of raw p-values) it would seem that drug response
# is heterogeneic but effects are small and compared to placebo non-significant.
# If their is higher variance than expecter, more samples would be required to detect differences.
#########
###########################


###########################
#########
# Normality
qqplot(as.matrix(by_group_GEx$placebo_baseline),
        as.matrix(by_group_GEx$placebo_final))

qqplot(as.matrix(by_group_GEx$baseline_2000), 
         as.matrix(by_group_GEx$final_2000))

qqplot(as.matrix(by_group_GEx$baseline_4000), 
         as.matrix(by_group_GEx$final_4000))

qqnorm(as.matrix(by_group_GEx$placebo_baseline))
qqnorm(as.matrix(by_group_GEx$placebo_baseline[1, ]))
shapiro.test(as.matrix(by_group_GEx$placebo_baseline[1, ]))
#########
###########################


###########################
#########
# Variance
# Variances for cytokines and vitamin D:
pheno <- read.csv('~/Desktop/BEST_D.DIR/mac_runs_to_upload/data.dir/BEST-D_phenotype_file_final_cytokines_and_transcripts.csv')

head(pheno)
dim(pheno)

summary(pheno$Ln_IFNgamma0)
summary(pheno$Ln_IFNgamma12)
IU_4000 <- pheno[which(pheno$arm2 == '4000_IU'), ]
summary(IU_4000$Ln_IFNgamma12)

IU_2000 <- pheno[which(pheno$arm2 == '2000_IU'), ]
placebo <- pheno[which(pheno$arm2 == 'Placebo'), ]

var.test(IU_4000$Ln_IFNgamma12, IU_4000$Ln_IFNgamma0)
var.test(IU_2000$Ln_IFNgamma12, IU_2000$Ln_IFNgamma0)
var.test(placebo$Ln_IFNgamma12, placebo$Ln_IFNgamma0)

leveneTest(IU_4000$Ln_IFNgamma12, pheno$arm)
bartlett.test(pheno$Ln_IFNgamma12, pheno$arm)

colnames(pheno)[grep('Ln', colnames(pheno))]

f <- 'Ln_IFNgamma12'
b <- 'Ln_IFNgamma0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]

f <- 'Ln_IL10_12'
b <- 'Ln_IL10_0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]

f <- 'Ln_IL6_12'
b <- 'Ln_IL6_0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]

f <- 'Ln_IL8_12'
b <- 'Ln_IL8_0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]

f <- 'Ln_TNFalpha12'
b <- 'Ln_TNFalpha0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]


f <- 'vitd12'
b <- 'vitd0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]

pheno$crp0
f <- 'crp12'
b <- 'crp0'
var.test(IU_4000[, f], IU_4000[, b])[3]
var.test(IU_2000[, f], IU_2000[, b])[3]
var.test(placebo[, f], placebo[, b])[3]
#########
###########################


###########################
#########
# TO DO: check
# Global:
# Check objects:
final_table$B.UI4000minusplacebo

by_group_GEx$placebo_baseline
x <- as.matrix(by_group_GEx$placebo_baseline)
y <- as.matrix(by_group_GEx$placebo_final)
x
y
class(x)
placebo_F <- var.test(x, y)
str(placebo_F)
placebo_F
#########

#########
# Between timepoints, but these are not independent groups:
var.test(as.matrix(by_group_GEx$placebo_baseline), 
         as.matrix(by_group_GEx$placebo_final))
# Not significant

var.test(as.matrix(by_group_GEx$baseline_2000), 
         as.matrix(by_group_GEx$final_2000))
# Significant

var.test(as.matrix(by_group_GEx$baseline_4000), 
         as.matrix(by_group_GEx$final_4000))
# Significant
#########

#########
# Between baseline groups:
var.test(as.matrix(by_group_GEx$baseline_4000),
         as.matrix(by_group_GEx$placebo_baseline))
# Not significant

var.test(as.matrix(by_group_GEx$baseline_2000),
         as.matrix(by_group_GEx$placebo_baseline))
# Significant

var.test(as.matrix(by_group_GEx$baseline_2000),
         as.matrix(by_group_GEx$baseline_4000))
# Significant
#########

#########
# Between 12 month groups:
var.test(as.matrix(by_group_GEx$final_4000),
         as.matrix(by_group_GEx$placebo_final))
# Significant

var.test(as.matrix(by_group_GEx$final_2000),
         as.matrix(by_group_GEx$placebo_final))
# Significant

var.test(as.matrix(by_group_GEx$final_2000),
         as.matrix(by_group_GEx$final_4000))
# Not significant
#########
###########################

###########################
##########
# Test variance for each probe in each group and timepoint
# Generate var for each probe at each arm and timepoint, already done in merge_all.. tables script

# Get data into one dataframe:
df_by_group_GEx <- as.data.frame(by_group_GEx)
dim(df_by_group_GEx)
head(df_by_group_GEx)[1:5, 1:5]
var_results <- data.frame('probe_ID' = rownames(df_by_group_GEx))
dim(var_results)

df <- df_by_group_GEx

# e.g. for one probe comparing placebo:
# var_test <- var.test(
#   as.numeric(df[1, grep('placebo_baseline', colnames(df))]),
#   as.numeric(df[1, grep('placebo_final', colnames(df))])
#   )
# 
# str(var_test)

# Run apply loop for F statistic
var_results$placebo_F_test <- apply(df, 1,
                                    function(x) var.test(
                                      as.numeric(x[grep('placebo_baseline', colnames(df))]),
                                      as.numeric(x[grep('placebo_final', colnames(df))])
                                      )[1]
                                    )

head(var_results)

# run for p-values
# Run apply loop for F statistic
var_results$placebo_p_value <- apply(df, 1,
                                     function(x) var.test(
                                       as.numeric(x[grep('placebo_baseline', colnames(df))]),
                                       as.numeric(x[grep('placebo_final', colnames(df))])
                                       )[3]
                                     )

head(var_results)
dim(var_results)
str(var_results)
var_results$placebo_p_value$ILMN_1651228$p.value
var_results2 <- as.data.frame(unlist(var_results))
str(var_results2)
placebo_p_value <- var_results[order(var_results$placebo_p_value), ]
qqnorm(var_results$placebo_p_value[order(var_results$placebo_p_value), ])
# TO DO:
# Calculate FDR

##########

##########
# Legacy:
# get_var_fvalue <- function(list1, df1, df2, col1) {
#   a <- as.data.frame(list1[df1])
#   a <- t(a)
#   b <- as.data.frame(list1[df2])
#   b <- t(b)
#   c <- var.test(a[, col1], b[, col1])
#   return(c)
# }

# list1 <- by_group_GEx
# df1 <- 'placebo_baseline'
# df2 <- 'placebo_final'
# head(a)[1:5, 1:5]
# colnames(a)[1:5]
# rownames(a)[1:5]
# dim(a)

# get_var_fvalue(by_group_GEx,
#                'placebo_baseline',
#                'placebo_final',
#                3
#                )

# TO DO:
# loop so that each pair of columns gets tested
##########

##########
# qqnorm for each (or -log10)

##########
###########################


###########################
sessionInfo()
q()
###########################