#########
# Exploration of variance and distributions

# F tests, qqplotting, etc.

# Re-load a previous R session, data and objects:
load('R_session_saved_image_diff_expression_merge_tables.RData', verbose=T)

#########

#########
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
# However, there is a consistent deviation from around 10-3 to lower p-values in the treated groups.

# TO DO: continue from here
qqplot(y = final_table$`t.joint_pairedtreated_2000+4000`,
       final_table$t.joint_pairedtreated_placebo)
abline(coef = c(0, 1))
# t-statistics 

# Together with SF7 plot (histograms of raw p-values) it would seem that drug response
# is heterogeneic but effects are small and compared to placebo non-significant.
#########


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


#########
# Variance
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

#########
sessionInfo()
q()
#########