#############################################
## Plot expression levels against genotypes
# Set-up ggplot for plotting:

genoLong <- tidyr::gather(geno2_matched, snp, genotype, -sample)
exprLong <- tidyr::gather(expr2_matched, gene, expression, -sample)
dataLong <- cbind(genoLong, exprLong['expression'])
dataLong$genotype <-as.factor(dataLong$genotype)

# Plot in ggplot2:
ggplot(dataLong, aes(genotype, expression)) + 
  geom_jitter(colour='darkgrey', position=position_jitter(width=0.25)) + 
  geom_boxplot(outlier.size=0, alpha=0.6, fill='grey') + facet_wrap(~snp) + 
  theme_bw()

#############################################

#############################################
## Fit a linear model through each SNP-gene pair:
fit <- mapply(function(e, g) lm(e ~ g), expr2_matched[-1], geno2_matched[-1], SIMPLIFY=F)

# View results for a particular gene/location:
summary(fit[[2]])

# Create confidence intervals:
ci <- sapply(fit, confint, 'g')
rownames(ci) <- c('lower', 'upper')
#############################################


#############################################
# Estimate model coefficients
betaHat <- sapply(fit, coef)[2,]
View(genoLong)
estimates <- data.frame(estimate=betaHat, t(ci), maf=maf)

# Plot coefficients/confidence intervals:
fig <- ggplot(estimates, aes(x=maf)) + geom_hline(yintercept=0, linetype='longdash') + 
  geom_errorbar(aes(ymin=lower, ymax=upper)) + geom_point(aes(y=estimate)) + theme_bw() 

fig <- fig + geom_hline(yintercept=1.5)

fig
#############################################


#############################################
## Test model coefficients for departure from 0 (ie SNPs with non-zero coefficients)
#Interpretation of p-values

#############################################


#############################################
## Repeat as above but generate a model with covariates

# Fit a linear model. For each SNP-gene pair, run a lm with covariates included:
covarFit <- mapply(function(e, g, var) lm(e ~ g + var), expr2_matched[-1], geno2_matched[-1], 
                   MoreArgs=list(as.matrix(covar[2:6])), SIMPLIFY=F)

# Beta coefficient calculation:
covarBetaHat <- sapply(covarFit, coef)[2,]

# CI (and add rownames)
covarCI <- sapply(covarFit, confint, 'g')
rownames(covarCI) <- c('lower', 'upper')

# View summary and plot
summary(covarFit[[2]])
ggplot(dataLong, aes(genotype, expression)) + 
  geom_jitter(colour='darkgrey', position=position_jitter(width=0.25)) + 
  geom_boxplot(outlier.size=0, alpha=0.6, fill='grey') + facet_wrap(~snp) + 
  theme_bw()
#############################################



#############################################
## Multi-colinearity analysis
# How to deal with covariates that are linearly correlated
# Check correlation between expalantory variables
# Check variance inflation factor (VIF)
# Choose only one from each group of correlated variables (those of interest)


#############################################


#############################################
## Analysis of real datasets where there are sources of variation that are not known
# Assumption: for non-genetic causes we don't need to know the sources necessarily
# Non-genetic variation may be the strongest component of gene expression variation
# Use principle components analysis (PCA): transform data into set of linearly uncorrelated variables (the PCs)
# Order PCs by the proportion of variance they explain
# Reduce dimensionality of a dataset by considering only the first k PCs.

#Considerations:
# Number of PCs to include is arbitrary/dataset dependent
# Need to test for PCs that do correlate with genotype
# Double check data formatting
# Data must be centred and scaled prior to PCA (ie transpose dataset)

# Compute PCs of gene expression data
# ?prcomp

# Read in the data:
# expr <- readr::read_tsv(file('/data/monocytes/expression/ifn_expression.tab.gz'))
# geno <- readr::read_tsv('/genotypes/genotypes/genotypes.tab')

# Compute the PCs, first transpose the expression values and ignore the first column, then run the PCs:
pca <- prcomp(t(expr2_matched[-1]), center=TRUE, scale=TRUE)

# Obtain values for all PCs output:
pc <- pca$x

# Explore dimensions and plot first 10 or so components:
dim(pc)
dim(expr2_matched)
str(pca)
plot(pca)

# Check how much of the variance in gene expression values is explained by the first x PCs:
sum(pca$sdev[1:10]^2)/length(expr2_matched[-1])

#############################################


#############################################
# Subset SNPs of interest:
rs_interest <- subset(x=geno2_matched, id == 'rs4077515')
probe_interest <- which(expr2_matched$Probe == '3710685')

# Create dataframe with variables of interest and PCs:
eQTL_interest <- data.frame(pc[,1:10], t(probe_interest[,-1]), t(rs_interest[,-1]))
View(eQTL_interest)
colnames(eQTL_interest)[11] = 'probe_interest'
colnames(eQTL_interest)[12] = 'rs_interest'
View(eQTL_interest)

#Fit linear models with and without PCs:
fit_mono_no_pc <- lm(probe_interest ~ rs_interest, data = eQTL_interest)
fit_mono_pc <- lm(probe_interest ~ ., data = eQTL_interest)

# Look at resuts of lm: 
summary(fit_mono_pc)

# Calculate beta coefficients (not required here) and CIs:
BetaHat_interest_pc <- coef(fit_mono_pc)
ci_pc <- confint(fit_mono_pc)


# rownames(ci_pc) <- c('lower', 'upper')

# Extract values of interest for PCs only and substract from effect of SNP of interest:
corrected <- eQTL_interest$probe_interest - rowSums(coef(fit_mono_pc)[(2:11)]*eQTL_interest[,11:12])
corrected <- data.frame(expression=corrected, genotype=factor(eQTL_interest$rs_interest))

# Boxplots of expression values after substracting effect of first 10 PCs:
ggplot(corrected, aes(genotype, expression)) + 
  geom_jitter(colour='darkgrey', position=position_jitter(width=0.25)) + 
  geom_boxplot(outlier.size=0, alpha=0.6, fill='grey') + 
  theme_bw()

#############################################