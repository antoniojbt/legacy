###################
# Produce pretty tables
# See:
# http://www.strengejacke.de/sjPlot/sjt.lm/
# http://www.strengejacke.de/sjPlot/sjtbasics/
# https://github.com/strengejacke/sjPlot

# Stargazer:
# https://cran.r-project.org/web/packages/stargazer/vignettes/stargazer.pdf
###################

###################
# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('sjPlot')
# biocLite('sjmisc')
# biocLite("stargazer")

library(sjPlot)
library(sjmisc)
library(stargazer)
###################

###################
# With sj libraries, HTML output
# Fit models, generate the table and specify labels:
summary(lm_Ln_IFNgamma12)
sjt.lm(lm_Ln_IFNgamma12,
       lm_Ln_IL10_12,
       lm_Ln_IL6_12,
       lm_Ln_IL8_12,
       lm_Ln_TNFalpha12,
       depvar.labels = c(  'IFNg 12 months',
                           'IL10 12 months',
                           'IL6 12 months',
                           'IL8 12 months',
                           'TNFa 12 months'),
       file = 'lm_table.html',
       group.pred = FALSE,
       pred.labels = c() # For predictor variables
       )

# sjt.lm(use.viewer = FALSE) # Suppress output
# separate.ci.col = FALSE) # ci in same cell as estimates
# pred.labels = c() # For predictor variables
###################


###################
# With stargazer, LATEX, HTML outputs
# Needs more work to specify output, more control though
# Example data:
stargazer(attitude)

# Run models, pass to stargazer:
stargazer(lm_Ln_IFNgamma12,
          lm_Ln_IL10_12,
          lm_Ln_IL6_12,
          lm_Ln_IL8_12,
          lm_Ln_TNFalpha12,
          dep.var.labels = c('IFNg 12 months',
                           'IL10 12 months',
                           'IL6 12 months',
                           'IL8 12 months',
                           'TNFa 12 months'),
          summary = TRUE,
          align = TRUE,
          out = 'lm_table.html')

###################

###################
# ANOVA tables with stargazer

# Run models, pass to stargazer:
stargazer(anova_Ln_IFNgamma12$anova.table,
          anova_Ln_IL10_12$anova.table,
          anova_Ln_IL6_12$anova.table,
          anova_Ln_IL8_12$anova.table,
          anova_Ln_TNFalpha12$anova.table,
          dep.var.labels = c('IFNg 12 months',
                             'IL10 12 months',
                             'IL6 12 months',
                             'IL8 12 months',
                             'TNFa 12 months'),
          # summary = TRUE,
          align = TRUE,
          out = 'anova_table.html')

stargazer(anova_Ln_IFNgamma12,
          anova_Ln_IL10_12,
          anova_Ln_IL6_12,
          anova_Ln_IL8_12,
          anova_Ln_TNFalpha12,
          dep.var.labels = c('IFNg 12 months',
                             'IL10 12 months',
                             'IL6 12 months',
                             'IL8 12 months',
                             'TNFa 12 months'),
          # summary = TRUE,
          align = TRUE,
          out = 'anova_table.html')
###################