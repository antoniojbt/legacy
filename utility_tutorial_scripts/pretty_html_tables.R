# Produce pretty html tables
# See:
# http://www.strengejacke.de/sjPlot/sjt.lm/
# http://www.strengejacke.de/sjPlot/sjtbasics/
# https://github.com/strengejacke/sjPlot

# source("https://bioconductor.org/biocLite.R")
# biocLite()
# biocLite('sjPlot')
# biocLite('sjmisc')

library(sjPlot)
library(sjmisc)

# sample data
data()

# fit first model
fit1 <- lm(barthtot ~ c160age + c12hour + c161sex + c172code, data = efc)
# fit second model 
fit2 <- lm(neg_c_7 ~ c160age + c12hour + c161sex + c172code, data = efc)
# Note that both models share the same predictors and only differ 
# in their dependent variable. See examples of stepwise models 
# later..

# Generate the table:
sjt.lm(fit1, fit2)

# You can specify the ‘model’ label via depvar.labels parameter:
sjt.lm(fit1, fit2, depvar.labels = c("Barthel-Index", "Negative Impact"))

