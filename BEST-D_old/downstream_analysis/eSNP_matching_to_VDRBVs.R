getwd()
setwd('Desktop/BEST_D_03_MAR.DIR/BEST-D_22_Apr_2016/')


eSNPs_BESTD <- read.csv('eSNPs_BEST-D.txt', header = FALSE)
head(eSNPs_BESTD)
class(eSNPs_BESTD)
dim(eSNPs_BESTD)

VDR_rBVs <- read.csv('VDR-recurrentBVs.txt', header = FALSE)
head(VDR_rBVs)
class(VDR_rBVs)
dim(VDR_rBVs)

VDR_singleBVs <- read.csv('VDR-singleBVs.txt', header = FALSE)
head(VDR_singleBVs)
class(VDR_singleBVs)
dim(VDR_singleBVs)

VDR_singleBVs_full <- read.csv('VDR-singleBVs_full.csv', header = TRUE, skip = 21)
head(VDR_singleBVs_full)
class(VDR_singleBVs_full)
dim(VDR_singleBVs_full)
View(VDR_singleBVs_full)

which(as.character(VDR_rBVs$V1) %in% as.character(eSNPs_BESTD$V1))
which(!is.na(match(x = eSNPs_BESTD$V1, VDR_rBVs$V1)))

matching_singleBVs_BESTD <- which(as.character(eSNPs_BESTD$V1) %in% as.character(VDR_singleBVs_full$ID))
length(matching_singleBVs_BESTD)

matching_singleBVs_BESTD <- as.character(eSNPs_BESTD[matching_singleBVs_BESTD, ])
matching_singleBVs_BESTD <- as.data.frame(matching_singleBVs_BESTD)
names(matching_singleBVs_BESTD)[1] <- 'ID'
matching_singleBVs_BESTD

singleBVs_BESTD <- merge(matching_singleBVs_BESTD, VDR_singleBVs_full)
head(singleBVs_BESTD)
dim(singleBVs_BESTD)
View(singleBVs_BESTD)


grep(pattern = 'POU2F2', ignore.case = TRUE, perl = TRUE, x = singleBVs_BESTD)


