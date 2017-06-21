#Workflow from https://sites.tufts.edu/cbi/files/2014/11/BioConductor-and-Rf.pdf for Illumina analysis


#Load libraries:
library(limma)
library(lumi)
library(genefilter)
library(lumiHumanAll.db)
library(gplots)


#Set working directory
#setwd("")
#getwd()


#File with names of files with data
lumi_file_list <- c('SampleProbeProfile_JK-P1-P140328.txt', 'Sample_Probe_Profile_P140328-Plate_2.txt', 'Sample_Probe_Profile_Plate_3.txt', 
                    'SampleProbeProfile-Plate_4.txt', 'SampleProbeProfile_Plate-5.txt', 'SampleProbeProfile_Plate-6.txt')

lumi_control_list <- c('ControlProbeProfile_JK-P1-P140328.txt', 'Control_Probe_Profile_P140328-Plate_2.txt', 'Control_Probe_Profile-Plate_4.txt', 
  'Control_Probe_Profile_Plate_3.txt', 'ControlProbeProfile_Plate5.txt', 'ControlProbeProfile_P6.txt')

lumi_file_list
lumi_control_list

#test_files <- c('SampleProbeProfile_JK-P1-P140328.txt', 'ControlProbeProfile_JK-P1-P140328.txt', 'SampleTableControl_JK-P1-P140328.txt')
#test_batch <- lumiR.batch(test_files)

#Read raw and make lumibatch object: ##This isn't returning a lumi batch object though, needs to be csv input?
# ?lumiR

lumi_read_data <- lumiR.batch(lumi_file_list, verbose=TRUE)
class(lumi_read_data)
lumi_read_data

lumi_read_controls <-lumiR.batch(lumi_control_list)
class(lumi_read_controls)
lumi_read_controls


lumi_one_file <- lumiR('SampleProbeProfile_JK-P1-P140328.txt', verbose=TRUE)
control_one <- lumiR('ControlProbeProfile_JK-P1-P140328.txt', verbose=TRUE)

lumi_one_file
class(lumi_one_file)

add_control_data <- addControlData2lumi('ControlProbeProfile_JK-P1-P140328.txt', x.lumi=lumi_one_file)
control_data_one <- getControlData(x=control_one, type='LumiBatch')


# control_data <- addControlData2lumi(controlData=lumi_control_list, x.lumi=lumi_read_data)
control_data <- getControlData(x=lumi_read_controls, type='LumiBatch')
control_data

control_data_dataframe <- getControlData(x=lumi_read_controls)
head(control_data_dataframe)
class(control_data_dataframe)

add_control_data <- addControlData2lumi(controlData=control_data_dataframe, x.lumi=lumi_read_data)


#?getControlData
#?addControlData2lumi


lumi_batch_object <- (log2(exprs(lumi_read_data)))
class(lumi_batch_object)
head(lumi_batch_object)[1:5,1:5]
dim(lumi_batch_object)
range(lumi_batch_object[,1:575])



#Pass samples that failed QC:
FAILED_QC <- failed_QC_input_file

#Get samples IDs and index numbers from EListRaw object:
array_sample_IDs <- colnames(read_files)
to_extract <- match(FAILED_QC, array_sample_IDs)

#Check indexes match ID:

to_extract

lumi_read_data[0,to_extract]
read_files$E[0,to_extract]

#Get clean EListRaw object:
read_files_cleaned_QC <- read_files[,-to_extract]

dim(read_files_cleaned_QC)
dim(read_files)

#Explore contents of file:
head(read_files_cleaned_QC)
class(read_files_cleaned_QC)
summary(read_files_cleaned_QC)
head(read_files_cleaned_QC$source)
head(read_files_cleaned_QC$E)[1:5,1:5]
range(read_files_cleaned_QC$E)
median(read_files_cleaned_QC$E)
head(read_files_cleaned_QC$genes)
summary(read_files_cleaned_QC$other)
head(read_files_cleaned_QC$other$'Detection Pval')[1:5,1:5]
range(read_files_cleaned_QC$other$'Detection Pval')
read_files_cleaned_QC$targets


#Quantile normalize:
lumi_batch_object_normalised <- lumiN(lumi_batch_object, method="quantile")
class(lumi_batch_object_normalised)
head(lumi_batch_object_normalised)
dim(lumi_batch_object_normalised)
####


lumi_normalised_vsn <- lumiN(x.lumi=lumi_read_data, method="vsn", verbose=T)
class(lumi_batch_object_normalised)
head(lumi_batch_object_normalised)
dim(lumi_batch_object_normalised)


#Expressions:
exprs(lumi_read_data)=lumi_batch_object_normalised


#Consider only those probes that are present in three arrays:
presentCount <- detectionCall(lumi_read_data, Th=0.05)
head(presentCount)
sele.N1<- lumi_batch_object_normalised[presentCount==3,]
head(sele.N1)
selprobeList <- rownames(sele.N1)
head(selprobeList)
probeList <- rownames(lumi_batch_object_normalised)
length(probeList)


#Annotation:
if (require(lumiHumanAll.db) & require(annotate)) {
  geneSymbol <- getSYMBOL(probeList, 'lumiHumanAll.db')
  selgeneSymbol <- getSYMBOL(selprobeList, 'lumiHumanAll.db')
  geneName <- sapply(lookUp(probeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
  EntrezID  <- sapply(lookUp(probeList, 'lumiHumanAll.db', 'ENTREZID'), function(x) x[1])
  selgeneName <- sapply(lookUp(selprobeList, 'lumiHumanAll.db', 'GENENAME'), function(x) x[1])
}

lumi_batch_object_normalised_sel = lumi_batch_object_normalised[rownames(lumi_batch_object_normalised) %in% selprobeList,]
dim(lumi_batch_object_normalised_sel)
head(lumi_batch_object_normalised_sel)

class(lumi_batch_object)
class(lumi_batch_object_normalised)
class(lumi_batch_object_normalised_sel)


#Design of Experiment:


#YES, otherwise order
e.Nsel columns by colnames.


#Analysis (Model Building)
TumorType=as.factor(des$TumorType)
Patient=as.factor(des$Patient)
tp=model.matrix(~-1+TumorType+Patient)
tp
fit1 <- lmFit(e.Nsel,design=tp)
boxplot(as.data.frame(fit1$coefficients))




#Comparison between groups
?makeContrasts

#Look at slide: contr1=makeContrasts(TumorTypeDL-TumorTypeL, levels=c("TumorTypeDL", "Patient2", "Patient3"));
fit2 <- contrasts.fit(fit1,contrasts=contr1)
fit3 <- eBayes(fit2)

#Annotation + Significance




#Heatmap


#The end:
#To save R workspace with all objects to use at a later time:
save.image("")

q()
