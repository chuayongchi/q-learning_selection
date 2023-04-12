library(dplyr)
library(data.table)
library(plyr)

data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/RAW_GSE53795.csv")
data[1:5,1:10]
# Renaming columns names
colnames(data)[1] <- "Gene_Entrez_ID"
data[1:5,1:10]

is.na(data)
na.fail(data)

sum(is.na(data))
rowSums(is.na(data))
colSums(is.na(data))

data_comp <- data[complete.cases(data),]
na.fail(data_comp)

##check duplicated rows
duplicated(data_comp$Gene_Entrez_ID)
sum(duplicated(data_comp$Gene_Entrez_ID))

## obtain the colnames for data_comp to take the average of the duplicate rows
colnames(data_comp)


#average the duplicate row
dat1 <- data_comp
dat2 <-data.table(dat1)
dat3 <- dat2[,list(GSM1300911_CRF_26_LST.CEL.gz=mean(GSM1300911_CRF_26_LST.CEL.gz),
                   GSM1300912_CRF_26_NST.CEL.gz=mean(GSM1300912_CRF_26_NST.CEL.gz), 
                   GSM1300913_CRF_27_LST.CEL.gz=mean(GSM1300913_CRF_27_LST.CEL.gz),
                   GSM1300914_CRF_27_NST.CEL.gz=mean(GSM1300914_CRF_27_NST.CEL.gz),
                   
                   GSM1300915_CRF_28_LST.CEL.gz=mean(GSM1300915_CRF_28_LST.CEL.gz),
                   GSM1300916_CRF_28_NST.CEL.gz=mean(GSM1300916_CRF_28_NST.CEL.gz),
                   GSM1300917_CRF_29_LST.CEL.gz=mean(GSM1300917_CRF_29_LST.CEL.gz),
                   GSM1300918_CRF_29_NST.CEL.gz=mean(GSM1300918_CRF_29_NST.CEL.gz),
                   
                   GSM1300919_CRF_30_LST.CEL.gz=mean(GSM1300919_CRF_30_LST.CEL.gz),
                   GSM1300920_CRF_30_NST.CEL.gz=mean(GSM1300920_CRF_30_NST.CEL.gz),
                   GSM1300921_CRF_31_LST.CEL.gz=mean(GSM1300921_CRF_31_LST.CEL.gz),
                   GSM1300922_CRF_31_NST.CEL.gz=mean(GSM1300922_CRF_31_NST.CEL.gz),
                   
                   GSM1300923_CRF_32_LST.CEL.gz=mean(GSM1300923_CRF_32_LST.CEL.gz),
                   GSM1300924_CRF_32_NST.CEL.gz=mean(GSM1300924_CRF_32_NST.CEL.gz),
                   GSM1300925_CRF_33_LST.CEL.gz=mean(GSM1300925_CRF_33_LST.CEL.gz),
                   GSM1300926_CRF_33_NST.CEL.gz=mean(GSM1300926_CRF_33_NST.CEL.gz),
                   
                   GSM1300927_CRF_34_LST.CEL.gz=mean(GSM1300927_CRF_34_LST.CEL.gz),
                   GSM1300928_CRF_34_NST.CEL.gz=mean(GSM1300928_CRF_34_NST.CEL.gz),
                   GSM1300929_CRF_35_LST.CEL.gz=mean(GSM1300929_CRF_35_LST.CEL.gz),
                   GSM1300930_CRF_35_NST.CEL.gz=mean(GSM1300930_CRF_35_NST.CEL.gz),
                   
                   GSM1300931_CRF_36_LST.CEL.gz=mean(GSM1300931_CRF_36_LST.CEL.gz),
                   GSM1300932_CRF_36_NST.CEL.gz=mean(GSM1300932_CRF_36_NST.CEL.gz),
                   GSM1300933_CRF_38_LST.CEL.gz=mean(GSM1300933_CRF_38_LST.CEL.gz),
                   GSM1300934_CRF_38_NST.CEL.gz=mean(GSM1300934_CRF_38_NST.CEL.gz)),list(Gene_Entrez_ID)]

dat4 <- with(dat3, aggregate(cbind(GSM1300911_CRF_26_LST.CEL.gz,
                                   GSM1300912_CRF_26_NST.CEL.gz,
                                   GSM1300913_CRF_27_LST.CEL.gz,
                                   GSM1300914_CRF_27_NST.CEL.gz,
                                   
                                   GSM1300915_CRF_28_LST.CEL.gz,
                                   GSM1300916_CRF_28_NST.CEL.gz,
                                   GSM1300917_CRF_29_LST.CEL.gz,
                                   GSM1300918_CRF_29_NST.CEL.gz,
                                   
                                   GSM1300919_CRF_30_LST.CEL.gz,
                                   GSM1300920_CRF_30_NST.CEL.gz,
                                   GSM1300921_CRF_31_LST.CEL.gz,
                                   GSM1300922_CRF_31_NST.CEL.gz,
                                   
                                   GSM1300923_CRF_32_LST.CEL.gz,
                                   GSM1300924_CRF_32_NST.CEL.gz,
                                   GSM1300925_CRF_33_LST.CEL.gz,
                                   GSM1300926_CRF_33_NST.CEL.gz,
                                   
                                   GSM1300927_CRF_34_LST.CEL.gz,
                                   GSM1300928_CRF_34_NST.CEL.gz,
                                   GSM1300929_CRF_35_LST.CEL.gz,
                                   GSM1300930_CRF_35_NST.CEL.gz,
                                   
                                   GSM1300931_CRF_36_LST.CEL.gz,
                                   GSM1300932_CRF_36_NST.CEL.gz,
                                   GSM1300933_CRF_38_LST.CEL.gz,
                                   GSM1300934_CRF_38_NST.CEL.gz), list(Gene_Entrez_ID), FUN=mean))
colnames(dat4)<-colnames(dat1)[c(1,2,3,4,5,6,7,8,9,10,
                                 11,12,13,14,15,16,17,18,19,20,
                                 21,22,23,24,25)]

dat5<-ddply(dat1,.(Gene_Entrez_ID),colwise(mean, c("GSM1300911_CRF_26_LST.CEL.gz",
                                                   "GSM1300912_CRF_26_NST.CEL.gz",
                                                   "GSM1300913_CRF_27_LST.CEL.gz",
                                                   "GSM1300914_CRF_27_NST.CEL.gz",
                                                   
                                                   "GSM1300915_CRF_28_LST.CEL.gz",
                                                   "GSM1300916_CRF_28_NST.CEL.gz",
                                                   "GSM1300917_CRF_29_LST.CEL.gz",
                                                   "GSM1300918_CRF_29_NST.CEL.gz",
                                                   
                                                   "GSM1300919_CRF_30_LST.CEL.gz",
                                                   "GSM1300920_CRF_30_NST.CEL.gz",
                                                   "GSM1300921_CRF_31_LST.CEL.gz",
                                                   "GSM1300922_CRF_31_NST.CEL.gz",
                                                   
                                                   "GSM1300923_CRF_32_LST.CEL.gz",
                                                   "GSM1300924_CRF_32_NST.CEL.gz",
                                                   "GSM1300925_CRF_33_LST.CEL.gz",
                                                   "GSM1300926_CRF_33_NST.CEL.gz",
                                                   
                                                   "GSM1300927_CRF_34_LST.CEL.gz",
                                                   "GSM1300928_CRF_34_NST.CEL.gz",
                                                   "GSM1300929_CRF_35_LST.CEL.gz",
                                                   "GSM1300930_CRF_35_NST.CEL.gz",
                                                   
                                                   "GSM1300931_CRF_36_LST.CEL.gz",
                                                   "GSM1300932_CRF_36_NST.CEL.gz",
                                                   "GSM1300933_CRF_38_LST.CEL.gz",
                                                   "GSM1300934_CRF_38_NST.CEL.gz")))

write.csv(dat5, "CLEAN_DATA_GSE53795.csv")
