library(dplyr)
library(data.table)
library(plyr)

data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE6475_RAW/RAW_GSE6475.csv")
data[1:5,1:10]
# Renaming columns names
colnames(data)[1] <- "Gene_Entrez_ID"
data[1:5,1:10]

is.na(data)
na.fail(data)

sum(is.na(data))
rowSums(is.na(data))
colSums(is.na(data))

#data <- read.csv("RAW_Final.csv")
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
dat3 <- dat2[,list(GSM148725.CEL.gz=mean(GSM148725.CEL.gz),
                   GSM148748.CEL.gz=mean(GSM148748.CEL.gz), 
                   GSM148762.CEL.gz=mean(GSM148762.CEL.gz),
                   GSM148763.CEL.gz=mean(GSM148763.CEL.gz),
                   GSM148764.CEL.gz=mean(GSM148764.CEL.gz),
                   
                   GSM148765.CEL.gz=mean(GSM148765.CEL.gz),
                   GSM148766.CEL.gz=mean(GSM148766.CEL.gz),
                   GSM148767.CEL.gz=mean(GSM148767.CEL.gz),
                   GSM148768.CEL.gz=mean(GSM148768.CEL.gz),
                   GSM148769.CEL.gz=mean(GSM148769.CEL.gz),
                   
                   GSM148770.CEL.gz=mean(GSM148770.CEL.gz),
                   GSM148771.CEL.gz=mean(GSM148771.CEL.gz),
                   GSM148887.CEL.gz=mean(GSM148887.CEL.gz),
                   GSM148888.CEL.gz=mean(GSM148888.CEL.gz),
                   GSM148889.CEL.gz=mean(GSM148889.CEL.gz),
                   
                   GSM148890.CEL.gz=mean(GSM148890.CEL.gz),
                   GSM148892.CEL.gz=mean(GSM148892.CEL.gz),
                   GSM148894.CEL.gz=mean(GSM148894.CEL.gz)),list(Gene_Entrez_ID)]

dat4 <- with(dat3, aggregate(cbind(GSM148725.CEL.gz,
                                   GSM148748.CEL.gz,
                                   GSM148762.CEL.gz,
                                   GSM148763.CEL.gz,
                                   GSM148764.CEL.gz,
                                   
                                   GSM148765.CEL.gz,
                                   GSM148766.CEL.gz,
                                   GSM148767.CEL.gz,
                                   GSM148768.CEL.gz,
                                   GSM148769.CEL.gz,
                                   
                                   GSM148770.CEL.gz,
                                   GSM148771.CEL.gz,
                                   GSM148887.CEL.gz,
                                   GSM148888.CEL.gz,
                                   GSM148889.CEL.gz,
                                   
                                   GSM148890.CEL.gz,
                                   GSM148892.CEL.gz,
                                   GSM148894.CEL.gz), list(Gene_Entrez_ID), FUN=mean))

colnames(dat4)<-colnames(dat1)[c(1,2,3,4,5,6,7,8,9,10,
                                 11,12,13,14,15,16,17,18,19)]

dat5<-ddply(dat1,.(Gene_Entrez_ID),colwise(mean, c("GSM148725.CEL.gz",
                                                   "GSM148748.CEL.gz",
                                                   "GSM148762.CEL.gz",
                                                   "GSM148763.CEL.gz",
                                                   "GSM148764.CEL.gz",
                                                   
                                                   "GSM148765.CEL.gz",
                                                   "GSM148766.CEL.gz",
                                                   "GSM148767.CEL.gz",
                                                   "GSM148768.CEL.gz",
                                                   "GSM148769.CEL.gz",
                                                   
                                                   "GSM148770.CEL.gz",
                                                   "GSM148771.CEL.gz",
                                                   "GSM148887.CEL.gz",
                                                   "GSM148888.CEL.gz",
                                                   "GSM148889.CEL.gz",
                                                   
                                                   "GSM148890.CEL.gz",
                                                   "GSM148892.CEL.gz",
                                                   "GSM148894.CEL.gz")))

##make sure no duplicated rows
duplicated(dat5$Gene_Entrez_ID)
sum(duplicated(dat5$Gene_Entrez_ID))

write.csv(dat5, "CLEAN_DATA_GSE6475.csv")
