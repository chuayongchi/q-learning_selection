library(dplyr)
library(data.table)
library(plyr)

data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE108110_RAW/RAW_GSE108110.csv")
View(data)

# Rename columns names
colnames(data)[1] <- "Gene_Entrez_ID"
View(data)

#check missing value
is.na(data)
na.fail(data)
sum(is.na(data))

data_comp <- data[complete.cases(data),]
na.fail(data_comp)

#check missing value
is.na(data_comp)
na.fail(data_comp)
sum(is.na(data_comp))

#average the duplicate row
dat1 <- data_comp
dat2 <-data.table(dat1)
dat3 <- dat2[,list(GSM2889960_9016_NPS02_L1.CEL.gz=mean(GSM2889960_9016_NPS02_L1.CEL.gz),
                   GSM2889961_9017_NPS01_L1.CEL.gz=mean(GSM2889961_9017_NPS01_L1.CEL.gz), 
                   GSM2889962_9019_NPS03_L1.CEL.gz=mean(GSM2889962_9019_NPS03_L1.CEL.gz),
                   GSM2889963_9020_NPS04_L1.CEL.gz=mean(GSM2889963_9020_NPS04_L1.CEL.gz),
                   GSM2889964_9021_NPS05_L1.CEL.gz=mean(GSM2889964_9021_NPS05_L1.CEL.gz),
                   GSM2889965_9025_NPS06_L1.CEL.gz=mean(GSM2889965_9025_NPS06_L1.CEL.gz),
                   GSM2889966_9026_NPS07_L1.CEL.gz=mean(GSM2889966_9026_NPS07_L1.CEL.gz),
                   GSM2889967_9028_NPS08_L1.CEL.gz=mean(GSM2889967_9028_NPS08_L1.CEL.gz),
                   GSM2889968_9032_NPS09_L1.CEL.gz=mean(GSM2889968_9032_NPS09_L1.CEL.gz),
                   GSM2889969_9016_NPS02_L3.CEL.gz=mean(GSM2889969_9016_NPS02_L3.CEL.gz),
                   GSM2889970_9017_NPS01_L3.CEL.gz=mean(GSM2889970_9017_NPS01_L3.CEL.gz),
                   GSM2889971_9019_NPS03_L3.CEL.gz=mean(GSM2889971_9019_NPS03_L3.CEL.gz),
                   GSM2889972_9020_NPS04_L3.CEL.gz=mean(GSM2889972_9020_NPS04_L3.CEL.gz),
                   GSM2889973_9021_NPS05_L3.CEL.gz=mean(GSM2889973_9021_NPS05_L3.CEL.gz),
                   GSM2889974_9025_NPS06_L3.CEL.gz=mean(GSM2889974_9025_NPS06_L3.CEL.gz),
                   GSM2889975_9026_NPS07_L3.CEL.gz=mean(GSM2889975_9026_NPS07_L3.CEL.gz),
                   GSM2889976_9028_NPS08_L3.CEL.gz=mean(GSM2889976_9028_NPS08_L3.CEL.gz),
                   GSM2889977_9032_NPS09_L3.CEL.gz=mean(GSM2889977_9032_NPS09_L3.CEL.gz),
                   GSM2889978_9016_NPS02_L.NI.CEL.gz=mean(GSM2889978_9016_NPS02_L.NI.CEL.gz),
                   GSM2889979_9017_NPS01_L.NI.CEL.gz=mean(GSM2889979_9017_NPS01_L.NI.CEL.gz),
                   GSM2889980_9019_NPS03_L.NI.CEL.gz=mean(GSM2889980_9019_NPS03_L.NI.CEL.gz),
                   GSM2889981_9020_NPS04_L.NI.CEL.gz=mean(GSM2889981_9020_NPS04_L.NI.CEL.gz),
                   GSM2889982_9021_NPS05_L.NI.CEL.gz=mean(GSM2889982_9021_NPS05_L.NI.CEL.gz),
                   GSM2889983_9025_NPS06_L.NI.CEL.gz=mean(GSM2889983_9025_NPS06_L.NI.CEL.gz),
                   GSM2889984_9026_NPS07_L.NI.CEL.gz=mean(GSM2889984_9026_NPS07_L.NI.CEL.gz),
                   GSM2889985_9028_NPS08_L.NI.CEL.gz=mean(GSM2889985_9028_NPS08_L.NI.CEL.gz),
                   GSM2889986_9032_NPS09_L.NI.CEL.gz=mean(GSM2889986_9032_NPS09_L.NI.CEL.gz),
                   GSM2889987_9018_PS01_L1.CEL.gz=mean(GSM2889987_9018_PS01_L1.CEL.gz),
                   GSM2889988_9022_PS02_L1.CEL.gz=mean(GSM2889988_9022_PS02_L1.CEL.gz),
                   GSM2889989_9023_PS03_L1.CEL.gz=mean(GSM2889989_9023_PS03_L1.CEL.gz),
                   GSM2889990_9024_PS04_L1.CEL.gz=mean(GSM2889990_9024_PS04_L1.CEL.gz),
                   GSM2889991_9027_PS05_L1.CEL.gz=mean(GSM2889991_9027_PS05_L1.CEL.gz),
                   GSM2889992_9029_PS06_L1.CEL.gz=mean(GSM2889992_9029_PS06_L1.CEL.gz),
                   GSM2889993_9033_PS08_L1.CEL.gz=mean(GSM2889993_9033_PS08_L1.CEL.gz),
                   GSM2889994_9034_PS09_L1.CEL.gz=mean(GSM2889994_9034_PS09_L1.CEL.gz),
                   GSM2889995_9035_PS10_L1.CEL.gz=mean(GSM2889995_9035_PS10_L1.CEL.gz),
                   GSM2889996_9018_PS01_L3.CEL.gz=mean(GSM2889996_9018_PS01_L3.CEL.gz),
                   GSM2889997_9022_PS02_L3.CEL.gz=mean(GSM2889997_9022_PS02_L3.CEL.gz),
                   GSM2889998_9023_PS03_L3.CEL.gz=mean(GSM2889998_9023_PS03_L3.CEL.gz),
                   GSM2889999_9024_PS04_L3.CEL.gz=mean(GSM2889999_9024_PS04_L3.CEL.gz),
                   GSM2890000_9027_PS05_L3.CEL.gz=mean(GSM2890000_9027_PS05_L3.CEL.gz),
                   GSM2890001_9029_PS06_L3.CEL.gz=mean(GSM2890001_9029_PS06_L3.CEL.gz),
                   GSM2890002_9033_PS08_L3.CEL.gz=mean(GSM2890002_9033_PS08_L3.CEL.gz),
                   GSM2890003_9034_PS09_L3.CEL.gz=mean(GSM2890003_9034_PS09_L3.CEL.gz),
                   GSM2890004_9035_PS10_L3.CEL.gz=mean(GSM2890004_9035_PS10_L3.CEL.gz),
                   GSM2890005_9018_PS01_L.NI.CEL.gz=mean(GSM2890005_9018_PS01_L.NI.CEL.gz),
                   GSM2890006_9022_PS02_L.NI.CEL.gz=mean(GSM2890006_9022_PS02_L.NI.CEL.gz),
                   GSM2890007_9023_PS03_L.NI.CEL.gz=mean(GSM2890007_9023_PS03_L.NI.CEL.gz),
                   GSM2890008_9024_PS04_L.NI.CEL.gz=mean(GSM2890008_9024_PS04_L.NI.CEL.gz),
                   GSM2890009_9027_PS05_L.NI.CEL.gz=mean(GSM2890009_9027_PS05_L.NI.CEL.gz),
                   GSM2890010_9029_PS06_L.NI.CEL.gz=mean(GSM2890010_9029_PS06_L.NI.CEL.gz),
                   GSM2890011_9033_PS08_L.NI.CEL.gz=mean(GSM2890011_9033_PS08_L.NI.CEL.gz),
                   GSM2890012_9034_PS09_L.NI.CEL.gz=mean(GSM2890012_9034_PS09_L.NI.CEL.gz),
                   GSM2890013_9035_PS10_L.NI.CEL.gz=mean(GSM2890013_9035_PS10_L.NI.CEL.gz)),list(Gene_Entrez_ID)]
dat4 <- with(dat3, aggregate(cbind(GSM2889960_9016_NPS02_L1.CEL.gz, GSM2889961_9017_NPS01_L1.CEL.gz, GSM2889962_9019_NPS03_L1.CEL.gz, 
                                   GSM2889963_9020_NPS04_L1.CEL.gz, GSM2889964_9021_NPS05_L1.CEL.gz, 
                                   GSM2889965_9025_NPS06_L1.CEL.gz, GSM2889966_9026_NPS07_L1.CEL.gz, 
                                   GSM2889967_9028_NPS08_L1.CEL.gz, GSM2889968_9032_NPS09_L1.CEL.gz, 
                                   GSM2889969_9016_NPS02_L3.CEL.gz, GSM2889970_9017_NPS01_L3.CEL.gz, 
                                   GSM2889971_9019_NPS03_L3.CEL.gz, GSM2889972_9020_NPS04_L3.CEL.gz, 
                                   GSM2889973_9021_NPS05_L3.CEL.gz, GSM2889974_9025_NPS06_L3.CEL.gz, 
                                   GSM2889975_9026_NPS07_L3.CEL.gz, GSM2889976_9028_NPS08_L3.CEL.gz, 
                                   GSM2889977_9032_NPS09_L3.CEL.gz, GSM2889978_9016_NPS02_L.NI.CEL.gz, 
                                   GSM2889979_9017_NPS01_L.NI.CEL.gz, GSM2889980_9019_NPS03_L.NI.CEL.gz, 
                                   GSM2889981_9020_NPS04_L.NI.CEL.gz, GSM2889982_9021_NPS05_L.NI.CEL.gz, 
                                   GSM2889983_9025_NPS06_L.NI.CEL.gz, GSM2889984_9026_NPS07_L.NI.CEL.gz, 
                                   GSM2889985_9028_NPS08_L.NI.CEL.gz, GSM2889986_9032_NPS09_L.NI.CEL.gz, 
                                   GSM2889987_9018_PS01_L1.CEL.gz, GSM2889988_9022_PS02_L1.CEL.gz, 
                                   GSM2889989_9023_PS03_L1.CEL.gz, GSM2889990_9024_PS04_L1.CEL.gz, 
                                   GSM2889991_9027_PS05_L1.CEL.gz, GSM2889992_9029_PS06_L1.CEL.gz, 
                                   GSM2889993_9033_PS08_L1.CEL.gz, GSM2889994_9034_PS09_L1.CEL.gz, 
                                   GSM2889995_9035_PS10_L1.CEL.gz, GSM2889996_9018_PS01_L3.CEL.gz, 
                                   GSM2889997_9022_PS02_L3.CEL.gz, GSM2889998_9023_PS03_L3.CEL.gz, 
                                   GSM2889999_9024_PS04_L3.CEL.gz, GSM2890000_9027_PS05_L3.CEL.gz, 
                                   GSM2890001_9029_PS06_L3.CEL.gz, GSM2890002_9033_PS08_L3.CEL.gz, 
                                   GSM2890003_9034_PS09_L3.CEL.gz, GSM2890004_9035_PS10_L3.CEL.gz, 
                                   GSM2890005_9018_PS01_L.NI.CEL.gz, GSM2890006_9022_PS02_L.NI.CEL.gz, 
                                   GSM2890007_9023_PS03_L.NI.CEL.gz, GSM2890008_9024_PS04_L.NI.CEL.gz, 
                                   GSM2890009_9027_PS05_L.NI.CEL.gz, GSM2890010_9029_PS06_L.NI.CEL.gz, 
                                   GSM2890011_9033_PS08_L.NI.CEL.gz, GSM2890012_9034_PS09_L.NI.CEL.gz, 
                                   GSM2890013_9035_PS10_L.NI.CEL.gz), list(Gene_Entrez_ID), FUN=mean))
colnames(dat4)<-colnames(dat1)[c(1,2,3,4,5,6,7,8,9,10,
                                 11,12,13,14,15,16,17,18,19,20,
                                 21,22,23,24,25,26,27,28,29,30,
                                 31,32,33,34,35,36,37,38,39,40,
                                 41,42,43,44,45,46,47,48,49,50,
                                 51,52,53,54,55)]

dat5<-ddply(dat1,.(Gene_Entrez_ID),colwise(mean, c("GSM2889960_9016_NPS02_L1.CEL.gz",
                                                   "GSM2889961_9017_NPS01_L1.CEL.gz","GSM2889962_9019_NPS03_L1.CEL.gz",
                                                   "GSM2889963_9020_NPS04_L1.CEL.gz","GSM2889964_9021_NPS05_L1.CEL.gz",
                                                   "GSM2889965_9025_NPS06_L1.CEL.gz","GSM2889966_9026_NPS07_L1.CEL.gz",
                                                   "GSM2889967_9028_NPS08_L1.CEL.gz","GSM2889968_9032_NPS09_L1.CEL.gz",
                                                   "GSM2889969_9016_NPS02_L3.CEL.gz","GSM2889970_9017_NPS01_L3.CEL.gz",
                                                   "GSM2889971_9019_NPS03_L3.CEL.gz","GSM2889972_9020_NPS04_L3.CEL.gz",
                                                   "GSM2889973_9021_NPS05_L3.CEL.gz","GSM2889974_9025_NPS06_L3.CEL.gz",
                                                   "GSM2889975_9026_NPS07_L3.CEL.gz","GSM2889976_9028_NPS08_L3.CEL.gz",
                                                   "GSM2889977_9032_NPS09_L3.CEL.gz","GSM2889978_9016_NPS02_L.NI.CEL.gz",
                                                   "GSM2889979_9017_NPS01_L.NI.CEL.gz","GSM2889980_9019_NPS03_L.NI.CEL.gz",
                                                   "GSM2889981_9020_NPS04_L.NI.CEL.gz","GSM2889982_9021_NPS05_L.NI.CEL.gz",
                                                   "GSM2889983_9025_NPS06_L.NI.CEL.gz","GSM2889984_9026_NPS07_L.NI.CEL.gz",
                                                   "GSM2889985_9028_NPS08_L.NI.CEL.gz","GSM2889986_9032_NPS09_L.NI.CEL.gz",
                                                   "GSM2889987_9018_PS01_L1.CEL.gz","GSM2889988_9022_PS02_L1.CEL.gz",
                                                   "GSM2889989_9023_PS03_L1.CEL.gz","GSM2889990_9024_PS04_L1.CEL.gz",
                                                   "GSM2889991_9027_PS05_L1.CEL.gz","GSM2889992_9029_PS06_L1.CEL.gz",
                                                   "GSM2889993_9033_PS08_L1.CEL.gz","GSM2889994_9034_PS09_L1.CEL.gz",
                                                   "GSM2889995_9035_PS10_L1.CEL.gz","GSM2889996_9018_PS01_L3.CEL.gz",
                                                   "GSM2889997_9022_PS02_L3.CEL.gz","GSM2889998_9023_PS03_L3.CEL.gz",
                                                   "GSM2889999_9024_PS04_L3.CEL.gz","GSM2890000_9027_PS05_L3.CEL.gz",
                                                   "GSM2890001_9029_PS06_L3.CEL.gz","GSM2890002_9033_PS08_L3.CEL.gz",
                                                   "GSM2890003_9034_PS09_L3.CEL.gz","GSM2890004_9035_PS10_L3.CEL.gz",
                                                   "GSM2890005_9018_PS01_L.NI.CEL.gz","GSM2890006_9022_PS02_L.NI.CEL.gz",
                                                   "GSM2890007_9023_PS03_L.NI.CEL.gz","GSM2890008_9024_PS04_L.NI.CEL.gz",
                                                   "GSM2890009_9027_PS05_L.NI.CEL.gz","GSM2890010_9029_PS06_L.NI.CEL.gz",
                                                   "GSM2890011_9033_PS08_L.NI.CEL.gz","GSM2890012_9034_PS09_L.NI.CEL.gz",
                                                   "GSM2890013_9035_PS10_L.NI.CEL.gz")))

write.csv(dat5, "CLEAN_DATA_GSE108110.csv")
