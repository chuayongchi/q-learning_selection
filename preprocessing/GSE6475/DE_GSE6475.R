install.packages("BiocManager")
install.packages("forcats")
install.packages("stringr")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("readr")
install.packages("tidyr")
install.packages("survminer")
install.packages("reshape2")
install.packages('outliers')
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

# obtain the gse object from the GEOquery
library(GEOquery)
my_id <- "GSE6475"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

# read the count data that have been processed
gse_norm <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE6475_RAW/CLEAN_DATA_GSE6475.csv")
gse_norm <- gse_norm[,2:ncol(gse_norm)]

# extract sample info
sampleInfo <- pData(gse)
sampleInfo

# making sure the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change colnames to match with sample info
## bcz the sample name in the sample info is GSM1300911, GSM1300912, etc
## rename in norm1
## rename the samples of data input and make sure the order of the sample is the same as the ori 
colnames(gse_norm[-1])

colnames(gse_norm)[2] <- "GSM148725"
colnames(gse_norm)[3] <- "GSM148748"
colnames(gse_norm)[4] <- "GSM148762"
colnames(gse_norm)[5] <- "GSM148763"
colnames(gse_norm)[6] <- "GSM148764"

colnames(gse_norm)[7] <- "GSM148765"
colnames(gse_norm)[8] <- "GSM148766"
colnames(gse_norm)[9] <- "GSM148767"
colnames(gse_norm)[10] <- "GSM148768"
colnames(gse_norm)[11] <- "GSM148769"

colnames(gse_norm)[12] <- "GSM148770"
colnames(gse_norm)[13] <- "GSM148771"
colnames(gse_norm)[14] <- "GSM148887"
colnames(gse_norm)[15] <- "GSM148888"
colnames(gse_norm)[16] <- "GSM148889"

colnames(gse_norm)[17] <- "GSM148890"
colnames(gse_norm)[18] <- "GSM148892"
colnames(gse_norm)[19] <- "GSM148894"

# check the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change the rownames to gene entrez_id
gse_norm2 <- gse_norm[,-1]
rownames(gse_norm2) <- gse_norm[,1]

# see the statistical info of dataset
summary(gse_norm2)

## get the distribution
b <- boxplot(gse_norm2)



library(dplyr)
## extract the characteristics_ch1 for obtaining the group of sample, such as lesional skin or normal skin
## also extract title to obtain skin sample by patient
sampleInfo <- select(sampleInfo, title, characteristics_ch1)

## rename the column to group
sampleInfo <- rename(sampleInfo, group = characteristics_ch1, patient = title)
sampleInfo
sampleInfo$patient
sampleInfo$patient[1] <- "Pt_1"
sampleInfo$patient[2] <- "Pt_1"
sampleInfo$patient[3] <- "Pt_2"
sampleInfo$patient[4] <- "Pt_2"
sampleInfo$patient[5] <- "Pt_3"
sampleInfo$patient[6] <- "Pt_3"
sampleInfo$patient[7] <- "Pt_4"
sampleInfo$patient[8] <- "Pt_4"
sampleInfo$patient[9] <- "Pt_5"
sampleInfo$patient[10] <- "Pt_5"
sampleInfo$patient[11] <- "Pt_6"
sampleInfo$patient[12] <- "Pt_6"
## non acne patient
sampleInfo$patient[13] <- "NPt_1"
sampleInfo$patient[14] <- "NPt_2"
sampleInfo$patient[15] <- "NPt_3"
sampleInfo$patient[16] <- "NPt_4"
sampleInfo$patient[17] <- "NPt_5"
sampleInfo$patient[18] <- "NPt_6"

sampleInfo
## there are correction need to be done 
## bcz the data extracted from GEOquery is wrong for the characteristics_ch1 
## the GSM148770 is acne skin but the extracted data is written as normal skin
sampleInfo$group[11] <- "acne skin sample"
sampleInfo

###
htree <- hclust(dist(t(gse_norm2)), method = "average")
plot(htree)

library(pheatmap)
## argument use="c" stops an error if there are any missing data points
corMatrix <- cor(gse_norm2,use="c")
corMatrix

### change colnames for matrix to match the sample info
pheatmap(corMatrix) 

## Print the rownames of the sample information and check it matches the correlation matrix
rownames(sampleInfo)
colnames(corMatrix)

# making sure the row names in colData matches to column names in counts_data
all(colnames(corMatrix) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(corMatrix) == rownames(sampleInfo))

pheatmap(corMatrix,
         annotation_col=sampleInfo) 



library(ggplot2)
library(ggrepel)

## To transpose the expression matrix for prcompt()
pca <- prcomp(t(gse_norm2))
pca$x

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste(" ", patient))) + geom_point() + geom_text_repel()

library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
## rename the column names
colnames(design) <- c("acne_skin","non_acne_patient_normal_skin", "normal_skin")

## calculate median expression level
cutoff <- median(unlist(gse_norm2), na.rm = TRUE)

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- gse_norm2 > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse_norm3 <- gse_norm2[keep,]
fit <- lmFit(gse_norm3, design)
head(fit$coefficients)

## view differential expressed by creating contrast
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)
contrasts <- makeContrasts( acne_skin - normal_skin,
                            normal_skin - non_acne_patient_normal_skin,
                            levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

#display the result of the contrast
topTable(fit2, coef=1)
topTable(fit2, coef=2)


##  to know how many genes are differentially-expressed overall we can use the decideTests function
decideTests(fit2)
table(decideTests(fit2))


### Further processing and visualization of the DE results
topTable(fit2)
full_results <- topTable(fit2, coef=1, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()
p <- ggplot(data=full_results, aes(x=logFC, y=-log10(`adj.P.Val`))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")



## set 5% of tolerance rate for false positive
p_cutoff <- 0.05 
fc_cutoff <- 0.13 # 0.58 ## Our experience suggests a minimal value, such as a 10% fold-change, corresponding to τ=log2(1.1)=0.13 on the log2-scale. 

summary(full_results)

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()


full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = -log10(`adj.P.Val`), col=Significant)) + geom_point()

library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, ID,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")












##Filtering and exporting the results table

## filter according to p-value (adjusted) and fold-change cut-offs
p_cutoff <- 0.05
fc_cutoff <- 0.13

filter(full_results, full_results$ID == 2357)
filter(full_results, adj.P.Val < 0.05, abs(logFC )> 0.13)
DEGenes <- filter(full_results, adj.P.Val < 0.05, abs(logFC )> 0.13)

DEGenes_ordered <- DEGenes[order(DEGenes$adj.P.Val), ]
top20_DEGenes <- DEGenes_ordered[1:20, ]
top20_DEG_list <- top20_DEGenes$ID
top20_features_DEG <- features[features$ENTREZ_GENE_ID %in% top20_DEG_list, ]


## check the existence of extracted gene ID to the input data
DEG_list <- DEGenes$ID
all(DEG_list %in% rownames(gse_norm2))

DEgenes_data <- gse_norm[gse_norm$Gene_Entrez_ID %in% DEG_list, ]
rownames(DEgenes_data) <- NULL
write.csv(DEgenes_data, "DEGenes_Matrix_GSE6475.csv")
write.csv(top20_DEG_list, "Top20_DEG_list_GSE6475.csv")


### obtain the features for the DEgenes
features <- fData(gse)
View(features)
features_DEG <- features[features$ENTREZ_GENE_ID %in% DEG_list, ]

temp_gene <- features[features$`Gene Symbol` == "FPR1", ]

write.csv(features_DEG, "Features_DEGenes_GSE6475.csv")
