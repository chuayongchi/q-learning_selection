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

library(GEOquery)
my_id <- "GSE53795"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

# read the count data that have been processed
gse_norm <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/CLEAN_DATA_GSE53795.csv")
gse_norm = gse_norm[,2:ncol(gse_norm)]

# extract sample info
sampleInfo <- pData(gse)
sampleInfo

# making sure the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change colnames to match with sample info
## bcz the sample name in the sample info is GSM1300911, GSM1300912, etc
## rename in norm and make sure the order is the same
colnames(gse_norm[-1])

colnames(gse_norm)[2] <- "GSM1300911"
colnames(gse_norm)[3] <- "GSM1300912"
colnames(gse_norm)[4] <- "GSM1300913"
colnames(gse_norm)[5] <- "GSM1300914"
colnames(gse_norm)[6] <- "GSM1300915"
colnames(gse_norm)[7] <- "GSM1300916"
colnames(gse_norm)[8] <- "GSM1300917"
colnames(gse_norm)[9] <- "GSM1300918"
colnames(gse_norm)[10] <- "GSM1300919"
colnames(gse_norm)[11] <- "GSM1300920"
colnames(gse_norm)[12] <- "GSM1300921"
colnames(gse_norm)[13] <- "GSM1300922"
colnames(gse_norm)[14] <- "GSM1300923"
colnames(gse_norm)[15] <- "GSM1300924"
colnames(gse_norm)[16] <- "GSM1300925"
colnames(gse_norm)[17] <- "GSM1300926"
colnames(gse_norm)[18] <- "GSM1300927"
colnames(gse_norm)[19] <- "GSM1300928"
colnames(gse_norm)[20] <- "GSM1300929"
colnames(gse_norm)[21] <- "GSM1300930"
colnames(gse_norm)[22] <- "GSM1300931"
colnames(gse_norm)[23] <- "GSM1300932"
colnames(gse_norm)[24] <- "GSM1300933"
colnames(gse_norm)[25] <- "GSM1300934"

# check the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change the rownames to gene enrez id
gse_norm2 <- gse_norm[,-1]
rownames(gse_norm2) <- gse_norm[,1]

# see the statistical info of dataset
summary(gse_norm2)

## get the distribution
b <- boxplot(gse_norm2)



library(dplyr)
#patient id#:ch1
## patient id#:ch1 and tissue subgroup:ch1 seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, 'tissue subgroup:ch1', 'patient id#:ch1')

## rename the column
sampleInfo <- rename(sampleInfo, group = 'tissue subgroup:ch1', patient = 'patient id#:ch1')
sampleInfo


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

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo) 



library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(gse_norm2))
pca$x

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Pt", patient))) + geom_point() + geom_text_repel()


library(limma)
design <- model.matrix(~0+sampleInfo$group)
design
## rename the column names
colnames(design) <- c("lesional_skin","non_lesional_skin")

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
contrasts <- makeContrasts( lesional_skin - non_lesional_skin,
                            levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
#display the result of the contrast
topTable(fit2)

##  to know how many genes are differentially-expressed overall we can use the decideTests function
decideTests(fit2)
table(decideTests(fit2))

## The arrayWeights function will assign a score to each sample; 
## with a value of 1 implying equal weight. 
## Samples with score less than 1 are down-weights, 
## and samples with scores greater than 1 are up-weighted. 
## calculate relative array weights
aw <- arrayWeights(gse_norm3,design)
aw

## The lmFit function can accept weights, 
## and the rest of the code proceeds as above.
fit3 <- lmFit(gse_norm3, design, weights = aw)
contrasts <- makeContrasts( lesional_skin - non_lesional_skin,
                            levels=design)

fit4 <- contrasts.fit(fit3, contrasts)
fit4 <- eBayes(fit4)

topTable(fit4)

##  to know how many genes are differentially-expressed overall we can use the decideTests function
decideTests(fit4)
table(decideTests(fit4))



### Further processing and visualization of the DE results
topTable(fit2)
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
## every coef is very different to each other; need to have a reference to decide what level
p_cutoff <- 0.05 ##or 0.001 #need to confrim back the log for the preprocessing ##0.05
fc_cutoff <- 0.13 ##Our experience suggests a minimal value, such as a 10% fold-change, corresponding to τ=log2(1.1)=0.13 on the log2-scale. 

summary(full_results)

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()



library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 0.13
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, ID,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")





##Filtering and exporting the results table
## filter according to p-value (adjusted) and fold-change cut-offs
p_cutoff <- 0.05
fc_cutoff <- 1

summary(full_results)
filter(full_results, adj.P.Val < 0.05)
filter(full_results, abs(logFC )> 1)
filter(full_results, adj.P.Val < 0.05, abs(logFC )> 0.13)

DEGenes <- filter(full_results, adj.P.Val < 0.05, abs(logFC )> 0.13)

DEG_list <- DEGenes$ID

all(DEG_list %in% rownames(gse_norm2))

DEgenes_data <- gse_norm[gse_norm$Gene_Entrez_ID %in% DEG_list, ]

rownames(DEgenes_data) <- NULL
write.csv(DEgenes_data, "DEGenes_Matrix_GSE53795.csv")


### obtain the features for the DEgenes
features <- fData(gse)
View(features)
features_DEG <- features[features$ENTREZ_GENE_ID %in% DEG_list, ]
write.csv(features_DEG, "Features_DEGenes_GSE53795.csv")
