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
install.packages('pca3d')
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("pheatmap")
BiocManager::install("org.Hs.eg.db")

# obtain the gse object from the GEOquery
library(GEOquery)
my_id <- "GSE108110"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

gse_norm <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE108110_RAW/CLEAN_DATA_GSE108110.csv")
gse_norm <- gse_norm[,2:ncol(gse_norm)]

# extract sample info
sampleInfo <- pData(gse)
sampleInfo

# checking the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change colnames to match with sample info
## bcz the sample name in the sample info is GSM2889960, GSM2889961, etc
## rename in norm
## rename the samples of data input and make sure the order of the sample is the same as the ori 
colnames(gse_norm[-1])

colnames(gse_norm)[2] <- "GSM2889960"
colnames(gse_norm)[3] <- "GSM2889961"
colnames(gse_norm)[4] <- "GSM2889962"
colnames(gse_norm)[5] <- "GSM2889963"
colnames(gse_norm)[6] <- "GSM2889964"
colnames(gse_norm)[7] <- "GSM2889965"
colnames(gse_norm)[8] <- "GSM2889966"
colnames(gse_norm)[9] <- "GSM2889967"
colnames(gse_norm)[10] <- "GSM2889968"
colnames(gse_norm)[11] <- "GSM2889969"
colnames(gse_norm)[12] <- "GSM2889970"
colnames(gse_norm)[13] <- "GSM2889971"
colnames(gse_norm)[14] <- "GSM2889972"
colnames(gse_norm)[15] <- "GSM2889973"
colnames(gse_norm)[16] <- "GSM2889974"
colnames(gse_norm)[17] <- "GSM2889975"
colnames(gse_norm)[18] <- "GSM2889976"
colnames(gse_norm)[19] <- "GSM2889977"
colnames(gse_norm)[20] <- "GSM2889978"
colnames(gse_norm)[21] <- "GSM2889979"
colnames(gse_norm)[22] <- "GSM2889980"
colnames(gse_norm)[23] <- "GSM2889981"
colnames(gse_norm)[24] <- "GSM2889982"
colnames(gse_norm)[25] <- "GSM2889983"
colnames(gse_norm)[26] <- "GSM2889984"
colnames(gse_norm)[27] <- "GSM2889985"
colnames(gse_norm)[28] <- "GSM2889986"
colnames(gse_norm)[29] <- "GSM2889987"
colnames(gse_norm)[30] <- "GSM2889988"
colnames(gse_norm)[31] <- "GSM2889989"
colnames(gse_norm)[32] <- "GSM2889990"
colnames(gse_norm)[33] <- "GSM2889991"
colnames(gse_norm)[34] <- "GSM2889992"
colnames(gse_norm)[35] <- "GSM2889993"
colnames(gse_norm)[36] <- "GSM2889994"
colnames(gse_norm)[37] <- "GSM2889995"
colnames(gse_norm)[38] <- "GSM2889996"
colnames(gse_norm)[39] <- "GSM2889997"
colnames(gse_norm)[40] <- "GSM2889998"
colnames(gse_norm)[41] <- "GSM2889999"
colnames(gse_norm)[42] <- "GSM2890000"
colnames(gse_norm)[43] <- "GSM2890001"
colnames(gse_norm)[44] <- "GSM2890002"
colnames(gse_norm)[45] <- "GSM2890003"
colnames(gse_norm)[46] <- "GSM2890004"
colnames(gse_norm)[47] <- "GSM2890005"
colnames(gse_norm)[48] <- "GSM2890006"
colnames(gse_norm)[49] <- "GSM2890007"
colnames(gse_norm)[50] <- "GSM2890008"
colnames(gse_norm)[51] <- "GSM2890009"
colnames(gse_norm)[52] <- "GSM2890010"
colnames(gse_norm)[53] <- "GSM2890011"
colnames(gse_norm)[54] <- "GSM2890012"
colnames(gse_norm)[55] <- "GSM2890013"

write.csv(gse_norm, "Data_Ready_GSE108110.csv")

# check the row names in colData matches to column names in counts_data
all(colnames(gse_norm[-1]) %in% rownames(sampleInfo))

# are they in the same order?
all(colnames(gse_norm[-1]) == rownames(sampleInfo))

## change the rownames to gene entrez_id
gse_norm2 <- gse_norm[,-1]
rownames(gse_norm2) <- gse_norm[,1]
#####################################
## extract the require sample
gse_norm3 <- select(gse_norm2, GSM2890005, GSM2890006, GSM2890007, GSM2890008, GSM2890009, GSM2890010, GSM2890011, GSM2890012, GSM2890013,
                    GSM2889996, GSM2889997, GSM2889998, GSM2889999, GSM2890000, GSM2890001, GSM2890002, GSM2890003, GSM2890004)
#####################################
# see the statistical info of dataset
summary(gse_norm3)


## get the distribution
b <- boxplot(gse_norm3)


library(dplyr)
## extract the source_name_ch1 for obtaining the group of sample, such as lesional skin or normal skin
## also extract title to obtain skin sample by patient
sampleInfo <- select(sampleInfo, title, source_name_ch1)

## rename the column to group
sampleInfo <- rename(sampleInfo, group = source_name_ch1, patient = title)
sampleInfo

## change patient to only contain series number 
sampleInfo$patient
sampleInfo$patient[1] <- "9016"
sampleInfo$patient[2] <- "9017"
sampleInfo$patient[3] <- "9019"
sampleInfo$patient[4] <- "9020"
sampleInfo$patient[5] <- "9021"
sampleInfo$patient[6] <- "9025"
sampleInfo$patient[7] <- "9026"
sampleInfo$patient[8] <- "9028"
sampleInfo$patient[9] <- "9032"

sampleInfo$patient[10] <- "9016"
sampleInfo$patient[11] <- "9017"
sampleInfo$patient[12] <- "9019"
sampleInfo$patient[13] <- "9020"
sampleInfo$patient[14] <- "9021"
sampleInfo$patient[15] <- "9025"
sampleInfo$patient[16] <- "9026"
sampleInfo$patient[17] <- "9028"
sampleInfo$patient[18] <- "9032"

sampleInfo$patient[19] <- "9016"
sampleInfo$patient[20] <- "9017"
sampleInfo$patient[21] <- "9019"
sampleInfo$patient[22] <- "9020"
sampleInfo$patient[23] <- "9021"
sampleInfo$patient[24] <- "9025"
sampleInfo$patient[25] <- "9026"
sampleInfo$patient[26] <- "9028"
sampleInfo$patient[27] <- "9032"

sampleInfo$patient[28] <- "9018"
sampleInfo$patient[29] <- "9022"
sampleInfo$patient[30] <- "9023"
sampleInfo$patient[31] <- "9024"
sampleInfo$patient[32] <- "9027"
sampleInfo$patient[33] <- "9029"
sampleInfo$patient[34] <- "9033"
sampleInfo$patient[35] <- "9034"
sampleInfo$patient[36] <- "9035"

sampleInfo$patient[37] <- "9018"
sampleInfo$patient[38] <- "9022"
sampleInfo$patient[39] <- "9023"
sampleInfo$patient[40] <- "9024"
sampleInfo$patient[41] <- "9027"
sampleInfo$patient[42] <- "9029"
sampleInfo$patient[43] <- "9033"
sampleInfo$patient[44] <- "9034"
sampleInfo$patient[45] <- "9035"

sampleInfo$patient[46] <- "9018"
sampleInfo$patient[47] <- "9022"
sampleInfo$patient[48] <- "9023"
sampleInfo$patient[49] <- "9024"
sampleInfo$patient[50] <- "9027"
sampleInfo$patient[51] <- "9029"
sampleInfo$patient[52] <- "9033"
sampleInfo$patient[53] <- "9034"
sampleInfo$patient[54] <- "9035"

sampleInfo2 <- sampleInfo[37:54, ]
library(pheatmap)
## argument use="c" stops an error if there are any missing data points
corMatrix <- cor(gse_norm3,use="c")
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
         annotation_col=sampleInfo2) 

###tree
htree <- hclust(dist(t(gse_norm3)), method = "average")
plot(htree)


library(ggplot2)
library(ggrepel)

## To transpose the expression matrix for prcompt()
pca <- prcomp(t(gse_norm3))
pca$x

## Join the PCs to the sample information
cbind(sampleInfo2, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste(" ", patient))) + geom_point() + geom_text_repel()

library(pca3d)
pca3d(pca, group=sampleInfo$group)
pca2d(pca, group=sampleInfo$group, legend="topleft")
pca3d(pca, group=sampleInfo$group, show.ellipses=TRUE,
      ellipse.ci=0.75, show.plane=FALSE)


## interpretation
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree Plot", xlab = "Principal Components", ylab = "Percent Variation")


library(limma)
design <- model.matrix(~0+sampleInfo2$group)
design
## rename the column names
### NSP_G1 <- "Non_scar_prone_non_lesional_skin"
### NSP_G2 <- "Non_scar_prone_papule_at_21_days"
### NSP_G3 <- "Non_scar_prone_papule_less_than_48_H_old"
### SP_G1 <- "Scar-prone non-lesional skin"
### SP_G2 <- "Scar-prone papule at 21 days"
### sP_G3 <- "Scar-prone papule less than 48 H old"
colnames(design) <- c(#"NSP_G1",
                      #"NSP_G2", 
                      #"NSP_G3",
                      "SP_G1",
                      "SP_G2")
                      #"SP_G3")

## calculate median expression level
cutoff <- median(unlist(gse_norm3), na.rm = TRUE)

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- gse_norm3 > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse_norm3 <- gse_norm3[keep,]
fit <- lmFit(gse_norm3, design)
head(fit$coefficients)

## view differential expressed by creating contrast
## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)
contrasts <- makeContrasts( SP_G2 - SP_G1,
                            #SP_G3 - SP_G1,
                            #SP_G3 - SP_G2,
                            #NSP_G2 - NSP_G1,
                            #NSP_G3 - NSP_G1,
                            levels=design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

#display the result of the contrast
topTable(fit2)
#topTable(fit2, coef=2)
#topTable(fit2, coef=3)
#topTable(fit2, coef=4)

##  to know how many genes are differentially-expressed overall we can use the decideTests function
decideTests(fit2)
table(decideTests(fit2))



### Further processing and visualization of the DE results
topTable(fit2)
full_results_co1 <- topTable(fit2, number=Inf)
full_results_co1 <- tibble::rownames_to_column(full_results_co1,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results_co1,aes(x = logFC, y=-log10(`adj.P.Val`))) + geom_point()

p <- ggplot(data=full_results_co1, aes(x=logFC, y=-log10(`adj.P.Val`))) + geom_point()
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
full_results_co1$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
full_results_co1$diffexpressed[full_results_co1$logFC > 0.58 & full_results_co1$adj.P.Val < 0.05] <- "UP"
full_results_co1$diffexpressed[full_results_co1$logFC < -0.58 & full_results_co1$adj.P.Val < 0.05] <- "DOWN"


# Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=full_results_co1, aes(x=logFC, y=-log10(`adj.P.Val`), col=diffexpressed)) + geom_point() + theme_minimal()

# Add lines as before...
p2 <- p + geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")

## Change point color 

# 1. by default, it is assigned to the categories in an alphabetical order):
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))

# 2. to automate a bit: ceate a named vector: the values are the colors to be used, the names are the categories they will be assigned to:
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)


# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
full_results_co1$delabel <- NA
full_results_co1$delabel[full_results_co1$diffexpressed != "NO"] <- full_results_co1$ID[full_results_co1$diffexpressed != "NO"]

ggplot(data=full_results_co1, aes(x=logFC, y=-log10(`adj.P.Val`), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text()


library(ggrepel)
# plot adding up all layers we have seen so far
ggplot(data=full_results_co1, aes(x=logFC, y=-log10(`adj.P.Val`), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "black", "red")) +
  geom_vline(xintercept=c(-0.58, 0.58), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


##Filtering and exporting the results table
DEGenes <- full_results_co1[full_results_co1$diffexpressed != "NO", ]
DEG_list <- DEGenes$ID


DEGenes_Up <- full_results_co1[full_results_co1$diffexpressed == "UP", ]
DEGenes_Up_ordered <- DEGenes_Up[order(DEGenes_Up$adj.P.Val), ]
top10_DEGenes_Up <- DEGenes_Up_ordered[1:10, ]
top10_DEG_Up_list <- top10_DEGenes_Up$ID
write.csv(top10_DEG_Up_list, "Top10_DEG_Up_list_GSE108110_coef1.csv")

DEGenes_Down <- full_results_co1[full_results_co1$diffexpressed == "DOWN", ]
DEGenes_Down_ordered <- DEGenes_Down[order(DEGenes_Down$adj.P.Val), ]
top10_DEGenes_Down <- DEGenes_Down_ordered[1:10, ]
top10_DEG_Down_list <- top10_DEGenes_Down$ID
write.csv(top10_DEG_Down_list, "Top10_DEG_Down_lis_GSE108110_coef1.csv")

features <- fData(gse)
View(features)

top20_DEG_list <- c(top10_DEG_Up_list, top10_DEG_Down_list)
features_DEG <- features[features$ENTREZ_GENE_ID %in% DEG_list, ]
write.csv(top20_features_DEG, "Top20_Features_DEGenes_GSE108110_coef1.csv")

DEgenes_data <- gse_norm[gse_norm$Gene_Entrez_ID %in% DEG_list, ]
DEgenes_data_final <- select(DEgenes_data, Gene_Entrez_ID, GSM2890005, GSM2890006, GSM2890007, GSM2890008, GSM2890009, GSM2890010, GSM2890011, GSM2890012, GSM2890013,
                    GSM2889996, GSM2889997, GSM2889998, GSM2889999, GSM2890000, GSM2890001, GSM2890002, GSM2890003, GSM2890004)
write.csv(DEgenes_data_final, "DEGenes_Matrix_GSE108110_Final.csv")



