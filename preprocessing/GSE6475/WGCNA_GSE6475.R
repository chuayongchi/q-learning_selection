BiocManager::install("impute")
BiocManager::install("WGCNA")
BiocManager::install("Biobase")
BiocManager::install("GO.db")
install.packages("devtools")
devtools::install_github("kevinblighe/CorLevelPlot")

library(devtools)
library(CorLevelPlot)
library(WGCNA)
library(tidyverse)
library(magrittr)

# obtain the gse object from the GEOquery
library(GEOquery)
my_id <- "GSE6475"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

# extract sample info
sampleInfo <- pData(gse)
sampleInfo


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
# tidyverse will pull in ggplot2, readr, other useful libraries
# magrittr provides the %>% operator


## load data
data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE6475_RAW/DEgenes/DEGenes_Matrix_GSE6475.csv")# <= path to the data file
data
sampleInfo <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE6475_RAW/DEGenes/sampleInfo_GSE6475.csv")
## extract the data with the gene_entrez_Id as rownames
data2 <- data[,3:ncol(data)]
rownames(data2) <- data[,2]
#colnames(data2) <- sampleInfo$group
sampleInfo2 <- sampleInfo[,2:ncol(sampleInfo)]
rownames(sampleInfo2) <- sampleInfo[,1]

colnames(data2)[1] <- "GSM148725"
colnames(data2)[2] <- "GSM148748"
colnames(data2)[3] <- "GSM148762"
colnames(data2)[4] <- "GSM148763"
colnames(data2)[5] <- "GSM148764"

colnames(data2)[6] <- "GSM148765"
colnames(data2)[7] <- "GSM148766"
colnames(data2)[8] <- "GSM148767"
colnames(data2)[9] <- "GSM148768"
colnames(data2)[10] <- "GSM148769"

colnames(data2)[11] <- "GSM148770"
colnames(data2)[12] <- "GSM148771"
colnames(data2)[13] <- "GSM148887"
colnames(data2)[14] <- "GSM148888"
colnames(data2)[15] <- "GSM148889"

colnames(data2)[16] <- "GSM148890"
colnames(data2)[17] <- "GSM148892"
colnames(data2)[18] <- "GSM148894"

dim(data2)
input_mat = t(data2)

input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns

#library(WGCNA)
allowWGCNAThreads()          # allow multi-threading (optional)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.8;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.80, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 20
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )
##############################################################
##### Relate Module (cluster) Assignments to Treatment Groups
# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)


module_df[1:5,]

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$sample = row.names(MEs0)
#################################################################
module_eigengenes <- MEs0



#####traits
traits <- sampleInfo2 %>%
  mutate(acne = ifelse(grepl('acne skin sample', group), 1, 0)) %>%
  select(3)

# Define numbers of genes and samples
nSamples <- nrow(input_mat)
nGenes <- ncol(input_mat)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[3],
             y = names(heatmap.data)[1:2],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

gene_list <- module_df %>% 
  filter(module_df$colors == 'turquoise' | module_df$colors == 'blue') %>% 
  rownames()

##top25 as the gene signature
# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, input_mat, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:2,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(input_mat, traits$acne, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene_sign_list <- gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(10)

#####100 not good
top_100_genes_sign_list <- gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(100)

top_100_genes_sign_list <- rownames(top_100_gene_sign_list)
write.csv(top_100_genes_sign_list,"top_100_genes_sign_list.csv")


write.csv(gene_sign_list,"all_gene_sign_list_GSE6475.csv")
# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-sample) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=sample, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")



##### Examine Expression Profiles
# pick out a few modules of interest here
modules_of_interest = c("turquoise", "blue")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

submod = module_df %>%
  subset(gene_id %in% top_100_genes_sign_list)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
data2[1:5,1:10]


subexpr = data2[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")



##### Generate and Export Networks
genes_of_interest = module_df %>%
  subset(gene_id %in% top_100_genes_sign_list)

expr_of_interest = data2[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

pheatmap(TOM) 

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

# Export Network file to be read into Cytoscape, VisANT, etc
write.csv(TOM,"WGCNA_Matrix_GSE6475_Final_top100.csv")
