BiocManager::install("impute")
BiocManager::install("WGCNA")
BiocManager::install("Biobase")
BiocManager::install("GO.db")
library(pheatmap)
library(devtools)
library(CorLevelPlot)
library(WGCNA)
# tidyverse will pull in ggplot2, readr, other useful libraries
library(tidyverse)
# provides the %>% operator
library(magrittr)

# obtain the gse object from the GEOquery
library(GEOquery)
my_id <- "GSE108110"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse

# extract sample info
sampleInfo <- pData(gse)
sampleInfo

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

## load data
sampleInfo <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE108110_RAW/DEGenes/sampleInfo2_GSE108110.csv")
data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE108110_RAW/DEGenes/DEGenes_Matrix_GSE108110_Final.csv")# <= path to the data file
data
sampleInfo2 <- sampleInfo[,2:ncol(sampleInfo)]
rownames(sampleInfo2) <- sampleInfo[,1]
## extract the data with the gene_entrez_Id as rownames
data2 <- data[,3:ncol(data)]
rownames(data2) <- data[,2]
#colnames(data2) <- sampleInfo$group

#################
## extract the require sample
data3 <- select(data2, GSM2890005, GSM2890006, GSM2890007, GSM2890008, GSM2890009, GSM2890010, GSM2890011, GSM2890012, GSM2890013,
                    GSM2889996, GSM2889997, GSM2889998, GSM2889999, GSM2890000, GSM2890001, GSM2890002, GSM2890003, GSM2890004)
###################

dim(data2)
input_mat = t(data2)
# Look at first 5 rows and 10 columns
input_mat[1:5,1:10]

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






##### Relate Module (cluster) Assignments to Treatment Groups
# netwk$colors[netwk$blockGenes[[1]]]
# table(netwk$colors)
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)
###########dont run
# Add treatment names
MEs0$treatment = row.names(MEs0)
###########dont run

module_eigengenes <- MEs0

sampleInfo2 <- sampleInfo[37:54, ]

acne_list <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

acne_list <- c(1,1,1,1,1,1,1,1,1,
               0,0,0,0,0,0,0,0,0)

#####traits
traits <- sampleInfo2 %>%
  mutate(acne = acne_list) %>%
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


##top25 as the gene signature
# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, input_mat, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:3,1:10]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(input_mat, traits$acne, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene_sign_list <- gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(10)
write.csv(gene_sign_list,"all_gene_sign_list_GSE108110.csv")


# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
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
modules_of_interest = c( "blue", "turquoise")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

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
  labs(x = "sample",
       y = "normalized expression")



##### Generate and Export Networks
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = data2[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

# range between 0 to 1
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

# Export Network file to be read into Cytoscape, VisANT, etc
write.csv(TOM,"WGCNA_Matrix_GSE108110_Final.csv")
