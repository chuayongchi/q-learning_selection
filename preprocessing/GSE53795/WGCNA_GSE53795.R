BiocManager::install("impute")
BiocManager::install("WGCNA")
BiocManager::install("Biobase")
BiocManager::install("GO.db")
library(pheatmap)
library(devtools)
library(CorLevelPlot)
library(WGCNA)
library(tidyverse)
library(magrittr)
# tidyverse will pull in ggplot2, readr, other useful libraries
# magrittr provides the %>% operator

library(GEOquery)
my_id <- "GSE53795"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]
gse


# extract sample info
sampleInfo <- pData(gse)
sampleInfo

library(dplyr)
#patient id#:ch1
## patient id#:ch1 and tissue subgroup:ch1 seem to contain factors we might need for the analysis. Let's pick just those columns
sampleInfo <- select(sampleInfo, 'tissue subgroup:ch1', 'patient id#:ch1')

## rename the column
sampleInfo <- rename(sampleInfo, group = 'tissue subgroup:ch1', patient = 'patient id#:ch1')
sampleInfo


## load data   DEGenes/DEGenes_Matrix_GSE53795.csv
gse_norm <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/CLEAN_DATA_GSE53795.csv")
gse_norm = gse_norm[,2:ncol(gse_norm)]


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

data <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/DEGenes/DEGenes_Matrix_GSE53795.csv")# <= path to the data file
data
sampleInfo <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/DEGenes/sampleInfo_GSE53795.csv")
## extract the data with the gene_entrez_Id as rownames
data2 <- data[,3:ncol(data)]
rownames(data2) <- data[,2]

sampleInfo2 <- sampleInfo[,2:ncol(sampleInfo)]
rownames(sampleInfo2) <- sampleInfo[,1]
#colnames(data2) <- sampleInfo$group
#data2 = data[,2:ncol(data)]
gse_norm2 <- gse_norm[,2:ncol(gse_norm)]
rownames(gse_norm2) <- gse_norm[,1]

dim(data2)
dim(gse_norm2)
input_mat = t(data2)
input_mat = t(gse_norm2)
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

# Add treatment names
MEs0$sample = row.names(MEs0)

module_eigengenes <- MEs0


#####traits
traits <- sampleInfo2 %>%
  mutate(acne = ifelse(grepl('non-lesional skin', group), 0, 1)) %>%
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
             x = names(heatmap.data)[4],
             y = names(heatmap.data)[1:3],
             col = c("blue1", "skyblue", "white", "pink", "red"))


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

write.csv(gene_sign_list,"all_gene_sign_list_GSE53795.csv")

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
  subset(colors %in% modules_of_interest)

expr_of_interest = data2[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]

TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

# Export Network file to be read into Cytoscape, VisANT, etc
write.csv(TOM,"WGCNA_Matrix_GSE53795_Final.csv")

GGC <- read.csv("C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE53795_RAW/DEGenes/WGCNA_Matrix_GSE53795_Final.csv")
pheatmap(GGC) 
