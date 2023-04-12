install.packages("BiocManager")
library("BiocManager")
BiocManager::install("affy")
BiocManager::install("annotate")
BiocManager::install("hgu133plus2.db")
# load packages
library("affy")
library("annotate")
library("hgu133plus2.db")

#set the working directory to C:/Users/USER/OneDrive/Documents/Reinforcement Learning/GSE108110_RAW
# Read CEL files from the working directory 
affydata = ReadAffy()

# normalizing data
# rma() - 'Robust Multi-Chip' average
# use rma to background correct and normalize probe levels
rma = rma(affydata)

# expression set
ed = exprs(rma)
write.csv(ed,"probe_GSE108110.csv")

# use hgu133plus2.db to convert to Entres ID
Entrez_IDs = unlist(mget(rownames(ed), hgu133plus2ENTREZID, ifnotfound=NA))
write.csv(Entrez_IDs,"probe_ID_GSE108110.csv")

# combine matrix with gene Entrez ID
mRNA_matrixTrain = cbind(Entrez_IDs,ed)
row.names(mRNA_matrixTrain) = NULL
row.names(mRNA_matrixTrain) = mRNA_matrixTrain[,1]

# remove probe ID
mRNA_matrixTrain = mRNA_matrixTrain[,2:ncol(mRNA_matrixTrain)]
mode(mRNA_matrixTrain) <- "numeric"
write.csv(mRNA_matrixTrain,"RAW_GSE108110.csv")
