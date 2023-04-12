install.packages("easyPubMed")
library(easyPubMed)
setwd("C:/Users/user/Desktop")
gene.acne.qlearning <- read.table("C:/Users/user/Desktop/gene acne qlearning.txt", quote="\"", comment.char="")
m=as.matrix(gene.acne.qlearning)

vdg <- matrix(NA, nrow = length(m), ncol = 1);
rownames(vdg) <- m;
colnames(vdg) <- c("acne PMIDs")

for ( i in 1 : length(m)){
    gname=m[i]
    
    # all cancers: text data mining
    dami_query_string=(paste(gname,"acne",sep=" AND "))
    dami_on_pubmed=get_pubmed_ids(dami_query_string)
    num_pmid=(dami_on_pubmed$Count)
    if(num_pmid>1){
        pmid1=(unlist(dami_on_pubmed$IdList))
        pmid=pmid1[1]
    }
    if(num_pmid==1){
        pmid=(unlist(dami_on_pubmed$IdList))
    }
    if(num_pmid==0){
        pmid=0
    }
    
    vdg[i] <- pmid
}

write.csv(vdg,"acne.csv")