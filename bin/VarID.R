#VarID
library(tidyverse)
library(viridis)
library(RaceID)
library(assertthat)
library(Matrix)
library(Seurat)

date = Sys.Date()
load("data/seurat-with-batch-effect-removal.RData")

prdata <- all@assays$RNA@counts
rm(all)

prdata <- as.matrix(prdata)
prdata <- prdata[rowSums(prdata)>0,]

sc <- SCseq(as.data.frame(prdata))

# filtering of expression data
#there are 2 batching variables: 
#1. is the date of sequencing, 
#2. the other is the question if the sample was digested with accumax or not

metadata <- data.frame(gender = factor(ifelse(grepl("^(c3|p3)", colnames(prdata)), "male", "female"), levels = c("female", "male")))
rownames(metadata) <- colnames(prdata)

#sanity check
assert_that(sum(is.na(metadata$gender))==0)

#build batch effect table
batch <- metadata[["gender"]]
names(batch) <- rownames(metadata)
  
# filtering of expression data
sc <- filterdata(sc, 
                 mintotal=500,
                 knn=10,
                 minexpr = 1,
                 CGenes=c("FOS",
                          "JUN"),
                 FGenes = c("^RPS*", 
                            "^RPL*")
)

d <- getExpData(sc)

#run pruneKnn and remove batch effects
res <- pruneKnn(d,large=TRUE, regNB=TRUE,pcaComp=20,batch=batch,metric="pearson",genes=NULL,knn=10,alpha=1,no_cores=3,FSelect=FALSE)

#plot background
plotBackVar(res)

y <- createKnnMatrix(res,pvalue=0.01)

#run clustering
cl <- graphCluster(res,pvalue=0.01)

#update raceid object
sc <- updateSC(sc,res=res,cl=cl)

#calculate tsne and umap
sc <- comptsne(sc)
sc <- compumap(sc, dimRed = T)

#save sc
save(sc, file = "sc-VarID.RData")

#plot transition probabilities
probs <-transitionProbs(res,cl,pvalue=0.01) 
plotTrProbs(sc,probs,tp=.5,prthr=0.03,cthr=0, um=T)

#plot marker genes
plotexpmap(sc,"MRC1", logsc=F,fr=F, um=T)
plotexpmap(sc,"LYVE1", logsc=F,fr=F, um=T)
plotexpmap(sc,"CD163", logsc=F,fr=F, um=T)
plotexpmap(sc,"TMEM119", logsc=F,fr=F, um=T)
plotexpmap(sc,"CX3CR1", logsc=F,fr=F, um=T)
plotexpmap(sc,"HIF1A", logsc=T,fr=F, um=T)
plotexpmap(sc,"PTPRC", logsc=F,fr=F, um=T)
plotexpmap(sc,"CD3E", logsc=F,fr=F, um=T)
plotexpmap(sc,"ITGAM", logsc=F,fr=F, um=T)
plotexpmap(sc,"CD8A", logsc=F,fr=F, um=T)
plotexpmap(sc,"CD4", logsc=F,fr=F, um=T)
plotexpmap(sc,"HLA-DRA", logsc=F,fr=F, um=T)
plotexpmap(sc,"ZBTB46", logsc=F,fr=F, um=T)
plotexpmap(sc,"IRF8", logsc=F,fr=F, um=T)
plotexpmap(sc,"FLT3", logsc=F,fr=F, um=T)
plotexpmap(sc,"IGKC", logsc=F,fr=F, um=T)
plotexpmap(sc,"VEGFA", logsc=F,fr=F, um=T)
plotexpmap(sc,"CST3", logsc=F,fr=F, um=T)
plotexpmap(sc,"NKG7", logsc=F,fr=F, um=T)
#plotexpmap(sc,"HEPACAM", logsc=F,fr=F, um=T)
plotexpmap(sc,"PLP1", logsc=F,fr=F, um=T)
plotexpmap(sc,"GFAP", logsc=F,fr=F, um=T)
plotexpmap(sc,"AQP4", logsc=F,fr=F, um=T)
plotexpmap(sc,"CALM1", logsc=F,fr=F, um=T)
plotexpmap(sc,"STAB1", logsc=F,fr=F, um=T)
plotexpmap(sc,"CD1C", logsc=F,fr=F, um=T)
plotexpmap(sc,"VIM", logsc=F,fr=F, um=T)
plotexpmap(sc,"S100A9", logsc=F,fr=F, um=T)
plotexpmap(sc,"SIGLEC8", logsc=F,fr=F, um=T)
plotexpmap(sc,"SPP1", logsc=F,fr=F, um=T)

