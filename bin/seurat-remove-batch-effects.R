#following this vignette url: https://satijalab.org/seurat/v3.0/sctransform_vignette.html
library(Seurat)
library(tidyverse)
library(clustree)
library(clusterProfiler)
library(org.Mm.eg.db)
library(Matrix)

all <- read.delim("data/counts/gbm_bcnorm.dat.txt", row.names = 1)

#rename genes
genes <- bitr(rownames(all), fromType = "ENSEMBL",
                   toType = c("ENSEMBL", "SYMBOL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db)

all <- all[genes$ENSEMBL,]
all <- Matrix(as.matrix(all), sparse=T)
rownames(all) <- genes$SYMBOL
all <- all[!duplicated(rownames(all)),]

#create seurat object
all <- CreateSeuratObject(all)

#mt genes
#all[["percent.mt"]] <-PercentageFeatureSet(all,pattern="^mt-")
#all[["percent.ribosomal"]] <-PercentageFeatureSet(all,pattern="RPS26")

#subset cells
all<-subset(all,subset=nFeature_RNA>350 & nFeature_RNA < 3000)

#assign batches and genders
all <- SCTransform(all) #, "percent.ribosomal", "gender", 

#run pca
all<-RunPCA(all,features=grep("^(Hsp|Jun|Fos|Zpf36|Neat1|Gm)", VariableFeatures(all), value=T, invert = T)) #features=grep("^(RPS|RPL|FOS|JUN|MTRNR|MT-|XIST|DUSP|ZFP36|RGS1|PMAIP1|HSP|NEAT1)", VariableFeatures(all), value=T, invert = T)

#
ElbowPlot(all)
all<-RunUMAP(all,dims=1:15)
all<-FindNeighbors(all,dims=1:15)

#find optimal resolution
all<-FindClusters(all,resolution=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,1))

#url https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html#seurat-objects
clustree(all)

pdf("plots/others/overview-cluster-resolutions.pdf")
clustree(all)
dev.off()

all<-FindClusters(all,resolution=.5)


DimPlot(all, label = TRUE) + NoLegend()

save(all, file = "data/seurat.RData")

#find cluster markers
all.markers<-FindAllMarkers(all,only.pos=T,min.pct=.25,logfc.threshold=.25,return.thresh = 0.05)

save(all.markers, file = "data/diffgenes.RData")
write_csv(all.markers, "data/diffgenes.csv")

#write_csv(data.frame(ID=colnames(all)[all@active.ident %in% c(7)]), "data/more-myeloid-cells.csv")
#write_csv(data.frame(ID=colnames(all)[which(all@assays$SCT@counts["ITGAM",]>0.2)]), "data/more-more-myeloid-cells.csv")
#write_csv(data.frame(ID=colnames(all)[which(all@assays$SCT@counts["CD79A",]>0.5)]), "data/b-cells.csv")
#write_csv(data.frame(ID=colnames(all)[which(all@assays$SCT@counts["GZMA",]>0.1 & all@assays$SCT@counts["CCL5",]>0.1)]), "data/nkt-cells.csv")
#write_csv(data.frame(ID=colnames(all)[which(all@assays$SCT@counts["CCL5",]>1)]), "data/more-nkt-cells.csv")
#write_csv(data.frame(ID=colnames(all)[all@active.ident %in% c(9)]), "data/outlier-cells.csv")
#write_csv(data.frame(ID=colnames(all)[all@active.ident %in% c(11)]), "data/more-outlier-cells.csv")
write_csv(data.frame(ID=colnames(all)[all@active.ident %in% c(8)]), "data/more-more-more-myeloid-cells.csv")
