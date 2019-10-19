library(tidyverse)
library(Seurat)

source('bin/functions.R')

load("data/seurat-with-batch-effect-removal-rps26-regressed-out.RData")
date = Sys.Date()

#build df
if (!file.exists("data/df.RData")) {
  df <- data.frame(all@reductions$umap@cell.embeddings,
                   Cluster=all@active.ident)
  df$ID <- rownames(df)
  df$Genotype <- factor(ifelse(grepl("^c", df$ID), "WT", "KO"), levels = c("WT", "KO"))
  df$Gender <- factor(ifelse(grepl("^(c3|p3)", df$ID), "Male", "Female"), levels = c("Female", "Male"))
  
  if (!file.exists("data/retain_cl.RData")) {
  retain_cl <- names(table(df$Cluster))[table(df$Cluster) > dim(all)[2]/100]
  save(retain_cl, file = "data/retain_cl.RData")
  } else {
    load("data/retain_cl.RData")
  }
  
  df$Group <- ifelse(grepl("^(c1|p1)", df$ID), "Group 1", 
                     ifelse(grepl("^(c2|p2)", df$ID), "Group 2", "Group 3"))
  
  df$Group <- factor(df$Group, levels = c("Group 1","Group 2","Group 3"))
  
  df$Sample <- ifelse(grepl("^c1", df$ID), "c1", 
                      ifelse(grepl("^p1", df$ID), "p1", 
                             ifelse(grepl("^c2", df$ID), "c2", 
                                    ifelse(grepl("^c3", df$ID), "c3", 
                                           ifelse(grepl("^p2", df$ID), "p2", "p3")))))
  
  df$Sample <- factor(df$Sample, levels = c("c2", "p2", "c3", "p3"))
  
  if (!file.exists('data/ord_clust.RData')) {
    clist <- split((as.data.frame(t(as.matrix(all@assays$SCT@counts)))), df[['Cluster']])
    clist <- lapply(clist, colMeans)
    clist <- bind_rows(clist)
    
    ord_clust <- hclust(dist(t(clist)))$order
    ord_clust <- c(0:13)[ord_clust]
    ord_clust <- ord_clust[ord_clust %in% retain_cl]
    
    save(ord_clust, file = 'data/ord_clust.RData')
  } else {
    load('data/ord_clust.RData')
  }
  ord_clust <- ord_clust[ord_clust %in% retain_cl]
  df$Cluster <- factor(df$Cluster, levels = ord_clust)
  
  df <- df[df$Cluster %in% retain_cl,]
  
  save(df, file = "data/df.RData")
} else {
  load("data/df.RData")
  load("data/retain_cl.RData")
  load("data/ord_clust.RData")
}

#plot clusters
tsne <- umap_plot_seurat(data = df, FILL = df$Cluster, fill_colors = c(colors_many, colors_pat)) +
  guides(colour = guide_legend(override.aes = list(size=5))) #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2

tsne

ggsave('plots/umap/clusters-umap.pdf')

svg('plots/umap/clusters-umap.svg', width = 8.48, height = 5.76)
tsne
dev.off()

#plot genotypes
        tsne <-umap_plot_seurat(data = df, FILL = df$Genotype, fill_colors = c(colors_pat)) +
          scale_color_brewer(palette = "Accent") +
          guides(colour = guide_legend(override.aes = list(size=5))) #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
        
        tsne
        ggsave('plots/umap/genotypes-umap.pdf')
        
        svg('plots/umap/genotypes-umap.svg', width = 8.48, height = 5.76)
        tsne
        dev.off()
        
        #marimekko genotypes
        mosaicGG2(df, X="Cluster", FILL="Genotype") +
          scale_fill_brewer(palette = "Accent")
        
        ggsave('plots/others/genotypes-marimekko.pdf')
        
        mosaicGG(df, X="Cluster", FILL="Genotype")
        
        ggsave('plots/others/genotypes-marimekko-stat.pdf')
        
#plot genders
        tsne <-umap_plot_seurat(data = df, FILL = df$Gender, fill_colors = c(colors_pat)) +
          scale_color_brewer(palette = "Set1") +
          guides(colour = guide_legend(override.aes = list(size=5))) #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
        
        tsne
        ggsave('plots/umap/Genders-umap.pdf')
        
        svg('plots/umap/Genders-umap.svg', width = 8.48, height = 5.76)
        tsne
        dev.off()
        
        #marimekko Genders
        mosaicGG2(df, X="Cluster", FILL="Gender") +
          scale_fill_brewer(palette = "Set1")
        
        ggsave('plots/others/Genders-marimekko.pdf')
        
        mosaicGG(df, X="Cluster", FILL="Gender")
        
        ggsave('plots/others/Genders-marimekko-stat.pdf')

#plot samples and groups
        tsne <-umap_plot_seurat(data = df, FILL = df$Sample, fill_colors = c(colors_pat)) +
          scale_color_brewer(palette = "Paired") +
          guides(colour = guide_legend(override.aes = list(size=5))) + #from url: https://stackoverflow.com/questions/20415963/how-to-increase-the-size-of-points-in-legend-of-ggplot2
          facet_wrap(~ Group, ncol=3)
        tsne
        
        ggsave('plots/umap/samples-by-groups-umap.pdf')
        
        svg('plots/umap/samples-by-groups-umap.svg', width = 8.48, height = 5.76)
        tsne
        dev.off()
        
        #marimekko samples
        mosaicGG2(df, X="Cluster", FILL="Sample") +
          scale_fill_brewer(palette = "Paired")
        
        ggsave('plots/others/Samples-marimekko.pdf')
        
        mosaicGG(df, X="Cluster", FILL="Sample")
        
        ggsave('plots/others/Samples-marimekko-stat.pdf')
        
#plot heatmap with markers
        all.markers <- read_csv("data/diffgenes.csv")
        
        top20 <- all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

        DoHeatmap(all, features = top20$gene) + NoLegend()
        
        ggsave("plots/heatmaps/top20-gene-heatmap.pdf", width = 30, height = 20)
        
# plot foxp3 expression
        plot_expmap_seurat("FOXP3" ,all, logsc = F)
        ggsave("plots/umap/foxp3-logsc.pdf")                
        
        plot_expmap_seurat("IL2RA" ,all, logsc = F)
        ggsave("plots/umap/cd25-logsc.pdf")                

#plot all marker genes
        genes <- all.markers$gene
        #genes <- genes[!grepl("MT-", genes)]
        for (i in genes) {
         tryCatch( {plot_expmap_seurat(i ,all, logsc = F)
          ggsave(paste0("plots/umap/single-genes/", i,".pdf")) }  )
          
        }
                