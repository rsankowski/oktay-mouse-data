library(tidyverse)
library(org.Mm.eg.db)

source('bin/functions.R')

load("data/seurat.RData")
date = Sys.Date()
    
    df <- read_csv('data/diffgenes.csv')
    
    enrich_up <- go_term_analysis_seurat(organism = 'org.Mm.eg.db')
    
    dir.create('GO-terms')
    dir.create('GO-terms/bp')
    write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms.csv')
    
    #mf terms
    enrich_up_mf <- go_term_analysis_seurat(ontogeny = 'MF', organism = 'org.Mm.eg.db')
    
    dir.create('GO-terms/mf')
    write.csv(enrich_up_mf, 'GO-terms/mf/mf_GO_terms.csv')

#plotting go terms
#bp
                enrich_up <- read_csv("GO-terms/bp/bp_GO_terms.csv")
                enrich <- enrich_up %>% na.omit 
                enrich <- enrich %>% filter(!duplicated(enrich$Description))
                
                for (i in 1:nrow(enrich)) {
                    tryCatch({
                        n=enrich[i,'Description']
                        svg(paste0('GO-terms/bp/umap/',as.character(enrich[i,'Description']),'.svg'), width = 8.57, height = 5.79) #, units = 'in', res = 300
                        try(
                            print(plot_expmap_seurat(features = c(stringr::str_split(enrich[i,'geneID'], pattern = '/')[[1]])))
                            
                        )
                        dev.off()
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    on.exit(dev.off())
                }
                
                #expanded antigen presentation signature
              #  genes <- c(unlist(stringr::str_split("HLA-DRB5/B2M/FCER1G/HLA-A/ARF1/CTSD/HLA-B/HLA-DRB1/HLA-C/HLA-DQB1/HLA-DPB1/CTSL/RAB7A/AP1B1/HLA-DPA1/HLA-DRA/CTSS/CANX/CD74", pattern = "/")),c("CD40", "CD80", "CD86"))
               # plot_expmap_seurat(features = name2id(genes, rownames(sc@ndata)), point_size = 5)
                
                #mf
                enrich_up <- read_csv("GO-terms/mf/mf_GO_terms.csv")
                enrich <- enrich_up %>% na.omit 
                enrich <- enrich %>% filter(!duplicated(enrich$Description))
                
                for (i in 1:nrow(enrich)) {
                    tryCatch({
                        n=enrich[i,'Description']
                        svg(paste0('GO-terms/mf/umap/',as.character(enrich[i,'Description']),'.svg'), width = 8.57, height = 5.79) #, units = 'in', res = 300
                        try(
                            print(plot_expmap_seurat(features = c(stringr::str_split(enrich[i,'geneID'], pattern = '/')[[1]])))
                            
                        )
                        dev.off()
                    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                    on.exit(dev.off())
                }
                
                #expanded antigen presentation signature
                genes <- c(unlist(stringr::str_split("HLA-DRB5/B2M/FCER1G/HLA-A/ARF1/CTSD/HLA-B/HLA-DRB1/HLA-C/HLA-DQB1/HLA-DPB1/CTSL/RAB7A/AP1B1/HLA-DPA1/HLA-DRA/CTSS/CANX/CD74", pattern = "/")),c("CD40", "CD80", "CD86"))
                plot_expmap_seurat(features = name2id(genes, rownames(sc@ndata)), point_size = 5)
                
                #dot plot
                #retain_cl
                cell_numbers <-as.numeric()
                for (i in 1:length(unique(sc@cpart)))
                {
                    cell_numbers[i] <- length(sc@cpart[sc@cpart==i])
                }
                names(cell_numbers) <- c(1:length(unique(sc@cpart)))
                retain_cl <- as.numeric(names(cell_numbers[cell_numbers > dim(sc@ndata)[2]/100]))
                
                enrich_up <- enrich_up[enrich_up$Cluster %in% retain_cl,]
                load('data/ord_clust.Robj')
                ord_clust <- ord_clust[ord_clust %in% retain_cl]
                enrich_up$Cluster <- factor(enrich_up$Cluster, levels = ord_clust)
                enrich_up$Description <- reorder(enrich_up$Description,enrich_up$Description,FUN=length)
                
                #re-order enrich_up based on the levels of Cluster
                enrich_up <- enrich_up[with(enrich_up, order(Cluster)),] #from url: https://stackoverflow.com/questions/1296646/how-to-sort-a-dataframe-by-columns
                enrich_up$Description <- factor(enrich_up$Description, levels = rev(enrich_up$Description[!duplicated(enrich_up$Description)]))
                colnames(enrich_up)[10] <- 'GeneCount'
                enrich_up <- enrich_up[!duplicated(enrich_up[,c("Description", "Cluster")]),]
                
                
                dot_plot <- ggplot(enrich_up, aes(Cluster, Description, size = GeneCount, fill= -log2(qvalue))) + #[enrich_up$GeneCount>4,]
                    geom_point(pch=21, stroke=0.25) +
                    scale_fill_viridis() +
                    theme_light() +
                    theme(text=element_text(size=10),
                          axis.title.y=element_blank())
                dot_plot
                
                ggsave(paste0('GO-terms/bp/',date,'-GBM_GoTerm_dot_plot.pdf'), height = 18, width = 9, units = 'in')
                
                svg(paste0('GO-terms/bp/',date,'-GBM_GoTerm_dot_plot.svg'), height = 18, width = 9) #, units = 'in', res = 300
                dot_plot #https://stackoverflow.com/questions/15678261/r-ggplot-does-not-work-if-it-is-inside-a-for-loop-although-it-works-outside-of
                dev.off()






#go term analysis of pseudo time
#wt
df <- read_csv('data/2019-05-21-nodes-stemid-vector-wt.csv', col_names = c("X1", "Cluster", "GENEID"), skip = 1)

enrich_up <- go_term_analysis(organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-trajectory-wt.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-trajectory-wt.csv')

#GO terms across all nodes
df$Cluster <- rep(1, nrow(df))

rm(enrich_up)
#analysis
enrich_up <- go_term_analysis(organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-across-all-nodes-trajectory-wt.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-across-all-nodes-trajectory-wt.csv')


#mut
df <- read_csv('data/2019-05-21-nodes-stemid-vector-mut.csv', col_names = c("X1", "Cluster", "GENEID"), skip=1)

enrich_up <- go_term_analysis(organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-trajectory-mut.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-trajectory-mut.csv')

#GO terms across all nodes
df$Cluster <- rep(1, nrow(df))

#analysis
enrich_up <- go_term_analysis(organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/bp/bp_GO_terms-across-all-nodes-trajectory-mut.csv')

#mf
enrich_up <- go_term_analysis(ontogeny = 'MF',organism = 'org.Hs.eg.db')

write.csv(enrich_up, 'GO-terms/mf/mf_GO_terms-across-all-nodes-trajectory-mut.csv')

#retain genes associated with certain go terms
#I downloaded te go term information from here: http://www.informatics.jax.org/vocab/gene_ontology/
#find genes in dataset
all_genes <- rownames(sc@ndata)[which(apply(as.matrix(sc@ndata) > 0,1,sum)>0)]
source('~/Documents/Single cell analysis/Advanced-plots/20190102_plot_expmap_seurat.R')

#before running loop please make sure to create the df object from plottin.R
for (i in list.files("data/go-terms-mirco")) {
        genes <- read_tsv(paste0("data/go-terms-mirco/", i)) %>%
          select(Symbol) %>%
          unique %>%
          unlist
        
        genes <- name2id(toupper(genes), rownames(sc@ndata))
        present_genes <- genes[genes %in% all_genes]
                 
        #plot_genes <-  as.character(unique(up_genes$GENEID))
            svg(paste0('plots/tsne/', i, '.svg'), width = 8.57, height = 5.79)
            pl <- plot_expmap_seurat(gene=present_genes, point_size = 5)
            print(pl)
            dev.off()
            
            svg(paste0('plots/tsne/', i, '-logsc.svg'), width = 8.57, height = 5.79)
            pl <- plot_expmap_seurat(gene=present_genes, point_size = 5, logsc=T)
            print(pl)
            dev.off()
            
    }
    
    