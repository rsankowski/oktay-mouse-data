#cluster-wise analysis
library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)
library(MASS)
library(rebus)
library(emmeans)

date = Sys.Date()
load("data/seurat-with-batch-effect-removal.RData")
source("bin/functions.R")


#build data frame
load("data/df.Robj")

#clusterwise gene expression

#analyze microglia cluster

data_t <- as.data.frame(t(FetchData(object = all,  slot = "data")))
data_t$ID <- rownames(data_t)
data_t <- data_t %>% left_join(df) %>% na.omit()

#plot genes
up_genes <- load_data(file.path('data/Cluster specific genes/Up'))
#genes <- as.character(unique(up_genes$GENEID))
genes <- c("Ctsb", "Ctsd", "P2ry12")

for (gene in genes) {
  plt <- ggplot(data_t, aes(Region, log2(data_t[[gene]]), fill = Region)) +
    geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
    geom_jitter(pch=21, width = 0.3, alpha=0.7, size=3) +
    scale_fill_brewer(palette = "Set3") +
    theme_minimal() +
    #scale_y_log10() +
    facet_wrap(~Cluster, ncol=13, scales = "free_x") +
    labs(y='log2(normalized counts)', title=gene) 
  print(plt)
  ggsave(paste0('plots/others/', gene,'-violin-plot.pdf'))
}


signature_genes <-  read_excel('/home/roman/Documents/Single cell analysis/EAE_Final/Cluster-information-dimRed.xlsx', 'Core signature',skip = 2)

terms <- colnames(signature_genes)
stats <- list()

for (i in terms) {
  data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[[i]])])
  data_t3 <- data_t[,c("Cluster", "Condition")] %>%
    bind_cols(data.frame("gene" = data_t2))
  plt <- ggplot(na.omit(data_t3), aes(Cluster, log2(gene), fill = Condition)) +
    geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
    scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    labs(y='log2(normalized counts)', title=i) 
  print(plt)
  ggsave(paste0('plots/others/', date,"-",i ,'-microglia-cluster-Condition.pdf'))
  
  #from urls: https://www.researchgate.net/post/Should_I_use_post-hoc_tukey_HSD_for_pairwise_comparisons_of_a_factor_on_a_zero-inflated_negative_binomial_mixed_models_ZINB
  #url:https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
  mod <- glm.nb(gene ~ Cluster * Condition, data = na.omit(data_t3))
  
  emms <- emmeans(mod, ~ Cluster : Condition)
  comps <- pairs(emms, adjust="tukey") 
  stats[[i]] <- as.data.frame(comps)
  
}

stats_all <- bind_rows(stats, .id = "Signature")

write_csv(stats_all, "data/statistical testing for APC,microglia and ahr MICROGLIA comparisons.csv")

stats2 <- stats_all %>% 
  filter(Signature %in% c("ahr", "APCs", "Microglia"), p.value < 0.05) %>%
  filter(grepl("RH", contrast) & grepl("WT", contrast))

stats2$comparison <- gsub(",| - ", " ", stats2$contrast)
  
pattern = START %R% capture(one_or_more(DGT)) %R% SPC %R% one_or_more(WRD) %R% SPC %R% capture(one_or_more(DGT))

ind <- str_match(stats2$comparison, pattern=pattern)
ind <- which(ind[,2] == ind[,3])

write_csv(stats2[ind,], "data/statistical testing for MICROGLIA APC, microglia and ahr signatures significant within cluster comparisons Condition.csv")

View(stats2[ind,])

stats1 <- stats_all %>% 
  filter(Signature %in% c("ahr", "APCs"), p.value < 0.05)

write_csv(stats1, "data/statistical testing for MICROLGIA APC, microglia and ahr signatures across cluster comparisons.csv")

#plot singnificant clusters
data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["APCs"]])])
data_t3 <- data_t[,c("Cluster", "Condition")] %>%
  bind_cols(data.frame("gene" = data_t2))
plt <- ggplot(na.omit(data_t3[data_t3$Cluster %in% c(13,16,32,4,14),]), aes(Cluster, log2(gene), fill = Condition)) +
  geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
  #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="APC signature") 
print(plt)
ggsave(paste0('plots/others/MICROGLIA-APC-signature-violin-plot-significant-cluster-Condition.pdf'))

data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["ahr"]])])
data_t3 <- data_t[,c("Cluster", "Condition")] %>%
  bind_cols(data.frame("gene" = data_t2))
plt <- ggplot(na.omit(data_t3[data_t3$Cluster %in% c(4),]), aes(Cluster, log2(gene), fill = Condition)) +
  geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
  #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="AHR signature") 
print(plt)
ggsave(paste0('plots/others/MICROGLIA-AHR-signature-violin-plot-significant-cluster-Condition.pdf'))

data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["Microglia"]])])
data_t3 <- data_t[,c("Cluster", "Condition")] %>%
  bind_cols(data.frame("gene" = data_t2))
plt <- ggplot(na.omit(data_t3[data_t3$Cluster %in% c(4),]), aes(Cluster, log2(gene), fill = Condition)) +
  geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
  #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(y='log2(normalized counts)', title="Homeostatic signature") 
print(plt)
ggsave(paste0('plots/others/MICROGLIA-Homeostatic-signature-violin-plot-significant-cluster-Condition.pdf'))

#analyze MONOCYTES cluster
          data_t <- as.data.frame(t(as.matrix(sc@ndata)*min(sc@counts)))
          data_t$ID <- rownames(data_t)
          data_t <- data_t %>% left_join(df) %>%
            filter(Cell_type_simpl == "Mo_der_Cells")
          
          signature_genes <-  read_excel('/home/roman/Documents/Single cell analysis/EAE_Final/Cluster-information-dimRed.xlsx', 'Core signature',skip = 2)
          signature_genes <- signature_genes %>% 
            bind_cols(data.frame("ahr"=c("Ahr", "Ahrr", "Cyp1a1", "Cyp1b1", "Entpd1", "Hif1a", "Il1b","Il27", "Klf4", "Pparg", "Stat1", "Stat3", "Tiparp", "Vegfa", rep(NA, 7)), stringsAsFactors = F))
          
          #terms <- colnames(signature_genes)
          terms <- c("APCs", "ahr")
          stats <- list()
          
          for (i in terms) {
            data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[[i]])])
            data_t3 <- data_t[,c("Cluster", "Condition")] %>%
              bind_cols(data.frame("gene" = data_t2))
            plt <- ggplot(na.omit(data_t3), aes(Cluster, log2(gene), fill = Condition)) +
              geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
              scale_fill_brewer(palette = "Set1") +
              theme_minimal() +
              labs(y='log2(normalized counts)', title=i) 
            print(plt)
            ggsave(paste0('plots/others/', date,"-",i ,'-MONOCYTES-cluster-Condition.pdf'))
            
            #from urls: https://www.researchgate.net/post/Should_I_use_post-hoc_tukey_HSD_for_pairwise_comparisons_of_a_factor_on_a_zero-inflated_negative_binomial_mixed_models_ZINB
            #url:https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
            mod <- glm.nb(gene ~ Cluster * Condition, data = na.omit(data_t3))
            
            emms <- emmeans(mod, ~ Cluster : Condition)
            comps <- pairs(emms, adjust="tukey") 
            stats[[i]] <- as.data.frame(comps)
          }
          
          stats_all <- bind_rows(stats, .id = "Signature")
          
          write_csv(stats_all, "data/MONOCYTES statistical testing for APC and ahr comparisons.csv")
          
          stats2 <- stats_all %>% 
            filter(Signature %in% c("ahr", "APCs", "Microglia"), p.value < 0.05) %>%
            filter(grepl("RH", contrast) & grepl("WT", contrast))
          
          stats2$comparison <- gsub(",| - ", " ", stats2$contrast)
          
          pattern = START %R% capture(one_or_more(DGT)) %R% SPC %R% one_or_more(WRD) %R% SPC %R% capture(one_or_more(DGT))
          
          ind <- str_match(stats2$comparison, pattern=pattern)
          ind <- which(ind[,2] == ind[,3])
          
          View(stats2[ind,])
          
          write_csv(stats2[ind,], "data/statistical testing for MONOCYTES APC and ahr signatures significant within cluster comparisons Condition.csv")
          
          stats1 <- stats_all %>% 
            filter(Signature %in% c("ahr", "APCs"),p.value < 0.05)
          
          write_csv(stats1, "data/MONOCYTES statistical testing for APC and ahr signatures across cluster comparisons.csv")
          
          #plot singnificant clusters
          data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["APCs"]])])
          data_t3 <- data_t[,c("Cluster", "Condition")] %>%
            bind_cols(data.frame("gene" = data_t2))
          plt <- ggplot(na.omit(data_t3[data_t3$Cluster %in% c(1,15),]), aes(Cluster, log2(gene), fill = Condition)) +
            geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
            #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
            scale_fill_brewer(palette = "Set1") +
            theme_minimal() +
            labs(y='log2(normalized counts)', title="APC signature") 
          print(plt)
          ggsave(paste0('plots/others/MONOCYTES-APC-signature-violin-plot-significant-cluster-Condition.pdf'))
          
          data_t2 <- rowSums(data_t[, colnames(data_t) %in% na.omit(signature_genes[["ahr"]])])
          data_t3 <- data_t[,c("Cluster", "Condition")] %>%
            bind_cols(data.frame("gene" = data_t2))
          plt <- ggplot(na.omit(data_t3[data_t3$Cluster %in% c(11,5,15,9,29),]), aes(Cluster, log2(gene), fill = Condition)) +
            geom_violin(scale = 'width', draw_quantiles = c(0.5)) +
            #geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(0.9))+
            scale_fill_brewer(palette = "Set1") +
            theme_minimal() +
            labs(y='log2(normalized counts)', title="AHR signature") 
          print(plt)
          ggsave(paste0('plots/others/MONOCYTES-AHR-signature-violin-plot-significant-cluster-Condition.pdf'))
          
          

#
ggplot(na.omit(data_t), aes(Cluster, VEGFA__chr6, fill = Condition)) +
  geom_violin(scale = 'width', lwd=0.25) +
  geom_boxplot(width=0.5, outlier.shape = NA, position=position_dodge(1))+
  scale_fill_manual(values = toupper(c('#f1a340','#998ec3'))) +
  theme_minimal() +
  labs(y='Gene expression', title='VEGFA') 

ggsave(paste0('plots/', date, 'VEGFA-cluster-Condition.pdf'))


mod <- glm.nb(VEGFA__chr6 ~ Cluster * Condition, data = na.omit(data_t))
summary(mod)
mod_aov <- aov(mod)
write.csv(broom::tidy(TukeyHSD(mod_aov)), paste0(date, '-stat-analysis-glm-nb-GM-WM-VEGFA.csv'))

#plot
