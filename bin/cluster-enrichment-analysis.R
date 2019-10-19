#cluster enrichment analysis
library(tidyverse)
library(viridis)
library(readxl)
library(RaceID)
library(RColorBrewer)
library(pheatmap)

date = Sys.Date()
load('data/sc.Robj')
source("bin/functions.R")

#df I created in the plotting function
hyper_test_seurat(df, var2 = "Genotype") %>%
  write_csv('data/hypergeometric-test-Genotypes.csv')

hyper_test_seurat(df, var2 = "Gender") %>%
  write_csv('data/hypergeometric-test-Genders.csv')
