# feature selection through rank-products algorithm

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RankProd")

library(RankProd)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggfortify)

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')
b1_all_f <- read.csv('RADicA_B1_NAfiltered_v2.csv')
b1_imp_L <- read.csv('RADicA_B1_NAfiltered_imputed_long.csv')

# change format of imputed dataset to wide
b1_imp <- b1_imp_L %>% dplyr::select(!c(logPeakAreaImp, logPeakArea, X)) %>%
  pivot_wider(names_from = comp, values_from = peakAreaImp)

rownames(b1_imp) <- b1_imp$Sample

# use RP to calculate binary contrasts of each class versus sample blanks 
# create separate datasets containing blank vs class of interest

# Blank vs BG
b_vs_bg <- b1_imp %>% filter(class %in% c('Blank', 'BG')) %>%
  dplyr::select(c(class, 4:230))
  
class1 <- b_vs_bg$class
class1 <- ifelse(class1 == 'Blank', 0, 1)

b_vs_bg <- t(b_vs_bg[, -1])

rp_b_bg <- RankProducts(data = b_vs_bg, 
                        cl = class1, 
                        logged = FALSE, 
                        na.rm = TRUE,
                        rand = 198)

### CANNOT DEAL WITH THIS SAMPLE SIZE
