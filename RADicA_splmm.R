# RADicA classification with elastic net LMM

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(forcats)
library(devtools)

install.packages('splmm')

library(splmm)

b1_corr_w1 <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b_corr_w1 <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]



# set columns in the same order
b1_corr_w1 <- b1_corr_w1[names(b_corr_w1)]

b_corr_w1 <- b_corr_w1 %>% left_join(clin_dyn %>% dplyr::select(RAD_ID, CoreVisit, FeNO)) %>%
  relocate(FeNO)

b_corr_w1$Diagnosis <- ifelse(b_corr_w1$Diagnosis == 'Asthma', 1, 0)
b_corr_w1$CoreVisit <- ifelse(b_corr_w1$CoreVisit == 'CV1', 0, 1)

ids <- data.frame(RAD_ID = unique(b_corr_w1$RAD_ID), ID = c(1:length(unique(b_corr_w1$RAD_ID))))

b_corr_w1 <- b_corr_w1 %>% left_join(ids) %>%
  relocate(ID) %>% dplyr::select(!RAD_ID)

b_corr_w1 <- drop_na(b_corr_w1)

assay <- cbind(Int = rep(1, 98), b_corr_w1[,c(6:ncol(b_corr_w1))])

model <- splmm(x = assay,
               y = b_corr_w1$FeNO,
               z = as.matrix(b_corr_w1$ID),
               penalty.b = 'lasso',
               penalty.L = 'lasso',
               grp = rep(1, 98),
               lam1 = 0.05,
               lam2 = 0)
plot(model$fitted.values, b_corr_w1$FeNO)

