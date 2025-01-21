## MULTI-GROUP PARTIAL LEAST SQUARES PERMUTATION TEST

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(forcats)
library(foreign)
library(stringr)

# load normalised bg-corrected data w/o multivariate outliers
b1_corr_w1 <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b_corr_w1 <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]


# 
out <- 'FeNO'

perm_data <- b_corr_w1 %>%
  left_join(clin_dyn_imp_b2 %>% 
              dplyr::select(Sample, out)) %>% relocate(out)

rownames(perm_data) <- perm_data$Sample

set.seed(1410)

cv_perm <- 
  replicate(n = 500, expr = {
    
    labs <- unique(perm_data$Sample)
    x <- sample(labs, length(labs), replace = FALSE)
    perm_data_shuff <- cbind(perm_data %>% dplyr::select(out),
                            data.frame(Sample = x) %>%
                              left_join(perm_data %>% dplyr::select(!out))) %>%
      mutate(study = ifelse(CoreVisit == 'CV1', 1 ,2)) %>%
      relocate(out, study) %>%
      rename(outcome = out) 
  
      
    model_perm <- mint.pls(X = perm_data_shuff[,-c(1:6)], 
                           Y = perm_data_shuff$outcome, 
                           study = as.factor(perm_data_shuff$study), 
                           ncomp = 1, scale = TRUE)
      
    model_perm_pred <- predict(model_perm,
                               newdata = perm_data_shuff[,-c(1:6)],
                               study.test = as.factor(perm_data_shuff$study),
                               scale = TRUE)
      
    prediction_perm <- prediction_fun(model_perm_pred, perm_data_shuff, 'outcome')
    
    mae_perm <- mae_fun(prediction_perm)
    mape_perm <- mape_fun(prediction_perm)
    me_perm <- me_fun(prediction_perm)
    cor_perm <- cor(prediction_perm$outcome, prediction_perm$Y.dim1, use = 'complete.obs')
    out <- c(mae_perm, mape_perm, me_perm, cor_perm) 
    
   })

par(mfrow = c(2,2))
    
hist(cv_perm[1,], main = 'mae')
hist(cv_perm[2,], main = 'mape')
hist(cv_perm[3,], main = 'me')
hist(cv_perm[4,], main = 'cor')

#

mae <- c(mean(cv_perf[1,], na.rm = TRUE), sd(cv_perf[1,], na.rm = TRUE))
mape <- c(mean(cv_perf[2,], na.rm = TRUE), sd(cv_perf[2,], na.rm = TRUE))
me <- c(mean(cv_perf[3,], na.rm = TRUE), sd(cv_perf[3,], na.rm = TRUE))
cor <- c(mean(cv_perf[4,], na.rm = TRUE), sd(cv_perf[4,], na.rm = TRUE))

mae
mape
me
cor