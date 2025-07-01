## Summarising breath VOC measurements

# author: Aggie Turlo
# project: RADicA
# date: 30/06/2025

#####################

library(tidyverse)
library(ggplot2)
library(tibble)
library(gridExtra)
library(stringr)
library(forcats)
library(grid)
library(RColorBrewer)
library(lme4)
library(lmerTest)
library(irr)
library(insight)

## load data
# load imputed normalised study data
b1_norm <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN.csv', check.names = FALSE)[,-1]
b2_norm <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN.csv', check.names = FALSE)[,-1]

# custom functions
'%ni%' <- Negate('%in%')

# keep only variables represented in both datasets
diff_voc2 <- setdiff(colnames(b2_norm), colnames(b1_norm))
diff_voc1 <- setdiff(colnames(b1_norm), colnames(b2_norm))

b1_norm <- b1_norm %>% subset(select = colnames(. ) %ni% c(diff_voc1, diff_voc2))
b2_norm <- b2_norm %>% subset(select = colnames(. ) %ni% c(diff_voc1, diff_voc2))

# remove internal standard
b1_norm <- b1_norm %>% dplyr::select(!Internal_Standard)
b2_norm <- b2_norm %>% dplyr::select(!Internal_Standard)

# load metadata
meta <- read.csv('RADicA_VOC_metadata.csv')

################################

## format datasets
# specify dataset for analysis
b_norm <- b1_norm # b2_norm
dat <- 'B1' # 'B2

#
#
#

# change format
b_norm_L <- b_norm %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea')
#b_norm_L$Analysis_date <- as.Date(b_norm_L$Analysis_date, format = '%d/%m/%Y')

b_imp_c <- b_norm_L %>% filter(class %in% c('BG','S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, ID),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  mutate(Sample = paste(ID, CoreVisit, sep = '_'))

b_imp_c1 <- b_imp_c %>%
  pivot_wider(id_cols = c(ID, comp, Analysis_date, CoreVisit, Batch, Sample), 
              names_from = class, values_from = peakArea) 
#
#
#

################################

## MEASURE DISPERSION BETWEEN TECHNICAL REPLICATES
# Intraclass correlation coefficient for each compound (ICC)
# log transform peak area data (assumption of normality)
b_imp_c2 <- b_imp_c1 %>% mutate(S1 = log(S1), S2 = log(S2))

# taking into account clustering by patient 
# using linear mixed effect model
gauss1_icc <- bind_rows(lapply(unique(b_imp_c2$comp), function(voc){
  test <- b_imp_c2 %>% filter(comp == voc) %>%
    pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')
  
  mod <- lmer(peakArea ~ Replicate + (1 | ID) + (1| ID:CoreVisit),
              data = test)
  mod1 <- lmer(peakArea ~ Replicate*CoreVisit + (1 | ID),
               data = test)
  
  vcov <- as.data.frame(VarCorr(mod), comp = 'Variance')[,c(1,4)]
  vcov1 <- as.data.frame(VarCorr(mod1), comp = 'Variance')[,c(1,4)]
  k <- 2 # number of measurements per sample
  varR1 <- vcov[1,2] # extracting random effect variance
  varR2 <- vcov[2,2]
  varR <- varR1 + varR2 # adding variance related to random effects
  varE <- vcov[3,2] # unexplained (error) variance
  varF <- get_variance_fixed(mod) # extracting fixed effect variance
  
  var1R1 <- vcov1[1,2]
  var1E <- vcov1[2,2]
  
  output <- data.frame(icc_consistency = ((k*varR)/((k*varR) + varE + ((k-1)*varE))),
                       icc_agreement = ((k*varR)/((k*varR) + varE + ((k-1)*varE) + 
                                                    (k/n_distinct(b_imp_c2$Sample))*
                                                    (n_distinct(b_imp_c2$Sample)*varF))),
                       icc_cvgen = ((varR2)/((varR2) + varE + varR1)),
                       icc_cvgen_consistency = ((k*var1R1)/((k*var1R1) + var1E + ((k-1)*var1E))),
                       comp = voc,
                       varR1 = varR1,
                       varR2 = varR2)
}))

#

# compare different ICC types
plot(gauss1_icc$icc_agreement, gauss1_icc$icc_consistency)
plot(gauss1_icc$icc_consistency, gauss1_icc$icc_cvgen_consistency)


#
#
#

# save dataset-specific results
b1_icc <- gauss1_icc
b1_imp_c1 <- b_imp_c1

# 

b2_icc <- gauss1_icc
b2_imp_c1 <- b_imp_c1

write.csv(b1_icc, 'ICC_breath_VOC_B1.csv')
write.csv(b2_icc, 'ICC_breath_VOC_B2.csv')

#
#
#

# plot distribution of ICC for all VOCs in the datset 
tiff(paste('Histogram_ICC_B1' , dat, '.tiff', sep = ''), res = 300, unit = 'mm', width = 110, height = 90)
hist(b1_icc$icc_agreement, xlab = 'ICC (agreement)', 
     main = 'Validation dataset', ylim = c(0,60), cex.main = 1)
abline(v = 0.5, col = 'red')
dev.off()

#

hist(b2_icc$icc_agreement, xlab = 'ICC (agreement)', 
     main = 'Training dataset', ylim = c(0,60), cex.main = 1)
abline(v = 0.5, col = 'red')

#################################

# FILTER VOCs ACCORDING TO ICC
# compare ICC values for the same VOCs across the datasets
icc <- b1_icc %>% left_join(b2_icc, by = 'comp')

icc_b1_f <- b1_icc %>% filter(icc_agreement > 0.5)
icc_b2_f <- b2_icc %>% filter(icc_agreement > 0.5)

nrow(icc_b1_f)
nrow(icc_b2_f)

both_pass <- intersect(icc_b1_f$comp, icc_b2_f$comp)
b1_pass <- setdiff(icc_b1_f$comp, icc_b2_f$comp)
b2_pass <- setdiff(icc_b2_f$comp, icc_b1_f$comp)
none_pass <- b1_icc$comp[b1_icc$comp %ni% c(both_pass, b1_pass, b2_pass)]
cor(icc$icc_agreement.x, icc$icc_agreement.y)

tiff('ICC_compare_B1_B2.tiff', res = 300, units = 'mm', width = 140, height = 120)
plot(
  icc$icc_agreement.x, icc$icc_agreement.y, 
  main = 'ICC agreement in Volatile Organic Compound breath replicates',
  xlab = 'Trainig dataset', ylab = 'Validation dataset', xlim = c(0,1), ylim = c(0,1),
  cex.main = 1) + 
  abline(0,1) + abline(v = 0.5, col = 'blue', lty = 2) + 
  abline(h = 0.5, col = 'blue', lty = 2) +
  text(y = 0.55, x = 0.99, label = length(both_pass), col = 'blue') +
  text(y = 0.55, x = 0.05, label = length(b1_pass), col = 'blue') +
  text(y = 0.01, x = 0.99, label = length(b2_pass), col = 'blue') +
  text(y = 0.01, x = 0.05, label = length(none_pass), col = 'blue')
dev.off()

#

# keep VOCs passing 0.5 ICC threshold in both datasets
b1_imp_c1 <- b1_imp_c1 %>% filter(comp %in% both_pass)
b2_imp_c1 <- b2_imp_c1 %>% filter(comp %in% both_pass)

#
#
#

################################

# SUMMARISE BREATH SAMPLES
b1_imp_c1 <- b1_imp_c1 %>% mutate(S = (S1+S2)/2) %>%
  dplyr::select(!c(S1, S2)) 

b2_imp_c1 <- b2_imp_c1 %>% mutate(S = (S1+S2)/2) %>%
  dplyr::select(!c(S1, S2)) 

# save output
write.csv(b1_imp_c1, 'RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised.csv')
write.csv(b2_imp_c1, 'RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised.csv')

#
#
#

###############################

# Figure S6
tiff('figS6.tiff', res = 300, unit = 'mm', width = 145, height = 100)

#
pdf('FigureS6.pdf', width = 5.7, height = 3.9)

#
s6 <- layout(matrix(c(1, 2, 3, 3), ncol = 2), widths = c(1,2))

hist(b2_icc$icc_agreement, xlab = 'ICC (agreement)', 
     main = 'Training dataset', ylim = c(0,60), cex.main = 1)
abline(v = 0.5, col = 'red')
text(xpd = NA, 'a)', x = 0, y = 135, pos = 2, offset = 2.5, cex = 1.5, font = 2)

hist(b1_icc$icc_agreement, xlab = 'ICC (agreement)', 
     main = 'Validation dataset', cex.main = 1)
abline(v = 0.5, col = 'red')

plot(
  icc$icc_agreement.x, icc$icc_agreement.y, 
  main = 'ICC agreement in VOC breath datasets',
  xlab = 'Trainig dataset', ylab = 'Validation dataset', xlim = c(0,1), ylim = c(0,1),
  cex.main = 1) + 
  abline(0,1) + abline(v = 0.5, col = 'blue', lty = 2) + 
  abline(h = 0.5, col = 'blue', lty = 2) +
  text(y = 0.55, x = 0.99, label = length(both_pass), col = 'blue') +
  text(y = 0.55, x = 0.05, label = length(b1_pass), col = 'blue') +
  text(y = 0.01, x = 0.99, label = length(b2_pass), col = 'blue') +
  text(y = 0.01, x = 0.05, label = length(none_pass), col = 'blue')
text(xpd = NA, 'b)', x = 0, y = 1.26, pos = 2, offset = 3, cex = 1.5, font = 2)
dev.off()

################################

#
#
#
