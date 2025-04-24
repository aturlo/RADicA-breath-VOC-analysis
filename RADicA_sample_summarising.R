## Summarising breath VOC measurements

# author: Aggie Turlo
# project: RADicA
# date: 21/01/2025

#####################

library(dplyr)
library(tidyr)
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
b1_norm <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN.csv')[,-1]
b2_norm <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN.csv')[,-1]

# keep only variables reoresented in both datasets
diff_voc2 <- setdiff(colnames(b2_norm), colnames(b1_norm))
diff_voc1 <- setdiff(colnames(b1_norm), colnames(b2_norm))

b1_norm <- b1_norm %>% subset(select = colnames(. ) %ni% c(diff_voc1, diff_voc2))
b2_norm <- b1_norm %>% subset(select = colnames(. ) %ni% c(diff_voc1, diff_voc2))

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

################################

## format datasets
# specify dataset for analysis
b_norm <- b1_norm
dat <- 'B1'

# change format
b_norm_L <- b_norm %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea')
b_norm_L$Analysis_date <- as.Date(b_norm_L$Analysis_date)

b_imp_c <- b_norm_L %>% filter(class %in% c('BG','S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(CoreVisit %ni% c('CV3', 'CV4')) %>% # keep CV1 and CV2 only!
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_'))

b_imp_c1 <- b_imp_c %>%
  pivot_wider(id_cols = c(RAD_ID, comp, Analysis_date, CoreVisit, Batch, Sample), 
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
  
  mod <- lmer(peakArea ~ Replicate + (1 | RAD_ID) + (1| RAD_ID:CoreVisit),
            data = test)
  mod1 <- lmer(peakArea ~ Replicate*CoreVisit + (1 | RAD_ID),
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

# model fit diagnostics
gauss_fit <- lapply(unique(b_imp_c2$comp), function(voc){
  test <- b_imp_c2 %>% filter(comp == voc) %>%
    pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')
  mod <- lmer(peakArea ~ Replicate + (1 | RAD_ID) + (1| RAD_ID:CoreVisit),
              data = test)
  plot(mod)})

pdf(paste('ICC_log_fit', dat, '.pdf', sep = ''))
marrangeGrob(gauss_fit, nrow = 2, ncol = 2)
dev.off()

# compare different ICC types
plot(gauss1_icc$icc_agreement, gauss1_icc$icc_consistency)
plot(gauss1_icc$icc_consistency, gauss1_icc$icc_cvgen_consistency)

# plot distribution of ICC for all VOCs in the datset 
tiff(paste('Histogram_ICC_' , dat, '.tiff', sep = ''), res = 300, unit = 'mm', width = 110, height = 90)
hist(join_icc$icc_consistency, xlab = 'ICC (agreement)', 
     main = 'ICC for breath sample replicates', ylim = c(0,55))
dev.off()

# filter VOCs with ICC > 0.5
View(gauss1_icc %>% filter(icc_agreement > 0.5))

#
set2 <- brewer.pal(3, 'Set2')

con <- gauss1_icc %>% ggplot() + 
  geom_histogram(aes(x = icc_consistency), fill = set2[1], bins = 12) +
  theme_bw() +
  ggtitle('Consistency') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'red') +
  geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'blue') +
  xlab('ICC')

con

agr <- gauss1_icc %>% ggplot() + 
  geom_histogram(aes(x = icc_agreement), fill = set2[2], bins = 12) +
  theme_bw() +
  ggtitle('Agreement')  +
  xlab('ICC') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'red') +
  geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'blue')

agr

comb <- gauss1_icc %>% ggplot() + 
  geom_histogram(aes(x = icc_consistency), alpha = 0.6, fill = set2[1], bins = 12) +
  geom_histogram(aes(x = icc_agreement), alpha = 0.6, fill = set2[2], bins = 12) +
  theme_bw() +
  xlab('ICC') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'red') +
  geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'blue') +
  ggtitle('Combined') 

hists <- arrangeGrob(con, agr, comb, nrow = 1)
plot(hists)
ggsave('ICC_CVrandom_comparison.tiff', hists, dpi = 300, units = 'mm', width = 150, height = 50)
ggsave('ICC_CVrandom_comparison_B2.tiff', hists, dpi = 300, units = 'mm', width = 150, height = 50)

#

# FILTER VOCs ACCORDING TO ICC
b_imp_c1 <- b_imp_c1 %>% filter(comp %in% gauss_icc_f$comp)

n_distinct(b_imp_c1$comp)


# VOCs passing 0.5 ICC threshold in both datasets
intersect(comp1_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp),
          comp_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp))

setdiff(comp_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp),
        comp1_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp))

write.csv(comp1_icc2, 'ICC_breath_VOC_B1.csv')

#
#
#

################################

# SUMMARISE BREATH SAMPLES
b_imp_c1 <- b_imp_c1 %>% mutate(S = (S1+S2)/2) %>%
  dplyr::select(!c(S1, S2)) 

b1_imp_c1 <- b1_imp_c1 %>% mutate(S = (S1+S2)/2) %>%
  dplyr::select(!c(S1, S2)) 

b_imp_c1$Analysis_date <- as.Date(b_imp_c1$Analysis_date, format = '%d/%m/%Y')

# filter VOCs with low ICC (based on B2 = train dataset)
excl <- comp_icc %>% filter(icc_consistency.x < 0.5 | icc_consistency.y < 0.5)

b_imp_c1 <- b_imp_c1 %>% filter(comp %ni% excl$comp)
b1_imp_c1 <- b1_imp_c1 %>% filter(comp %ni% excl$comp)

n_distinct(b1_imp_c1$comp)

# summarised peak areas - OUTPUT
write.csv(b_imp_c1, 'RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised.csv')
write.csv(b1_imp_c1, 'RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised.csv')

comp_icc <- comp_icc2 %>% left_join(comp1_icc2, by = 'comp')

icc_b1_f <- comp1_icc2 %>% filter(icc_consistency > 0.5)
icc_b_f <- comp_icc2 %>% filter(icc_consistency > 0.5)

nrow(icc_b1_f)
nrow(icc_b_f)
both_pass <- intersect(icc_b1_f$comp, icc_b_f$comp)
b1_pass <- setdiff(icc_b1_f$comp, icc_b_f$comp)
b_pass <- setdiff(icc_b_f$comp, icc_b1_f$comp)
cor(comp_icc$icc_consistency.x, comp_icc$icc_consistency.y)

tiff('ICC_compare_b1_b2.tiff', res = 300, units = 'mm', width = 140, height = 120)
plot(
  comp_icc$icc_consistency.x, comp_icc$icc_consistency.y, 
     main = 'ICC consistency',
     xlab = 'Post-covid (B2)', ylab = 'Pre-covid (B1)') + 
  abline(0,1) + abline(v = 0.5, col = 'blue', lty = 2) + 
  abline(h = 0.5, col = 'blue', lty = 2) +
  text(x = 0.95, y = 0.55, label = '163', col = 'blue') +
  text(y = 0.55, x = 0.2, label = '19', col = 'blue') +
  text(x = 0.95, y = 0.15, label = '30', col = 'blue') +
  text(x = 0.2, y = 0.15, label = '23', col = 'blue')
dev.off()




#
#
#

# BELOW POSSIBLY REDUNDANT
# change acquisition date in es and blank classes to Analysis Date
blanks <- b1_pqn1_L %>% filter(class %in% c('Blank', 'ES')) #%>%
  #dplyr::select(!c(Time)) %>%
  mutate(Analysis_date = substr(Sample, start = 1, stop = 6)) %>%
  mutate(Analysis_date = paste(str_sub(Analysis_date, start = 5),
                               str_sub(Analysis_date, start = 3, end = -3),
                               '20',
                               str_sub(Analysis_date, end = - 5), sep = ''))

blanks$Analysis_date <- as.Date(blanks$Analysis_date, format = '%d%m%Y')

# pivot longer and merge summarised sample values with other classes
b_imp_c <- b_imp_c %>% mutate(class = ifelse(class %in% c('S1', 'S2'), 'S', class)) %>%
  dplyr::select(Sample, class, Batch) %>% distinct()

b1_imp_c_L <- b1_imp_c1 %>% pivot_longer(cols = c(BG, S), 
                                         names_to = 'class', 
                                         values_to = 'peakArea') %>%
  left_join(b1_imp_c) #%>%
  #dplyr::select(!RAD_ID)

b1_imp_L_sum <- rbind(b1_imp_c_L, blanks)


# save the imputed summarised dataset in long format
write.csv(b1_imp_L_sum, 'RADicA_B1_NAfiltered_imputed_PQN_summarised_long.csv')

# filter out VOCs based on the ICC threshold
comp_icc_f <- gauss_icc %>% 
  filter(icc_consistency >= 0.5)


b1_imp_L_sum <- b1_imp_c_L %>% filter(comp %in% comp_icc_f$comp)
  
unique(b1_imp_L_sum$comp)

write.csv(b1_imp_L_sum, 'RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC.csv')

# visualise change in median peak area with time in summarised samples
dev.new()

b1_imp_L_sum %>% group_by(Analysis_date, class) %>%
  summarise(median = median(log(peakArea), na.rm = TRUE)) %>%
  ggplot(aes(x = Analysis_date, y = median)) +
  geom_line() +
  facet_wrap(~class, scale = 'free')

#
#
#

# evaluate instrument drift within each MS run using blanks/ES
drift <- b1_imp_L %>%
  filter(class == 'Blank') %>%
  group_by(Date, Time) %>%
  summarise(median = median(peakArea)) %>%
  group_by(Date) %>%
  do(plots = ggplot(data = ., aes(x = Time, y = median)) +
       geom_point() + ggtitle(.$Date))

marrangeGrob(drift$plots, nrow = 4, ncol = 4)

# save icc analysis results for datasets 1 and 2
write.csv(comp1_icc2, 'B1_ICC_results.csv')
write.csv(comp_icc2, 'B2_ICC_results.csv')

#
#
#

# explore data properties associated with good/bad ICC
expl_icc <- b1_imp_c2 %>%
  group_by(comp) %>%
  summarise(meanS1 = mean(S1),
            meanS2 = mean(S2),
            medS1 = median(S1),
            medS2 = median(S2)) %>%
  mutate(dataset = 'Dataset_1') %>%
  left_join(comp1_icc2) %>%
  rbind(b_imp_c2 %>% 
          group_by(comp) %>%
          summarise(meanS1 = mean(S1),
                    meanS2 = mean(S2),
                    medS1 = median(S1),
                    medS2 = median(S2)) %>%
          mutate(dataset = 'Dataset_2') %>%
          left_join(comp_icc2)) %>%
  mutate(pass = ifelse(comp %in% b1_pass, 'Dataset_1',
                       ifelse(comp %in% b_pass, 'Dataset_2',
                              ifelse(comp %in% both_pass, 'Both', 'None')))) 

expl_icc %>%
  ggplot(aes(x = icc_consistency, y = meanS1, colour = as.factor(pass))) +
  geom_point(size = 1, alpha = 0.8) +
  facet_grid(~ dataset) +
  theme_bw() +
  labs(colour = 'ICC > 0.05') +
  ylab('mean log(S1)') +
  xlab('ICC consistency') +
  ggtitle('Relationship between VOC abundance in breath and reproducibility') 

ggsave('ICC_vs_abundance.tiff', unit = 'mm', dpi = 300, width = 170, height = 80)

table(expl_icc$pass)/2

