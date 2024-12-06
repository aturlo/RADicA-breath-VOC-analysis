## Summarising breath VOC measurements

# author: Aggie Turlo
# project: RADicA
# date: 22/07/2024

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


'%ni%' <- Negate('%in%')

# load data
# imputed and normalised
b1_imp <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_noB5.csv')[,-1]
b_imp <- read.csv('RADicA_B2_NAfiltered_imputed_PQN.csv')[,-1]

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')
data <- read.csv('240508_Radica_VOC_peak_areas_v5.1.csv', check.names = FALSE)


b_imp_L <- b_imp_L %>% left_join(meta %>% dplyr::select(Sample_ID, Analysis_date),
                                   by = c('Sample' = 'Sample_ID'))
b1_imp_L <- b1_imp_L %>% left_join(meta %>% dplyr::select(Sample_ID, Analysis_date),
                                 by = c('Sample' = 'Sample_ID'))

# imputed and normalised (CC + PQN)
b1_norm <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_noB5.csv')[,-1]
b_norm <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN.csv')[,-1]

# change to long format
b_norm_L <- b_norm %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea')
b1_norm_L <- b1_norm %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea')

b_norm_L$Analysis_date <- as.Date(b_norm_L$Analysis_date)
b1_norm_L$Analysis_date <- as.Date(b1_norm_L$Analysis_date)

#
#
#
x <- b1_norm_L


# identify sample duplicates and remove from the database
# Batch 1
b1_imp_s <- x %>% pivot_wider(names_from = comp, values_from = peakArea) %>%
  filter(class %in% c('BG', 'S1', 'S2')) %>%
  mutate(Sample1 = str_sub(Sample, end = -13)) %>%
  relocate(Sample1)

View(b1_imp_s %>% 
  group_by(Sample1, class) %>% 
  mutate(dupe = n()>1) %>%
    filter(dupe == TRUE))

dupes <- c('200320_RaDICA_RAD094_B1_434919_1', '200320_RaDICA_RAD074_B1_609366_1', '200417_RaDICA_RAD114_S2_367150_1')

x <- x %>% filter(Sample %ni% dupes)

write.csv(b1_imp_L, 'RADicA_B1_NAfiltered_imputed_long.csv')

# Batch 2
b_imp_s <- b_pqn1_L %>%
  filter(class %in% c('BG', 'S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, RAD_ID),
                        by = c('Sample' = 'Sample_ID'))


b_imp_s %>% group_by(RAD_ID, CoreVisit, class, comp) %>%
  summarise(n = n()) %>%
  filter(n > 1)

#
#
#

# change format to wide (columns = classes)
x <- b1_norm_L

b1_imp_c <- x %>% filter(class %in% c('BG','S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) #%>%
  #mutate(Sample = str_sub(Sample, end = -13)) # Sample name formatting for B1


which(is.na(b1_imp_c$peakArea) == TRUE, arr.ind = TRUE)

b1_imp_c1 <- b1_imp_c %>%
  pivot_wider(id_cols = c(RAD_ID, comp, Analysis_date, CoreVisit, Batch), 
              names_from = class, values_from = peakArea) %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  filter(CoreVisit %ni% c('CV3', 'CV4'))

#

#
#
# identify breath samples with missing BG sample 
incomplete_s <- b1_imp_c1[which(is.na(b1_imp_c1$BG), arr.ind = TRUE),] #%>% 
  group_by(Sample) %>% summarise(n = n())

# B1 - exclude from the dataset (BG analysed in different Batch)
b_imp_c1 <- b_imp_c1 %>% filter(Sample %ni% incomplete_s$Sample)

#

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

#
#
#

# B2
# check if VOCs with missing BG data need to be BG corrected
model1_out <- read.csv('Model1_results.csv')[,-1]
model1_out %>% filter(predictor == 'logBG') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(comp %in% unique(incomplete_s$comp)) %>%
  filter(adj.p.value < 0.05)

# Pentanal

# find the BG VOC most correlated with Pentanal BG levels
b_all <- read.csv('RADicA_B2.csv')[,-1] # repeat for B1 and B2
b_all_bg <- b_all %>% filter(class == 'BG')

b1_all <- read.csv('RADicA_B1.csv')[,-1] # repeat for B1 and B2
b1_all_bg <- b1_all %>% filter(class == 'BG')

# missing rate for BG data
table(is.na(b_all_bg$Pentanal) == TRUE)
104/(104+82)

cor_mat <- cor(b_all_bg[,-c(1:5)], method = 'pearson', use = 'pairwise.complete.obs') %>%
  as.data.frame()
View(cor_mat %>% dplyr::select(Pentanal) %>% arrange(desc(Pentanal)) %>%
       mutate(comp = rownames(.)))

View(b_all_bg %>% dplyr::select(Pentanal, Heptane._2.2.4_trimethyl_) %>% drop_na()) 
# 47/186 paired obs Propanoic_acid
# 79/186 Heptane._2.2.4_trimethyl_

cor_mat1 <- cor(b1_all_bg[,-c(1:5)], method = 'pearson', use = 'pairwise.complete.obs') %>%
  as.data.frame()
View(cor_mat1 %>% dplyr::select(Pentanal) %>% arrange(desc(Pentanal)) %>%
       mutate(comp = rownames(.)))






#
#
#

# EVALUATE RELATIONSHIP BETWEEN TECHNICAL REPLICATES
# scatter plot of log peak areas
tr <- b_imp_c1 %>% group_by(comp) %>%
  do(plots = ggplot(data = .) + (aes(y = log(S1), x = log(S2))) + geom_point() +
       ggtitle(.$comp) + coord_fixed() +
       theme_bw(base_size = 6) + geom_abline(intercept = 0, slope = 1, colour = 'red'))

pdf('Radica_S1vsS2_PQN.pdf')
marrangeGrob(grobs = tr$plots, nrow = 4, ncol = 4)
dev.off()

# correlation matrix
corrs <- b_imp_c1 %>% dplyr::select(!BG) %>% 
  drop_na() %>%
  group_by(comp) %>% 
  summarise(cor = cor(S1, S2, method = 'pearson'))


#
#
#

# measure of dispersion between technical replicates
# coefficient of variation
b1_imp_c2 <- b1_imp_c %>% 
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  filter(class != 'BG') %>%
  group_by(Sample, comp, Analysis_date) %>%
  summarise(mean = mean(peakArea, na.rm = TRUE),
            sd = sd(peakArea, na.rm = TRUE),
            cv = sd/mean)

b1_imp_c2 <- b1_imp_c2 %>%
  left_join(b1_imp_c2 %>% 
              group_by(comp) %>%
              summarise(median = median(cv))) %>%
  arrange(median)

dev.new()

set2 <- brewer.pal(3, 'Set2')

cv_comp <- b1_imp_c2 %>% ggplot(aes(x = fct_inorder(comp), y = cv)) + 
  geom_boxplot(outliers = FALSE, lwd = 0.2, fill = set2[1]) +
  geom_hline(yintercept = 0.3, colour = 'red') +
  theme(axis.text.x = element_blank()) +
  xlab('Volatile Organic Compounds (VOCs)') +
  ylab('CV')


b1_imp_c3 <- b1_imp_c %>% filter(class != 'BG') %>%
  group_by(Sample, comp, Analysis_date) %>%
  summarise(mean = mean(peakArea, na.rm = TRUE),
            sd = sd(peakArea, na.rm = TRUE),
            cv = sd/mean) %>%
  left_join(b1_imp_c2 %>% 
              group_by(Sample) %>%
              summarise(median = median(cv))) %>%
  arrange(median)

dev.new()

cv_sample <- b1_imp_c3 %>% ggplot(aes(x = fct_inorder(Sample), y = cv)) + 
  geom_boxplot(outliers = FALSE, lwd = 0.2, fill = set2[2]) +
  geom_hline(yintercept = 0.3, colour = 'red') +
  theme(axis.text.x = element_blank()) +
  xlab('Breath samples') +
  ylab('CV')

cv_plots <- arrangeGrob(cv_comp, cv_sample, nrow = 2, top = textGrob('Coefficient of variation in breath sample replicates'))
plot(cv_plots)

ggsave('CV_plots_breath_replicates_PQN.tiff', cv_plots, dpi = 300, unit = 'mm', width = 200, height = 120)


#
#
#


# intraclass correlation coefficient for each compound
library(irr)
voc <- 'Benzonitrile'

# assumption of normality - data log transformation 
b_imp_c2 <- b_imp_c1 %>% mutate(S1 = log(S1), S2 = log(S2))

b_imp_c2$Analysis_date <- as.Date(b_imp_c2$Analysis_date)

b1_imp_c2 <- b1_imp_c1 %>% mutate(S1 = log(S1), S2 = log(S2))

b1_imp_c2$Analysis_date <- as.Date(b1_imp_c2$Analysis_date)

# IMPORTANT
# KEEP ONLY CV1 and CV2
b1_imp_c2 <- b1_imp_c2 %>% 
  left_join(meta %>% dplyr::select(RAD_ID, Analysis_date, CoreVisit)) %>%
  distinct() %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))


# works only with normally distributed data
# using R function icc
comp1_icc2 <- bind_rows(lapply(unique(b1_imp_c2$comp), function(voc){
   input <-  b1_imp_c2 %>% 
     filter(comp == voc) %>%
     dplyr::select(S1, S2, Sample) %>%
     as.data.frame()
   
  rownames(input) <- input$Sample
  input <- input[,-3]
   
   icc_res <- icc(input, model = 'twoway', type = 'agreement', unit = 'single')
   output <- data.frame(comp = voc, 
                        ICC = icc_res$value, 
                        p.value = icc_res$p.value,
                        CI_lwr = icc_res$lbound, CI_upr = icc_res$ubound)
   
 })) %>%
   mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

hist(comp1_icc2$icc_consistency, breaks = 12)

View(comp1_icc2 %>% filter(icc_consistency > 0.5))

# VOCs passing 0.5 ICC threshold in both datasets
intersect(comp1_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp),
          comp_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp))

setdiff(comp_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp),
          comp1_icc2 %>% filter(icc_consistency > 0.5) %>% dplyr::select(comp))

write.csv(comp1_icc2, 'ICC_breath_VOC_B1.csv')



# try to replicate ICC by hand
# using linear regression on transformed data
library(insight)

gauss1_icc <- bind_rows(lapply(unique(b1_imp_c2$comp), function(voc){
  test <- b1_imp_c2 %>% filter(comp == voc) %>%
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

pdf('ICC_log_fit.pdf')
marrangeGrob(gauss_fit, nrow = 2, ncol = 2)
dev.off()

library(insight)

plot(gauss_icc$icc_agreement, gauss_icc$icc_consistency)
plot(gauss_icc$icc_consistency, gauss_icc$icc_cvgen_consistency)

# compare with results not considering donor effect
comp_icc2 <- comp_icc2 %>% left_join(gauss_icc)
comp1_icc2 <- comp1_icc2 %>% left_join(gauss1_icc)
plot(comp_icc2$ICC, comp_icc2$icc_agreement)

tiff('Histogram_ICC_B2.tiff', res = 300, unit = 'mm', width = 110, height = 90)
hist(comp_icc2$icc_consistency, xlab = 'ICC (consistency)', 
     main = 'ICC for breath sample replicates in B2', ylim = c(0,55))
dev.off()

View(comp1_icc2 %>% filter(icc_agreement < 0.5))

#
set2 <- brewer.pal(3, 'Set2')

con <- gauss_icc %>% ggplot() + 
  geom_histogram(aes(x = icc_consistency), fill = set2[1], bins = 12) +
  theme_bw() +
  ggtitle('Consistency') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'red') +
  geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'blue') +
  xlab('ICC')

con

agr <- gauss_icc %>% ggplot() + 
  geom_histogram(aes(x = icc_agreement), fill = set2[2], bins = 12) +
  theme_bw() +
  ggtitle('Agreement')  +
  xlab('ICC') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.5, linetype = 'dashed', colour = 'red') +
  geom_vline(xintercept = 0.8, linetype = 'dashed', colour = 'blue')

comb <- gauss_icc %>% ggplot() + 
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


par(mfrow = c(1,2))

View(comp_icc2 %>% filter(adj.p.value < 0.05) %>% filter(icc >= 0.6))

gauss_icc_f <- gauss_icc %>% filter(icc_agreement > 0.5)


##
hists <- lapply(unique(b1_imp_c2$comp), function(voc){
  subset <- b1_imp_c2 %>%
    filter(comp == voc) %>%
    pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')
  
  plot <- subset %>% ggplot(aes(x = peakArea, fill = Replicate)) +
    geom_histogram(alpha = 0.6) + ggtitle(voc) +
    theme_bw(base_size = 8) +
    theme(legend.position = 'none')
  })

pdf('Histograms_S1S2_B2.pdf')
marrangeGrob(hists, nrow = 4, ncol = 4)
dev.off()

# FILTER VOCs ACCORDING TO ICC
b_imp_c1 <- b_imp_c1 %>% filter(comp %in% gauss_icc_f$comp)

n_distinct(b_imp_c1$comp)

#
#
#


# using generalizability theory
library(gtheory)

test <- b1_imp_c2 %>% filter(comp == voc) %>%
  pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')

g <- gstudy(data = test, 
            formula = peakArea ~ Replicate + (1 | RAD_ID),
            colname.objects = 'RAD_ID', colname.strata = 'CoreVisit')

#
test <- b1_imp_c2 %>% filter(comp == voc) %>%
  pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')

mod <- lmer(peakArea ~ Replicate*CoreVisit + (1 | RAD_ID),
            data = test)

summary(mod)
plot(mod)
var <- as.data.frame(VarCorr(mod), comp = 'Variance')[,c(1,4)]

# overall reliability (test-retest) generalising over replicate and core visit
icc <- var[1,2]/(var[1,2] + var[2,2])

# 
mod <- lmer(peakArea ~ Replicate + (1 | RAD_ID/CoreVisit),
            data = test)

summary(mod)
plot(mod)
var <- as.data.frame(VarCorr(mod), comp = 'Variance')[,c(1,4)]

# generalising over replicate and core visit
icc1 <- (var[1,2] + var[2,2])/(var[1,2] + var[2,2] + var[3,2])


#
#
#


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

