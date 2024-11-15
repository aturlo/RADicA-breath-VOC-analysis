## Background correction

# author: Aggie Turlo
# project: RADicA
# date: 23/07/2024

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(lme4)
library(lmerTest)
library(nlme)
library(MASS)
library(foreign)
library(sfsmisc)
library(robustlmm)
library(stringr)
library(insight)
library(performance)


# load normalised dataset
b1_imp_L_sum <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC.csv')[,-1] 
b1_imp_L_sum$Analysis_date <- as.Date(b1_imp_L_sum$Analysis_date)

b1_imp_sum <- b1_imp_L_sum %>% pivot_wider(names_from = comp, values_from = peakArea) %>%
  dplyr::select(!Internal_Standard)

# load normalised dataset filtered for VOCs similar in breath and BG
b1_imp_sum_c <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC_BGfiltered.csv')[,-1]

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

# format dataset
b1_imp_sum_c <- b1_imp_L_sum %>% filter(class %in% c('BG', 'S')) %>%
  pivot_wider(names_from = class, values_from = peakArea) %>%
  mutate(RAD_ID = str_sub(Sample, start = 15)) %>%
  left_join(meta %>% dplyr::select(RAD_ID, Analysis_date, Diagnosis, CoreVisit)) %>%
  distinct()
         

# plot S vs BG log abundances
bg_s_plot <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      filter(Sample %ni% c('190903_RaDICA_RAD032')) %>%
      filter(CoreVisit %in% c('CV1', 'CV2'))
    #filter(logS < mean(logS, na.rm = TRUE) + 3*sd(logS, na.rm = TRUE)) %>%
    #filter(logS > (mean(logS, na.rm = TRUE) - 3*sd(logS, na.rm = TRUE)))
  
  plot <- ggplot(data = subset, 
                 aes(x = logBG, y = logS)) + #colour = as.factor(Diagnosis))) +
       geom_point(alpha = 0.6, size = 1, aes(colour = as.factor(Diagnosis))) + 
    ggtitle(voc) +
    theme_bw() +
    theme(plot.title = element_text(size = 8),
          legend.position = 'none') +
    geom_smooth(method = 'lm', se = FALSE)
  })
  }


plots <- bg_s_plot(b1_imp_sum_c)

# only sign S-BG associations
b1_imp_sum_f <- b1_imp_sum_c %>% filter(comp %in% reg_results_pqn_sign$comp)

plots1 <- bg_s_plot(b1_imp_sum_f)
plots2 <- bg_s_plot(b1_imp_c1)

pdf('BGvsS_scatterPlots_PQN_CV12_BGfiltered_sign.pdf')
marrangeGrob(plots1, nrow = 4, ncol = 4)
dev.off()

# histograms of BG and S log abundances
bg_s_hist <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
    #filter(logS < mean(logS, na.rm = TRUE) + 3*sd(logS, na.rm = TRUE)) %>%
    #filter(logS > (mean(logS, na.rm = TRUE) - 3*sd(logS, na.rm = TRUE)))
      pivot_longer(cols = c(logBG, logS), names_to = 'class', values_to = 'peakArea')
    
    plot <- ggplot(data = subset, 
                   aes(x = peakArea, fill = as.factor(class))) + # change to logS
      geom_histogram(alpha = 0.6, position = 'identity') + ggtitle(voc) +
      theme_bw() +
      theme(plot.title = element_text(size = 8),
            legend.position = 'none') +
      xlab('log(peakArea)')
  })
}


plots1 <- bg_s_hist(b1_imp_sum_c)

pdf('BGvsS_histograms_PQN.pdf')
marrangeGrob(plots1, nrow = 4, ncol = 4)
dev.off()

#

# histograms of BG and S abundances
bg_s_hist1 <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      #filter(logS < mean(logS, na.rm = TRUE) + 3*sd(logS, na.rm = TRUE)) %>%
      #filter(logS > (mean(logS, na.rm = TRUE) - 3*sd(logS, na.rm = TRUE)))
      pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea')
    
    plot <- ggplot(data = subset, 
                   aes(x = peakArea, fill = as.factor(class))) + # change to logS
      geom_histogram(alpha = 0.6, position = 'identity') + ggtitle(voc) +
      theme_bw() +
      theme(plot.title = element_text(size = 8),
            legend.position = 'none') +
      xlab('peakArea')
  })
}


plots2 <- bg_s_hist1(b1_imp_sum_c)

pdf('BGvsS_histograms_PQN_untransf.pdf')
marrangeGrob(plots2, nrow = 4, ncol = 4)
dev.off()

#
#
#

# difference between BG and S distributions
df <- b1_imp_L_sum %>% 
  left_join(b1_imp_sum_c %>% 
              dplyr::select(Sample, Analysis_date, CoreVisit)) %>%
  distinct()

# identifying outliers
bg_s_outliers <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc){
    subset <- df %>%
      mutate(RAD_ID = str_sub(Sample, start = 15)) %>%
      left_join(meta %>% dplyr::select(RAD_ID, Analysis_date, CoreVisit)) %>%
      distinct() %>%
      #filter(CoreVisit %in% c('CV1','CV2')) %>% # optional
      filter(class %in% c('S', 'BG')) %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logPeakArea = log(peakArea))
    
    x <- subset$logPeakArea
    m <- median(x, na.rm = TRUE)
    abs.dev   <- abs(x - m)
    left.mad  <- median(abs.dev[x<=m])
    right.mad <- median(abs.dev[x>=m])
    
    subset <- subset %>% 
      mutate(mad = ifelse(peakArea > m, right.mad, left.mad),
                                mad.dist = abs(x-m)/mad) %>%
      mutate(outlier = ifelse(mad.dist > 3, 'yes', 'no'))}))}
    
outliers <- bg_s_outliers(b1_imp_L_sum)
out_sum <- outliers %>% filter(outlier == 'yes') %>% group_by(class, Sample) %>%
  summarise(n = n())
View(table(out_sum$Sample))


# linear modelling
input <- b1_imp_sum_c %>%
  pivot_longer(cols = c('BG', 'S'), names_to = 'class', values_to = 'peakArea')


bg_s_gauss <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logPeakArea = log(peakArea))
    
    #x <- subset$logPeakArea
    #m <- median(x, na.rm = TRUE)
    #abs.dev   <- abs(x - m)
    #left.mad  <- median(abs.dev[x<=m])
    #right.mad <- median(abs.dev[x>=m])
    
    #subset <- subset %>% 
      #mutate(mad = ifelse(peakArea > m, right.mad, left.mad),
             #mad.dist = abs(x-m)/mad) %>%
      #filter(mad.dist <= 3)
    
    model <- lmer(logPeakArea ~ class + (1 | RAD_ID), 
      data = subset)
    
    coefs <- as.data.frame(coef(summary(model))) %>%
      mutate(comp = rep(voc, 2),
             coef = rownames(.),
             CI_low = confint(model)[4,1],
             CI_upr = confint(model)[4,2])
    
  }))
}


# THE OUTPUT
bg_s_output_f <- bg_s_output %>% filter(coef == 'classS') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(adj.p.value > 0.05)

bg_s_output

setdiff(wilcox_bg_f$comp, bg_s_output_f$comp)
setdiff(bg_s_output_f$comp, wilcox_bg_f$comp)
intersect(bg_s_output_f$comp, wilcox_bg_f$comp)

# diagnostic plots
bg_s_gauss_plot <- function(df) {
  lapply(unique(df$comp), function(voc){
    subset <- df %>%
      mutate(RAD_ID = str_sub(Sample, start = 15)) %>%
      left_join(meta %>% dplyr::select(RAD_ID, Analysis_date, CoreVisit)) %>%
      distinct() %>%
      filter(class %in% c('S', 'BG')) %>%
      filter(comp == voc) %>%
      drop_na() %>% 
      mutate(logPeakArea = log(peakArea))
    
    #x <- subset$logPeakArea
    #m <- median(x, na.rm = TRUE)
    #abs.dev   <- abs(x - m)
    #left.mad  <- median(abs.dev[x<=m])
    #right.mad <- median(abs.dev[x>=m])
    
    #subset <- subset %>% 
      #mutate(mad = ifelse(peakArea > m, right.mad, left.mad),
             #mad.dist = abs(x-m)/mad) %>%
      #filter(mad.dist <= 3)
    
    model <- lmer(logPeakArea ~ class + (1 | RAD_ID) + (1 | RAD_ID:CoreVisit) 
                   + (1 | Analysis_date), 
                   data = subset)
    
    #plot(model, main = voc)
    lattice::qqmath(model, main = voc)
    
  })
}


diag_plots <- bg_s_gauss_plot(b1_imp_L_sum)

pdf('BGvsS_lmer_diagPlots_QQ.pdf')
marrangeGrob(diag_plots, nrow = 2, ncol = 2)
dev.off()

# remove VOCs with no difference in sample mean between S and BG
b1_imp_sum_c <- b1_imp_sum_c %>% filter(comp %ni% bg_s_output_f$comp)

b1_imp_sum_c <- b1_imp_sum_c %>% filter(comp %ni% wilcox_agreement)

write.csv(b1_imp_sum_c, 'RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC_BGfiltered.csv')

b1_imp_sum_c <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC_BGfiltered.csv')[,-1]

#
#
#

# APPROACH 1 - SUBTRACTION
# correct for blank/background values by subtraction
b1_imp_sum_c <- b1_imp_sum_c %>% mutate(Scorr = S - BG)


# plot global distribution in corrected/uncorrected data
b1_imp_sum_c %>% drop_na() %>% ggplot() + 
  geom_density(aes(x = S), alpha = 0.6) + 
  geom_density(aes(x = Scorr), colour = 'blue', alpha = 0.6) +
  theme_bw() +
  xlim(-18737.66, 75847.08)


quantile(b1_imp_sum_c$Scorr, na.rm = TRUE, c(0.1, 0.9))
min(b1_imp_sum_c$Scorr, na.rm = TRUE)
max(b1_imp_sum_c$Scorr, na.rm = TRUE)

quantile(b1_imp_sum_c$S, na.rm = TRUE, c(0.1, 0.9))
min(b1_imp_sum_c$S, na.rm = TRUE)
max(b1_imp_sum_c$S, na.rm = TRUE)


# effect on variance
sdS <- sd(b1_imp_sum_c$S, na.rm = TRUE)
sdScorr <- sd(b1_imp_sum_c$Scorr, na.rm = TRUE)

meanS <- mean(b1_imp_sum_c$S, na.rm = TRUE)
meanScorr <- mean(b1_imp_sum_c$Scorr, na.rm = TRUE)

sdS/abs(meanS)
sdScorr/abs(meanScorr)

# of individual compounds
pqn_corr_cvs <- b1_imp_sum_c %>% 
  group_by(comp) %>%
  summarise(cvS = sd(S, na.rm = TRUE)/abs(mean(S, na.rm = TRUE)),
            cvScorr = sd(Scorr, na.rm = TRUE)/abs(mean(Scorr, na.rm = TRUE)),
            madS = mad(S, na.rm = TRUE),
            madScorr = mad(Scorr, na.rm = TRUE))

  
pqn_corr_cvs %>%  pivot_longer(cols = c(madS, madScorr), names_to = 'Sample', values_to = 'cv') %>%
  ggplot(aes(x = Sample, y = cv, group = Sample)) + geom_boxplot(outliers = FALSE) 



# heteroscedastic?
b1_imp_sum_c %>% group_by(comp) %>%
  summarise(mean = mean(S, na.rm = TRUE),
            sd = sd(S, na.rm = TRUE)) %>%
  arrange(mean) %>%
  mutate(meanRank = row_number()) %>%
  ggplot(aes(x = meanRank, y = sd)) + geom_point()


b1_imp_sum_c %>% group_by(comp) %>%
  summarise(mean = abs(mean(Scorr, na.rm = TRUE)),
            sd = sd(Scorr, na.rm = TRUE)) %>%
  arrange(mean) %>%
  mutate(meanRank = row_number()) %>%
  ggplot(aes(x = meanRank, y = sd)) + geom_point()

# 

# distrubution of individual compounds after correction
voc_hist_Scorr <- b1_imp_sum_c %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = Scorr, fill = Diagnosis)) +
       geom_histogram(alpha = 0.7) +
       ggtitle(.$comp) +
       theme_bw(base_size = 6) +
       geom_vline(xintercept = 0, colour = 'red', linetype = 'dashed', size = 1) +
       theme(legend.position = 'none'))

pdf('Histograms_diff_SvsBG.pdf')
marrangeGrob(voc_hist_Scorr$plots, nrow = 4, ncol = 4)
dev.off()

#
iqr <- b1_imp_sum_c %>% 
  group_by(comp) %>%
  summarise(Q1 = quantile(Scorr, 0.25, na.rm = TRUE),
            Q3 = quantile(Scorr, 0.75, na.rm = TRUE),
            median = median(Scorr, na.rm = TRUE),
            mean = mean(Scorr, na.rm = TRUE),
            sd = sd(Scorr, na.rm = TRUE),
            mad = mad(Scorr, na.rm = TRUE))


iqr %>% filter(Q1 < 0 & Q3 < 0) %>% nrow()


#
#
#


# Wilcoxon paired test for distribution of differences
wilcox_bg <- bind_rows(lapply(unique(b1_imp_sum_c$comp), function(voc) {
    subset1 <- b1_imp_sum_c %>% 
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV2')) %>%
      drop_na()
    test <- wilcox.test(subset1$S, subset1$BG, paired = TRUE,)
    df <- data.frame(comp = voc, 
                     p.value = test[["p.value"]],
                     min = min(subset1$Scorr), 
                     max = max(subset1$Scorr),
                     Q1 = quantile(subset1$Scorr, 0.25), 
                     Q3 = quantile(subset1$Scorr, 0.75))
  }))


wilcox_bg <- wilcox_bg %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

wilcox_bg_f_cv2 <- wilcox_bg %>% filter(adj.p.value > 0.05)

wilcox_agreement <- intersect(wilcox_bg_f_cv1$comp, wilcox_bg_f_cv2$comp)

# remove VOCs where null hypothesis accepted in CV1 and CV2
b1_imp_sum_c <- b1_imp_sum_c %>% filter(comp %ni% wilcox_agreement)


intersect(wilcox_agreement, bg_s_output_f$comp)

setdiff(wilcox_agreement, bg_s_output_f$comp)
setdiff(bg_s_output_f$comp, wilcox_agreement)


View(wilcox_bg %>% filter(comp %ni% wilcox_bg_f$comp) %>% 
       filter(Q1 <0 & Q3 <0))

intersect(wilcox_agreement, wilcox_bg_f$comp)
setdiff(wilcox_bg_f_cv4$comp, wilcox_agreement)

#

# visualise histograms of differences
diff_hist <- b1_imp_sum_c %>% 
  filter(comp %ni% wilcox_bg_f$comp) %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = Scorr, fill = Diagnosis)) + 
                      geom_histogram(position = 'identity', alpha = 0.6) +
       #xlim(quantile(.$Scorr1, c(0.1, 0.9))) +
       theme_bw(base_size = 8) +
       ggtitle(.$comp) +
       theme(legend.position = 'none'))

pdf('Histograms_diff_logFC_SvsBG.pdf')
marrangeGrob(diff_hist$plots, nrow = 4, ncol = 4)
dev.off()

intersect(wilcox_bg_f$comp, reg_results_pqn2_f$comp)

#

# relationship between difference magnitude and peak area values
scorr_vs_s <- b1_imp_corr_L %>% 
  pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = peakArea, y = abs(Scorr), colour = class)) +
       geom_point() + ggtitle(.$comp))

marrangeGrob(scorr_vs_s$plots, nrow = 2, ncol = 2)


b1_imp_corr_L %>% 
  pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  filter(comp == 'X2_methylbutanal') %>%
  ggplot(aes(x = Sample, y = peakArea, colour = class)) +
       geom_point() + geom_line(aes(group = Sample), colour = 'black')

b1_imp_corr_L %>% 
  pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  filter(comp == 'Pentane') %>%
  ggplot(aes(x = peakArea, y = Scorr, colour = class)) +
  geom_point()

# calculate cvs for each pair and plot against mean
dist_cv <- b1_imp_corr_L %>% 
  pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  group_by(comp, Sample) %>%
  summarise(mean = mean(peakArea, na.rm = TRUE),
            sd = sd(peakArea, na.rm = TRUE),
            cv = sd/abs(mean)) %>%
  left_join(b1_imp_corr_L %>% dplyr::select(Sample, comp, Scorr))
  
cv_dist <- dist_cv %>% 
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = Scorr, y = cv)) + 
       geom_point() + ggtitle(.$comp) + theme_bw(base_size = 8))

marrangeGrob(cv_dist$plots, nrow = 4, ncol = 4)


#
#
#

# PRINCIPAL COMPONENT ANALYSIS

# APPROACH 0
b1_sub_w <- b1_imp_sum_c %>% dplyr::select(!c(Scorr, BG)) %>%
  pivot_wider(names_from = comp, values_from = S) %>% drop_na() %>%
  as.data.frame() %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

rownames(b1_sub_w) <- b1_sub_w$Sample

# keep only sample with both observations present
c_reps <- b1_sub_w %>% group_by(RAD_ID) %>% summarise(n = n()) %>%
  filter(n == 2)

b1_sub_w <- b1_sub_w %>% filter(RAD_ID %in% c_reps$RAD_ID)

#
assay0 <- b1_sub_w[,-c(1:6)]

# variance stabilise
vsn <- vsn::vsnMatrix(as.matrix(assay0))
vsnout <- vsn@hx
assay0v <- as.data.frame(vsnout)

pc0 <- pca(log(assay0), center = TRUE, scale = TRUE, ncomp = 4) 
           #multilevel = b1_sub_w$RAD_ID)
p0 <- plotIndiv(pc0,
          #pch = 1,
          #ind.names = b1_sub_w$RAD_ID,
          pch = as.factor(b1_sub_w$CoreVisit),
          cex = 1,
          group = b1_sub_w$Diagnosis,
          legend = TRUE,
          title = 'No correction',
          comp = c(1,2),
          size.title = 10)

ggsave('PCA_S_noCorr.tiff', p0$graph, unit = 'mm', dpi = 300, width = 120, height = 80 )

pc0$prop_expl_var$X
  

#
#
#

# APPROACH 1 cd
b1_sub_corr <- b1_imp_sum_c

# keep only compounds with positive IQR
pos_comps <- b1_imp_sum_c %>% group_by(comp) %>%
  summarise(Q1 = quantile(Scorr, 0.25, na.rm = TRUE),
            Q3 = quantile(Scorr, 0.75, na.rm = TRUE)) %>%
  filter(Q1 > 0 & Q3 > 0)

neg_comps <- b1_imp_sum_c %>% group_by(comp) %>%
  summarise(Q1 = quantile(Scorr, 0.25, na.rm = TRUE),
            Q3 = quantile(Scorr, 0.75, na.rm = TRUE)) %>%
  filter(Q1 < 0 & Q3 < 0)

pos_med <- b1_imp_sum_c %>% group_by(comp) %>%
  summarise(median = median(Scorr, na.rm = TRUE)) %>%
  filter(median >= 0)

b1_sub_corr1 <- b1_imp_sum_c %>% filter(comp %ni% neg_comps$comp)
b1_sub_corr1 <- b1_imp_sum_c %>% filter(comp %in% pos_comps$comp)
b1_sub_corr1 <- b1_imp_sum_c %>% filter(comp %in% pos_med$comp)

test <- b1_sub_corr %>% filter(comp %ni% c(pos_med$comp, neg_comps$comp))
unique(test$comp)


#

b1_sub_corr_w <- b1_sub_corr1 %>% dplyr::select(!c(S, BG)) %>%
  pivot_wider(names_from = comp, values_from = Scorr) %>% drop_na() %>%
  as.data.frame() #%>%
  #filter(CoreVisit %in% c('CV1', 'CV2'))

b1_sub_w <- b1_sub_corr1 %>% dplyr::select(!c(Scorr, BG)) %>%
  pivot_wider(names_from = comp, values_from = S) %>% drop_na() %>%
  as.data.frame() %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

rownames(b1_sub_corr_w) <- b1_sub_corr_w$Sample
rownames(b1_sub_w) <- b1_sub_w$Sample

write.csv(b1_sub_corr_w, 'RADiCA_breath_BGcorrected.csv')

#

b1_sub_corr_w1 <- b1_sub_corr1 %>% dplyr::select(!c(S, BG)) %>%
  pivot_wider(names_from = comp, values_from = Scorr) %>% drop_na() %>%
  as.data.frame() %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

rownames(b1_sub_corr_w1) <- b1_sub_corr_w1$Sample

# keep only sample with both observations present
c_reps <- b1_sub_corr_w %>% group_by(RAD_ID) %>% summarise(n = n()) %>%
  filter(n >= 2)

b1_sub_corr_w <- b1_sub_corr_w %>% filter(RAD_ID %in% c_reps$RAD_ID)
b1_sub_corr_w1 <- b1_sub_corr_w1 %>% filter(RAD_ID %in% c_reps$RAD_ID)

#

assay <- b1_sub_corr_w[,-c(1:6)]
assay1 <- b1_sub_w[,-c(1:6)]

# variance stabilise
library(vsn)

vsn <- vsn::vsnMatrix(as.matrix(assay))
vsnout <- vsn@hx
assay <- as.data.frame(vsnout)

vsn <- vsn::vsnMatrix(as.matrix(assay1))
vsnout <- vsn@hx
assay1 <- as.data.frame(vsnout)

# check results of variance stabilisation
assay %>% 
  pivot_longer(cols = colnames(vsnout), names_to = 'comp', values_to = 'peakArea') %>%
  group_by(comp) %>%
  summarise(mean = abs(mean(peakArea)),
            sd = sd(peakArea)) %>%
  arrange(mean) %>%
  mutate(meanRank = row_number()) %>%
  ggplot(aes(x = meanRank, y = sd)) + geom_point()

#
library(mixOmics)

pc1 <- pca(assay, center = TRUE, scale = TRUE, ncomp = 4)
           #multilevel = b1_sub_corr_w$RAD_ID)

pc1a <- pca(assay1, center = TRUE, scale = TRUE, ncomp = 4)
           #multilevel = b1_sub_corr_w1$RAD_ID)

barplot(pc1$prop_expl_var$X)
barplot(pc1a$prop_expl_var$X)

p1 <- plotIndiv(pc1,
          #pch = 1,
          cex = 3,
          group = b1_sub_corr_w$Diagnosis,
          ind.names = b1_sub_corr_w$RAD_ID,
          #pch = as.factor(b1_sub_corr_w$CoreVisit),
          legend = TRUE,
          title = 'Breath sample VOCs (corrected)',
          comp = c(1,2),
          size.title = 10)
          #ellipse = TRUE)

p1a <- plotIndiv(pc1a,
                #pch = 1,
                cex = 1,
                group = b1_sub_corr_w$Diagnosis,
                #ind.names = b1_sub_corr_w$RAD_ID,
                pch = as.factor(b1_sub_corr_w$CoreVisit),
                legend = TRUE,
                title = 'Breath sample VOCs (uncorrected)',
                comp = c(1,2),
                size.title = 10)
#ellipse = TRUE)


ggsave('PCA_S_PosMed.tiff', p1a$graph, unit = 'mm', dpi = 300, width = 130, height = 90 )

# CV1 and CV2 separately
cv1 <- b1_sub_corr_w %>% filter(CoreVisit == 'CV1')
cv2 <- b1_sub_corr_w %>% filter(CoreVisit == 'CV2')

assay_c1 <- assay %>% filter(rownames(.) %in% cv1$Sample)
assay_c2 <- assay %>% filter(rownames(.) %in% cv2$Sample)

assay1_c1 <- assay1 %>% filter(rownames(.) %in% cv1$Sample)
assay1_c2 <- assay1 %>% filter(rownames(.) %in% cv2$Sample)
#

pc1_c1 <- pca(assay_c1, center = TRUE, scale = TRUE, ncomp = 4)
pc1a_c1 <- pca(assay1_c1, center = TRUE, scale = TRUE, ncomp = 4)

pc1_c2 <- pca(assay_c2, center = TRUE, scale = TRUE, ncomp = 4)
pc1a_c2 <- pca(assay1_c2, center = TRUE, scale = TRUE, ncomp = 4)


barplot(pc1a_c2$prop_expl_var$X)
#

p1_c2 <- plotIndiv(pc1_c2,
                #pch = 1,
                cex = 1,
                group = cv2$Diagnosis,
                #ind.names = b1_sub_corr_w$RAD_ID,
                pch = as.factor(cv2$CoreVisit),
                legend = FALSE,
                title = 'CV2 - Breath sample VOCs (corrected)',
                comp = c(1,2),
                size.title = 8)
#ellipse = TRUE)

p1_cv <- arrangeGrob(p1_c1$graph, p1_c2$graph, nrow = 1)
plot(p1_cv)

ggsave('PCA_S_SubPosMed_CV_sep.tiff', p1_cv, unit = 'mm', dpi = 300, width = 160, height = 75 )

#

p1a_c2 <- plotIndiv(pc1a_c2,
                 #pch = 1,
                 cex = 1,
                 group = cv2$Diagnosis,
                 #ind.names = b1_sub_corr_w$RAD_ID,
                 pch = as.factor(cv2$CoreVisit),
                 legend = FALSE,
                 title = 'CV2 - Breath sample VOCs (uncorrected)',
                 comp = c(1,2),
                 size.title = 8)

p1a_cv <- arrangeGrob(p1a_c1$graph, p1a_c2$graph, nrow = 1)
plot(p1a_cv)

ggsave('PCA_S_PosMed_CV_sep.tiff', p1a_cv, unit = 'mm', dpi = 300, width = 160, height = 75 )

#
#
#


# APPROACH 2

# identify outliers
# robust PCA
library(pcaPP)

# for breath data
s_assay <- b1_imp_sum %>% filter(class == 'S') %>% 
  left_join(b1_imp_sum_c %>% dplyr::select(Sample, CoreVisit)) %>%
  distinct() %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  as.data.frame()
rownames(s_assay) <- s_assay$Sample 
s_assay <- s_assay %>% dplyr::select(!c('Sample', 'Batch', 'CoreVisit', 'class',
                                        'Analysis_date'))

Rpc <- PCAgrid(log(s_assay), scale = 'sd', center = 'mean')

summary(Rpc)

scores <- as.data.frame(Rpc$scores) %>%
  mutate(Sample = rownames(s_assay)) 
scores %>% ggplot(aes(x = Comp.1, y = Comp.2)) + geom_point()

sdod <- PCdiagplot(s_assay, Rpc, plotbw = FALSE)

# identify outliers based on the critical OD and SD values (0.999 quantile)
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(s_assay)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,1] | V2 > sdod$critSD[2,1])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(s_assay)
odOut <- od %>% filter(V1 > sdod$critOD[1,1] | V2 > sdod$critOD[2,1])

# for BG data
bg_assay <- b1_imp_sum %>% filter(class == 'BG') %>% 
  left_join(b1_imp_sum_c %>% dplyr::select(Sample, CoreVisit)) %>%
  distinct() %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  as.data.frame()
rownames(bg_assay) <- bg_assay$Sample 
bg_assay <- bg_assay %>% dplyr::select(!c('Sample', 'Batch', 'CoreVisit', 'class',
                                        'Analysis_date'))

Rpc_bg <- PCAgrid(log(bg_assay), scale = 'sd', center = 'mean')

summary(Rpc_bg)

scores <- as.data.frame(Rpc_bg$scores) %>%
  mutate(Sample = rownames(bg_assay)) 
scores %>% ggplot(aes(x = Comp.1, y = Comp.2)) + geom_point()

bg_sdod <- PCdiagplot(bg_assay, pc_bg, plotbw = FALSE)

# identify outliers based on the critical OD and SD values (0.999 quantile)
sd_bg <- as.data.frame(bg_sdod$SDist)
rownames(sd_bg) <- rownames(bg_assay)
sdOut_bg <- sd_bg %>% filter(V1 > bg_sdod$critSD[1,1] | V2 > bg_sdod$critSD[2,1])

od_bg <- as.data.frame(bg_sdod$ODist)
rownames(od_bg) <- rownames(bg_assay)
odOut_bg <- od_bg %>% filter(V1 > bg_sdod$critOD[1,1] | V2 > bg_sdod$critOD[2,1])

#

intersect(rownames(odOut), rownames(odOut_bg))

# 14 samples flagged as outliers according to stringent criteria
# 6 samples according to less stringent criteria


# linear regression
reg_coefs <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
      #filter(Sample %ni% c('190903_RaDICA_RAD032'))
    
    
    model <- lmer(logS ~ logBG*Diagnosis + CoreVisit + (1 | RAD_ID),
                 data = subset)
    
    coefs <- cbind(as.data.frame(coef(summary(model))),
                   as.data.frame(confint(model)[3:7,])) %>%
      mutate(comp = rep(voc, 5),
             n = rep(nrow(subset), 5),
             r2 = rep(r2(model)[[2]], 5),
             predictor = rownames(.),
             AIC = rep(AIC(model), 5),
             BIC = rep(BIC(model), 5))
    
    
  }))
  
}

#

reg_results0_pqn <- reg_coefs(b1_imp_sum_c)
colnames(reg_results0_pqn)[5] <- 'p.value'
reg_results0_pqn_itc <- reg_results0_pqn %>%
  filter(predictor == 'logBG:DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'fdr')) %>%
           filter(p.value < 0.05)

View(reg_results0_pqn %>%
       filter(predictor == 'DiagnosisNot Asthma') %>%
       mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
       filter(p.value < 0.05))

#

infl <- dfbetas(model) %>% as.data.frame()
hist(infl$`logBG:DiagnosisNot Asthma`)
t <- 2/sqrt(94)
View(infl %>% filter(abs(`logBG:DiagnosisNot Asthma`) > t))

#

plot_voc <- function(voc){
  subset <- df %>%
    filter(comp == voc) %>%
    filter(CoreVisit %in% c('CV1', 'CV2')) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) 

  
  subset %>%
    #filter(Diagnosis == 'Not Asthma') %>%
    ggplot(aes(x = logBG, y = logS, colour = Diagnosis)) + 
    geom_point(size = 1, alpha = 0.8) +
    ylab('logS') + xlab('logBG') +
    theme_bw(base_size = 10) +
    #theme(legend.position = 'none') +
    #geom_label(aes(label = Sample)) +
    geom_smooth(method = 'lm', se = FALSE) +
    ggtitle(voc)
#}

dev.new()
plot_voc('Azetidine')
ggsave('BGvsS_Ethyl_acetate.tiff', unit = 'mm', dpi = 300, width = 120, height = 80)

subset %>%
  mutate(logS1 = logS - 0.41*logBG) %>%
  ggplot(aes(x = Diagnosis, y = logS, fill = Diagnosis)) + 
  geom_violin() +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  ggtitle('Ethyl acetate uncorrected') +
  ylab('logS') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('Boxplot_uncorrected_Ethyl_acetate.tiff', unit = 'mm', dpi = 300, width = 120, height = 80)


t <- subset %>%
  mutate(logS1 = logS - 0.41*logBG) 

t %>% group_by(Diagnosis) %>%
  summarise(mean = mean(logS),
            corr_mean = mean(logS1),
            sd = sd(logS),
            corr_sd = sd(logS1))

summary(lm(logS1 ~ Diagnosis,
   data = t))

#
#
#

reg_coefs1 <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
      #filter(Sample %ni% c('190903_RaDICA_RAD032'))#
    #))
    
    #x <- subset$logS
    #m <- median(x, na.rm = TRUE)
    #abs.dev   <- abs(x - m)
    #left.mad  <- median(abs.dev[x<=m])
    #right.mad <- median(abs.dev[x>=m])
    
    #subset <- subset %>% 
    #mutate(mad = ifelse(x > m, right.mad, left.mad),
    #mad.dist = abs(x-m)/mad) %>%
    #filter(mad.dist <= 3)
    
    model <- lmer(logS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID),
                  data = subset)
    
    coefs <- cbind(as.data.frame(coef(summary(model))),
                   as.data.frame(confint(model)[3:6,])) %>%
      mutate(comp = rep(voc, 4),
             n = rep(nrow(subset), 4),
             r2 = rep(r2(model)[[2]], 4), 
             predictor = rownames(.),
             AIC = rep(AIC(model), 4),
             BIC = rep(BIC(model), 4))
    
    
  }))
  
}

reg_results1_pqn <- reg_coefs1(b1_imp_sum_c)

colnames(reg_results1_pqn)[5] <- 'p.value'

write.csv(reg_results1_pqn, 'Model1_results.csv')

asthma_pred_pqn <- reg_results1_pqn %>% filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

asthma_pred_pqn_sign <- asthma_pred_pqn %>% filter(p.value < 0.05)

hist(asthma_pred_pqn$p.value, breaks = 30)

# log FC asthma and not asthma
fc <- b1_imp_sum_c %>% 
  left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
  mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
  filter(obs %ni% c('infl')) %>%
  mutate(logS = log2(S)) 

fc_sum <- fc %>% group_by(comp, Diagnosis) %>%
  summarise(mean = mean(logS, na.rm = TRUE)) %>%
  pivot_wider(names_from = Diagnosis, values_from = mean) %>%
  rename(Not_Asthma = 'Not Asthma') %>%
  mutate(logFC = Asthma - Not_Asthma) %>%
  left_join(asthma_pred_pqn %>%
              dplyr::select(comp, p.value, adj.p.value)) %>%
  mutate(lab = ifelse(p.value < 0.05 & logFC > 0.6, comp, ifelse(
    p.value < 0.05 & logFC < -0.6, comp, NA))) %>%
  mutate(col = ifelse(is.na(lab) == TRUE, 'no', 'yes'))

vol_plot <- fc_sum %>%
  ggplot(aes(y = -log10(p.value), x = logFC)) +
  geom_point(aes(colour = col)) +
  geom_vline(xintercept = c(-1, 1), linetype = 'dashed', colour = 'darkgrey') +
  geom_vline(xintercept = c(-0.6, 0.6), linetype = 'dashed', colour = 'darkgrey') +
  geom_hline(yintercept = 1.3, linetype = 'dashed', colour = 'darkgrey') +
  theme_bw() +
  geom_text_repel(aes(label = lab), size = 2) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c('yes' = 'red',
                                 'no' = 'black')) +
  ggtitle('Volcano plot showing effect of asthma diagnosis on breath VOCs') +
  xlim(-2.06, 2.06) +
  xlab('log2 Fold Change')

vol_plot

ggsave('Volcano_asthma.tiff', dpi = 300, unit = 'mm', width = 150, height = 120)

#
#
#

# compare model fit with and without interaction term
# with and without random effect
plot(unique(reg_results0_pqn$BIC), unique(reg_results1_pqn$BIC)) +
  abline(0, 1)

plot(unique(reg_results_pqn$BIC), unique(reg_results1_pqn$BIC)) +
  abline(0, 1)

mod_fit <- rbind(
  reg_results0_pqn %>% mutate(model = 'Fixed_itc') %>%
  filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2),
  reg_results1_pqn %>% mutate(model = 'Fixed') %>%
    filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2),
  reg_results_pqn %>% mutate(model = 'Mixed') %>%
    filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2))
  
comp <- mod_fit %>% 
  mutate(model = ifelse(model == 'Fixed', 'Model_1', ifelse(
    model == 'Fixed_itc', 'Model_2', 'Baseline_model'))) %>%
  pivot_longer(cols = c(AIC, BIC), names_to = 'Parameter', values_to = 'value') %>%
  ggplot(aes(x = model, y = value, fill = model)) + geom_violin() +
  geom_boxplot(width = 0.5) +
  facet_wrap(~Parameter, scale = 'free') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom')

comp

ggsave('BIC_AIC_violin.tiff', dpi = 300, unit = 'mm', width = 120, height = 70)

comp1 <- mod_fit %>% 
  mutate(model = ifelse(model == 'Fixed', 'Model_1', ifelse(
    model == 'Fixed_itc', 'Model_2', 'Baseline_model'))) %>%
  pivot_longer(cols = c(AIC, BIC), names_to = 'Parameter', values_to = 'value') %>%
  ggplot(aes(x = value, colour = model)) + geom_density(alpha = 0.6) +
  facet_wrap(~Parameter) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = 'bottom')

comp1

ggsave('BIC_AIC_density.tiff', dpi = 300, unit = 'mm', width = 120, height = 70)

mod_fit %>% group_by(model) %>% summarise(mean_AIC = mean(AIC),
                                         mean_BIC = mean(BIC),
                                         mean_r2 = mean(r2),
                                         median_AIC = median(AIC),
                                         median_BIC = median(BIC),
                                         median_r2 = median(r2))
#

p1 <- mod_fit %>% 
  mutate(model = ifelse(model == 'Fixed', 'Model_1', ifelse(
    model == 'Fixed_itc', 'Model_2', 'Baseline_model'
  ))) %>%
  pivot_wider(id_cols = comp, names_from = model, values_from = AIC) %>%
  ggplot(aes(y = Baseline_model, x = Model_1)) + 
  #geom_point(alpha = 0.3, size = 2) +
  stat_binhex() +
  theme_bw() +
  theme(legend.key.size = unit(5, 'mm')) +
  coord_fixed() + geom_abline(intercept = 0, slope = 1, colour = 'red') +
  ggtitle('AIC')
 

p1

p2 <- mod_fit %>% 
  mutate(model = ifelse(model == 'Fixed', 'Model_1', ifelse(
    model == 'Fixed_itc', 'Model_2', 'Baseline_model'
  ))) %>%
  pivot_wider(id_cols = comp, names_from = model, values_from = AIC) %>%
  ggplot(aes(y = Model_1, x = Model_2)) + 
  #geom_point(alpha = 0.3, size = 2) +
  theme_bw() +
  stat_binhex() +
  theme(legend.key.size = unit(5, 'mm')) +
  coord_fixed() + geom_abline(intercept = 0, slope = 1, colour = 'red') +
  ggtitle('AIC')


p2


ggsave('AIC_model_comp2.tiff', p2, dpi = 300, unit = 'mm', width = 100, height = 70)

#

reg_coefs2 <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      drop_na() %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
      #filter(Sample %ni% c('190903_RaDICA_RAD032'))#
      #))
    
    #x <- subset$logS
    #m <- median(x, na.rm = TRUE)
    #abs.dev   <- abs(x - m)
    #left.mad  <- median(abs.dev[x<=m])
    #right.mad <- median(abs.dev[x>=m])
    
    #subset <- subset %>% 
      #mutate(mad = ifelse(x > m, right.mad, left.mad),
            #mad.dist = abs(x-m)/mad) %>%
      #filter(mad.dist <= 3)
    
    model <- lmer(logS ~ logBG + (1 | RAD_ID),
                 data = subset) 
    
    coefs <- as.data.frame(coef(summary(model))) %>%
      mutate(predictor = rownames(.),
             comp = voc,
             n = nrow(subset),
             CI_lwr = ifelse(predictor == 'logBG', 
                             confint(model)[4,1],
                             confint(model)[3,1]),
             CI_upr = ifelse(predictor == 'logBG',
                             confint(model)[4,2],
                             confint(model)[3,2]),
             r2 = r2(model)[[2]],
             AIC = AIC(model),
             BIC = BIC(model)) 

  
  
  }))
  
}


reg_results_pqn <- reg_coefs2(b1_imp_sum_c)

colnames(reg_results_pqn)[5] <- 'p.value'

reg_results_pqn <- reg_results_pqn %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

test <- reg_results_pqn %>% 
       arrange(p.value) %>%
       mutate(rank = c(1:152),
              fdr = 0.05*(rank/nrow(.)))
View(test %>% filter(p.value < fdr))

reg_results_pqn_sign <- reg_results_pqn %>% filter(adj.p.value < 0.05)

write.csv(reg_results_pqn, 'PQN_regression_results_BG_filtered.csv')

#

reg_results_pqn_infl <- reg_coefs2(b1_imp_sum_c)

colnames(reg_results_pqn_infl)[5] <- 'p.value'

reg_results_pqn_infl <- reg_results_pqn_infl %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

reg_results_pqn_inflsign <- reg_results_pqn_infl_BG %>% filter(adj.p.value < 0.05)

# assess FDR in regression results
reg_results_pqn %>% ggplot(aes(x = p.value)) + 
  geom_histogram(bins = 20) + 
  scale_y_continuous(breaks = seq(0, 105, 5)) +
  geom_hline(yintercept = 5) +
  geom_vline(xintercept = 0.05)

reg_results_pqn_infl %>% ggplot(aes(x = p.value)) + 
  geom_histogram(bins = 20) + 
  scale_y_continuous(breaks = seq(0, 105, 5)) +
  geom_hline(yintercept = 5) +
  geom_vline(xintercept = 0.05)

reg_results1_pqn %>%
  filter(predictor == 'logBG') %>%
  ggplot(aes(x = p.value)) + 
  geom_histogram(bins = 20) + 
  scale_y_continuous(breaks = seq(0, 105, 5)) +
  geom_hline(yintercept = 5) +
  geom_vline(xintercept = 0.05)

# compare the effect of excluding influential observations on regression estimates
# number of VOCs with significant logBG effect
library(VennDiagram)

setdiff(reg_results_pqn_sign$comp, reg_results_pqn_infl_sign$comp)
setdiff(reg_results_pqn_infl_sign$comp, reg_results_pqn_sign$comp)

venn.diagram(
  x = list(reg_results_pqn_sign$comp, reg_results_pqn_infl_sign$comp),
  category.names = c('All \n data', 'Excl. \n outliers'),
  filename = 'Venn_mixed_reg_sign.png',
  output = TRUE,
  cex = 2,
  fontface = 'bold',
  cat.cex = 2,
  cat.fontface = 'bold',
  cat.default.pos = "outer",
  fill = c('pink', 'grey'),
  fontfamily = "sans",
  cat.fontfamily = 'sans',
  margin = 0.08,
  cat.dist = c(0.05, 0.05),
  main = 'Number of VOCs with \n significant logBG effect',
  main.cex = 2,
  main.fontface = 'bold',
  main.fontfamily = 'sans',
  main.just = c(0.5, -0.5)
)

# slopes
reg_results_comp <- reg_results_pqn %>%
  mutate(Data = 'All_data') %>%
  rbind(reg_results_pqn_infl %>%
          mutate(Data = 'Excl_outliers')) %>%
  mutate(Sign = ifelse(adj.p.value < 0.05, 'Yes', 'No')) %>%
  mutate(CI_width = CI_upr - CI_lwr) 


reg_results_comp1 <- reg_results_comp %>%
  pivot_wider(id_cols = comp, names_from = Data, values_from = c(Estimate, Sign, CI_width, r2)) %>%
  mutate(Significant = ifelse(Sign_All_data == 'Yes' & Sign_Excl_outliers == 'Yes',
                       'Both', ifelse(
                         Sign_All_data == 'Yes' & Sign_Excl_outliers == 'No',
                         'All_data', ifelse(
                           Sign_All_data == 'No' & Sign_Excl_outliers == 'No',
                           'None', 'Excl_outliers'
                         )
                       ))) 

 
estimate_comp <- 
  reg_results_comp1 %>%
  ggplot(aes(y = Estimate_All_data, x = Estimate_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.7) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 10) +
  xlab('Excl. outliers') +
  ylab('All data') +
  ggtitle('Regression logBG coefficient') +
  theme(legend.position = 'none')

estimate_comp

CI_comp <- 
  reg_results_comp1 %>%
  ggplot(aes(y = CI_width_All_data, x = CI_width_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.7) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 10) +
  xlab('Excl. outliers') +
  ylab('All data') +
  ggtitle('Regression logBG 95% CI width') +
  theme(legend.position = 'none')

r2_comp <- 
  reg_results_comp1 %>%
  ggplot(aes(y = r2_All_data, x = r2_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.7) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 10) +
  xlab('Excl. outliers') +
  ylab('All data') +
  ggtitle('Regression logBG R2')


library(cowplot)

comp_plots <- plot_grid(estimate_comp, CI_comp, r2_comp, nrow = 1, align = 'h',
          rel_widths = c(0.3, 0.3, 0.4))

comp_plots

ggsave('Regression_comparison_del_diag.tiff' , comp_plots, dpi = 300, unit = 'mm',
       width = 260, height = 80)

#
#
#

# compare the effect of including covariates on logBG slope estimation
reg_results <- reg_results_pqn_infl %>%
  filter(predictor == 'logBG') %>%
  mutate(model = 'Mixed') %>%
  rbind(reg_results1_pqn %>% 
          filter(predictor == 'logBG') %>%
          mutate(model = 'Fixed',
                 CI_lwr = `2.5 %`,
                 CI_upr = `97.5 %`,
                 adj.p.value = p.adjust(p.value, method = 'BH')) %>%
          dplyr::select(!c('2.5 %', '97.5 %')))
  
# number of VOCs with significant effect of background
m_sbg <- reg_results %>% filter(model == 'Mixed') %>% filter(adj.p.value < 0.05)
f_sbg <- reg_results %>% filter(model == 'Fixed') %>% filter(adj.p.value < 0.05)

m_sbg_only <- setdiff(m_sbg$comp, f_sbg$comp)
f_sbg_only <- setdiff(f_sbg$comp, m_sbg$comp)
both <- intersect(f_sbg$comp, m_sbg$comp)

reg_results <- reg_results %>% 
  mutate(Significant = ifelse(comp %in% m_sbg_only, 'Baseline', ifelse(
  comp %in% f_sbg_only, 'Model_1', ifelse(
    comp %in% both, 'Both', 'Neither'
  )
)))

table(reg_results$Significant)

#
est <- reg_results %>% 
  pivot_wider(id_cols = c(comp, Significant), names_from = model, values_from = Estimate) %>%
  ggplot(aes(x = Fixed, y = Mixed)) + 
  #geom_point() +
  geom_hex() +
  xlab('Model_1') + ylab('Baseline_model') +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme_bw() +
  ggtitle('Regression logBG coefficient') #+
  #geom_text(aes(label = comp))

reg_results %>% ggplot(aes(x = model, y = Estimate, group = model)) + 
  geom_violin() + geom_boxplot()

reg_results <- reg_results %>% 
  mutate(CI_width = CI_upr - CI_lwr)

ci <- reg_results %>%
  pivot_wider(id_cols = c(comp, Significant), names_from = model, values_from = CI_width) %>%
  ggplot(aes(x = Fixed, y = Mixed)) + 
  #geom_point() +
  geom_hex() +
  xlab('Model_1') + ylab('Baseline_model') +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  theme_bw() +
  ggtitle('Regression logBG 95% CI width')

reg_results %>% 
  ggplot(aes(x = model, y = CI_width, group = model)) + 
  geom_violin() + geom_boxplot()

reg_results %>% 
  ggplot(aes(x = adj.p.value, fill = model)) + 
  geom_histogram(position = 'identity', alpha = 0.6)



baseline_vs_m1 <- arrangeGrob(est, ci, nrow = 1)
plot(baseline_vs_m1)
ggsave('Baseline_vs_M1_Est_CI.tiff', baseline_vs_m1, dpi = 300, unit = 'mm', width = 160, height = 80)


#
#
#

tiff(filename = 'PQN_reg_results_BGfiltered.tiff', res = 300, units = 'mm',
     width = 180, height = 80)
par(mfrow = c(1,2))
hist(reg_results_pqn_sign$Estimate, main = 'BG regression coefficients',
     xlab = '')
hist(reg_results_pqn_sign$r2, breaks = 12, main = 'R2',
     xlab = '')
dev.off()

#
#
#


#pent_pqn <- 
  subset %>%
  ggplot(aes(x = log(BG), y = log(S), colour = Diagnosis)) + 
  geom_point(size = 1) +
  ylab('logS') + xlab('logBG') +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none') +
  #geom_label(aes(label = Sample), size = 3)
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = FALSE) 

is_pqn
is_pqn_wo

is_nn_plot <- arrangeGrob(is_nn, is_nn_wo, nrow = 2)
ggsave('IS_regression.tiff', is_nn_plot, dpi = 300, unit = 'mm', width = 50, height = 100)

is_pqn_plot <- arrangeGrob(is_pqn, is_pqn_wo, nrow = 2)
ggsave('IS_regression_pqn.tiff', is_pqn_plot, dpi = 300, unit = 'mm', width = 50, height = 100)

#
acetone_plot <- arrangeGrob(acet_pqn, acet_nn, nrow = 1,
                            top = textGrob('Acetone'))

ggsave('Acetone_regression.tiff', acetone_plot, dpi = 300, unit = 'mm', width = 100, height = 50)

#

isoprene_plot <- arrangeGrob(iso_pqn, iso_nn, nrow = 1,
                            top = textGrob('Isoprene'))

ggsave('Isoprene_regression.tiff', isoprene_plot, dpi = 300, unit = 'mm', width = 100, height = 50)

#

pentane_plot <- arrangeGrob(pent_pqn, pent_nn, nrow = 1,
                             top = textGrob('Pentane'))

ggsave('Pentane_regression.tiff', pentane_plot, dpi = 300, unit = 'mm', width = 100, height = 50)

#

hist(subset$logBG)
sd(subset1$S/mean(subset1$S))
sd(subset$S/mean(subset$S))
sd(subset$logBG)
range(subset1$logS)
mad(subset1$logS)
IQR(subset$logBG)

subset <- subset %>% filter(Sample %ni% c('190903_RaDICA_RAD032'))
subset1 <- subset1 %>% filter(Sample %ni% c('191105_RaDICA_RAD019'))

# deletion diagnostics
del_diag <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
      subset <- df %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) 
   
    
    model <- lmer(logS ~ logBG + (1 | RAD_ID),
                  data = subset)
    
    infl <- influence(model)
    dfB <- dfbetas(infl)[,2] %>% as.data.frame() %>%
      mutate(comp = voc) %>%
      cbind(cooks.distance(model)) %>%
      rename(dfB = '.',
             cooksD = 'cooks.distance(model)')})
  )
}

pqn_diag <- del_diag(b1_imp_sum_c)

# visualise deletion diagnostics results
library(forcats)
library(ggrepel)

pqn_diag$obs <- gsub("\\..*","", rownames(pqn_diag))

del_diag_cooks <- 
  lapply(unique(pqn_diag$comp), function(voc) {
  data <- pqn_diag %>% filter(comp == voc) %>%
  mutate(label = ifelse(cooksD > 0.5, obs, NA))
  
  data %>% ggplot(aes(x = fct_inorder(as.factor(obs)), y = abs(cooksD))) + 
  geom_col(fill = 'seagreen4', width = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 0.5, colour = 'seagreen4', linetype = 'dashed') +
  geom_text_repel(aes(label = label), size = 2) +
  ylab('Cooks distance') +
  ggtitle(voc)
    })

pdf('Deletion_diagnostic_cooksD.pdf')
marrangeGrob(del_diag_cooks, nrow = 4, ncol = 1)
dev.off()

#

del_diag_dfbetas <- 
  lapply(unique(pqn_diag$comp), function(voc) {
    data <- pqn_diag %>% filter(comp == voc) %>%
      mutate(label = ifelse(abs(dfB) > (2/sqrt(95)), obs, NA))
    
    data %>% ggplot(aes(x = fct_inorder(as.factor(obs)), y = abs(dfB))) + 
      geom_col(fill = 'darkorange3', width = 0.6) +
      theme_bw() +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank()) +
      geom_hline(yintercept = 2/sqrt(95), colour = 'darkorange3', linetype = 'dashed') +
      geom_text_repel(aes(label = label), size = 2) +
      ylab('Dfbetas') +
      ggtitle(voc)
  })

pdf('Deletion_diagnostic_dfbetas.pdf')
marrangeGrob(del_diag_dfbetas, nrow = 4, ncol = 1)
dev.off()

#

infl_obs <- read.csv('Deletion_diagnostics.csv')
infl_obs <- infl_obs %>% 
  rename(comp = VOC) %>%
  separate(Potential_outlier, 
  into = c('o1', 'o2', 'o3', 'o4', 'o5', 'o6', 'o7', 'o8'), 
  sep = ',') %>%
  pivot_longer(cols = !comp, names_to = 'outlier_no', values_to = 'obs') %>%
  drop_na()

infl_obs$obs <- gsub(' ', '', infl_obs$obs)

table(infl_obs$outlier_no) %>% as.data.frame() %>% 
  ggplot(aes(x = gsub('o', '', Var1), y = Freq)) +
  geom_col() +
  xlab('Number of influential observations') +
  ggtitle('Distribution of VOCs according to \n the influential observation count') +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10)) +
  ylab('VOC number') +
  geom_vline(xintercept = 5, colour = 'red') +
  annotate('text',label = '5% of data', x = 7, y = 140, colour = 'red')
  
ggsave('Histogram_infl_obs_VOCs.tiff', unit = 'mm', dpi = 300, width = 65, height = 50)  
  
table(infl_obs$obs) %>% as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  ggplot(aes(x = fct_inorder(Var1), y = Freq)) + 
  geom_col() +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank()) +
  xlab('observation ID') +
  ylab('VOC number') +
  theme(plot.title = element_text(size = 10)) +
  ggtitle('Number of VOCs with observation flagged as influential')

ggsave('Histogram_infl_obs.tiff', unit = 'mm', dpi = 300, width = 100, height = 60)

# add Sample IDs to influential observations database
infl_obs <- infl_obs %>%
  left_join(b1_imp_sum_c %>% filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  dplyr::select(Sample) %>% distinct() %>%
  mutate(obs = rownames(.)))

# highlight influential observations on scatter plots
bg_s_plot_infl <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      left_join(infl_obs %>% dplyr::select(!outlier_no)) %>%
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl'))
    
    plot <- ggplot(data = subset, 
                   aes(x = logBG, y = logS)) +
      geom_point(alpha = 0.6, size = 1, aes(colour = obs)) + 
      ggtitle(voc) +
      theme_bw() +
      theme(plot.title = element_text(size = 10),
            legend.position = 'none') +
      scale_colour_manual(values = c('infl' = 'red')) +
      geom_smooth(method = 'lm', se = FALSE, colour = 'grey50') +
      geom_smooth(data = subset %>% filter(obs %ni% c('infl')), 
                  method = 'lm', se = FALSE)
  })
}

cooks <- plot
dfbetas <- plot

del_diag_ex <- arrangeGrob(cooks, dfbetas, nrow = 2)
plot(del_diag_ex)
ggsave('Del_diag_examples.tiff', del_diag_ex, unit = 'mm', dpi = 300, width = 60, height = 120)

plots_infl <- bg_s_plot_infl(b1_imp_sum_c)

pdf('BGvsS_scatterPlots_PQN_CV12_BGfiltered_DelDiag.pdf')
marrangeGrob(plots_infl, nrow = 4, ncol = 4)
dev.off()

#

OD <- rbind(odOut %>% mutate(Sample = rownames(.), odS = 'Yes', odBG = 'No'),
            odOut %>% mutate(Sample = rownames(.), odS = 'No', odBG = 'Yes')) %>%
  dplyr::select(Sample, odS, odBG)

top_out <- table(pqn_diag$obs) %>% as.data.frame() %>%
  left_join(subset %>% mutate(Var1 = rownames(.)) %>% dplyr::select(Var1, Sample)) %>%
  left_join(OD, by = 'Sample')



#    
    
del_coefs <- infl[["fixed.effects[-case]"]] %>% as.data.frame()

del_coefs %>%
  ggplot(aes(x = logBG)) + geom_histogram()


max(del_coefs$logBG)
min(del_coefs$logBG)

cooks.distance(infl)
hist(dfbeta(infl)[,2])
hist(dfbetas(infl)[,2])

# dfbetas threshold
t <- 2/sqrt(142)
which(abs(dfbetas(infl)[,2]) > t)

#

infl1 <- influence(model)
coef_change <- infl1[["coefficients"]] %>% as.data.frame()
coef_loo <- data.frame(logBG = coef_change$logBG + coef(summary(model))[2,1])
hist(coef_loo$logBG)
confint(model)

PRESS <- sum(rstandard(model, type="pred")^2)

#

colnames(reg_results)[5] <- 'p.value'
reg_results$adj.p.value <- p.adjust(reg_results$p.value, method = 'BH')

reg_results_f <- reg_results %>% filter(adj.p.value < 0.05)

View(reg_results_f %>% left_join(wilcox_bg %>% dplyr::select(comp, Q1, Q3)))

#

#
#
#

## use linear regression coefficients to correct the background signal
b1_sub_corr1 <- b1_imp_sum_c %>% left_join(reg_results_f %>%
                                            dplyr::select(comp, Estimate)) %>%
  mutate(logS_adj = ifelse(is.na(Estimate) == TRUE, log(S),
                           (log(S) - Estimate*log(BG))))

# 

b1_sub_corr1_w <- b1_sub_corr1 %>% 
  dplyr::select(!c(S, BG, Estimate, Scorr)) %>%
  pivot_wider(names_from = comp, values_from = logS_adj) %>% 
  drop_na() %>%
  as.data.frame() %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

rownames(b1_sub_corr1_w) <- b1_sub_corr1_w$Sample

# keep only sample with both observations present
#c_reps <- b1_sub_corr1_w %>% group_by(RAD_ID) %>% summarise(n = n()) %>%
  #filter(n == 2)
#b1_sub_corr1_w <- b1_sub_corr1_w %>% filter(RAD_ID %in% c_reps$RAD_ID)

#

assay2 <- b1_sub_corr1_w[,-c(1:6)]


pc2 <- pca(assay2, center = TRUE, scale = TRUE, ncomp = 4)
           #multilevel = b1_sub_corr1_w$RAD_ID)
p2 <- plotIndiv(pc2,
          #pch = 1,
          #ind.names = b1_sub_corr1_w$RAD_ID,
          group = b1_sub_corr1_w$Diagnosis,
          legend = TRUE,
          title = 'Regression',
          comp = c(1,2),
          pch = as.factor(b1_sub_corr_w$CoreVisit),
          size.title = 10,
          cex = 1)


ggsave('PCA_S_Reg.tiff', p2$graph, unit = 'mm', dpi = 300, width = 120, height = 80 )

grid.arrange(p0$graph, p1$graph, p2$graph, ncol = 3, nrow = 1)

#
#
#

bg_s_hist <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(class = 'S')
    
    plot <- ggplot(data = subset, 
                   aes(x = logS_adj, fill = class)) +
      geom_histogram(alpha = 0.6, position = 'identity') + 
      ggtitle(voc) +
      theme_bw() +
      theme(plot.title = element_text(size = 8),
            legend.position = 'none') +
      xlab('log(peakArea)')
  })
}

plots_after <- bg_s_hist(b1_sub_corr1)

pdf('S_hist_after_correction.pdf')
marrangeGrob(plots_after, nrow = 4, ncol = 4)
dev.off()

#
#
#

# effect of BG correction method on on variance
# coefficient of variation calculated for each compound
library(PKNCA)

pqn_corr_cvs <- b1_imp_sum_c %>% 
  group_by(comp) %>%
  summarise(sdS = sd(S, na.rm = TRUE),
            cvS = sdS/abs(mean(S, na.rm = TRUE)),
            sdScorr = sd(Scorr, na.rm = TRUE),
            cvScorr = sdScorr/abs(mean(Scorr, na.rm = TRUE)),
            madS = mad(S, na.rm = TRUE),
            madScorr = mad(Scorr, na.rm = TRUE))


corr1_cvs <- b1_sub_corr1 %>% group_by(comp) %>%
  summarise(sd = sd(logS_adj, na.rm = TRUE),
            var = sd^2,
            cvGeo1 = sqrt(exp( (log(10))^2 * var ) - 1))

corr2_cvs <- b1_sub_corr1 %>% 
  mutate(S_adj = exp(logS_adj)) %>%
  group_by(comp) %>%
  summarise(sd = sd(S_adj, na.rm = TRUE),
            cvS1 = sd/abs(mean(S_adj, na.rm = TRUE)),
            madS1 = mad(S_adj, na.rm = TRUE))

all_cvs <- pqn_corr_cvs %>% 
  #left_join(corr1_cvs %>% dplyr::select(comp, cvGeo1)) %>%
  left_join(corr2_cvs %>% dplyr::select(comp, cvS1, madS1)) %>%
  left_join(S_log)
  

all_cvs %>% pivot_longer(cols = c(cvS1, cvS, cvScorr), 
                         names_to = 'Correction', 
                         values_to = 'CV') %>%
  ggplot(aes(x = Correction, y = CV, group = Correction)) + geom_boxplot(outliers = FALSE)

all_cvs %>% pivot_longer(cols = c(madS1, madS, madScorr), 
                         names_to = 'Correction', 
                         values_to = 'MAD') %>%
  ggplot(aes(x = Correction, y = MAD, group = Correction)) + geom_boxplot(outliers = FALSE)


# test for geometric CV
var1 <- 0.18^2
cvGeo <- function(var) {sqrt(exp( (log(10))^2 * var ) - 1)*100}
cvGeo(var1)

# geometric CV of log distribution of uncorrected data
S_log <- b1_imp_sum_c %>%
  mutate(logS = log(S)) %>%
  group_by(comp) %>%
  summarise(sd = sd(logS, na.rm = TRUE),
            var = sd^2,
            cvGeo = sqrt(exp( (log(10))^2 * var ) - 1))
  
sdLogS <- sd(az$logS, na.rm = TRUE)
varLogS <- sdLogS^2
cvGeo <- sqrt(exp( (log(10))^2 * varLogS ) - 1)




#

pqn_corr_cvs %>%  pivot_longer(cols = c(cvS, cvScorr), names_to = 'Sample', values_to = 'cv') %>%
  ggplot(aes(x = Sample, y = cv, group = Sample)) + geom_boxplot(outliers = FALSE) 

