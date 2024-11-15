## Background correction

# author: Aggie Turlo
# project: RADicA
# date: 30/09/2024

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
library(pcaPP)

# load normalised dataset with summarised breath samples
b1_imp_sum_c <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised.csv')[,-1]
b_imp_sum_c <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised.csv')[,-1]

pqn_imp_sum_c <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC_BGfiltered.csv')[,-1] %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

infl_obs <- read.csv()

#
#
#

# plot S vs BG log abundances
bg_s_plot <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      #filter(Sample %ni% c('190903_RaDICA_RAD032')) %>%
      filter(CoreVisit %in% c('CV1', 'CV2'))
    
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


plots <- bg_s_plot(b_imp_sum_c)

pdf('BGvsS_scatterPlots_PQN_CV12_BGfiltered_sign.pdf')
marrangeGrob(plots1, nrow = 4, ncol = 4)
dev.off()

#############################

# Filtering VOCs with similar BG and S distributions
# B2 - exclude VOCs with missing BG data
incomplete_s <- b_imp_sum_c[which(is.na(b_imp_sum_c$BG), arr.ind = TRUE),] 
bg_mis <- unique(incomplete_s$comp)

# Wilcoxon paired test for distribution of differences
vocs1 <- unique(b1_imp_sum_c$comp)
vocs <- vocs1[vocs1 %ni% bg_mis] # B2

wilcox_bg_fun <- function(data, voc_list, CV) {
  bind_rows(lapply(voc_list, function(voc) {
  subset1 <- data %>% 
    filter(comp == voc) %>%
    filter(CoreVisit == CV) %>%
    drop_na()
  test <- wilcox.test(subset1$S, subset1$BG, paired = TRUE,)
  df <- data.frame(comp = voc, 
                   p.value = test[["p.value"]],
                   min = min(subset1$S), 
                   max = max(subset1$S),
                   Q1 = quantile(subset1$S, 0.25), 
                   Q3 = quantile(subset1$S, 0.75),
                   logFC = mean(log(subset1$S) - mean(log(subset1$BG))))
  }))
  }


wilcox_bg_bcv1 <- wilcox_bg_fun(b_imp_sum_c, vocs, 'CV1') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

wilcox_bg_bcv2 <- wilcox_bg_fun(b_imp_sum_c, vocs, 'CV2') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

wilcox_bg_bcv1_f <- wilcox_bg_bcv1 %>% filter(adj.p.value > 0.05)
wilcox_bg_bcv2_f <- wilcox_bg_bcv2 %>% filter(adj.p.value > 0.05)

wilcox_b_agreement <- intersect(wilcox_bg_bcv1_f$comp, wilcox_bg_bcv2_f$comp)

View(wilcox_bg_bcv1 %>% filter(comp %ni% wilcox_b_agreement))

endo_b_cv1 <- wilcox_bg_bcv1 %>% filter(comp %ni% wilcox_b_agreement) %>%
  filter(logFC > 0.4)
endo_b_cv2 <- wilcox_bg_bcv2 %>% filter(comp %ni% wilcox_b_agreement) %>%
  filter(logFC > 0.4)

cons_b <- intersect(endo_b_cv1$comp, endo_b_cv2$comp)
setdiff(endo_b_cv1$comp, endo_b_cv2$comp)

#

wilcox_bg_b1cv1 <- wilcox_bg_fun(b1_imp_sum_c, vocs, 'CV1') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
wilcox_bg_b1cv2 <- wilcox_bg_fun(b1_imp_sum_c, vocs, 'CV2') %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

wilcox_bg_b1cv1_f <- wilcox_bg_b1cv1 %>% filter(adj.p.value > 0.05)
wilcox_bg_b1cv2_f <- wilcox_bg_b1cv2 %>% filter(adj.p.value > 0.05)

wilcox_b1_agreement <- intersect(wilcox_bg_b1cv1_f$comp, wilcox_bg_b1cv2_f$comp)

endo_b1_cv1 <- wilcox_bg_b1cv1 %>% filter(comp %ni% wilcox_b1_agreement) %>%
  filter(logFC > 0.4)
endo_b1_cv2 <- wilcox_bg_b1cv2 %>% filter(comp %ni% wilcox_b1_agreement) %>%
  filter(logFC > 0.4)

cons_b1 <- intersect(endo_b1_cv1$comp, endo_b1_cv2$comp)
setdiff(endo_b1_cv1$comp, endo_b1_cv2$comp)

intersect(cons_b, cons_b1)

# define VOCs that need removing from the dataset due to similarity with BG

all_vocs <- wilcox_bg_bcv1$comp
intersect(wilcox_b_agreement, wilcox_b1_agreement)

keep_wilcox_min <- intersect(all_vocs[all_vocs %ni% wilcox_b_agreement],
                             all_vocs[all_vocs %ni% wilcox_b1_agreement]) %>%
  as.data.frame() %>%
  rename(comp = '.') 

keep_wilcox_max <- all_vocs[all_vocs %ni% 
                              intersect(wilcox_b_agreement, wilcox_b1_agreement)] %>%
  as.data.frame() %>%
  rename(comp = '.') %>%
  mutate(wilcox_filter = ifelse(comp %in% keep_wilcox_min$comp, 'Both_datasets', 'One_dataset')) %>%
  mutate(FC_filter = ifelse(comp %in% intersect(cons_b, cons_b1), 'Endo_both_datasets',
                            ifelse(comp %in% setdiff(cons_b, cons_b1), 'Endo_one_dataset',
                                   ifelse(comp %in% setdiff(cons_b1, cons_b), 'Endo_one_dataset', 'Exo')))) 

table(keep_wilcox_max$FC_filter)

# remove VOCs where null hypothesis accepted in CV1 and CV2 (in B1 and B2)
b_imp_sum_c <- b_imp_sum_c %>% filter(comp %in% c(keep_wilcox_max$comp, bg_mis)) %>% 
  filter(comp != 'Internal_Standard') 

b1_imp_sum_c <- b1_imp_sum_c %>% filter(comp %in%  c(keep_wilcox_max$comp, bg_mis)) %>%
  filter(comp != 'Internal_Standard') 

n_distinct(b1_imp_sum_c$comp)

write.csv(b_imp_sum_c, 'RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')
write.csv(b1_imp_sum_c, 'RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')

#

b_imp_sum_c %>%
  filter(comp == 'Camphor') %>%
  ggplot() +
  geom_histogram(aes(x = log(BG)), alpha = 0.6, fill = 'blue') +
  geom_histogram(aes(x = log(S)), alpha = 0.6, fill = 'orange') +
  theme_bw() +
  xlab('log(peak area)') +
  ggtitle('B1')


#

View(as.data.frame(intersect(b_imp_sum_c$comp, unique(model1_out$comp))))
setdiff(b_imp_sum_c$comp, unique(model1_out$comp))
setdiff(unique(model1_out$comp), b_imp_sum_c$comp)


View(model1_out %>% filter(predictor == 'logBG') %>% filter(predictor == 'logBG') %>% 
       filter(comp %in% wilcox_agreement))

dot_plot <- function(voc) {
  b_imp_sum_c %>% filter(comp == voc) %>% 
  ggplot(aes(x = log(BG), y = log(S))) + 
  geom_point(aes(colour = as.integer(Analysis_date))) +
  theme_bw() +
  scale_colour_gradientn(name = 'Date', colours = c('yellow', 'blue'), labels = as.Date) +
  ggtitle(voc) +
    theme(plot.title = element_text(hjust = 0.5))}

cam <- dot_plot("X2_ethylfuran")
cam


p <- arrangeGrob(benzene, sevo, benzoic, ncol = 1)
plot(p)

ggsave('B2_wilcox_drift.tiff', p, dpi = 300, unit = 'mm', width = 90, height = 170)

setdiff(wilcox_agreement, bg_sign$comp)

dev.new()
b_all %>% dplyr::select(Internal_Standard, Sample, class) %>%
  left_join(meta %>% dplyr::select(Sample_ID, RAD_ID, CoreVisit),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(class %in% c('BG', 'S1', 'S2')) %>%
  pivot_wider(id_cols = c(RAD_ID, CoreVisit), 
              names_from = class, values_from = Internal_Standard) %>%
  ggplot(aes(x = log(BG))) + geom_point(aes(y = log(S1)), colour = 'blue') +
  geom_point(aes(y = log(S2))) +
  xlim(15, NA) + ylim(15, NA)

b_pqn1 %>% group_by(class) %>% summarise(mean = mean(Internal_Standard, na.rm = TRUE),
                                        sd = sd(Internal_Standard, na.rm = TRUE),
                                        cv = sd/mean)

# PCA for BG and S VOC levels following different normalisation methods
b_imp_c1 <- b_imp_c1 %>% mutate(S = (S1+S2)/2) %>%
  dplyr::select(!c(S1, S2)) 

b_imp_c1$Analysis_date <- as.Date(b_imp_c1$Analysis_date, format = '%d/%m/%Y')

b_imp_sum_c <- b_imp_c1

b_imp_sum_c_w <- b_imp_sum_c %>% pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  pivot_wider(names_from = comp, values_from = peakArea)

assay_sum_c <- b_imp_sum_c_w[,-c(1:6)]

pc_sum_cc <- pca(log(assay_sum_c), scale = TRUE, center = TRUE)

pc_sum_cc2$variates$X %>% as.data.frame() %>% 
  cbind(b_imp_sum_c_w$Analysis_date) %>%
  rename(Analysis_date = 'b_imp_sum_c_w$Analysis_date') %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(colour = as.integer(as.Date(Analysis_date)))) +
  scale_color_gradientn(name = 'Date', colours = c('yellow', 'blue'), labels = as.Date) +
  theme_bw() +
  ggtitle('PCA of breath and background B2 samples')

dev.new()
plotLoadings(pc_sum_c, comp = 2, ndisplay = 20)

ggsave('B2_PCA_CC2_drift.tiff', dpi = 300, width = 105, height = 60, unit = 'mm')

pc_sum_cc 
pc_sum_pqn
pc_sum_cc2

pqn_load <- pc_sum_pqn$loadings$X %>% as.data.frame()
cc_load <-  pc_sum_cc$loadings$X %>% as.data.frame()
cc2_load <- pc_sum_cc2$loadings$X %>% as.data.frame()

top_cc2 <- cc2_load %>% arrange(desc(abs(PC2))) %>% slice(1:20)
top_cc <- cc_load %>% arrange(desc(abs(PC2))) %>% slice(1:20)
top_pqn1 <- pqn_load %>% arrange(desc(abs(PC1))) %>% slice(1:20)

intersect(rownames(top_cc), rownames(top_cc2))
intersect(rownames(top_cc), rownames(top_pqn1))
intersect(rownames(top_cc2), rownames(top_pqn1))


#########################################################truehist()#
#############################

# BG correction based on regression modelling
# linear regression

# Baseline model
library(performance)

reg_coefs2 <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      drop_na() %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      mutate(logS = log(S), logBG = log(BG)) #%>%
      #left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      #mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      #filter(obs %ni% c('infl'))

    
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

# in test dataset (B2)

reg_results_b <- reg_coefs2(b_imp_sum_c)
reg_results_b1 <- reg_coefs2(b1_imp_sum_c)

colnames(reg_results_b)[5] <- 'p.value'

reg_results_b <- reg_results_b %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

reg_results_b_sign <- reg_results_b %>% filter(adj.p.value < 0.05)

# compare with results from B1 dataset normalised using PQN only
pqn_reg_b1 <- read.csv('PQN_regression_results_BG_filtered.csv')[,-1]
colnames(pqn_reg_b1)[5] <- 'p.value'

pqn_reg_b1 <- pqn_reg_b1 %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

pqn_reg_b1_sign <- pqn_reg_b1 %>% filter(adj.p.value < 0.05)

reg_comp_b1 <- pqn_reg_b1 %>% full_join(reg_results_b1, by = 'comp') %>%
  drop_na()

pqn_s <- reg_comp_b1 %>% filter(adj.p.value.x < 0.05)
ccpqn_s <- reg_comp_b1 %>% filter(adj.p.value.y < 0.05)

View(reg_comp_b1 %>% filter(comp %in% setdiff(ccpqn_s$comp, pqn_s$comp)))

plot(reg_comp_b1$adj.p.value.x, reg_comp_b1$adj.p.value.y) + abline(h = 0.05)  + abline(v = 0.05)

#

reg_comp_m1 <- reg_results1_b1 %>% 
  full_join(reg_results1_b %>% filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')), by = 'comp') %>%
  drop_na()

b1_s <- reg_comp_m1 %>% filter(adj.p.value.x < 0.05)
b_s <- reg_comp_m1 %>% filter(adj.p.value.y < 0.05)

View(reg_comp_m1 %>% filter(comp %in% setdiff(b_s$comp, b1_s$comp)))
intersect(b1_s$comp, b_s$comp)
  
intersect(reg_results_b1_sign$comp, pqn_reg_b1_sign$comp)
intersect(reg_results_b_sign$comp, pqn_reg_b1_sign$comp)



#pqn_imp_sum_c %>% 
b1_imp_sum_c %>%  
  filter(comp == 'X1_Propene._1_.methylthio._._.Z._') %>%
  ggplot(aes(x = log(BG), y = log(S))) +
  geom_point(aes(colour = as.integer(as.Date(Analysis_date)))) +
  scale_color_gradientn(name = 'Date', colours = c('yellow', 'blue'), labels = as.Date)  +
  #xlim(NA, 12.5) +
  theme_bw()

comp_regs_b1 <- pqn_reg_b1_sign %>% full_join(reg_results_b1_sign, by = 'comp')
comp_regs <- reg_results1_b_sign %>% full_join(reg_results_b1_sign, by = 'comp')

plot(comp_regs$Estimate.x, comp_regs$Estimate.y) + abline(0,1)

#

test <- reg_results_pqn %>% 
  arrange(p.value) %>%
  mutate(rank = c(1:152),
         fdr = 0.05*(rank/nrow(.)))
View(test %>% filter(p.value < fdr))

reg_results_pqn_sign <- reg_results_pqn %>% filter(adj.p.value < 0.05)

write.csv(reg_results_b, 'B2_CC_PQN_regression_results_BG_filtered.csv')

# Baseline model on data without influential observations
# infl observations detected based on the code below
reg_results_pqn_infl <- reg_coefs2(b1_imp_sum_c)

colnames(reg_results_pqn_infl)[5] <- 'p.value'

reg_results_pqn_infl <- reg_results_pqn_infl %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

reg_results_pqn_inflsign <- reg_results_pqn_infl_BG %>% filter(adj.p.value < 0.05)

#
#
#

##############################
##############################

# DELETION DIAGNOSTICS OF INFLUENTIAL OBSERVATIONS
del_diag <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) 
    
    
    model <- lmer(logS ~ logBG + (1 | RAD_ID),
                  data = subset)
    
    infl <- stats::influence(model)
    dfB <- stats::dfbetas(infl)[,2] %>% 
      as.data.frame() %>%
      mutate(comp = voc) %>%
      cbind(cooks.distance(model)) %>%
      rename(dfB = '.',
             cooksD = 'cooks.distance(model)')
    
    }))
}

b1_diag <- del_diag(b1_imp_sum_c)
write.csv(b1_diag, 'B1_CCPQN_deletion_diagnostics.csv')

# visualise deletion diagnostics results
library(forcats)
library(ggrepel)

b1_diag$obs <- gsub("\\..*","", rownames(b1_diag))

del_diag_cooks_fun <- function(data){ 
  lapply(unique(data$comp), function(voc) {
    data <- data %>% filter(comp == voc) %>%
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
  })}

del_diag_cooks1 <- del_diag_cooks_fun(b1_diag)

pdf('B1_Deletion_diagnostic_cooksD.pdf')
marrangeGrob(del_diag_cooks1, nrow = 4, ncol = 1)
dev.off()

#

del_diag_dfbetas_fun <- function(data){ 
  lapply(unique(data$comp), function(voc) {
    data <- data %>% filter(comp == voc) %>%
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
  })}

del_diag_dfbetas1 <- del_diag_dfbetas_fun(b1_diag)

pdf('B1_Deletion_diagnostic_dfbetas.pdf')
marrangeGrob(del_diag_dfbetas1, nrow = 4, ncol = 1)
dev.off()

#

# READ CSV FILE WITH INFLUENTIAL POINTS FLAGGED BY HAND 
infl_obs <- read.csv('B2_Deletion_diagnostics.csv')
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
  geom_vline(xintercept = 6, colour = 'red') +
  annotate('text',label = '5% of data', x = 7, y = 140, colour = 'red')

ggsave('Histogram_infl_obs_VOCs_B1_PQNCC.tiff', unit = 'mm', dpi = 300, width = 65, height = 50)  

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

ggsave('Histogram_infl_obs_B2.tiff', unit = 'mm', dpi = 300, width = 100, height = 60)

# add Sample IDs to influential observations database
infl_obs <- infl_obs %>%
  left_join(b_imp_sum_c %>% filter(CoreVisit %in% c('CV1', 'CV2')) %>%
              dplyr::select(Sample) %>% distinct() %>%
              mutate(obs = rownames(.)))

View(infl_obs %>% group_by(Sample) %>% summarise(n = n()))

write.csv(infl_obs, 'Deletion_diagnostics_formatted_B2.csv')
infl_obs1 <- read.csv('Deletion_diagnostics_formatted_B1_CCPQN.csv')[,-1]
infl_obs <- read.csv('Deletion_diagnostics_formatted_B2.csv')[,-1]

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

plots_infl <- bg_s_plot_infl(b_imp_sum_c)

pdf('B2_BGvsS_scatterPlots_PQN_CV12_BGfiltered_DelDiag.pdf')
marrangeGrob(plots_infl, nrow = 4, ncol = 4)
dev.off()

#
#
#

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

# confidence intervals
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

# save figures
library(cowplot)

comp_plots <- plot_grid(estimate_comp, CI_comp, r2_comp, nrow = 1, align = 'h',
                        rel_widths = c(0.3, 0.3, 0.4))

comp_plots

ggsave('Regression_comparison_del_diag.tiff' , comp_plots, dpi = 300, unit = 'mm',
       width = 260, height = 80)

reg_results_pqn_sign <- reg_results_pqn_sign %>% mutate(model = 'All data')
reg_results_pqn_infl_BG_sign <- reg_results_pqn_infl_BG_sign %>% mutate(model = 'Excl. outlier')

BGresults <- rbind(reg_results_pqn_infl_BG_sign, reg_results_pqn_sign)
p1 <- BGresults %>% ggplot(aes(x = Estimate, fill = model)) + 
  geom_histogram(position = 'identity', alpha = 0.5, bins = 12) +
  theme_bw() 
p2 <- BGresults %>% 
  mutate(CI_width = CI_upr - CI_lwr) %>%
  ggplot(aes(x = CI_width, fill = model)) + 
  geom_histogram(position = 'identity', alpha = 0.5, bins = 12) +
  theme_bw() 
p3 <- BGresults %>% 
  ggplot(aes(x = r2, fill = model)) + 
  geom_histogram(position = 'identity', alpha = 0.5, bins = 12) +
  theme_bw() 

p <- arrangeGrob(p1, p2, p3, nrow = 1)
plot(p)
ggsave('Histograms_reg_results_outliers.tiff', p, dpi = 300, unit = 'mm', width = 220, height = 40)


#
#
#

# MODEL 2 (with interaction)
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

# function for plotting individual VOCs
plot_voc <- function(voc){
  subset <- df %>%
    filter(comp == voc) %>%
    filter(CoreVisit %in% c('CV1', 'CV2')) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) 
  
  
  subset %>%
    #filter(Diagnosis == 'Not Asthma') %>%
    ggplot(aes(x = logBG, y = logS)) + 
    geom_point(size = 1, alpha = 0.8, aes(colour = Diagnosis)) +
    ylab('logS') + xlab('logBG') +
    theme_bw(base_size = 10) +
    #theme(legend.position = 'none') +
    #geom_label(aes(label = Sample)) +
    geom_smooth(method = 'lm', se = FALSE) +
    ggtitle(voc) +
    geom_abline(intercept = 0, slope = 1)
  }
  
dev.new()
plot_voc('Azetidine')
ggsave('Ethyl_Acetate.tiff', unit = 'mm', dpi = 300, width = 120, height = 80)
 
# 
#
#

# MODEL 1 (Final)
b_imp_sum_c <- b_imp_sum_c %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct()

reg_coefs1 <- function(df, voc_list) {
  bind_rows(lapply(voc_list, function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      #mutate(SlogS = logS/sd(logS, na.rm = TRUE), # code for standarised regression coefficients
             #SlogBG = logBG/sd(logBG, na.rm = TRUE)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
    
    model <- lmer(logS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID),
                  data = subset)
    
    #model <- lmer(SlogS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID), # code for standarised regression coefficients
                  #data = subset)
    
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

vocs <- unique(b_imp_sum_c$comp)
vocs <- vocs[vocs %ni% c('Pentanal', 'Methylthioacetate')]

reg_results1_b <- reg_coefs1(b_imp_sum_c, vocs)

colnames(reg_results1_b)[5] <- 'p.value'

#

reg_results1_b1_ccpqn_sign <- reg_results1_b1_ccpqn %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(adj.p.value < 0.05)

write.csv(reg_results1_b1_ccpqn, 'Model1_results_B1_CCPQN.csv')
reg_results1_b1 <- read.csv('Model1_results.csv')
write.csv(reg_results1_b, 'Model1_results_B2.csv')

asthma_pred_b1 <- reg_results1_b1_ccpqn %>% 
  filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

asthma_pred_b_sign <- asthma_pred_b %>% filter(p.value < 0.05)

hist(asthma_pred_b1$p.value, breaks = 20)

# FINAL MODEL
# compare B1 regression results depending on normalisation method
b1_reg <- reg_results1_b1_ccpqn %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  left_join(reg_results1_b1 %>% 
              filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')), 
            by = 'comp') %>%
  drop_na() %>%
  mutate(sign = ifelse(adj.p.value.x < 0.05 & adj.p.value.y < 0.05, 'both',
                       ifelse(adj.p.value.x < 0.05 & adj.p.value.y > 0.05, 'CC PQN', 
                              ifelse(adj.p.value.x > 0.05 & adj.p.value.y < 0.05, 'PQN', 'none'))))

table(b1_reg$sign)

View(b1_reg %>% filter(adj.p.value.y < 0.05))

tiff('Scatter_plot_compare_b1_model1_pvals.tiff', res = 300, units = 'mm', width = 110, height = 110)
plot(-log10(b1_reg$adj.p.value.x), -log10(b1_reg$adj.p.value.y), 
     col = as.factor(b1_reg$sign), cex.main = 0.9,
     xlab = 'CC2 PQN',  ylab= 'PQN',
     main = 'Adjusted -log10(p-value) related to the logBG effect') +
  abline(h = -log10(0.05), lty = 2) +
  abline(v = -log10(0.05), lty = 2) +
  abline(0, 1) +
  text(x = 10, y = 0.1, label = 6, col = 'red') +
  text(x = 0.1, y = 10, label = 16, col = 'blue') +
  text(x = 10, y = 3, label = 62)
dev.off() 

#
tiff('Scatter_plot_compare_b1_model1_coefs.tiff', res = 300, units = 'mm', width = 110, height = 110)
plot(b1_reg$Estimate.x, b1_reg$Estimate.y, 
     col = as.factor(b1_reg$sign),
     main = 'Regression coefficients related to logBG effect',
     xlab = 'CC2 PQN',  ylab= 'PQN', cex.main = 0.9) +
  abline(0, 1) +
  abline(v = 0, lty = 2) +
  abline(h = 0, lty = 2)
dev.off()

# compare B1 and B2
intersect(reg_results1_b_sign$comp, reg_results1_b1_ccpqn_sign$comp)

comp_b1b <- reg_results1_b1_ccpqn %>% filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  left_join(reg_results1_b%>% filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')), by = 'comp') %>%
  mutate(sign = ifelse(adj.p.value.x < 0.05 & adj.p.value.y < 0.05, 'both',
                       ifelse(adj.p.value.x < 0.05 & adj.p.value.y > 0.05, 'B1',
                              ifelse(adj.p.value.x > 0.05 & adj.p.value.y < 0.05, 'B2', 'none'))))

table(comp_b1b$sign)

tiff('Scatter_plot_compare_b1VSb_model1_pvals.tiff', res = 300, units = 'mm', width = 110, height = 110)
plot(-log10(comp_b1b$adj.p.value.x), -log10(comp_b1b$adj.p.value.y), 
     col = as.factor(comp_b1b$sign), cex.main = 0.9,
     xlab = 'B1',  ylab= 'B2',
     main = 'Adjusted -log10(p-value) related to the logBG effect') +
  abline(h = -log10(0.05), lty = 2) +
  abline(v = -log10(0.05), lty = 2) +
  abline(0, 1) +
  text(x = 10, y = 0.1, label = 34) +
  text(x = 0.1, y = 10, label = 29, col = 'red') +
  text(x = 10, y = 4, label = 54, col = 'green')
dev.off() 

#
tiff('Scatter_plot_compare_b1VSb_model1_coefs.tiff', res = 300, units = 'mm', width = 110, height = 110)
plot(comp_b1b$Estimate.x, comp_b1b$Estimate.y, 
     col = as.factor(comp_b1b$sign),
     main = 'Regression coefficients related to logBG effect',
     xlab = 'B1',  ylab= 'B2', cex.main = 0.9) +
  abline(0, 1) +
  abline(v = 0, lty = 2) +
  abline(h = 0, lty = 2)

dev.off()


#
#
#

# standarised regression
reg_results1_pqn_stand <- reg_coefs1(b1_imp_sum_c)

colnames(reg_results1_pqn_stand)[5] <- 'p.value'

asthma_pred_pqn_stand <- reg_results1_pqn_stand %>% filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

#
#
#

##############################
##############################

# MODEL SELECTION
# compare model fit with and without interaction term
# with and without random effect
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
#
#

# compare the effect of including covariates on logBG slope estimation
# Model 1 vs Baseline
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

# number of VOCs with significant effect of background across datasets
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

# slopes
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

# confidence intervals
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

# save figures
baseline_vs_m1 <- arrangeGrob(est, ci, nrow = 1)
plot(baseline_vs_m1)
ggsave('Baseline_vs_M1_Est_CI.tiff', baseline_vs_m1, dpi = 300, unit = 'mm', width = 160, height = 80)


# visualise effect of asthma based on confidence and effect size
# calculation of fold change (not useful?)
hist(asthma_pred_pqn$p.value, breaks = 30)
sd(asthma_pred_pqn$p.value)

asthma_pred_pqn <- asthma_pred_pqn %>%
  mutate(lab = ifelse(p.value < 0.05 & Estimate > 0.4, comp, ifelse(
    p.value < 0.05 & Estimate < -0.4, comp, NA))) %>%
  mutate(col = ifelse(is.na(lab) == TRUE, 'no', 'yes'))

vol_plot <- asthma_pred_pqn %>%
  ggplot(aes(y = -log10(p.value), x = Estimate)) +
  geom_point(aes(colour = col)) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = 'dashed', colour = 'darkgrey') +
  geom_hline(yintercept = 1.3, linetype = 'dashed', colour = 'darkgrey') +
  theme_bw() +
  geom_text_repel(aes(label = lab), size = 2) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c('yes' = 'red',
                                 'no' = 'black')) +
  ggtitle('Volcano plot showing effect of asthma diagnosis on breath VOCs') +
  xlim(-1, 1) +
  xlab('partial regression coefficient')

vol_plot

ggsave('Volcano_asthma.tiff', dpi = 300, unit = 'mm', width = 150, height = 120)

##############################
##############################

#
#
#

# APPLYING BG CORRECTION
## use linear regression coefficients to correct the background signal
b1_imp_sum_corr <- b1_imp_sum_c %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  left_join(reg_results1_b %>%
              filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
              dplyr::select(comp, Estimate, '2.5 %', '97.5 %', adj.p.value)) %>%
  mutate(logS_adj = ifelse(adj.p.value > 0.05 | is.na(adj.p.value), log(S), # account for VOCs with missing BG info
                           log(S) - Estimate*log(BG)))

colnames(b1_imp_sum_corr)[10] <- 'CI_lwr'
colnames(b1_imp_sum_corr)[11] <- 'CI_upr'

#
b_imp_sum_corr <- b_imp_sum_c %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  left_join(reg_results1_b %>%
              filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
              dplyr::select(comp, Estimate, '2.5 %', '97.5 %', adj.p.value)) %>%
  mutate(logS_adj = ifelse(adj.p.value > 0.05 | is.na(adj.p.value), log(S),
                           (log(S) - Estimate*log(BG))))

colnames(b_imp_sum_corr)[10] <- 'CI_lwr'
colnames(b_imp_sum_corr)[11] <- 'CI_upr'

# exclude VOCs where 95% CI for logBG coefficient includes 1 
excl <- b1_imp_sum_corr %>%
  filter(is.na(Estimate) == FALSE & CI_upr > 1 ) 

b1_imp_sum_corr <- b1_imp_sum_corr %>%
  filter(comp %ni% excl$comp)
b_imp_sum_corr <- b_imp_sum_corr %>%
  filter(comp %ni% excl$comp)

# BG-CORRECTED OUTPUT
write.csv(b1_imp_sum_corr, 'RADicA_B1_BG_adjusted.csv')
write.csv(b_imp_sum_corr, 'RADicA_B2_BG_adjusted.csv')

#
#
#

# Principal Component Analysis of corrected and uncorrected breath data
# transform to wide format
pivot_wide <- function(data, column){
  data_w <- data %>% 
  left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>%
  distinct() %>%
  dplyr::select(c(Sample, comp, logS_adj, Diagnosis, CoreVisit, RAD_ID)) %>%
  pivot_wider(names_from = comp, values_from = logS_adj) %>% 
  as.data.frame()
}

b1_imp_sum_corr_w <- pivot_wide(data = b1_imp_sum_corr, column = 'logS_adj')
rownames(b1_imp_sum_corr_w) <- b1_imp_sum_corr_w$Sample

b_imp_sum_corr_w <- pivot_wide(data = b_imp_sum_corr, column = 'logS_adj')
rownames(b_imp_sum_corr_w) <- b_imp_sum_corr_w$Sample

# exclude VOCs considered exogenous based on logFC
exo <- reg_valid_join %>% filter(FC_filter == 'Exo')
exo <- unique(exo$comp)

b1_endo <- b1_imp_sum_corr_w %>% dplyr::select(!exo)
b2_endo <- b_imp_sum_corr_w %>% dplyr::select(!exo)

# 
pca_cv12 <- function(data){
  assay <- data[,-c(1:4)]
  pc <- pca(assay, center = TRUE, scale = TRUE, ncomp = 4)
  #multilevel = b1_sub_corr1_w$RAD_ID)
  pc_plot <- plotIndiv(pc,
                #pch = 1,
                #ind.names = data$RAD_ID,
                group = data$Diagnosis,
                legend = TRUE,
                comp = c(1,2),
                pch = as.factor(data$CoreVisit),
                size.title = 10,
                title = 'Breath sample VOCs',
                cex = 3)
}

pc_plot1 <- pca_cv12(b1_imp_sum_corr_w1)
pc_plot_unc1 <- pca_cv12(b1_imp_sum_w)
pc_plot1_endo <- pca_cv12(b1_endo)

pc_plot <- pca_cv12(b_imp_sum_corr_w1)
pc_plot_unc <- pca_cv12(b_imp_sum_w)
pc_plot_endo <- pca_cv12(b2_endo)

# Each CV separately
pca_cv1 <- function(data){
  cv1 <- data %>% filter(CoreVisit == 'CV1')
  assay1 <- cv1[,-c(1:4)]
  pc1 <- pca(assay1, center = TRUE, scale = TRUE, ncomp = 4)
  pc_plot1 <- plotIndiv(pc1,
                     pch = 1,
                     #ind.names = cv1$Sample,
                     group = cv1$Diagnosis,
                     #legend = TRUE,
                     comp = c(1,2),
                     size.title = 10,
                     title = 'Core Visit 1',
                     cex = 1.5)
}

pc_plot1_b1 <- pca_cv1(b1_imp_sum_corr_w)
pc_plot1_unc_b1 <- pca_cv1(b1_imp_sum_w)

pc_plot1 <- pca_cv1(b_imp_sum_corr_w)
pc_plot1_unc <- pca_cv1(b_imp_sum_w)

#

pca_cv2 <- function(data) {
  cv2 <- data %>% filter(CoreVisit == 'CV2')
  assay2 <- cv2[,-c(1:4)]
  pc2 <- pca(assay2, center = TRUE, scale = TRUE, ncomp = 4)
  pc_plot2 <- plotIndiv(pc2,
                     pch = 1,
                     #ind.names = cv2$Sample,
                     group = cv2$Diagnosis,
                     #legend = TRUE,
                     comp = c(1,2),
                     size.title = 10,
                     title = 'Core Visit 2',
                     cex = 1.5)
}

#

pc_plot2_b1 <- pca_cv2(b1_imp_sum_corr_w)
pc_plot2_unc_b1 <- pca_cv2(b1_imp_sum_w)

pc_plot2_b2 <- pca_cv2(b_imp_sum_corr_w)
pc_plot2_unc_b2 <- pca_cv2(b_imp_sum_w)

#

pca_plots <- plot_grid(pc_plot$graph, pc_plot1$graph, pc_plot2$graph, nrow = 1, align = 'h',
                       rel_widths = c(0.41, 0.295, 0.295))
pca_plots

ggsave('PCA_regression_corrected.tiff', pca_plots, dpi = 300, unit = 'mm', width = 300, height = 80)

#

pca_plots_unc <- plot_grid(pc_plot_unc$graph, pc_plot1_unc$graph, pc_plot2_unc$graph, nrow = 1, align = 'h',
                       rel_widths = c(0.41, 0.295, 0.295))
pca_plots_unc

ggsave('PCA_regression_uncorrected.tiff', pca_plots_unc, dpi = 300, unit = 'mm', width = 300, height = 80)


#

out <- infl_obs %>% group_by(Sample) %>% summarise(n = n())

out_s <- data.frame(Sample = c('191206_RaDICA_RAD078', '191018_RaDICA_RAD056', '190823_RaDICA_RAD022', 
           '190828_RaDICA_RAD016', '190828_RaDICA_RAD022', '190813_RaDICA_RAD012'),
           obs = c(61, 44, 20, 21, 23, 15))

out_s <- out_s %>% left_join(out)

#
#
#

# OUTLIER DETECTION
# detect multivariate outliers using robust PCA
data <- b1_endo

assay <- data[,-c(1:4)]
Rpc <- PCAgrid(assay, scale = 'sd', center = 'mean')
sdod <- PCdiagplot(assay, Rpc, plotbw = FALSE, crit = 0.98)
# identify outliers based on the critical OD and SD values 
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(assay)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,1] & V2 > sdod$critSD[2,1])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(assay)
odOut <- od %>% filter(V1 > sdod$critOD[1,1] & V2 > sdod$critOD[2,1])

outliers <- c(rownames(odOut), rownames(sdOut)) %>% unique()

distances <- od %>% 
  mutate(Sample = rownames(.)) %>%
  pivot_longer(cols = c(V1, V2), names_to = 'k', values_to = 'OD') %>%
  left_join(sd %>% 
              mutate(Sample = rownames(.)) %>%
              pivot_longer(cols = c(V1, V2), names_to = 'k', values_to = 'SD')) %>%
  left_join(as.data.frame(data$Sample) %>%
            rename(Sample = 'data$Sample') %>%
              mutate(obs = rownames(.))) %>%
  mutate(outliers = ifelse(Sample %in% outliers, 'yes', 'no'),
         outlier_obs = ifelse(outliers == 'yes', Sample, NA))

distances %>%
  filter(k == 'V1') %>%
  ggplot(aes(x = OD, y = SD, colour = outliers)) +
  geom_point() +
  geom_hline(yintercept = sdod$critSD[1,1], linetype = 'dashed') +
  geom_vline(xintercept = sdod$critOD[2,1], linetype = 'dashed') +
  theme_bw() +
  geom_text_repel(aes(label = outlier_obs), colour = 'black') +
  ylab('score distance') + xlab('orthogonal distance') +
  ggtitle('Diagnostic plot to the PCA \n on corrected breath sample VOCs \n (Dataset 1)') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c('no' = 'black',
                                 'yes' = 'red'))


ggsave('PCA_diagnostic_plot_B1.tiff', dpi = 300, unit = 'mm', width = 100, height = 80)

View(b2_endo %>% filter(Sample %in% outliers))

# PCA following outlier removal
b1_endo1 <- b1_endo %>% filter(Sample %ni% outliers)
b2_endo1 <- b2_endo %>% filter(Sample %ni% outliers)

b1_imp_sum_corr_w1 <- b1_imp_sum_corr_w %>% filter(Sample %ni% outliers)
b_imp_sum_corr_w1 <- b_imp_sum_corr_w %>% filter(Sample %ni% outliers)

b1_endo1 <- b1_imp_sum_corr_w1 %>% dplyr::select(!exo)
b2_endo1 <- b_imp_sum_corr_w1 %>% dplyr::select(!exo)

write.csv(b1_imp_sum_corr_w1, 'RADicA_BG_adjusted_B1_outl_removed.csv')
write.csv(b_imp_sum_corr_w1, 'RADicA_BG_adjusted_B2_outl_removed.csv')

pca_cv12 <- function(data){
  assay <- data[,-c(1:4)]
  pc <- pca(assay, center = TRUE, scale = TRUE, ncomp = 4)
  #multilevel = b1_sub_corr1_w$RAD_ID)
  pc_plot <- plotIndiv(pc,
                       #pch = 1,
                       ind.names = data$RAD_ID,
                       group = data$Diagnosis,
                       legend = TRUE,
                       comp = c(1,2),
                       pch = as.factor(data$CoreVisit),
                       size.title = 10,
                       title = 'Dataset 2',
                       cex = 1.5)
}



dev.new()
pc_plot_out_b1 <- pca_cv12(b1_imp_sum_corr_w1)
pc_plot_out_b1 <- pca_cv12(b1_endo1)

pc_plot_out <- pca_cv12(b_imp_sum_corr_w1)
pc_plot_out <- pca_cv12(b2_endo1)

pc_plots_out <- arrangeGrob(pc_plot_out_b1$graph, pc_plot_out$graph, ncol = 1,
                            top = textGrob('PCA of breath sample VOCs'))
plot(pc_plots_out)
ggsave('PCA_corr_out_B1B2.tiff', pc_plots_out, dpi = 300, unit = 'mm', width = 120, height = 160)



#
pc_corr_b1 <- pca(b1_imp_sum_corr_w1[,-c(1:4)], scale = TRUE, center = TRUE)
pc_corr_b2 <- pca(b_imp_sum_corr_w1[,-c(1:4)], scale = TRUE, center = TRUE) 

pc_corr_b1_endo <- pca(b1_endo1[,-c(1:4)], scale = TRUE, center = TRUE)
pc_corr_b2_endo <- pca(b2_endo1[,-c(1:4)], scale = TRUE, center = TRUE) 

score_boxplot <- function(pc, data, PC){
  scores <- {{pc}}$variates$X %>% as.data.frame() %>% mutate(Sample = rownames(.)) %>%
    left_join(data %>% dplyr::select(Sample, Diagnosis)) %>%
    dplyr::select(Sample, Diagnosis, PC)
  colnames(scores)[3] <- 'PC'
  
  box <- scores %>% ggplot(aes(x = Diagnosis, y = PC)) + geom_boxplot() +
    ylab(PC)
  plot(box)
  
  scores %>% mutate(PCpos = ifelse(PC > 0, 'yes', 'no')) %>%
  group_by(PCpos, Diagnosis) %>%
  summarise(n = n())

  }

score_boxplot(pc_corr_b2, b_imp_sum_corr_w1, 'PC2')
score_boxplot(pc_corr_b2_endo, b2_endo1, 'PC1')

loads_fun <- function(pc, data, PC) {
  loads <- {{pc}}$loadings$X %>% as.data.frame() %>% mutate(comp = rownames(.)) %>%
  dplyr::select(PC, comp) %>%
  left_join(reg_valid_join %>% dplyr::select(comp, FC_filter, wilcox_filter)) %>%
  distinct()
  colnames(loads)[1] <- 'PC'
  
  loads %>% arrange(desc(abs(PC))) %>%
  slice(1:15) %>%
  ggplot(aes(x = PC, y = fct_inorder(as.factor(comp)), fill = FC_filter)) + geom_col() +
  theme_bw() + scale_fill_brewer(palette = 'Set2') +
  xlab(PC)}

loads_fun(pc_corr_b1, b1_imp_sum_corr_w1, 'PC1') + ggtitle('B1 all')
dev.new()
loads_fun(pc_corr_b1_endo, b1_endo1, 'PC1') + ggtitle('B1 endo')

dev.new()
loads_fun(pc_corr_b2, b_imp_sum_corr_w1, 'PC2') + ggtitle('B2 all')
loads_fun(pc_corr_b2_endo, b2_endo1, 'PC1') + ggtitle('B2 endo')

ggsave('PCA_loadings_B2_corr.tiff', b2_pc2, dpi = 300, unit = 'mm', width = 130, height = 80)

#
#
#

# hierarchical clustering of corrected breath samples
library(stats)

names <- paste(b1_imp_sum_corr_w1$RAD_ID, b1_imp_sum_corr_w1$CoreVisit, b1_imp_sum_corr_w1$Diagnosis,
               sep = '_')

rownames(b1_imp_sum_corr_w1) <- names

d_matrix <- dist(b1_imp_sum_corr_w1[,-c(1:4)], method = 'euclidean')
clust <- hclust(d_matrix)
dev.new()
plot(clust)


#
#
#

# MODEL VALIDATION
# validation of regression model prediction on B1
b12_voc <- unique(b_imp_sum_c$comp)

b_imp_sum_c <- b_imp_sum_c %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>%
  distinct()

b1_imp_sum_c <- b1_imp_sum_c %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>%
  distinct()
                              
reg_valid_out <- bind_rows(lapply(b12_voc, function(voc) {
  sub_train <- b_imp_sum_c %>%
      filter(comp == voc) %>%
      filter(CoreVisit %in% c('CV1', 'CV2')) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      #mutate(SlogS = logS/sd(logS, na.rm = TRUE), # code for standarised regression coefficients
      #SlogBG = logBG/sd(logBG, na.rm = TRUE)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
  
  sub_test <- b1_imp_sum_c %>%
    filter(comp == voc) %>%
    filter(CoreVisit %in% c('CV1', 'CV2')) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) %>%
    distinct() 
    
    model <- lmer(logS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID),
                  data = sub_train)
    
    #model <- lmer(SlogS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID), # code for standarised regression coefficients
    #data = subset)
    
    pred_test <- sub_test %>% 
      mutate(logSpred = predict(model, newdata = sub_test, allow.new.levels = TRUE),
             error = logSpred - logS)
    
    perf <- data.frame(mse = mean((pred_test$error)^2),
                       mae = mean(abs(pred_test$error)),
                       me = mean(pred_test$error),
                       mape = mean(abs(pred_test$error/pred_test$logSpred)*100),
                       comp = voc)
    
}))

#

reg_valid_out1 <- reg_valid_out %>% 
  left_join(reg_results1_b %>% filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')))

reg_valid_out_sign <- reg_valid_out1 %>% filter(adj.p.value < 0.05)

# error for the baseline model
reg_valid_base1 <- reg_valid_base %>% 
  left_join(reg_results1_b %>% filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')))

reg_valid_base_sign <- reg_valid_base1 %>% filter(adj.p.value < 0.05)

plot(reg_valid_base_sign$mae, reg_valid_out_sign$mae) + abline(0,1)

hist(reg_valid_base_sign$mae - reg_valid_out_sign$mae)

median(reg_valid_base_sign$mae)
median(reg_valid_out_sign$mae)
median(reg_valid_null_sign$mae)


# error for the null model
reg_valid_null <- bind_rows(lapply(b12_voc, function(voc) {
  sub_train <- b_imp_sum_c %>%
    filter(comp == voc) %>%
    filter(CoreVisit %in% c('CV1', 'CV2')) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) %>%
    #mutate(SlogS = logS/sd(logS, na.rm = TRUE), # code for standarised regression coefficients
    #SlogBG = logBG/sd(logBG, na.rm = TRUE)) %>%
    left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
    mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
    filter(obs %ni% c('infl'))
  
  sub_test <- b1_imp_sum_c %>%
    filter(comp == voc) %>%
    filter(CoreVisit %in% c('CV1', 'CV2')) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) %>%
    distinct() 
  
  model <- lmer(logS ~ 1 + (1 | RAD_ID),
                data = sub_train)
  
  #model <- lmer(SlogS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID), # code for standarised regression coefficients
  #data = subset)
  
  pred_test <- sub_test %>% 
    mutate(logSpred = predict(model, newdata = sub_test, allow.new.levels = TRUE),
           error = logSpred - logS)
  
  perf <- data.frame(mse = mean((pred_test$error)^2),
                     mae = mean(abs(pred_test$error)),
                     me = mean(pred_test$error),
                     mape = mean(abs(pred_test$error/pred_test$logSpred)*100),
                     comp = voc)
  
}))


reg_valid_null1 <- reg_valid_null %>% 
  left_join(reg_results1_b %>% filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')))

reg_valid_null_sign <- reg_valid_null1 %>% filter(adj.p.value < 0.05)

plot(reg_valid_out$mape, reg_valid_null$mape) + abline(0,1)


hist(reg_valid_null$mae - reg_valid_out$mae)
hist(reg_valid_null$mape - reg_valid_out$mape)

#

# compare error to p value and regression estimate
# effect of different filtering strategies for endogenous VOCs on regression results
reg_valid_out <- reg_valid_out %>% mutate(model = 'Model1')
reg_valid_base <- reg_valid_base %>% mutate(model = 'Baseline')
reg_valid_null <- reg_valid_null %>% mutate(model = 'Null')

reg_valid_join <- rbind(reg_valid_out, #reg_valid_base, 
                        reg_valid_null) %>%
  full_join(reg_results1_b %>% 
              filter(predictor == 'logBG') %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH'))) %>%
  mutate(sign = ifelse(adj.p.value < 0.05, 'yes', 'no')) %>%
  left_join(keep_wilcox_max)

write.csv(reg_valid_join, 'RegValid_results.csv')

# effect on regression results (dataset 2)
# p-value
reg_valid_join %>% group_by(wilcox_filter, sign) %>% summarise(n = n())

  reg_valid_join %>%
  filter(model == 'Model1') %>%
  mutate(adj.p.value = -log10(adj.p.value)) %>%
  pivot_longer(cols = c(r2, adj.p.value), names_to = 'parameter', values_to = 'value') %>%
  ggplot(aes(y = value, x = wilcox_filter, colour = FC_filter)) +
  facet_wrap(~ parameter, scale = 'free') +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(position = position_dodge(width = 0.9), width = 0.15, outliers = FALSE) +
  theme_bw(base_size = 10) +
  geom_hline(data = data.frame(yint = -log10(0.05), 
                               parameter = "adj.p.value"), 
             aes(yintercept = yint), linetype = "dashed") +
  ggtitle('Model results in Dataset 2 (training)') +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_colour_brewer(palette = 'Set2') +
    ylab('-log10(adj.p.value)')

ggsave('RegRes_vs_filtering.tiff', dpi = 300, unit = 'mm', width = 170, height = 60)


# effect on prediction error
reg_valid_join %>%
  filter(adj.p.value < 0.05) %>%
  filter(model == 'Model1') %>%
  filter(mae < 4) %>%
  pivot_longer(cols = c(mae, mape), names_to = 'error', values_to = 'value') %>%
  ggplot(aes(y = value, x = wilcox_filter, colour = FC_filter)) +
  facet_wrap(~error, scale = 'free') +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw(base_size = 10) +
  ggtitle('Model prediction error') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))  +
  scale_colour_brewer(palette = 'Set2')

ggsave('RegPred_vs_filtering.tiff', dpi = 300, unit = 'mm', width = 170, height = 60)

# effect on difference in prediction error with null model
reg_valid_join %>%
  filter(adj.p.value < 0.05) %>%
  filter(mae < 4) %>%
  pivot_longer(cols = c(mae, mape), names_to = 'error', values_to = 'value') %>%
  dplyr::select(error, value, wilcox_filter, FC_filter, comp, model) %>%
  pivot_wider(names_from = model, values_from = value) %>%
  ggplot(aes(y = Null-Model1, x = wilcox_filter, colour = FC_filter)) +
  facet_wrap(~error, scale = 'free') +
  geom_violin(position = position_dodge(width = 0.9)) +
  geom_boxplot(outliers = FALSE, position = position_dodge(width = 0.9), width = 0.2) +
  theme_bw(base_size = 10) +
  geom_hline(yintercept = 0, linetype = 'dashed')  +
  ggtitle('Prediction error difference with null model') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))  +
  scale_colour_brewer(palette = 'Set2')

ggsave('RegPredDiff_vs_filtering.tiff', dpi = 300, unit = 'mm', width = 170, height = 60)

#
#
#

b_imp_sum_w <- b_imp_sum_c %>% 
  pivot_longer(cols = c(BG, S), names_to = 'class', values_to = 'peakArea') %>%
  pivot_wider(names_from = comp, values_from = peakArea)

b_pc <- pca(log(b_imp_sum_w[, -c(1:7)]), scale = TRUE, center = TRUE)
scores <- as.data.frame(b_pc$variates$X) %>% cbind(b_imp_sum_w[, c(1:7)])

b_pc$prop_expl_var

scores %>% ggplot(aes(x = PC1, y = PC2, colour = class)) + geom_point()

scores %>% ggplot(aes(x = PC1, y = PC2, colour = as.integer(as.Date(Analysis_date)),
                      shape = class)) + 
  geom_point() +
  scale_colour_gradientn(name = 'Date', 
                         colours = c('yellow', 'blue'), labels = as.Date) +
  theme_bw() +
  ggtitle('CC (Blank) + PQN') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('PC1 (19.6%)') +
  ylab('PC2 (8.6%)')

ggsave('B2_PCA_time_BG_S.tiff', unit = 'mm', dpi = 300, width = 120, height = 80)

dev.new()
plotLoadings(b_pc, comp = 2, ndisplay = 10)

scatp <- function(voc){
  b_imp_sum_c %>%
  filter(comp == voc) %>%
  ggplot(aes(x = log(BG), y = log(S), colour = as.integer(as.Date(Analysis_date)))) +
  geom_point() +
  scale_colour_gradientn(name = 'Date', colours = c('yellow', 'blue'), labels = as.Date) +
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle(voc) +
    coord_fixed()
  }

scatp('Undecane')

scats <- lapply(rownames(top_loads), scatp)
marrangeGrob(scats, nrow = 3, ncol = 3)

hist_voc <- function(voc) {
  b_imp_sum_c %>%
  filter(comp == voc) %>%
  ggplot() +
  geom_histogram(aes(x = log(S)), bins = 15, alpha = 0.6, fill = 'yellow3')
  geom_histogram(aes(x = log(BG)), alpha = 0.6, bins = 15, fill = 'blue') +
  theme_bw() +
  ggtitle(voc) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = 'none')
  }

top_loads <- as.data.frame(b_pc$loadings$X[,2]) %>%
  rename(load = 'b_pc$loadings$X[, 2]') %>%
  arrange(desc(abs(load))) %>%
  slice(1:9)

hists <- lapply(rownames(top_loads) , hist_voc)
marrangeGrob(hists, nrow = 3, ncol = 3)

hist_voc('Cyclohexane')

timep <- function(voc){
  b_imp_sum_c %>%
  filter(comp == voc) %>%
  group_by(Analysis_date) %>%
  summarise(BG = median(log(BG)),
            S = median(log(S))) %>%
  pivot_longer(cols = c(BG, S), names_to = 'sample', values_to = 'med_peak') %>%
  ggplot(aes(x = as.Date(Analysis_date), y = med_peak, colour = sample)) +
  geom_point() +
  #geom_line(aes(group = sample)) +
  theme_bw() +
  xlab('Analysis_date') +
    ylab('median log(peakArea)') +
  ggtitle(voc)
  }

timep('X3_Carene')
timeps <- lapply(rownames(top_loads) , timep)
marrangeGrob(timeps, nrow = 3, ncol = 1)

View(as.data.frame(b_pass))
