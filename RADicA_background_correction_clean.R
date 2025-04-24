## Background correction

# author: Aggie Turlo
# project: RADicA
# date: 27/01/2025

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
library(clusrank)
library(performance)
library(forcats)
library(ggrepel)
library(VennDiagram)
library(shades)

## load data
# load normalised dataset with summarised breath samples
b1_imp_sum <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]
b2_imp_sum <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]

b1_imp_sum$Analysis_date <- as.Date(b1_imp_sum$Analysis_date)
b2_imp_sum$Analysis_date <- as.Date(b2_imp_sum$Analysis_date)

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

# load list of outliers in background-breath linear model (dataset 2)
infl_obs <- read.csv('Deletion_diagnostics_formatted_B2.csv')[,-1]

# load VOC annotation based on logFC filtering
endo_exo <- read.csv('Endo_Exo_filters.csv')[,-1]

# custom functions
'%ni%' <- Negate('%in%')


#
#
#


#############################

## BACKGROUND CORRECTION
# Baseline model of the relationship between breath and background
reg_base <- function(df) {
  bind_rows(test <- lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))

    
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


base_results_b2 <- reg_base(b2_imp_sum)

colnames(base_results_b2)[5] <- 'p.value'

write.csv(base_results_b2, 'Baseline_BG_model_results_B2.csv')

base_results_b2 <- read.csv('Baseline_BG_model_results_B2.csv')[,-1]

base_results_b2 <- base_results_b2 %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))


#
#
#

###############################

# DELETION DIAGNOSTICS OF INFLUENTIAL OBSERVATIONS
del_diag <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
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

b2_diag <- del_diag(b2_imp_sum)

write.csv(b2_diag, 'B2_deletion_diagnostics.csv')

# visualise deletion diagnostics results
b2_diag$obs <- gsub("\\..*","", rownames(b2_diag))

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

del_diag_cooks1 <- del_diag_cooks_fun(b2_diag)

pdf('B2_Deletion_diagnostic_cooksD.pdf')
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

del_diag_dfbetas1 <- del_diag_dfbetas_fun(b2_diag)

pdf('B2_Deletion_diagnostic_dfbetas.pdf')
marrangeGrob(del_diag_dfbetas1, nrow = 4, ncol = 1)
dev.off()

#

# load file of influential observations selected based on inspection of deletion diagnostics graphs 
infl_obs_raw <- read.csv('B2_deletion_diagnostics_raw.csv')

# format the file to annotate observations with sample IDs
infl_obs_raw1 <- cbind(infl_obs_raw[,1], str_split(infl_obs_raw$Potential_outlier, pattern = ',', 
                                                   n = 6, simplify = TRUE)) %>% as.data.frame()

colnames(infl_obs_raw1) <- c('comp', 'o1', 'o2', 'o3', 'o4', 'o5', 'o6')

infl_obs_raw_l <- infl_obs_raw1 %>% pivot_longer(cols =! comp, names_to = 'outlier_no', values_to = 'obs') %>%
  filter(obs != '') %>%
  mutate(obs = as.numeric(obs)) 

obs_id <- b2_imp_sum %>% filter(comp == '.beta._Myrcene') %>%
  mutate(obs = as.numeric(rownames(.))) %>% dplyr::select(obs, Sample)

infl_obs <- infl_obs_raw_l %>% left_join(obs_id)

write.csv(infl_obs, 'Deletion_diagnostics_formatted_B2.csv')


#

# visualise % of data removed across VOCs
infl_obs %>% group_by(comp) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = n)) +
  geom_histogram(binwidth = 0.5) +
  xlab('Number of influential observations') +
  ggtitle('Distribution of VOCs according to \n the influential observation count') +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(size = 10)) +
  ylab('VOC number') +
  geom_vline(xintercept = 6, colour = 'red') +
  annotate('text',label = '5% of data', x = 4.8, y = 40, colour = 'red', size = 3)

ggsave('Histogram_infl_obs_VOCs_B2.tiff', unit = 'mm', dpi = 300, width = 65, height = 50)  

infl_sum <- infl_obs %>% group_by(comp) %>% summarise(n = n())
table(infl_sum$n)

#

# visualise number of VOCs where observation was flagged as outlier
table(infl_obs$obs) %>% as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  ggplot(aes(x = fct_inorder(Var1), y = Freq)) + 
  geom_col() +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank()) +
  xlab('observation ID') +
  ylab('VOC number') +
  theme(plot.title = element_text(size = 9)) +
  ggtitle('Number of VOC where observation was flagged as influential')

ggsave('Histogram_infl_obs_B2.tiff', unit = 'mm', dpi = 300, width = 100, height = 60)

#

# highlight influential observations on scatter plots
bg_s_plot_infl <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
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


plots_infl <- bg_s_plot_infl(b2_imp_sum)

pdf('B2_BGvsS_scatterPlots_outliers.pdf')
marrangeGrob(plots_infl, nrow = 3, ncol = 3)
dev.off()

#
#
#

# Effect of removing influential points on baseline model results
# re-run reg_base function with activated code removing outliers
base_results_b2_infl <- reg_base(b2_imp_sum) 

colnames(base_results_b2_infl)[5] <- 'p.value'

write.csv(base_results_b2_infl, 'Baseline_no_infl_BG_model_results_B2.csv')

base_results_b2_infl <- read.csv('Baseline_no_infl_BG_model_results_B2.csv')[,-1]

base_results_b2_infl <- base_results_b2_infl %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

#

# Compare the effect of excluding influential observations on regression estimates
# slopes (coefficients)
base_results_comp <- base_results_b2 %>%
  mutate(Data = 'All_data') %>%
  rbind(base_results_b2_infl %>%
          mutate(Data = 'Excl_outliers')) %>%
  mutate(Sign = ifelse(adj.p.value < 0.05, 'Yes', 'No')) %>%
  mutate(CI_width = CI_upr - CI_lwr) 


base_results_comp1 <- base_results_comp %>%
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
  base_results_comp1 %>%
  ggplot(aes(y = Estimate_All_data, x = Estimate_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.6) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 10) +
  xlab('Excl. outliers') +
  ylab('All data') +
  ggtitle('Regression logBG coefficient') +
  theme(legend.position = 'none')


# confidence intervals
CI_comp <- 
  base_results_comp1 %>%
  ggplot(aes(y = CI_width_All_data, x = CI_width_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.6) +
  coord_fixed() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw(base_size = 10) +
  xlab('Excl. outliers') +
  ylab('All data') +
  ggtitle('Regression logBG 95% CI width') +
  theme(legend.position = 'none')

# R2 (varaince explained)
r2_comp <- 
  base_results_comp1 %>%
  ggplot(aes(y = r2_All_data, x = r2_Excl_outliers)) + 
  geom_point(aes(colour = Significant), size = 0.8, alpha = 0.6) +
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
dev.new()
comp_plots

ggsave('Regression_comparison_del_diag.tiff' , comp_plots, dpi = 300, unit = 'mm',
       width = 260, height = 80)


#
#
#

###############################

#  Final model of the relationship between breath and background
# including covariates 

# join with clinical metadata
b2_imp_sum_annot <- b2_imp_sum %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() 

reg_fin <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc){
    subset <- df %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
    
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

fin_results_b2 <- reg_fin(b2_imp_sum_annot)

colnames(fin_results_b2)[5] <- 'p.value'

write.csv(fin_results_b2, 'Final_BG_model_results_B2.csv')

fin_results_b2 <- read.csv('Final_BG_model_results_B2.csv')[-1]

# number of VOCs with significant effect of background
fin_results_b2_bg <- fin_results_b2 %>% 
  filter(predictor == 'logBG')

fin_results_b2_bg_sign <- fin_results_b2 %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(adj.p.value < 0.05)

# number of VOCs with significant effect of diagnosis
asthma_pred_b2 <- fin_results_b2 %>% 
  filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

asthma_pred_b_sign <- asthma_pred_b2 %>% filter(p.value < 0.05)

hist(asthma_pred_b2$p.value, breaks = 12)

#
#
#

##############################

# MODEL SELECTION
# compare model fit after exclusion of influential observations
# after inclusion of covariates
mod_fit <- rbind(
  base_results_b2 %>% mutate(model = 'Baseline') %>%
    filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2, Estimate, p.value, adj.p.value),
  base_results_b2_infl %>% mutate(model = 'Baseline_no_infl') %>%
    filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2, Estimate, p.value, adj.p.value),
  fin_results_b2 %>% mutate(model = 'Final') %>%
    filter(predictor == 'logBG') %>% dplyr::select(comp, model, AIC, BIC, r2, Estimate, p.value) %>%
    mutate(adj.p.value = p.adjust(p.value, method = 'BH'))) 

mod_fit_L <- mod_fit %>% pivot_longer(cols = c(AIC, BIC, r2), names_to = 'parameter', values_to = 'value')

plot_pairwise <- function(Measure){
  anno_plot <- compare_means(value ~ model,
                             paired = TRUE,
                             method = 'wilcox.test',
                             data = mod_fit_L %>% filter(parameter == Measure),
                             p.adjust.method = 'BH') %>%
    mutate(y = c(570, 535, 500))  %>% 
    as.data.frame() %>%
    filter(p.adj < 0.05) %>%
    mutate(p.adj = '< 0.001')
  
  mod_fit_L %>%
    filter(parameter == Measure) %>%
  ggplot(aes(x = model, y = value, fill = model)) + geom_violin() +
  geom_boxplot(width = 0.5) +
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8),
        legend.position = 'bottom') +
  geom_signif(data = anno_plot, 
              aes(annotations = p.adj, xmin = group1, xmax = group2, y_position = y), manual = TRUE,
              inherit.aes = FALSE, tip_length = 0.01, textsize = 2, size = 0.2) +
    ylim(NA, 600) 
  }


aic <- plot_pairwise('AIC') + ggtitle('BIC') + theme(plot.title = element_text(hjust = 0.5),
                                                     legend.position = 'none') 
bic <- plot_pairwise('BIC') + ggtitle('AIC') + theme(plot.title = element_text(hjust = 0.5),
                                                     legend.text = element_text(size = 7), 
                                                     legend.key.size = unit(3, 'mm'))

comp_mod <- arrangeGrob(aic, bic, nrow = 2, heights = c(0.45, 0.55))
plot(comp_mod)

ggsave('BIC_AIC_violin.tiff', comp_mod, dpi = 300, unit = 'mm', width = 75, height = 120)

# summarise AIC and BIC value distributions for each model
mod_fit %>% group_by(model) %>% summarise(mean_AIC = mean(AIC),
                                          sd_AIC = sd(AIC),
                                          mean_BIC = mean(BIC),
                                          sd_BIC = sd(BIC),
                                          mean_r2 = mean(r2),
                                          median_AIC = median(AIC),
                                          median_BIC = median(BIC),
                                          median_r2 = median(r2),
                                          mean_b = mean(Estimate),
                                          med_b = median(Estimate))

mod_fit_w <- mod_fit %>% dplyr::select(comp, model, AIC, BIC) %>%
  pivot_wider(names_from = model, values_from = c(AIC, BIC)) %>%
  mutate(AIC_diff_BF = (AIC_Baseline - AIC_Final)/AIC_Baseline,
         BIC_diff_BF = (BIC_Baseline - BIC_Final)/BIC_Baseline,
         AIC_diff_BBni = (AIC_Baseline - AIC_Baseline_no_infl)/AIC_Baseline,
         BIC_diff_BBni = (BIC_Baseline - BIC_Baseline_no_infl)/BIC_Baseline,
         AIC_diff_BniF = (AIC_Baseline_no_infl - AIC_Final)/AIC_Baseline_no_infl,
         BIC_diff_BniF = (BIC_Baseline_no_infl - BIC_Final)/BIC_Baseline_no_infl)

hist(mod_fit_w$AIC_diff_BBni)

median(mod_fit_w$AIC_diff_BF)
median(mod_fit_w$BIC_diff_BF)

median(mod_fit_w$AIC_diff_BBni)
median(mod_fit_w$BIC_diff_BBni)

median(mod_fit_w$AIC_diff_BniF)
median(mod_fit_w$BIC_diff_BniF)

wilcox.test(mod_fit_w$AIC_Baseline, mod_fit_w$AIC_Baseline_no_infl, paired = TRUE)

#

# number of VOCs with significant effect of background
base_results_b2_sign <- base_results_b2 %>% filter(adj.p.value < 0.05)
base_results_b2_infl_sign <- base_results_b2_infl %>% filter(adj.p.value < 0.05)


venn <- venn.diagram(
  x = list(base_results_b2_sign$comp, 
           base_results_b2_infl_sign$comp,
           fin_results_b2_bg_sign$comp),
  category.names = c('Baseline', 'Baseline \n w/o outliers', 'Final'),
  filename = 'Venn_mixed_reg_sign.png',
  output = TRUE,
  cex = 1.2,
  fontface = 'bold',
  cat.cex = 1.2,
  cat.fontface = 'bold',
  cat.default.pos = "outer",
  fill = c('pink', 'grey', 'lightblue'),
  fontfamily = "sans",
  cat.fontfamily = 'sans',
  margin = 0.07,
  cat.dist = c(0.08, 0.12, 0.06),
  main = 'Number of VOCs with significant logBG effect',
  main.cex = 1.2,
  main.fontface = 'bold',
  main.fontfamily = 'sans',
  main.just = c(0.5, -0.5),
  resolution = 300,
  height = 120,
  width = 120,
  units = 'mm'
)

# plot of regression coefficients and r2
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

fin_res_plot <- fin_results_b2_bg %>%
  dplyr::select(comp, r2, Estimate, p.value) %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'),
         sign = ifelse(adj.p.value < 0.05, '< 0.05', '> 0.05')) %>%
  ggplot(aes(x = r2, y = Estimate, colour = sign)) + 
  geom_point(size = 0.8, alpha = 0.6) +
  xlim(0, NA) + ylim(0, NA) +
  theme_bw(base_size = 10) +
  scale_colour_manual(name = 'adj.p value',
                      values = c('> 0.05' = my_pal[[7]],
                                 '< 0.05' = 'black')) +
  ylab('Regression coefficient') +
  ggtitle('Results of the final mixed-effect model') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.key.size = unit(3, 'mm'))
  
ggsave('Final_model_r2VSbeta.tiff', units = 'mm', width = 85, height = 60)

quantile(fin_results_b2_bg_sign$Estimate)
quantile(fin_results_b2_bg_sign$r2)

# relationship between linear model outputs and VOC origin annotation
fin_results_b2_bg_class <- fin_results_b2_bg %>% 
  left_join(endo_exo %>% filter(wilcox_filter != 'Neither_dataset')) %>%
  mutate(bg_sign = ifelse(comp %in% fin_results_b2_bg_sign$comp, 'adj.p < 0.05', 'adj.p > 0.05'))


fin_class_sum <- fin_results_b2_bg_class %>%
  group_by(bg_sign, FC_filter) %>%
  summarise(n = n()/n_distinct(fin_results_b2_bg$comp)*100) 

fin_class_sum
  
fin_class_sum %>%
  ggplot(aes(x = bg_sign, y = n/n_distinct(fin_results_b2_bg$comp)*100, fill = FC_filter)) +
  geom_col() +
  theme_bw() +
  ylab('% VOCs') +
  xlab('Background effect')

#
#
#

################################

# MODEL VALIDATION
# Validation of regression models prediction in test dataset
# Baseline model validation
base_valid <- function(train, test){
  bind_rows(lapply(unique(train$comp), function(voc) {
  sub_train <- train %>%
    filter(comp == voc) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) #%>%
    #left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
    #mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
    #filter(obs %ni% c('infl'))
  
  sub_test <- test %>%
    filter(comp == voc) %>%
    drop_na() %>%
    mutate(logS = log(S), logBG = log(BG)) %>%
    distinct() 
  
  model <- lmer(logS ~ logBG + (1 | RAD_ID),
                data = sub_train)
  
  pred_test <- sub_test %>% 
    mutate(logSpred = predict(model, newdata = sub_test, allow.new.levels = TRUE),
           error = logSpred - logS)
  
  perf <- data.frame(mse = mean((pred_test$error)^2),
                     mae = mean(abs(pred_test$error)),
                     me = mean(pred_test$error),
                     mape = mean(abs(pred_test$error/pred_test$logSpred)*100),
                     comp = voc)
  }))
  }

# Fit baseline model to test dataset
base_results_b1 <- base_valid(b2_imp_sum, b1_imp_sum)

# Fit baseline model trained without influential observations to test dataset
# activate code excluding outliers in the function body above
base_results_b1_infl <- base_valid(b2_imp_sum, b1_imp_sum)

#

# Final model validation
fin_valid <- function(train, test){
  bind_rows(lapply(unique(train$comp), function(voc) {
    sub_train <- train %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
    
    sub_test <- test %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      distinct() 
    
    model_fin <- lmer(logS ~ logBG + Diagnosis + CoreVisit + (1 | RAD_ID),
                  data = sub_train)
    
    pred_test_fin <- sub_test %>% 
      mutate(logSpred = predict(model_fin, newdata = sub_test, allow.new.levels = TRUE),
             error = logSpred - logS)
    
    perf <- data.frame(mse = mean((pred_test_fin$error)^2),
                       mae = mean(abs(pred_test_fin$error)),
                       me = mean(pred_test_fin$error),
                       mape = mean(abs(pred_test_fin$error/pred_test_fin$logSpred)*100),
                       comp = voc)
  }))
}

# join with clinical metadata
b2_imp_sum_annot <- b2_imp_sum %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() 
b1_imp_sum_annot <- b1_imp_sum %>% left_join(meta %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() 

# Fit final model to test dataset
fin_results_b1 <- fin_valid(b2_imp_sum_annot, b1_imp_sum_annot)

#

# Fit null model (intercept only)
null_valid <- function(train, test) {
  bind_rows(lapply(unique(train$comp), function(voc) {
    sub_train <- train %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs) %>% # code to exclude influential points based on deletion diagnostics
      mutate(obs = ifelse(is.na(obs) == TRUE, NA, 'infl')) %>%
      filter(obs %ni% c('infl'))
    
    sub_test <- test %>%
      filter(comp == voc) %>%
      drop_na() %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      distinct() 
    
    model_null <- lm(logS ~ 1,
                data = sub_train)
    
    
    pred_test_null <- sub_test %>% 
      mutate(logSpred = predict(model_null, newdata = sub_test, allow.new.levels = TRUE),
             error = logSpred - logS)
    
    perf <- data.frame(mse = mean((pred_test_null$error)^2),
                       mae = mean(abs(pred_test_null$error)),
                       me = mean(pred_test_null$error),
                       mape = mean(abs(pred_test_null$error/pred_test_null$logSpred)*100),
                       comp = voc)
  }))
}

#
null_results_b1 <- null_valid(b2_imp_sum, b1_imp_sum)

# compare prediction error in test dataset across three models
# examine only VOCs where effect of background considered significant
valid_results <- base_results_b1 %>%
  pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
  mutate(model = 'Baseline') %>%
  rbind(base_results_b1_infl %>%
          pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
          mutate(model = 'Baseline_no_infl')) %>%
  rbind(fin_results_b1 %>%
          pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
          mutate(model = 'Final')) %>%
  rbind(null_results_b1 %>%
          pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
          mutate(model = 'Null')) %>%
  filter(error != 'mse') %>%
  filter(comp %in% fin_results_b2_bg_sign$comp)

valid_results %>% pivot_wider(names_from = 'error', values_from = 'value') %>%
  ggplot(aes(x = model, y = mae, group = model)) + geom_violin() + 
  geom_boxplot(outliers = FALSE, width = 0.1)  

valid_results %>% pivot_wider(names_from = 'error', values_from = 'value') %>%
  group_by(model) %>% 
  summarise(med_mae = median(mae),
            iqr_mae = IQR(mae),
            med_mape = median(mape),
            iqr_mape = IQR(mape),
            med_me = median(me),
            iqr_me = IQR(me))

valid_results %>% 
  filter(model %in% c('Final', 'Null')) %>%
  filter(error %ni% c('mse', 'me')) %>%
  ggplot(aes(x = model, y = value)) + 
  geom_violin() +
  geom_boxplot(width = 0.4, outliers = FALSE) +
  facet_wrap(~ error, scale = 'free') +
  theme_bw() 



View(compare_means(value ~ model,
              paired = TRUE,
              data = valid_results,
              group.by = c('error'),
              p.adjust.method = "BH"))

valid_results_w <- 
  valid_results %>% 
  dplyr::select(comp, error, value, model) %>%
  #filter(model %in% c('Final', 'Null')) %>%
  filter(error != 'mse') %>%
  pivot_wider(names_from = error, values_from = value) %>%
  pivot_wider(names_from = model, values_from = c(mae, mape, me)) %>%
  mutate(mae_diff_NF = mae_Null - mae_Final,
         mape_diff_NF = mape_Null - mape_Final)

median(valid_results_w$mae_diff_NF)/median(valid_results_w$mae_Null)
median(valid_results_w$mape_diff_NF)/median(valid_results_w$mape_Null)

#

valid_results1 <- fin_results_b1 %>%
  pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
  mutate(model = 'Final') %>%
  rbind(null_results_b1 %>%
          pivot_longer(cols =! comp, names_to = 'error', values_to = 'value') %>%
          mutate(model = 'Null'))

write.csv(valid_results1, 'Final_BG_model_errors_B1.csv')

valid_results1 <- read.csv('Final_BG_model_errors_B1.csv')[,-1]

#
#
#

#################################

# visualise effect of asthma based on confidence and effect size
hist(asthma_pred_b2$p.value, breaks = 30)
sd(asthma_pred_b2$p.value)

#
fin_results_b1 <- read.csv('Final_BG_model_errors_B1.csv')[,-1]
fin_results_b1a <- fin_results_b1 %>% filter(model == 'Final') %>%
  pivot_wider(names_from = error, values_from = value)

asthma_pred_b2b1 <- asthma_pred_b2 %>%
  left_join(fin_results_b1a) %>%
  mutate(lab = ifelse(p.value < 0.05 & abs(Estimate) > 0.4, comp, NA))

colvalquant <- base::seq(from = 0, to = 1, length.out = 5)
quants <- quantile(asthma_pred_b2b1$mape, colvalquant)
asthma_pred_b2b1$ptile_var <- ecdf(asthma_pred_b2b1$mape)(asthma_pred_b2b1$mape)

#

vol_plot1 <- asthma_pred_b2b1 %>%
  ggplot(aes(y = -log10(p.value), x = Estimate)) +
  geom_point(aes(fill = ptile_var), pch = 21, size = 1.3) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = 'dashed', colour = 'darkgrey', lwd = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed', colour = 'darkgrey', lwd = 0.3) +
  theme_bw(base_size = 8) +
  geom_text_repel(aes(label = str_sub(lab, end = 26)), size = 1.8,
                  point.padding = 0.1, nudge_y = 0.01, lwd = 0.3) +
  scale_fill_viridis_c(name = 'Validation dataset error (MAPE)',
                         labels = c(round(quants, 1)),
                       option = 'plasma') +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
        legend.key.height = unit(1, 'mm'),
        legend.position = 'bottom',
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0, 0.1), 'cm')) +
  ggtitle('Volcano plot showing effect of asthma diagnosis on breath VOC') +
  xlim(-0.7, 1) +
  xlab('Partial regression coefficient') +
  ylab('-Log10(p.value)') +
  annotate(geom = 'text', label = 'Asthma', x = -0.6, y = -0.5, size = 3)  +
  annotate(geom = 'text', label = 'Not Asthma', x = 0.75, y = -0.5, size = 3)


dev.new()
vol_plot1

ggsave('Volcano_plot_univariate_diagnosis.tiff', unit = 'mm', width = 140, height = 90)

################################

# Figure 2a
plot2a <- vol_plot1

#
#
#

##############################

# APPLYING BG CORRECTION
## use linear regression coefficients to correct the background signal
colnames(fin_results_b2_bg)[6:7] <- c('CI_lwr', 'CI_upr')

b1_imp_sum_corr <- b1_imp_sum %>%
  left_join(fin_results_b2_bg %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
              dplyr::select(comp, Estimate, CI_lwr, CI_upr, adj.p.value)) %>%
  mutate(logS_adj = log(S) - Estimate*log(BG))

#

b2_imp_sum_corr <- b2_imp_sum %>%
  left_join(fin_results_b2_bg %>%
              mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
              dplyr::select(comp, Estimate, CI_lwr, CI_upr, adj.p.value)) %>%
  mutate(logS_adj = log(S) - Estimate*log(BG))

# exclude VOCs where 95% CI for logBG coefficient includes 1 
excl <- fin_results_b2_bg %>%
  filter(CI_upr > 1) 

b1_imp_sum_corr <- b1_imp_sum_corr %>%
  filter(comp %ni% excl$comp)

b2_imp_sum_corr <- b2_imp_sum_corr %>%
  filter(comp %ni% excl$comp)

# BG-CORRECTED OUTPUT
write.csv(b1_imp_sum_corr, 'RADicA_B1_BG_adjusted.csv')
write.csv(b2_imp_sum_corr, 'RADicA_B2_BG_adjusted.csv')

#
#
#

###############################

# FIGURE 1b
plot_1b <- function(df, voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      left_join(infl_obs %>% dplyr::select(Sample, obs, comp)) %>%
      filter(is.na(obs) == TRUE)
    
    plot <- ggplot(data = subset, 
                   aes(x = logBG, y = logS)) +
      geom_point(alpha = 0.5, size = 0.8) + 
      ggtitle(voc) +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 8),
            legend.position = 'none',
            axis.title = element_blank(),
            panel.grid = element_blank()) +
      geom_smooth(method = 'lm', se = FALSE, colour = my_pal[[7]]) +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
    
    print(plot)
  }

ace <- plot_1b(b2_imp_sum, 'Acetophenone') +
  annotate(geom = 'text', label = 'b = 0.62', x = 10.5, y = 14, size = 2.2, colour = 'purple') +
  annotate(geom = 'text', label = 'r2 = 0.58', x = 10.5, y = 13.6, size = 2.2, colour = 'purple') +
  annotate(geom = 'text', label = 'p < 0.01', x = 10.5, y = 13.2, size = 2.2, colour = 'purple')


cyc <- plot_1b(b2_imp_sum, 'Cyclopentane') +
  annotate(geom = 'text', label = 'b = -0.09', x = 14.05, y = 14.9, size = 2.2, colour = 'purple') +
  annotate(geom = 'text', label = 'r2 = 0.06', x = 14.05, y = 14.65, size = 2.2, colour = 'purple') +
  annotate(geom = 'text', label = 'p = 0.26', x = 14.05, y = 14.4, size = 2.2, colour = 'purple')

plot1b <- plot_grid(ace, cyc, ncol = 1) +
  draw_label(label = 'Background (Log peak area)', size = 8, x = 0.5, y = 0, vjust = 1) +
  draw_label(label = 'Breath (Log peak area)', x = 0, y = 0.5, angle = 90, size = 8, vjust = 0) +
  theme(plot.margin = unit(c(0, 0, 0.5, 0.5), "cm"))
  
plot(plot1b)

# merge top figure panels
plot1ab <- plot_grid(plot1a, plot1b, nrow = 1, rel_widths = c(0.73, 0.27), labels = c('a)', 'b)'), label_size = 12)
dev.new()
plot(plot1ab)



#
#
#

# Figure 1c
plot1c <- fin_results_b2_bg %>%
  dplyr::select(comp, r2, Estimate, p.value) %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'),
         sign = ifelse(adj.p.value < 0.05, '< 0.05', '> 0.05')) %>%
  ggplot(aes(x = r2, y = Estimate, colour = sign)) + 
  geom_point(size = 0.9, shape = 1, alpha = 0.9) +
  xlim(0, NA) + ylim(0, NA) +
  theme_bw(base_size = 8) +
  scale_colour_manual(name = 'Adj.p value',
                      values = c('> 0.05' = my_pal[[7]],
                                 '< 0.05' = 'black')) +
  ylab('Regression coefficient') +
  ggtitle('Estimates of background effect \n in mixed-effect models') +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 8),
        legend.position = 'bottom',
        plot.margin = unit(c(0.6, 0.5, 1, 0.5), 'cm'),
        panel.grid = element_blank(),
        legend.direction = 'vertical',
        legend.title = element_text(size = 6))

plot1c

