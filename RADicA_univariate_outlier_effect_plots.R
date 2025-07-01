# RADicA background correction

# effect of outlier removal on model estimates

# requires running RADicA_background_correction_clean.R script first

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

# Validation of regression models prediction in test dataset
# Baseline model validation
base_valid <- function(train, test){
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
    
    model <- lmer(logS ~ logBG + (1 | ID),
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

