## Comparing effect of background correction on univariate analysis

# author: Aggie Turlo
# project: RADicA
# date: 10/03/2025

#####################
# load normalised dataset with summarised breath samples
b1_imp_sum <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]
b2_imp_sum <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]

b1_imp_sum$Analysis_date <- as.Date(b1_imp_sum$Analysis_date)
b2_imp_sum$Analysis_date <- as.Date(b2_imp_sum$Analysis_date)

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

# load mixed-effect model results
bg <- read.csv('Final_BG_model_results_B2.csv')
no_bg <- read.csv('Final_BG_model_results_B2_noBG.csv')

asthma_pred_bg <- bg %>% 
  filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

asthma_pred_bg_sign <- asthma_pred_bg %>% filter(p.value < 0.05)

bg_pred <- bg %>% 
  filter(predictor == 'logBG') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

bg_pred_sign <- bg_pred %>% filter(adj.p.value < 0.05)

#

asthma_pred_no_bg <- no_bg %>% 
  filter(predictor == 'DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

asthma_pred_no_bg_sign <- asthma_pred_no_bg %>% filter(p.value < 0.05)

#

intersect(asthma_pred_bg_sign$comp, asthma_pred_no_bg_sign$comp)

setdiff(asthma_pred_bg_sign$comp, asthma_pred_no_bg_sign$comp)
setdiff(asthma_pred_no_bg_sign$comp, asthma_pred_bg_sign$comp)

#
asthma_pred <- rbind(asthma_pred_bg %>% mutate(BGcorrection = 'yes'),
                     asthma_pred_no_bg %>% mutate(BGcorrection = 'no'))

#
vol_plot1 <- asthma_pred %>%
  #mutate(BG_sign = ifelse(comp %in% bg_pred_sign$comp, 'yes', 'no')) %>%
  filter(comp %in% bg_pred_sign$comp) %>%
  ggplot(aes(y = -log10(p.value), x = Estimate)) +
  geom_point(aes(fill = BGcorrection), pch = 21, size = 1) +
  geom_vline(xintercept = c(-0.4, 0.4), linetype = 'dashed', colour = 'tomato') +
  geom_hline(yintercept = 1.3, linetype = 'dashed', colour = 'tomato') +
  theme_bw(base_size = 8) +
  geom_path(aes(group = comp), 
            arrow = arrow(length = unit(0.2, "cm"), ends = 'first'),
            colour = 'grey5', size = 0.1) +
  #scale_colour_manual(values = c('yes' = 'gray20', 'no' = 'grey')) +
  #geom_text_repel(aes(label = str_sub(lab, end = 26)), size = 1.8) +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.key.size = unit(3, 'mm')) +
  ggtitle('Volcano plot showing effect of asthma diagnosis on breath VOC levels') +
  xlim(-0.7, 1) +
  xlab('partial regression coefficient') +
  annotate(geom = 'text', label = 'Asthma', x = -0.55, y = -0.5, size = 3.5)  +
  annotate(geom = 'text', label = 'Not Asthma', x = 0.75, y = -0.5, size = 3.5) +
  scale_fill_manual(values = c('yes' = 'grey18', 'no' = 'grey75')) +
  annotate(geom = 'text', label = 'p = 0.05', x = -0.65, y = 1.45, colour = 'tomato', size = 3.5)

dev.new()
vol_plot1

ggsave('VOlcano_plot_BG_effect.tiff', unit = 'mm', dpi = 300, width = 140, height = 83)

#
asthma_pred_w <- asthma_pred %>% 
  filter(comp %in% bg_pred_sign$comp) %>% 
  dplyr::select(comp, model, Estimate, p.value, AIC, BIC) %>%
  pivot_wider(names_from = model, values_from = c(Estimate, p.value, AIC, BIC)) %>%
  mutate(Estimate_diff = Estimate_noBG - Estimate_BG,
         p.value_diff = p.value_noBG - p.value_BG,
         Estimate_diff_p = (abs(Estimate_noBG - Estimate_BG)/abs(Estimate_noBG))*100,
         AIC_diff = AIC_noBG - AIC_BG,
         BIC_diff = BIC_noBG - BIC_BG)

hist(asthma_pred_w$Estimate_diff)
hist(abs(asthma_pred_w$Estimate_diff))
hist(abs(asthma_pred_w$p.value_diff))              
hist(asthma_pred_w$Estimate_diff_p)

quantile(abs(asthma_pred_w$Estimate_diff))
quantile(asthma_pred_w$Estimate_diff)
quantile(abs(asthma_pred_w$p.value_diff))
quantile(asthma_pred_w$Estimate_diff_p)

wilcox.test(asthma_pred_w$AIC_BG, asthma_pred_w$AIC_noBG, paired = TRUE)

hist(asthma_pred_w$AIC_diff)
hist(asthma_pred_w$BIC_diff)
