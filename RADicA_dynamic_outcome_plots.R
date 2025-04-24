library(dplyr)
library(ggplot2)
library(forcats)

dyn <- read.csv('Dynamic_outcomes_model_perf.csv')

dyn$test_set <- factor(test_set, levels = c('train','CV','test'))

mae_plot <- dyn %>% filter(measure == 'MAE') %>%
  mutate(outcome  = ifelse(outcome == 'FEVPPred', 'FEVPPred (%)',
                       ifelse(outcome == 'FVC', 'FVC (l)', 
                              ifelse(outcome == 'FeNO', 'FeNO (ppb)', 'FEV/FVC')))) %>%
  ggplot(aes(x = fct_inorder(test_set), y = mean, fill = fct_inorder(test_set))) +
  geom_bar(stat = 'identity', position = 'dodge', colour = 'black') +
  facet_wrap(~outcome, scale = 'free', nrow = 1) +
  xlab('') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
                width =.2) +
  theme_bw() +
  ylab('Mean absolute error') +
  scale_fill_brewer(palette = 'Greys') +
  theme(legend.position = 'none')

cor_plot <- dyn %>% filter(measure == 'r') %>%
  mutate(outcome  = ifelse(outcome == 'FEVPPred', 'FEVPPred (%)',
                       ifelse(outcome == 'FVC', 'FVC (l)', 
                              ifelse(outcome == 'FeNO', 'FeNO (ppb)', 'FEV/FVC')))) %>%
  ggplot(aes(x = fct_inorder(test_set), y = mean, fill = fct_inorder(test_set))) +
  geom_bar(stat = 'identity', position = 'dodge', colour = 'black') +
  facet_wrap(~ outcome, scale = 'free', nrow = 1) +
  xlab('') +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd),
                width =.2) +
  theme_bw() +
  ylab('Correlation') +
  scale_fill_brewer(palette = 'Greys') +
  theme(legend.position = 'none')

dyn_plots <- arrangeGrob(mae_plot, cor_plot, nrow = 2)
plot(dyn_plots)

ggsave('Dyn_outcomes_results.tiff', dyn_plots, unit = 'mm', dpi = 300, width = 160, height = 100)
