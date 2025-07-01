# RADicA backgroun filtering

# comparison of measures of spread across VOC classes

# reguires running RADicA_background_filtering.R scrip first

# compare the peak area variance in difference VOC classes
# calculate relative measures of data spread for each VOC
spread_fc <- function(Data) {
  b_class_sum <- Data %>%
    left_join(annot_voc) %>%
    filter(wilcox_filter != 'Neither_dataset') %>%
    group_by(comp) %>%
    summarise(rIQR_BG = IQR(BG)/median(BG),
              rIQR_S = IQR(S)/median(S),
              rMAD_BG = mad(BG)/median(BG),
              rMAD_S = mad(S)/median(S)) %>%
    left_join(annot_voc) %>%
    pivot_longer(cols = c(rIQR_BG, rIQR_S, rMAD_BG, rMAD_S),
                 names_to = 'measure', values_to = 'value') %>%
    mutate(sample = str_sub(measure, start = 6),
           measure = str_sub(measure, end = 4),
           sample = ifelse(sample == 'BG', 'BG (Background)', 'S (Breath)'))
}

b1_class_sum <- spread_fc(b1_imp_sum)
b2_class_sum <- spread_fc(b2_imp_sum)

View(b1_class_sum %>% group_by(FC_filter, measure, sample) %>%
       summarise(median_iqr  = median(value)))

View(b2_class_sum %>% group_by(FC_filter, measure, sample) %>%
       summarise(median_iqr  = median(value)))


# plot comparison in spread measures between VOC classses
# add results of pairwise wilcox tests
plot_pairwise <- function(Data, Measure){
  anno_plot <- compare_means(value ~ FC_filter,
                             group.by = c('sample'),
                             method = 'wilcox.test',
                             data = Data %>% filter(measure == Measure)) %>%
    as.data.frame() %>%
    filter(p.adj < 0.05) %>%
    mutate(y = seq(from = 3.25, to = (3.25+(0.25*(nrow(.)-1))), by = 0.25)) # parameters for rIQR
  #mutate(y = seq(from = 1.5, to = (1.5+(0.2*(nrow(.)-1))), by = 0.2)) # parameters for rMAD
  
  fc_plot <- Data %>%
    filter(measure == Measure) %>%
    ggplot(aes(x = FC_filter, y = value, fill = FC_filter)) +
    geom_boxplot(outliers = FALSE, 
                 width = 0.6,
                 position = position_dodge(0.75)) +
    theme_bw(base_size = 8) +
    facet_wrap(~ sample) +
    geom_signif(data = anno_plot, 
                aes(annotations = p.adj, xmin = group1, xmax = group2, y_position = y), manual = TRUE,
                inherit.aes = FALSE, tip_length = 0.01, textsize = 2, size = 0.2) +
    theme(axis.text.x = element_blank()) +
    ylab(Measure) +
    #ylim(NA, 2.25) # rMAD
    ylim(NA, 4) # rIQR
  
  fc_plot
}

# change plot titles accordingly
b1_iqr <- plot_pairwise(b1_class_sum, 'rIQR') + theme(legend.position = 'none', 
                                                      plot.title = element_text(hjust = 0.5, size = 8)) + ggtitle('Dataset 1')

b2_iqr <- plot_pairwise(b2_class_sum, 'rIQR') + theme(plot.title = element_text(hjust = 0.5, size = 8),
                                                      legend.text = element_text(size = 6), 
                                                      legend.key.size = unit(3, 'mm')) + ggtitle('Dataset 2') 

#

b1_rmad <- plot_pairwise(b1_class_sum, 'rMAD') + theme(legend.position = 'none', 
                                                       plot.title = element_text(hjust = 0.5, size = 8)) + ggtitle('Dataset 1')

b2_rmad <- plot_pairwise(b2_class_sum, 'rMAD') + theme(plot.title = element_text(hjust = 0.5, size = 8),
                                                       legend.text = element_text(size = 6), 
                                                       legend.key.size = unit(3, 'mm')) + ggtitle('Dataset 2')

iqr_plots <- arrangeGrob(b1_iqr, b2_iqr, ncol = 2, widths = c(0.4, 0.6))
rmad_plots <- arrangeGrob(b1_rmad, b2_rmad, ncol = 2, widths = c(0.4, 0.6))

fc_plots <- arrangeGrob(iqr_plots, rmad_plots, nrow = 2)
plot(fc_plots)

ggsave('FCfilter_vs_spread_boxplots.tiff', fc_plots, unit = 'mm', dpi = 300, width = 155, height = 85)

#
#
#
