## Background correction (filtering VOCs)

# author: Aggie Turlo
# project: RADicA
# date: 28/01/2025

#####################

library(dplyr)
library(tidyr)
library(clusrank)
library(ggpubr)
library(ggplot2)
library(shades)
library(grafify)
library(cowplot)

## load data
# load normalised dataset with summarised breath samples
b1_imp_sum <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised.csv')[,-1]
b2_imp_sum <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised.csv')[,-1]

b1_imp_sum$Analysis_date <- as.Date(b1_imp_sum$Analysis_date)
b2_imp_sum$Analysis_date <- as.Date(b2_imp_sum$Analysis_date)

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

# custom functions
'%ni%' <- Negate('%in%')

#
#
#

# remove samples from the same patients from B2
rep <- intersect(b1_imp_sum$RAD_ID, b2_imp_sum$RAD_ID)

b2_imp_sum <- b2_imp_sum %>% filter(RAD_ID %ni% rep)

#

## format study datasets
# identify VOCs without matching background information (dataset 2)
no_bg_voc <- b2_imp_sum %>% group_by(comp) %>% summarise(na_bg = sum(is.na(BG) == TRUE)) %>%
  filter(na_bg == n_distinct(b2_imp_sum$Sample)) %>% pull(comp)

# remove from both datasets
b1_imp_sum <- b1_imp_sum %>% filter(comp %ni% no_bg_voc)
b2_imp_sum <- b2_imp_sum %>% filter(comp %ni% no_bg_voc)

# identify breath samples without matching background sample (dataset 1)
no_bg_s <- b1_imp_sum[which(is.na(b1_imp_sum$BG), arr.ind = TRUE),] %>% 
  group_by(Sample) %>% summarise(n = n()) %>% pull(Sample)

# remove from dataset 1
b1_imp_sum <- b1_imp_sum %>% filter(Sample %ni% no_bg_s)

#
#
#

################################

## Filtering VOCs with similar breath and background distributions

# plot breath vs background log abundance distributions
# specify dataset for analysis
data <- b1_imp_sum
dat <- 'B1'

#

bg_s_plot <- function(df) {
  lapply(unique(df$comp), function(voc) {
    subset <- df %>% filter(comp == voc) %>%
      mutate(logS = log(S), logBG = log(BG)) %>%
      dplyr::select(!c(S, BG)) %>%
      pivot_longer(cols = c(logS, logBG), names_to = 'Sample_type', values_to = 'logPeakArea')
    
    plot <- ggplot(data = subset, 
                   aes(x = logPeakArea, fill = Sample_type)) +
      geom_histogram(alpha = 0.6, position = 'identity', bins = 12) + 
      ggtitle(voc) +
      theme_bw(base_size = 8) +
      theme(plot.title = element_text(size = 8),
            legend.position = 'none')
  })
}

#

plots <- bg_s_plot(data)

pdf(paste('BGvsS_histograms_', dat, '.pdf', sep = ''))
marrangeGrob(plots, nrow = 4, ncol = 4)
dev.off()


#


# clustered Wilcoxon rank-sum test
# function applying clustered Wilcoxon test to each voc
wilcox_bg_fun <- function(df) {
  bind_rows(lapply(unique(df$comp), function(voc) {
    
    subset <- df %>% 
      filter(comp == voc) 
    
    subset1 <- subset %>%
      pivot_longer(cols = c(BG, S), names_to = 'Sample_type', values_to = 'peakArea') %>%
      drop_na() 
    
    test <- clusWilcox.test(peakArea ~  Sample_type + cluster(Sample) + stratum(RAD_ID),
                            data = subset1, paired = FALSE, method = 'rgl', alternative = 'two')
    
    
    df <- data.frame(comp = voc, 
                     p.value = test$p.value,
                     logFC = mean(log(subset$S)) - mean(log(subset$BG)))
  }))
}

# apply to each dataset
wilcox_bg_b1 <- wilcox_bg_fun(b1_imp_sum) %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

wilcox_bg_b2 <- wilcox_bg_fun(b2_imp_sum) %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

# filter according to statistical significance
wilcox_b1_f <- wilcox_bg_b1 %>% filter(adj.p.value > 0.05)
wilcox_b2_f <- wilcox_bg_b2 %>% filter(adj.p.value > 0.05)

# identify VOCs with no difference between distributions in both datasets
wilcox_agreement <- intersect(wilcox_b1_f$comp, wilcox_b2_f$comp)

# filter according to log fold change
logfc_b1_f <- wilcox_bg_b1 %>% filter(logFC > 0.4)
logfc_b2_f <- wilcox_bg_b2 %>% filter(logFC > 0.4)

logfc_agreement <- intersect(logfc_b1_f$comp, logfc_b2_f$comp)

#

# annotate VOCs depending on results of filtering based on Wilcoxon test and logFC 
wilcox_b1 <- setdiff(wilcox_b1_f$comp, wilcox_b2_f$comp)
wilcox_b2 <- setdiff(wilcox_b2_f$comp, wilcox_b1_f$comp)

logfc_b1 <- setdiff(logfc_b1_f$comp, logfc_b2_f$comp)
logfc_b2 <- setdiff(logfc_b2_f$comp, logfc_b1_f$comp)

all_vocs <- unique(b1_imp_sum$comp)

annot_voc <- #all_vocs[all_vocs %ni% # remove VOCs if Wilcoxon null hypothesis accepted in both datasets
                        #wilcox_agreement] %>%
  data.frame(comp = all_vocs) %>%
  mutate(wilcox_filter = ifelse(comp %in% wilcox_agreement, 'Neither_dataset',
                                ifelse(comp %ni% c(wilcox_b1, wilcox_b2), 'Both_datasets', 
                                ifelse(comp %in% wilcox_b1, 'Dataset_2', 'Dataset_1')))) %>%
  mutate(FC_filter = ifelse(comp %in% logfc_agreement, 'Endo_both_datasets',
                            ifelse(comp %in% logfc_b1, 'Endo_dataset_1',
                                   ifelse(comp %in% logfc_b2, 'Endo_dataset_2', 'Exo')))) %>%
  left_join(wilcox_bg_b1) %>%
  mutate(p.value_Dataset1 = p.value,
         logFC_Dataset1 = logFC,
         adj.p.value_Dataset1 = adj.p.value) %>%
  dplyr::select(!c(p.value, adj.p.value, logFC)) %>%
  left_join(wilcox_bg_b2) %>%
  mutate(p.value_Dataset2 = p.value,
         logFC_Dataset2 = logFC,
         adj.p.value_Dataset2 = adj.p.value)  %>%
  dplyr::select(!c(p.value, adj.p.value, logFC))

write.csv(annot_voc, 'Endo_Exo_filters.csv')

# remove VOCs where null hypothesis accepted in both datasets
b1_imp_sum <- b1_imp_sum %>% filter(comp %ni% wilcox_agreement)
b2_imp_sum <- b2_imp_sum %>% filter(comp %ni% wilcox_agreement) 

# save the outputs
write.csv(b1_imp_sum, 'RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')
write.csv(b2_imp_sum, 'RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')


#
#
#

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

################################

# FIGURE 1a
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

plot_1a <- function(df, voc) {
  subset <- df %>% filter(comp == voc) %>%
  mutate(logS = log(S), logBG = log(BG)) %>%
  dplyr::select(!c(S, BG)) %>%
  pivot_longer(cols = c(logS, logBG), names_to = 'Sample_type', values_to = 'logPeakArea') %>%
    mutate(Sample_type = ifelse(Sample_type == 'logS', 'Breath', 'Background'))
  
  plot <- ggplot(data = subset, 
               aes(x = logPeakArea, fill = Sample_type, col = Sample_type)) +
    geom_histogram(alpha = 0.3, position = 'identity', bins = 10, colour = NA) +
    geom_step(stat = "bin", bins = 10, direction = "mid", linewidth = 1) +
    theme_bw(base_size = 8) +
    theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 8),
        panel.grid = element_blank(),
        legend.key.width = unit(3, 'mm'),
        legend.title = element_text(size = 6)) +
    scale_fill_manual(values = c(Background = my_pal[[4]],
                                 Breath = my_pal[[3]]),
                      name = 'Sample type') +
    scale_colour_manual(values = c(Background = my_pal[[4]],
                                 Breath = my_pal[[3]]),
                      name = 'Sample type') +
    geom_step(stat = "bin", bins = 10, direction = "mid", size = 0.1)

}

cd1 <- plot_1a(b1_imp_sum, 'D_Carvone')
cd2 <- plot_1a(b2_imp_sum, 'D_Carvone') + ggtitle('D-Carvone')
cd3 <- plot_1a(b1_imp_sum, 'D_Carvone')

ba1 <- plot_1a(b1_imp_sum, 'Benzaldehyde')
ba2 <- plot_1a(b2_imp_sum, 'Benzaldehyde') + ggtitle('Benzaldehyde')

fu1 <- plot_1a(b1_imp_sum, 'Furan')
fu2 <- plot_1a(b2_imp_sum, 'Furan') + ggtitle('Furan')

# legend
leg <- plot_1a(b1_imp_sum, 'D_Carvone') 
leg1 <- get_legend(leg)
dev.new()
plot(leg1)

#
plot <- plot_grid(cd2, fu2, ba2, cd1, fu1, ba1, align = 'vh', vjust = 1, scale = 1) +
  draw_label(label = 'Log peak area', size = 8, x = 0.5, y = 0, vjust = 1) +
  draw_label(label = 'Count', x = 0, y = 0.5, angle = 90, size = 8, vjust = 0) +
  draw_label(label = 'Training dataset', x = 0, y = 0.75, angle = 90, size = 8, vjust = -1.5) +
  draw_label(label = 'Validation dataset', x = 0, y = 0.25, angle = 90, size = 8, vjust = -1.5) +
  draw_label(label = 'Contaminants \n (p > 0.05)', x = 0.15, y = 1, size = 8, vjust = -0.5, fontface = 'bold') +
  draw_label(label = 'Breath-enriched VOCs \n (p < 0.05, logFC > 0.4)', x = 0.5, y = 1, size = 8, vjust = -0.5, fontface = 'bold') +
  draw_label(label = 'Ambiguous VOCs \n (p < 0.05, logFC < 0.4)', x = 0.85, y = 1, size = 8, vjust = -0.5, fontface = 'bold') +
  theme(plot.margin = unit(c(0.9, 0.1, 0.5, 1), "cm")) 

plot1a <- plot_grid(plot, leg1, rel_widths = c(0.85, 0.15))

# FINAL FIGURE PANEL
plot(plot1a)

