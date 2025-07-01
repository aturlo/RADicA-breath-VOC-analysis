## Exploratory multivariate analysis of corrected breath data

# author: Aggie Turlo
# project: RADicA
# date: 01/07/2025

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(mixOmics)
library(gridExtra)
library(forcats)
library(ggpattern)
library(scales)
library(cowplot)
library(grid)
library(pcaPP)
library(ggrepel)

## load data
# load the file with corrected VOC peak areas
b1_corr <- read.csv('RADicA_B1_BG_adjusted.csv')[,-1] 
b2_corr <- read.csv('RADicA_B2_BG_adjusted.csv')[,-1]

# load metadata
meta <- read.csv('RADicA_VOC_metadata.csv')

# load VOC origin annotations
endo_exo <- read.csv('Endo_Exo_filters.csv')[,-1] %>%
  left_join(b1_corr %>% dplyr::select(comp, adj.p.value) %>% distinct()) %>%
  mutate(adj = ifelse(adj.p.value < 0.05, 'yes', 'no'))

# load the list of univaroate outliers
infl_obs <- read.csv('Deletion_diagnostics_formatted_B2.csv')

# custom functions
'%ni%' <- Negate('%in%')

#
#
#

# remove samples from the same patients from B2 (information leakage)
rep <- intersect(b1_corr$ID, b2_corr$ID)

b2_corr <- b2_corr %>% filter(ID %ni% rep)


## Principal Component Analysis of corrected and uncorrected breath data
# transform data to wide format
pivot_wide <- function(data, values){
  data_w <- data %>%
    dplyr::select(c(Sample, comp, values, CoreVisit, ID)) %>%
    pivot_wider(names_from = comp, values_from = values) %>% 
    left_join(meta %>% dplyr::select(ID, Diagnosis), relationship = 'many-to-many') %>%
    distinct() %>%
    relocate(Diagnosis) %>% 
    as.data.frame()
}

b1_corr_w <- pivot_wide(data = b1_corr, values = 'logS_adj')
rownames(b1_corr_w) <- b1_corr_w$Sample

b1_w <- pivot_wide(data = b1_corr, values = 'S')
rownames(b1_w) <- b1_w$Sample

#

b2_corr_w <- pivot_wide(data = b2_corr, values = 'logS_adj')
rownames(b2_corr_w) <- b2_corr_w$Sample

b2_w <- pivot_wide(data = b2_corr, values = 'S')
rownames(b2_w) <- b2_w$Sample

# log transform uncorrected peak area data to match the corrected data format
b1_w <- cbind(b1_w[,1:4], log(b1_w[,5:ncol(b1_w)]))
b2_w <- cbind(b2_w[,1:4], log(b2_w[,5:ncol(b2_w)]))

# PCA 
pca_b1 <- pca(b1_w[,-c(1:4)], scale = TRUE, center = TRUE)
pca_b1_corr <- pca(b1_corr_w[,-c(1:4)], scale = TRUE, center = TRUE)

pca_b2 <- pca(b2_w[,-c(1:4)], scale = TRUE, center = TRUE)
pca_b2_corr <- pca(b2_corr_w[,-c(1:4)], scale = TRUE, center = TRUE)

# scores plots
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

plot_pca <- function(pca_b, data, Title){
  p <- plotIndiv(pca_b,
            #ind.names = data$ID,
            group = data$Diagnosis,
            legend = TRUE, 
            comp = c(1,2),
            pch = as.factor(data$CoreVisit),
            size.title = 8,
            cex = 0.9,
            title = Title,
            legend.title.pch = 'Collection visit',
            size.xlabel = 8,
            size.ylabel = 8,
            size.legend.title = 9,
            point.lwd = 0.8
            )
  
  p$graph + scale_colour_manual(name = 'Diagnosis',
                                values = c('Asthma' = my_pal[[2]],
                                           'Not Asthma' = my_pal[[1]]))  +
    theme(text = element_text(size = 6),
          plot.margin = unit(c(0, 0.1, 0.1, 0.1), 'cm'),
          legend.key.width = unit(2, 'mm'),
          panel.grid = element_blank(),
          strip.background = element_blank())
}

b1_plot <- plot_pca(pca_b1, b1_w, '') 
  
b1_corr_plot <- plot_pca(pca_b1_corr, b1_corr_w, '')

b2_plot <- plot_pca(pca_b2, b2_w, 'Before correction')

b2_corr_plot <- plot_pca(pca_b2_corr, b2_corr_w, 'After correction')

##################################

# Figure 1d
legd <- get_legend(b1_plot)
plot(legd)

#

plotd <- plot_grid(b2_plot, b2_corr_plot, b1_plot, b1_corr_plot,
                   rel_heights = c(0.52, 0.48)) 
dev.new()
plot(plotd)

plot1d <- plot_grid(plotd, legd, rel_widths = c(0.82, 0.18)) +
  draw_label(label = 'PCA score plots of breath VOC data \n before and after background correction', 
             size = 8, x = 0.5, y = 0.97, vjust = -0.5, fontface = 'bold') +
  theme(plot.margin = unit(c(0.7, 0, 0, 0.5), "cm")) +
  draw_label(label = 'Training dataset', x = 0, y = 0.75, angle = 90, size = 8, vjust = -0.5) +
  draw_label(label = 'Validation dataset', x = 0, y = 0.25, angle = 90, size = 8, vjust = -0.5)

plot(plot1d)

# merge bottom figure panels
plot1cd <- plot_grid(plot1c, plot1d, rel_widths = c(0.35, 0.65), labels = c('c)', 'd)'), label_size = 12)

dev.new()
plot(plot1cd)

# Assemble the final figure - FIGURE 1
fig1 <- plot_grid(plot1ab, plot1cd, ncol = 1, rel_heights = c(0.47, 0.53))
dev.new()
plot(fig1)

pdf('Figure1.pdf', width = 6.18, height = 6.3)
plot(fig1)
dev.off()

ggsave('fig1.tiff', unit = 'mm', dpi = 300, width = 157, height = 160)

#
#
#

#################################

## Multivariate outlier detection in BG corrected datasets
# detect multivariate outliers using robust PCA
# specify dataset for analysis
data <- b2_corr_w
dat <- 'B2'

#
assay <- data[,-c(1:4)]
Rpc <- PCAgrid(assay, scale = 'sd', center = 'mean')
sdod <- PCdiagplot(assay, Rpc, plotbw = FALSE, crit = 0.999)

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
  geom_point(size = 0.8) +
  geom_hline(yintercept = sdod$critSD[1,1], linetype = 'dashed') +
  geom_vline(xintercept = sdod$critOD[1,1], linetype = 'dashed') +
  theme_bw(base_size = 8) +
  ylab('score distance') + xlab('orthogonal distance') +
  theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  scale_colour_manual(values = c('no' = 'black',
                                 'yes' = 'red')) + 
  ggtitle('Diagnostic plot to the PCA on breath sample VOCs \n (Dataset 2)') # change title accordingly

ggsave(paste('rPCA_diagnostic_plot_', dat, '.tiff', sep =''), dpi = 300, unit = 'mm', width = 80, height = 60)

#

distances %>%
  filter(k == 'V1') %>%
  ggplot(aes(x = OD)) + geom_histogram(bins = 20) +
  geom_vline(xintercept = sdod$critOD[2,1], linetype = 'dashed', colour = 'red') +
  theme_bw(base_size = 8) +
  ggtitle('Diagnostic plot to the PCA on breath sample VOCs \n (Dataset 1)') +
  theme(plot.title =  element_text(hjust = 0.5, size = 8))

ggsave(paste('rPCA_OD_histogram_', dat, '.tiff', sep =''), dpi = 300, unit = 'mm', width = 80, height = 60)


# examine profile of outlying observations
data[,1:4] %>% filter(Sample %in% outliers)

# match the multivariate outliers with univariate outliers (dataset 2)
infl_obs_sum <- infl_obs %>% group_by(Sample) %>% summarise(n = n())

infl_obs_sum %>% filter(Sample %in% outliers)

# remove the outlying observations from the dataset
b1_corr_out <- b1_corr_w %>% filter(Sample %ni% outliers)
b2_corr_out <- b2_corr_w %>% filter(Sample %ni% outliers)

# save VOC datasets after multivariate outlier removal

write.csv(b1_corr_out, 'RADicA_BG_adjusted_B1_outl_removed.csv')
write.csv(b2_corr_out, 'RADicA_BG_adjusted_B2_outl_removed.csv')


#
#
#

# PCA after outlier removal
pca_b1_out <- pca(b1_corr_out[,-c(1:4)], scale = TRUE, center = TRUE)
pca_b2_out <- pca(b2_corr_out[,-c(1:4)], scale = TRUE, center = TRUE)

#

b1_out_plot <- plot_pca(pca_b1_out, b1_corr_out, 'Validation dataset') + # change plot name depending on input
  theme(text = element_text(size = 8),
        legend.key.width = unit(2, 'mm'),
        legend.key.height = unit(5, 'mm')) 

b2_out_plot <- plot_pca(pca_b2_out, b2_corr_out, 'Training dataset') + # change plot name depending on input
  theme(text = element_text(size = 8),
        legend.key.width = unit(2, 'mm'),
        legend.key.height = unit(5, 'mm')) 

#

b1_out_plot
b2_out_plot

#

########################

# Figure S7
s7 <- arrangeGrob(b2_out_plot, b1_out_plot, nrow = 1,
                  widths = c(0.42, 0.58),
                  top = textGrob('PCA score plots of breath VOC data following outlier removal',
                                 gp = gpar(fontsize = 10, fontface = 'bold')))

plot(s7)

ggsave('figS7.tiff',s7, unit = 'mm', dpi = 300, width = 157, height = 75)

pdf('FigureS7.pdf', width = 6.18, height = 2.95)
plot(s7)
dev.off()

#########################

#

pca_plots <- arrangeGrob(b1_plot, b1_corr_plot, b1_out_plot,
                         b2_plot, b2_corr_plot, b2_out_plot,
                         nrow = 2, widths = c(0.28, 0.28, 0.44),
                         top = textGrob('PCA score plots of breath VOC datasets'))

dev.new()
plot(pca_plots)

ggsave('PCA_BGcorr_effect_scores.tiff', 
       pca_plots, dpi = 300, unit = 'mm', width = 157, height = 95)


#

# loadings plots
cols <- hue_pal()(4)
swatch(cols)

plot_loadings <- function(pca_b, PC){
  loadings <- as.data.frame(pca_b$loadings$X) %>%
    mutate(comp = rownames(.)) %>%
    left_join(endo_exo %>% dplyr::select(comp, FC_filter, adj)) %>%
    rename(pc = PC) %>%
    arrange(desc(abs(pc))) %>%
    slice(1:10)
  
  plot <- loadings %>% 
    ggplot(aes(x = pc, y = fct_inorder(as.factor(comp)))) +#, fill = FC_filter)) +
    geom_col_pattern(aes(pattern = adj), pattern_fill = 'white',
                     pattern_size = 0.5, alpha = 1,
                     pattern_colour = NA) +
    #geom_col() +
    theme_bw(base_size = 10) +
    scale_pattern_manual(name = 'BG correction', 
                         values = c('no' = 'none',
                                    'yes' = 'stripe')) +
    #scale_fill_manual(values = c('Endo_both_datasets' = cols[1],
    #                             'Endo_dataset_1' = cols[2],
    #                             'Exo' = cols[4],
    #                             'Endo_dataset_2' = cols[3])) +
    theme(#legend.position = 'none',
          plot.title = element_text(hjust = 0.5)) +
    ylab('') +
    xlab(PC)
  
} 

l1_corr <- plot_loadings(pca_b1_corr, 'PC1')
l1_corr2 <- plot_loadings(pca_b1_corr, 'PC2')

l1_corr
l1_corr2

dev.new()
lplots <- arrangeGrob(l1_corr, l1_corr2, nrow = 1, widths = c(0.45, 0.55),
                      top = textGrob('PCA loadings in background-corrected dataset'))
plot(lplots)

#

l1_corr_out <- plot_loadings(pca_b1_out, 'PC1') + ggtitle('Dataset 1')
l2_corr_out <- plot_loadings(pca_b2_out, 'PC1') + ggtitle('Dataset 2')

# create figure legends
fill_plot <- plot_loadings(pca_b2_out, 'PC1') + 
  theme(text = element_text(size = 8),
        legend.key.width = unit(3, 'mm'),
        legend.key.height = unit(5, 'mm'))

fill_legend <- get_legend(fill_plot)

pattern_plot <- plot_loadings(pca_b1_out, 'PC1') + 
  theme(text = element_text(size = 8),
        legend.key.width = unit(3, 'mm'),
        legend.key.height = unit(5, 'mm'))

pattern_legend <- get_legend(pattern_plot)

#

l_plots <- arrangeGrob(l1_corr_out, l2_corr_out, 
                       arrangeGrob(pattern_legend, fill_legend, nrow = 2),
                       nrow = 1, widths = c(0.37, 0.43, 0.2))
dev.new()
plot(l_plots)

ggsave('PCA_BGcorr_effect_loadings_PC1.tiff', 
       l_plots, dpi = 300, unit = 'mm', width = 157, height = 70)


#
#
#
