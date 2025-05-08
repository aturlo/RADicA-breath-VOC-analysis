## Missing value imputation

# author: Aggie Turlo
# project: RADicA
# date: 17/01/2025

###############################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(GMSimpute)
library(ggplotify)
library(stringr)
library(gridExtra)
library(lme4)
library(mixOmics)
library(cowplot)
library(grafify)

#

## load data
# load filtered study data with missing values
b1_all_f <- read.csv('RADicA_B1_NAfiltered.csv')[,-1]
b2_all_f <- read.csv('RADicA_B2_NAfiltered.csv')[,-1]

rownames(b1_all_f) <- b1_all_f$Sample
rownames(b2_all_f) <- b2_all_f$Sample

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

#
#
#

## EXPLORATION OF MISSING DATA TYPES
# specify dataset for analysis
b_all_f <- b1_all_f
dat <- 'B1'

#

# change dataset to long format
b_all_fL <- b_all_f %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                     names_to = 'comp', values_to = 'peakArea')
## NAs according to peak intensity
sumint <- b_all_fL %>% 
  group_by(comp, class) %>%
  summarise(med = median(log(peakArea), na.rm = TRUE),
            nonas = sum((is.na(peakArea) == TRUE)/length(peakArea))*100)

# visualise relationship between % of missing values and feature median abundance
pdf(paste('NAvsMedianIntensity_', dat, '.pdf', sep =''), width = 7, height = 4.5)

medVSna <- sumint %>% group_by(class) %>%
  arrange(med) %>%
  do(plots = ggplot(data = ., aes(x = med, y = nonas)) +
       geom_point(alpha = 0.6) + ggtitle(.$class) + theme_bw() +
       theme(axis.title.x = element_blank(),
             axis.title.y = element_blank()))

marrangeGrob(medVSna$plots, nrow = 2, ncol = 3)

dev.off()

################################

# Figure S2B

medVSna_b1 <- sumint %>% filter(class %in% c('S1', 'S2')) %>%
  ggplot(aes(x = med, y = nonas)) +
  geom_point(alpha = 0.6, size = 0.8) + theme_bw(base_size = 10) +
  ylab('% of missing values') + xlab('Median peak intensity') +
  facet_wrap(~ class) +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  ggtitle('Validation dataset')


medVSna_b2 <- sumint %>% filter(class %in% c('S1', 'S2')) %>%
  ggplot(aes(x = med, y = nonas)) +
  geom_point(alpha = 0.6, size = 0.8) + theme_bw(base_size = 10) +
  ylab('% of missing values') + xlab('Median peak intensity') +
  facet_wrap(~ class) +
  theme(strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  ggtitle('Training dataset')

s2_b <- arrangeGrob(medVSna_b2, medVSna_b1, nrow = 1, ncol = 2, 
                    top = textGrob('Missing value ratio in breath VOCs',
                                   gp = gpar(fontsize = 10, fontface = 'bold')))
plot(s2_b)

##################################

#

# number of complete variables across groups
sumint %>% group_by(class) %>% summarise(complete_var = sum(nonas == 0)/n_distinct(comp))

#

# Fraction of different types of missingness across breath samples
s_tech <- b_all_fL %>% filter(class %in% c('S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Analysis_date, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) %>%
  #mutate(Sample = str_sub(Sample, end = -13)) %>% # Batch 1
  pivot_wider(id_cols = c(RAD_ID, comp, Batch, Analysis_date, CoreVisit), 
              names_from = class, values_from = peakArea) %>%
  mutate(techNA = ifelse(is.na(S1) & is.na(S2), 'Missing2', 
                         ifelse(is.na(S1) | is.na(S2), 'Missing1', 
                                'Non-missing'))) %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_'))

s_sum <- s_tech %>% group_by(techNA, Batch) %>% summarise(n = n()) 

s_sum1 <- s_sum %>% left_join(s_sum %>% group_by(Batch) %>% summarise(BatchSum = sum(n))) %>%
  mutate(fraqBatch = n/BatchSum *100)

s_sum %>% group_by(techNA) %>% summarise(fraqAll = sum(n)/sum(s_sum$n)*100)

s_sum1 %>%
  ggplot(aes(x = as.factor(Batch), y = fraqBatch, fill = techNA)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_brewer(palette = 'Set2') +
  ylab('% of all observations') +
  xlab('Batch') +
  ggtitle('Missingness in breath sample replicates') +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.title = element_blank())

ggsave(paste('Missingness_S1S2_barchart_', dat, '.tiff', sep =''), dpi = 300, unit = 'mm', width = 100, height = 80)

#

# Effect of missingness on abundance of all VOCs within sample
Mis2 <- s_tech %>% filter(techNA == 'Missing2') %>%
  group_by(comp, Batch) %>% summarise(n = n())

violin_plots <- s_tech %>% filter(techNA != 'Missing2') %>%
  pivot_longer(cols = c(S1, S2) , names_to = 'Replicate', values_to = 'peakArea') %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  group_by(Sample) %>%
  do(plots = ggplot(data = ., aes(x = techNA, y = log(peakArea), colour = Replicate)) +
       geom_violin(aes(colour = Replicate), position = position_dodge(0.9)) + 
       geom_boxplot(aes(colour = Replicate), position = position_dodge(0.9), width = 0.1) + 
       theme_bw(base_size = 8) +
       ggtitle(str_sub(.$Sample)) +
       theme(legend.position = 'none') +
       scale_color_brewer(palette = 'Set2'))

pdf(paste('Missing_pattern_sample_', dat, '.pdf'))
marrangeGrob(violin_plots$plots, nrow = 4, ncol = 4)
dev.off()

#################################

# Figure S2C
pals <- grafify::graf_palettes
my_pal = pals$fishy

s3_b1 <- s_tech %>% filter(techNA != 'Missing2') %>%
  pivot_longer(cols = c(S1, S2) , names_to = 'Replicate', values_to = 'peakArea') %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  filter(Sample %in% c('RAD002_CV1', 'RAD002_CV2')) %>%
  ggplot(aes(x = techNA, y = log(peakArea), colour = Replicate)) +
  geom_violin(aes(colour = Replicate), position = position_dodge(0.9)) + 
  geom_boxplot(aes(colour = Replicate), position = position_dodge(0.9), width = 0.1) + 
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  facet_wrap(~ Sample) +
  scale_color_manual(values = c('S1' = my_pal[[4]], 'S2' = my_pal[[5]])) +
  ggtitle('Validation dataset')

s3_b1 

s3_b2 <- s_tech %>% filter(techNA != 'Missing2') %>%
  pivot_longer(cols = c(S1, S2) , names_to = 'Replicate', values_to = 'peakArea') %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  filter(Sample %in% c('RAD110_CV1', 'RAD110_CV2')) %>%
  ggplot(aes(x = techNA, y = log(peakArea), colour = Replicate)) +
  geom_violin(aes(colour = Replicate), position = position_dodge(0.9)) + 
  geom_boxplot(aes(colour = Replicate), position = position_dodge(0.9), width = 0.1) + 
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  facet_wrap(~ Sample) +
  scale_color_manual(values = c('S1' = my_pal[[4]], 'S2' = my_pal[[5]])) +
  ggtitle('Training dataset')

s3_b2 

s3_c <- arrangeGrob(s3_b2, s3_b1, nrow = 1, top = textGrob('VOC abundance in breath sample replicates',
                                                           gp = gpar(fontsize = 10, fontface = 'bold')))
plot(s3_c)

##################################

#
#
#

## EFFECT OF EXPERIMENTAL FACTORS ON MISSINGNESS
# heatmap of NAs according to experimental condition and sample type
# specify dataset 
data <- b1_all_f
dat <- 'B1'

#

data_heat <- data[,6:ncol(data)]
data_heat[is.na(data_heat)] <- 0
data_heat[data_heat > 0] <- 1
rownames(data_heat) <- data$Sample

data$class <- factor(data$class, levels = c('ES', 'Blank', 'BG', 'S1', 'S2'))



cols1 <- my_pal[1:(n_distinct(data$class))]
cols2 <- my_pal[(10-(n_distinct(data$Batch))):9] 
names(cols1) <- c('ES', 'Blank', 'BG', 'S1', 'S2')
names(cols2) <- unique(data$Batch)
annotCol <- list(class = cols1,
                 Batch = cols2)

data$Batch <- as.factor(data$Batch)
rownames(data) <- data$Sample


tiff(paste('Heatmap_NA_', dat, '_class_batch.tiff', sep =''), res = 300, width = 250, height = 200, unit = 'mm')
heatmap <- pheatmap(as.matrix(t(data_heat)), 
                    annotation_col = data %>% 
                      dplyr::select(Batch, class), 
                    show_colnames = FALSE, cluster_cols =  T, 
                    show_rownames = FALSE, cluster_rows = FALSE,
                    fontsize = 10, main = 'Missing values across samples',
                    #'Training dataset', # title for supplementary figure
                    color = c('Tomato1', 'Black'),
                    treehight_col = 100,
                    annotation_colors = annotCol,
                    treeheight_col = 70,
                    legend = FALSE) # use for analysis date mapping
heatmap
dev.off()

#

# in breath samples only, including diagnosis and core visit
b_s <- data %>% filter(class %in% c('S1','S2')) %>% 
  dplyr::select(!1:5)

data_heat1 <- b_s
data_heat1[is.na(data_heat1)] <- 0
data_heat1[data_heat1 > 0] <- 1
rownames(data_heat1) <- rownames(b_s)

b_all_an <- data %>% filter(class %in% c('S1','S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

data_heat1 <- data_heat1 %>% filter(rownames(.) %in% b_all_an$Sample)

b_all_an$Batch <- as.factor(b_all_an$Batch)

rownames(b_all_an) <- b_all_an$Sample 

cols3 <- my_pal[1:(n_distinct(b_all_an$Batch))]
cols4 <- my_pal[5:(4 + (n_distinct(b_all_an$Diagnosis)))]
cols5 <- my_pal[7:(6 + (n_distinct(b_all_an$CoreVisit)))]
names(cols3) <- na.omit(unique(b_all_an$Batch))
names(cols4) <- na.omit(unique(b_all_an$Diagnosis))
names(cols5) <- unique(b_all_an$CoreVisit)
annotCol1 <- list(Batch = cols3,
                  Diagnosis = cols4,
                  CoreVisit = cols5)


tiff(paste('Heatmap_NAs_', dat, '_Breath.tiff', sep = ''), res = 300, width = 250, height = 200, unit = 'mm')
heatmap1 <- pheatmap(as.matrix(t(data_heat1)), 
                     annotation_col = b_all_an %>% dplyr::select(Batch, Diagnosis, CoreVisit), 
                     show_colnames = FALSE, cluster_cols =  T, 
                     show_rownames = FALSE, cluster_rows = FALSE,
                     fontsize = 10, main = #'Missing values across breath samples',
                       'Validation dataset',
                     color = c('Tomato1', 'Black'),
                     annotation_colors = annotCol1,
                     treeheight_col = 70,
                     legend = FALSE)

heatmap1
dev.off()

##################################

# Figure S3

s3_a <- heatmap1
s3_b <- heatmap1

s3 <- plot_grid(s3_a$gtable, s3_b$gtable, nrow = 2)
plot(s3)

ggsave('figS3.tiff', dpi = 300, unit = 'mm', width = 120, height = 210)

pdf('FigureS3.pdf', width = 4.72, height = 7.87)
plot(s3)
dev.off()

#
#
#


###############################

# IMPUTATION
# specify dataset for analysis
df <- b2_all_f
rownames(df) <- df$Sample

#
#
#

# examine colinearity of VOC log peak areas
cor_vals <- cor(log(df[,-c(1:5)]), use = 'pairwise.complete.obs')
sum(table(cor_vals[cor_vals > 0.5]) %>% as.data.frame() %>% dplyr::select(Freq))/(243^2)


# create assay dataset for each class
assay <- function(x){
  df %>% filter(class == x) %>% dplyr::select(6:ncol(.))
}

# 
assay_blank <- assay('Blank') %>% filter(rownames(.) %ni% '200131_RaDICA_Batch372_tracker_497813_1') # B1 poor quality blank sample
assay_es <- assay('ES') 
assay_bg <- assay('BG')
assay_s1 <- assay('S1')
assay_s2 <- assay('S2')


# use GMSimpute (two-step Lasso) 
impute <- function(x){
  t(GMS.Lasso(t(x), log.scale = TRUE, TS.Lasso = TRUE)) %>%
    as.data.frame()
}

assay_blank_imp <- impute(assay_blank)
assay_es_imp <- impute(assay_es)
assay_bg_imp <- impute(assay_bg)
assay_s1_imp <- impute(assay_s1)
assay_s2_imp <- impute(assay_s2)

# merge data class subsets into one dataset
imp_long <- function(assay) {
  assay %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols =! Sample, names_to = 'comp', values_to = 'peakArea') %>%
    left_join(df %>% dplyr::select(Sample, Date, Time, Batch, class))

}

blank_impL <- imp_long(assay_blank_imp)
es_impL <- imp_long(assay_es_imp)
bg_impL <- imp_long(assay_bg_imp)
s1_impL <- imp_long(assay_s1_imp)
s2_impL <- imp_long(assay_s2_imp)

b_imp_L <- rbind(es_impL, blank_impL, bg_impL, s1_impL, s2_impL) 

#
#
#

# change new object name accordingly
b1_imp <- b_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea)
b2_imp <- b_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea)

write.csv(b1_imp, 'RADicA_B1_NAfiltered_imputed.csv')
write.csv(b2_imp, 'RADicA_B2_NAfiltered_imputed.csv')

#
#
#

################################

# EFFECT OF IMPUTATION ON DATA
# specify dataset for analysis
b_imp <- b1_imp %>% as.data.frame()
b_all_f <- b1_all_f  %>% filter(Sample %ni% '200131_RaDICA_Batch372_tracker_497813_1')
dat <- 'B2'

# change data format to long
b_imp_L <- b_imp %>% pivot_longer(cols = c(6:ncol(b_imp)), names_to = 'comp', values_to = 'peakArea')
b_all_fL <- b_all_f %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                     names_to = 'comp', values_to = 'peakArea')

#
#
#

# density plots
dens_plot <-
  b_imp_L %>% 
  filter(is.na(class) == FALSE) %>%
  ggplot() + 
  geom_density(aes(x = log(peakArea), colour = 'Imputed dataset')) +
  geom_density(data = b_all_fL %>% drop_na(),
                 aes(x = log(peakArea), colour = 'Observed dataset')) +
  theme_bw(base_size = 10) +
  xlab('log( Peak area)') +
  facet_wrap(~ class, scale = 'free') +
  xlim(0, NA) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

dens_plot

ggsave(paste('Density_imputed_', dat, '.tiff', sep = ''), dens_plot, dpi = 300, units = 'mm', width = 240, height = 120)

#

#################################

# Figure S4A

s4a_b2 <- b_imp_L %>% 
  filter(class %in% c('S1', 'S2')) %>%
  ggplot() + 
  geom_density(aes(x = log(peakArea), colour = 'Imputed dataset'), lwd = 0.3) +
  geom_density(data = b_all_fL %>% drop_na() %>% 
                 filter(class %in% c('S1', 'S2')),
               aes(x = log(peakArea), colour = 'Observed dataset'), lwd = 0.3) +
  theme_bw(base_size = 10) +
  xlab('log( Peak area)') +
  facet_wrap(~ class, scale = 'free') +
  xlim(0, NA) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        legend.position = 'none',
        axis.title = element_text(size = 8)) +
  ggtitle('Training dataset')


s4a_b2

s4a_b1 <- b_imp_L %>% 
  filter(class %in% c('S1', 'S2')) %>%
  ggplot() + 
  geom_density(aes(x = log(peakArea), colour = 'Imputed dataset'), lwd = 0.3) +
  geom_density(data = b_all_fL %>% drop_na() %>%
                 filter(class %in% c('S1', 'S2')),
               aes(x = log(peakArea), colour = 'Observed dataset'), lwd = 0.3) +
  theme_bw(base_size = 10) +
  xlab('log( Peak area)') +
  facet_wrap(~ class, scale = 'free') +
  xlim(0, NA) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 8)) +
  ggtitle('Validation dataset')

s4a_b1

s4a <- arrangeGrob(s4a_b2 + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm')), 
                   s4a_b1 + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm')), 
                   nrow = 1, widths = c(0.37, 0.63),
                   top = textGrob('Breath sample VOC peak area distributions',
                   gp = gpar(fontsize = 10, fontface = 'bold')))
plot(s4a)

##################################

#

# time-course change in dataset median
b_imp_L$Date <- as.Date(b_imp_L$Date)
b_all_fL$Date <- as.Date(b_all_fL$Date)

time_plot <- 
  b_imp_L %>% 
  drop_na(class) %>%
  group_by(class, Date) %>%
  summarise(medianImp = median(log(peakArea), na.rm = TRUE)) %>%
  ggplot(aes(x = Date)) + 
  #geom_line(aes(y = medianObs, colour = 'Observed values')) +
  geom_line(aes(y = medianImp, group = 'Imputed dataset', colour = 'Imputed dataset')) +
  geom_line(data = b_all_fL %>% 
              drop_na(class) %>%
              group_by(class, Date) %>%
              summarise(medianObs = median(log(peakArea), na.rm = TRUE)),
            aes(y = medianObs, group = 'Observed dataset', colour = 'Observed dataset')) +
  theme_bw() +
  facet_wrap(~ class, scale = 'free') +
  ylab('median log( Peak area)') 

time_plot

ggsave(paste('Time_imputed_', dat, '.tiff', sep =''), time_plot, dpi = 300, units = 'mm', width = 240, height = 120)

#

## PCA complete vs imputed variables
# create data class subsets with complete variables only
complete_set <- function(x, y, CL){ 
  df <- x %>% filter(class == CL)
  rownames(df) <- df$Sample
  sumint <- y %>% 
    group_by(comp, class) %>%
    summarise(nonas = sum(is.na(peakArea) == TRUE))
  comp_vars <- sumint %>% filter(class == CL) %>% filter(nonas == 0)
  comp_df <- df[,comp_vars$comp]
}

std <- c('Acetone', 'Isoprene', 'Benzene', 'X3_Pentanone', 'X1.4_dioxane',
          'Pyridine', 'Toluene', 'Octane', 'P_Xylene', 'Nonane', 'Benzaldehyde',
          'X1_Heptanol', 'Decane', 'X3_Carene', 'Limonene', 'Undecane', 'Nonanal',
          'X1.2.3.4_tetrahydronaphthalene', 'Dodecane', 'X1_methylindole', 'Tridecane',
          'Pentadecane')

b_blank_f_c <- complete_set(b_all_f, b_all_fL, 'Blank')
b_es_f_c <- b_all_f %>% filter(class == 'ES') %>% subset(select = colnames(.) %in% std) # for ES use spike-in data only
b_es_f_c <- complete_set(b_all_f, b_all_fL, 'ES')
b_bg_f_c <- complete_set(b_all_f, b_all_fL, 'BG')
b_s1_f_c <- complete_set(b_all_f, b_all_fL, 'S1')
b_s2_f_c <- complete_set(b_all_f, b_all_fL, 'S2')

# create data class subsets with imputed observations
b_imp_blank <- complete_set(b_imp, b_imp_L, 'Blank')
b_imp_es <- complete_set(b_imp, b_imp_L, 'ES')
b_imp_bg <- complete_set(b_imp, b_imp_L, 'BG')
b_imp_s1 <- complete_set(b_imp, b_imp_L, 'S1')
b_imp_s2 <- complete_set(b_imp, b_imp_L, 'S2')

# create data annotation subset for the experimental variables (group, class)
b_es_f <- b_all_f %>% filter(class == 'ES') 
b_blank_f <- b_all_f %>% filter(class == 'Blank') 
b_bg_f <- b_all_f %>% filter(class == 'BG') 
b_s1_f <- b_all_f %>% filter(class == 'S1') 
b_s2_f <- b_all_f %>% filter(class == 'S2') 

#

plot_pca <- function(x, y, z){
  pca <- mixOmics::pca(log(x), scale = TRUE, center = TRUE)
  plot_pca <- plotIndiv(pca,
                        pch = 1,
                        cex = 0.6,
                        group = y$Batch,
                        title = z,
                        size.title = 8,
                        size.xlabel = 8,
                        size.ylabel = 8,
                        size.legend.title = 6,
                        size.axis = 6)
  plot_pca <- plot_pca$graph + theme(strip.background = element_blank(),
                                     plot.title = element_text(face = 'plain'),
                                     plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm')) + 
    scale_colour_grafify(palette = 'fishy', name = 'Batch')
    
}

#

pc_blank_c <- plot_pca(b_blank_f_c, b_blank_f, 'Blank')
pc_blank_i <- plot_pca(b_imp_blank, b_blank_f, 'Blank')

pc_es_c <- plot_pca(b_es_f_c, b_es_f, 'ES')
pc_es_i <- plot_pca(b_imp_es, b_es_f, 'ES')

pc_bg_c <- plot_pca(b_bg_f_c, b_bg_f, 'BG')
pc_bg_i <- plot_pca(b_imp_bg, b_bg_f, 'BG')

pc_s1_c <- plot_pca(b_s1_f_c, b_s1_f, 'S1')
pc_s1_i <- plot_pca(b_imp_s1, b_s1_f, 'S1')

pc_s2_c <- plot_pca(b_s2_f_c, b_s2_f, 'S2')
pc_s2_i <- plot_pca(b_imp_s2, b_s2_f, 'S2')

pc_plots <- arrangeGrob(pc_blank_c$graph, pc_es_c$graph, pc_bg_c$graph, pc_s1_c$graph, pc_s2_c$graph,
                        pc_blank_i$graph, pc_es_i$graph, pc_bg_i$graph, pc_s1_i$graph, pc_s2_i$graph,
                        ncol = 5, nrow = 2)

dev.new()
plot(pc_plots)

ggsave(paste('PCA_complete_vs_imputed_' , dat, '.tiff', sep = ''),
       pc_plots, dpi = 300, unit = 'mm', width = 310, height = 120)

#################################

# Figure S4B

s4b_b2 <- arrangeGrob(pc_s1_c, 
                      pc_s2_c, 
                      pc_s1_i, 
                      pc_s2_i, 
                      nrow = 2,
                      top = textGrob('Training dataset',
                                     gp = gpar(fontsize = 10)))
plot(s4b_b2)

s4b_b1 <- arrangeGrob(pc_s1_c, 
                      pc_s2_c, 
                      pc_s1_i, 
                      pc_s2_i, 
                      nrow = 2,
                      top = textGrob('Validation dataset',
                                     gp = gpar(fontsize = 10)))
plot(s4b_b1)

#

get_legend_35 <- function(plot) {
  # return all legend candidates
  legends <- get_plot_component(plot, "guide-box", return_all = TRUE)
  # find non-zero legends
  nonzero <- vapply(legends, \(x) !inherits(x, "zeroGrob"), TRUE)
  idx <- which(nonzero)
  # return first non-zero legend if exists, and otherwise first element (which will be a zeroGrob) 
  if (length(idx) > 0) {
    return(legends[[idx[1]]])
  } else {
    return(legends[[1]])
  }
}

legb2 <- pc_s1_c + theme(legend.position = 'bottom') +
  scale_colour_grafify(palette = 'fishy', name = 'Batch \n (training dataset)') +
  theme(legend.title = element_text(size = 8))
legb2a <- get_legend_35(legb2)
plot(legb2a)


legb1 <- pc_s1_c + theme(legend.position = 'bottom') +
  scale_colour_grafify(palette = 'fishy', name = 'Batch \n (validation dataset)') +
  theme(legend.title = element_text(size = 8))
legb1a <- get_legend_35(legb1)
plot(legb1a)

legs <- arrangeGrob(legb2a, legb1a, nrow = 1)
plot(legs)

#

s4b <- arrangeGrob(s4b_b2, s4b_b1, nrow = 1, widths = c(0.44, 0.44),
                   top = textGrob('PCA score plots of breath sample VOC profiles',
                                  gp = gpar(fontface = 'bold', fontsize = 10)))
plot(s4b)



#



################################

# Assemble figure S4
s4 <- plot_grid(s4a, s4b, legs, nrow = 3, rel_heights = c(0.3, 0.65, 0.05),
                labels = c('a)', 'b)'))

dev.new()
plot(s4)

ggsave('figS4.tiff', unit = 'mm', dpi = 300, width = 157, height = 145)

pdf('FigureS4.pdf', width = 6.18, height = 5.7)
plot(s4)
dev.off()

################################

#
#
#

## Effect of imputation on technical replicates in breath samples
s_tech <- b_all_fL %>% filter(class %in% c('S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Analysis_date, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) %>%
  #mutate(Sample = str_sub(Sample, end = -13)) %>% # Batch 1
  pivot_wider(id_cols = c(RAD_ID, comp, Batch, Analysis_date, CoreVisit), 
              names_from = class, values_from = peakArea) %>%
  mutate(techNA = ifelse(is.na(S1) & is.na(S2), 'Missing2', 
                         ifelse(is.na(S1) | is.na(S2), 'Missing1', 
                                'Non-missing'))) %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_'))

s_tech_imp <- b_imp_L %>% filter(class %in% c('S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Analysis_date, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  dplyr::select(!Time) %>%
  pivot_wider(id_cols = c(Sample, comp), names_from = class, values_from = peakArea) %>%
  left_join(s_tech %>% dplyr::select(Sample, comp, Batch, techNA))

#

violin_plots_imp <- s_tech_imp %>% filter(techNA != 'Missing2') %>%
  pivot_longer(cols = c(S1, S2) ,names_to = 'Replicate', values_to = 'peakArea') %>%
  group_by(Sample) %>%
  do(plots = ggplot(data = ., aes(x = techNA, y = log(peakArea), colour = Replicate)) +
       geom_violin(aes(colour = Replicate), position = position_dodge(0.9)) + 
       geom_boxplot(aes(colour = Replicate), position = position_dodge(0.9), width = 0.1) + 
       theme_bw(base_size = 8) +
       ggtitle(str_sub(.$Sample)) +
       theme(legend.position = 'none') +
       scale_color_brewer(palette = 'Set2'))

pdf(paste('Missing_pattern_sample_imp_', dat, '.pdf'))
marrangeGrob(violin_plots_imp$plots, nrow = 4, ncol = 4)
dev.off()

#

# compare distributions of breath data with and without imputation
s_tech1 <- s_tech %>% filter(techNA == 'Missing1')
s_tech_imp1 <- s_tech_imp %>% filter(techNA == 'Missing1')

nas <- s_tech1 %>% group_by(Sample) %>% summarise(n = n())

# exclude comparisons with < 10 observations 
ncomp <- n_distinct(s_tech1$comp)-10
nas1 <- nas %>% filter(nas$n < ncomp & nas$n> 10)

# Wilcoxon paired test
imp_test <- 
  sapply(nas1$Sample, function(id){ 
  subset <- s_tech1 %>%
    filter(Sample == id) %>%
    pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea')
  
  subset_imp <- s_tech_imp1 %>%
    filter(Sample == id) %>%
    pivot_longer(cols = c(S1, S2), names_to = 'Replicate', values_to = 'peakArea') %>%
    left_join(subset, by = c('Sample', 'comp', 'Replicate')) %>%
    mutate(Imp = ifelse(is.na(peakArea.y) == TRUE, 'Yes', 'No')) %>%
    pivot_wider(id_cols = c('Sample', 'comp'), names_from = Imp, values_from = peakArea.x)
  
  t <- wilcox.test(subset_imp$No, subset_imp$Yes)
  out <- data.frame(p = t$p.value, n = nrow(subset_imp))
})

imp_test_df <- imp_test %>% t() %>% as.data.frame() %>%
  mutate(p = unlist(p),
         n = unlist(n),
         adj.p = p.adjust(p, method = 'BH')) 

# edit file name accordingly
write.csv(imp_test_df, 'Wilcox_paired_breath_imputedVSobserved_B1.csv')

#
#
#

#################################

# Figure S2 assembly
s2 <- plot_grid(s2_a, s2_b, s3_c, nrow = 3, labels = c('a)', 'b)', 'c)'),
                label_size = 12,
                rel_heights = c(0.325, 0.35, 0.325))
plot(s2)

ggsave('figS2.tiff', s2, unit = 'mm', dpi = 300, width = 150, height = 150)

pdf('FigureS2.pdf', width = 6.18, height = 6.18)
plot(s2)
dev.off()

#
#
#