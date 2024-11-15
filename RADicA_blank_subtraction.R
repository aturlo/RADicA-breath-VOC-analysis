# RADicA blank subtraction 

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggfortify)
library(gridExtra)

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')
b1_imp_L <- read.csv('RADicA_B1_NAfiltered_imputed_long.csv')
b1_all_f <- read.csv('RADicA_B1_NAfiltered_v2.csv')
data <- read.csv('240508_Radica_VOC_peak_areas_v5.1.csv', check.names = FALSE)

# format the imputed dataset 
b1_imp_L <- b1_imp_L[,-1]
b1_imp <- b1_imp_L %>% 
  dplyr::select(!c(logPeakAreaImp, logPeakArea, Date)) %>%
  pivot_wider(names_from = comp, values_from = peakAreaImp)

b1_imp <- b1_imp[ order(match(b1_imp$Sample, data$Sample)), ] %>%
  as.data.frame()

# format the filtered (non-imputed) dataset
b1_all_fL <- b1_all_f[,3:ncol(b1_all_f)] %>% 
  pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
               names_to = 'comp', values_to = 'peakArea')

#

# Wilcoxon test blank vs study groups for each variable
# compunds with missing blank values
na_comps <- c("Anethole","Cyclohexane._2_propyl_1.1.3_trimethyl_", "Geranyl_acetone",
              "Dodecane._2.6.11_trimethyl_", 'Eucalyptol')

wilcox_blank <- function(CL){
  blank <- b1_imp %>% filter(class == 'Blank') %>% dplyr::select(!all_of(na_comps))
  sample <- b1_imp %>% filter(class == CL) %>% dplyr::select(!all_of(na_comps))
  
  results <- bind_rows(lapply(4:ncol(blank),function(x) {
  blank_comp <- blank[,x]
  sample_comp <- sample[,x]
  test <- wilcox.test(blank_comp,
                      sample_comp)
  df <- data.frame(comp = colnames(b1_imp)[x],
                   wilcox.p = test$p.value,
                   Mblank = median(blank_comp),
                   Msample = median(sample_comp))
  }))
  
  results <- results %>% mutate(adj.wilcox.p = p.adjust(wilcox.p, method = 'BH'))
  }

blank_vs_bg <- wilcox_blank('BG')
b_bg <- blank_vs_bg %>% filter(adj.wilcox.p > 0.05)
View(blank_vs_bg %>% filter(adj.wilcox.p < 0.05) %>% filter(Mblank > Msample))

blank_vs_s1 <- wilcox_blank('S1')
b_s1 <- blank_vs_s1 %>% filter(adj.wilcox.p > 0.05)
View(blank_vs_s1 %>% filter(adj.wilcox.p < 0.05) %>% filter(Mblank > Msample))

blank_vs_s2 <- wilcox_blank('S2')
b_s2 <- blank_vs_s2 %>% filter(adj.wilcox.p > 0.05)
View(blank_vs_s2 %>% filter(adj.wilcox.p < 0.05) %>% filter(Mblank > Msample))

setdiff(b_s1$comp, b_s2$comp)

#

# using gamma regression instead
gamma_blank <- function(CL){
  subset <- b1_imp_L %>% filter(class %in% c('Blank', CL)) #%>% dplyr::select(!all_of(na_comps))
  results <- bind_rows(lapply(colnames(blank[4:ncol(blank)]),function(x) {
    subset1 <- subset %>% filter(comp == x) %>% drop_na() 
    test <- glm(peakArea ~ class,
                family = 'Gamma'(link = 'identity'),
                data = subset1)
    df <- c(x, coef(summary(test))[2,]) %>% as.data.frame() %>% t() %>% as.data.frame()
  }))
  colnames(results)[5] <- 'p.value'
  results <- results %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
}

b_bg_G <- gamma_blank('BG') %>% filter(adj.p.value > 0.05) 
b_s1_G <- gamma_blank('S1') %>% filter(adj.p.value > 0.05)
b_s2_G <- gamma_blank('S2') %>% filter(adj.p.value > 0.05)

remove <- intersect(b_s2_G$comp, b_s1_G$comp)
remove <- remove[remove != 'Internal_Standard']

# remove variables with equal distribution in Blank/S1/S2
b1_imp <- b1_imp %>% dplyr::select(!all_of(remove))

#
#
#

# BLANK CORRECTION
# add analysis time
b1_imp <- b1_imp %>% left_join(data %>% dplyr::select(Sample, Acq_Date_Time)) %>%
  relocate(Acq_Date_Time) %>%
  separate(Acq_Date_Time, c('Date', 'Time'), sep = ' ')

b1_imp$Date <- as.Date(b1_imp$Date, format = "%d/%m/%Y")

b1_imp <- b1_imp %>% mutate(Date_Time = paste(Date, Time)) %>%
  relocate(Date_Time) %>%
  dplyr::select(!c(Date, Time))

#

b1_all_f <- b1_all_f %>% mutate(Date_Time = paste(Date, Time)) %>%
  relocate(Date_Time) %>%
  dplyr::select(!c(Date, Time))

b1_all_f1 <- b1_all_f %>% dplyr::select(!c(Analysis_date, RAD_ID, X, Batch))

# subtract peak areas in the last blank measured before the sample from sample peak areas
df = b1_all_f1
CL = 'S1'
sample = '200121_RaDICA_RAD096_S1_373830_1'

corr_blank <- function(df, CL) {
  b1_cl <- df %>% filter(class == CL)
  names <- b1_cl$Sample
  corr_df <- bind_rows(lapply(names, function(sample){
  sample <- b1_imp %>% filter(Sample == sample)
  date_time <- sample$Date_Time
  sample_assay <- sample %>% dplyr::select(!c(Date_Time, Sample, class, Internal_Standard))
  near_blank <- df %>% filter(class == 'Blank' & Date_Time < date_time) %>%
  tail(n = 1) %>% dplyr::select(!c(Date_Time, Sample, class, Internal_Standard))
  near_blank[is.na(near_blank) == TRUE] <- 0
  corr_sample <- sample_assay-near_blank
  corr_sample1 <- cbind(sample %>% dplyr::select(c(Date_Time, Sample, class)),
                     corr_sample, sample %>% dplyr::select(Internal_Standard))
  }))
  corr_assay <- corr_df[,4:ncol(corr_df)]
  corr_assay[corr_assay < 0] <- min(corr_assay[corr_assay > 0]/2)
  corr_assay[corr_assay == 0] <- min(corr_assay[corr_assay > 0]/2)
  corr_df1 <- cbind(corr_df[,1:3], corr_assay) 
}


bg_corr <- corr_blank(b1_imp, 'BG')
s1_corr <- corr_blank(b1_imp, 'S1')
s2_corr <- corr_blank(b1_imp, 'S2')

# concatenate corrected sample values with QCs
b1_imp_corr <- rbind(b1_imp %>% filter(class %in% c('ES', 'Blank')),
                                       bg_corr, 
                                       s1_corr,
                                       s2_corr) 

write.csv(b1_imp_corr, 'RADicA_B1_NAfiltered_imputed_blankCorr.csv')

b1_imp_corrL <- b1_imp_corr %>% pivot_longer(cols =! c(1:3),
                                             names_to = 'comp',
                                             values_to = 'peakArea')

#

# check distribution of 0s
b1_imp_corrS <- b1_imp_corr %>% filter(class %in% c('BG', 'S1', 'S2')) %>% as.data.frame() %>%
  left_join(data %>% dplyr::select(Sample, Batch)) %>% relocate(Batch)
b1_imp_corrS$Batch <- as.factor(b1_imp_corrS$Batch)
rownames(b1_imp_corrS) <- b1_imp_corrS$Sample
heat <- b1_imp_corrS[,5:ncol(b1_imp_corrS)]
heat[heat > 0] <- 1
heat <- as.matrix(heat)
rownames(heat) <- b1_imp_corrS$Sample

heatmap_0s <- pheatmap(t(heat), 
                       show_rownames = FALSE,
                       show_colnames = FALSE,
                       annotation_col = b1_imp_corrS %>% dplyr::select(Batch, class))
dev.new()
heatmap_0s

# check % of 0s across variables in data subsets
barplot(colSums(bg_corr[,4:ncol(bg_corr)] == 0)/162*100)
View(table(colSums(bg_corr[,4:ncol(bg_corr)] == 0)/162*100))

bg_sparsity <- as.data.frame(colSums(bg_corr[,4:ncol(bg_corr)] == 0)/162*100) %>%
  mutate(comp = rownames(.))
s1_sparsity <- as.data.frame(colSums(s1_corr[,4:ncol(s1_corr)] == 0)/160*100) %>%
  mutate(comp = rownames(.))
s2_sparsity <- as.data.frame(colSums(s2_corr[,4:ncol(s2_corr)] == 0)/161*100) %>%
  mutate(comp = rownames(.))

# compare with results of wilcoxon test
b_bg <- b_bg %>% left_join(bg_sparsity)
b_s1 <- b_s1 %>% left_join(s1_sparsity)
b_s2 <- b_s2 %>% left_join(s2_sparsity)

# compare with results of Gamma regression
colnames(b_bg_G)[1] <- 'comp'
b_bg_G <- b_bg_G %>% left_join(bg_sparsity)

colnames(b_s1_G)[1] <- 'comp'
b_s1_G <- b_s1_G %>% left_join(s1_sparsity)

colnames(b_s2_G)[1] <- 'comp'
b_s2_G <- b_s2_G %>% left_join(s2_sparsity)

# visualise distribution of blank and sample values in comps classified for removal
# wilcoxon test
dev.new()

b1_imp_L %>% filter(comp %in% intersect(b_s1$comp, b_s2$comp)) %>%
  filter(class %in% c('Blank', 'S1', 'S2')) %>%
  ggplot(aes(x = log(peakArea), fill = class, alpha = 0.6)) + 
  geom_histogram(position = 'identity') +
  facet_wrap(~ comp, scale = 'free')

# gamma regression
intersect(b_s1_G$comp, b_s2_G$comp)

pdf('Hist_noise_features.pdf')
hists <- b1_imp_L %>% filter(comp %in% intersect(b_s1_G$comp, b_s2_G$comp)) %>%
  filter(class %in% c('Blank', 'S1', 'S2')) %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = log(peakArea), fill = class, alpha = 0.6)) + 
                      geom_histogram(position = 'identity') + 
                      ggtitle(.$comp) + theme(legend.position = 'none',
                                              plot.title = element_text(size = 10)))
  
marrangeGrob(hists$plots, nrow = 4, ncol = 4)
dev.off()

#

pdf('Hist_noise_features_BG.pdf')
hists <- b1_imp_L %>% filter(comp %in% b_bg_G$comp) %>%
  filter(class %in% c('Blank', 'BG')) %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = log(peakArea), fill = class, alpha = 0.6)) + 
       geom_histogram(position = 'identity') + 
       ggtitle(.$comp) + theme(legend.position = 'none',
                               plot.title = element_text(size = 10)))

marrangeGrob(hists$plots, nrow = 4, ncol = 4)
dev.off()

#
#
#

# IS normalisation following blank correction
assay_imp_corr_IS <- b1_imp_corr[,4:ncol(b1_imp_corr)]/b1_imp_corr$Internal_Standard
b1_imp_corr_IS <- cbind(b1_imp_corr[,1:3], assay_imp_corr_IS) %>%
  dplyr::select(!Internal_Standard)

b1_imp_corr_ISL <- b1_imp_corr_IS %>% pivot_longer(cols =! c(1:3), names_to = 'comp', values_to = 'peakArea')

# calculate summary statistic of the ES groups
b1_imp_L <- b1_imp_L %>% mutate(peakArea = peakAreaImp)

rsd <- function(dataL){
  #assay <- dataL %>% filter(class == CL)
  rsd <- dataL %>% group_by(comp, class) %>%
    summarise(mean = mean(peakArea, na.rm = TRUE),
              sd = sd(peakArea, na.rm = TRUE),
              rsd = sd/abs(mean)*100,
              mad = mad(peakArea))}

rsd_imp <- rsd(b1_imp_L)
rsd_imp_corr <- rsd(b1_imp_corrL)
rsd_imp_corr_IS <- rsd(b1_imp_corr_ISL)

imp <- rsd_imp %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()
imp_corr <- rsd_imp_corr %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()
imp_corrIS <- rsd_imp_corr_IS %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()

dev.new()
grid.arrange(imp_corr, imp_corrIS, ncol = 2)  
  

# density plot
'%ni%' <- Negate('%in%')

b1_imp_corrL %>% ggplot(aes(x = log2(peakArea), 
                       colour = factor(as.factor(class), levels = c('Blank', 'ES', 'BG', 'S1', 'S2')))) + 
  geom_density() +
  theme_bw() + ggtitle('Density plot of blank corrected log peak areas in pre-covid dataset') +
  xlab('log(Peak area)') +
  theme(plot.title = element_text(hjust = 0.3, size = 11)) +
  labs(colour = 'Class')

ggsave('Density_class_blankCorr_B1.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

# change of median in time
pa_time <- b1_imp_L %>% filter(class %ni% c('ES', 'Blank','BG')) %>%
  group_by(class, Date) %>%
  summarise(median = median(log(peakAreaImp))) %>%
  ggplot(aes(x = Date, y = median, group = class)) + 
  geom_line(aes(colour = class)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ylab('median log( Peak area)') +
  ggtitle('Change in median peak area in time') +
  ylim(6, 12)

bc_pa_time <- b1_imp_corrL %>% filter(class %ni% c('ES', 'Blank','BG')) %>%
  separate(Date_Time, c('Date', 'Time'), sep = ' ') %>%
  group_by(class, Date) %>%
  summarise(median = median(log(peakArea))) %>%
  ggplot(aes(x = Date, y = median, group = class)) +
  geom_line(aes(colour = class)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  ylab('median log( Peak area)') +
  ggtitle('Change in median peak area in time (blank corrected)') +
  ylim(6,12)

bc_pa_time

pa_time1 <- arrangeGrob(pa_time, bc_pa_time, ncol = 2, nrow = 1)
plot(pa_time1)

ggsave('Median_PA_time_blankCorr.tiff', pa_time1, dpi = 300, units = 'mm', width = 210, height = 70)

# pca of blank corrected values
library(mixOmics)
b1_imp_corr1 <- b1_imp_corr %>% left_join(data %>% dplyr::select(Sample, Batch)) %>%
  relocate(Batch)

d <- b1_imp_corr1 %>% filter(class != 'Blank')
pc_blank_corr <- pca(log(d[,5:ncol(d)]), scale = TRUE, center = TRUE, ncomp = 4)
corr_pc23 <- plotIndiv(pc_blank_corr,
          pch = 1,
          cex = 2,
          group = d$Batch,
          legend = TRUE,
          comp = c(1,2),
          title = 'Peak Area (blank corrected)',
          size.title = 10)


corr_pc <- arrangeGrob(corr_pc12$graph, corr_pc23$graph, nrow = 2, ncol = 1)
ggsave('PCA_raw_peak_area_imputed_correctedBatch.tiff', corr_pc, dpi = 300, unit = 'mm', width = 100, height = 150)

