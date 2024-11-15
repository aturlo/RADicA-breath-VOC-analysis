## Missing value imputation

# author: Aggie Turlo
# project: RADicA
# date: 17/07/2024

#####################

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

# load filtered study data with missing values
b_all_f <- read.csv('RADicA_B2_NAfiltered.csv')[,-1]
b_all_f <- read.csv('RADicA_Train_NAfiltered.csv')[,-1]
rownames(b_all_f) <- b_all_f$Sample

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

b_all_fL <- b_all_f %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                       names_to = 'comp', 
                                       values_to = 'peakArea')

table(duplicated(b_all_f$Sample))

dupes <- c('200320_RaDICA_RAD094_B1_434919_1', '200320_RaDICA_RAD074_B1_609366_1', '200417_RaDICA_RAD114_S2_367150_1')

'%ni%' <- Negate('%in%')

b_all_fL <- b_all_fL %>% filter(Sample %ni% dupes)

#

# blank filtered
b1_all_f_corr <- read.csv('RADicA_B1_NAfiltered_blank_filtered.csv')
b1_all_f_corr <- b1_all_f_corr[,-1]
rownames(b1_all_f_corr) <- b1_all_f_corr$Sample

b1_all_f_corrL <- b1_all_f_corr %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                                 names_to = 'comp', 
                                                 values_to = 'peakArea')
'%ni%' <- Negate('%in%')
b1_all_f_corrL <- b1_all_f_corrL %>% filter(Sample %ni% dupes)

#

# IMPUTATION
# create assay dataset for each class
df <- b_all_f # wide format of a data frame

assay <- function(x){
  df %>% filter(class == x) %>% dplyr::select(6:ncol(.))
}

assay_blank <- assay('Blank')
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
    left_join(b_all_f %>% dplyr::select(Sample, Date, Time, Batch, class))

}

blank_impL <- imp_long(assay_blank_imp)
es_impL <- imp_long(assay_es_imp)
bg_impL <- imp_long(assay_bg_imp)
s1_impL <- imp_long(assay_s1_imp)
s2_impL <- imp_long(assay_s2_imp)

b_imp_L <- rbind(es_impL, blank_impL, bg_impL, s1_impL, s2_impL) 

b_imp_L <- b_imp_L %>% filter(Sample %ni% dupes)
b_imp <- b_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea)

is <- b_all_f %>% dplyr::select(Sample, Internal_Standard) %>%
  left_join(b_imp %>% dplyr::select(Sample, Internal_Standard),
            by = 'Sample')

plot(is$Internal_Standard.x, is$Internal_Standard.y)

write.csv(b_imp_L, 'RADicA_B2_NAfiltered_imputed_long.csv')
write.csv(b_imp, 'RADicA_B2_NAfiltered_imputed.csv')

write.csv(b_imp_L, 'RADicA_Train_NAfiltered_blank_filtered_imputed_long.csv')
write.csv(b_imp, 'RADicA_Train_NAfiltered_blank_filtered_imputed.csv')

b_imp <- read.csv('RADicA_B2_NAfiltered_blank_filtered_imputed.csv')[,-1]
b_imp <- read.csv('RADicA_Train_NAfiltered_blank_filtered_imputed.csv')[,-1]

b_imp_L <- b_imp %>% pivot_longer(cols = c(6:ncol(b_imp)), names_to = 'comp', values_to = 'peakArea')

#
#
#

# COMPARE IMPUTED AND OBSERVED DATA
# density plots
dens_plot <-
  b_imp_L %>% 
  filter(is.na(class) == FALSE) %>%
  ggplot() + 
  geom_density(aes(x = log(peakArea), colour = 'Imputed values')) +
  geom_density(data = b_all_fL %>% drop_na(),
                 aes(x = log(peakArea), colour = 'Observed values')) +
  theme_bw() +
  xlab('log( Peak area)') +
  facet_wrap(~ class, scale = 'free')

dens_plot

ggsave('Density_imputed_Train.tiff', dens_plot, dpi = 300, units = 'mm', width = 240, height = 120)

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
  geom_line(aes(y = medianImp, group = 'Imputed values', colour = 'Imputed values')) +
  geom_line(data = b_all_fL %>% 
              drop_na(class) %>%
              group_by(class, Date) %>%
              summarise(medianObs = median(log(peakArea), na.rm = TRUE)),
            aes(y = medianObs, group = 'Observed values', colour = 'Observed values')) +
  theme_bw() +
  facet_wrap(~ class, scale = 'free') +
  ylab('median log( Peak area)') 

dev.new()
time_plot

ggsave('Time_imputed_Train.tiff', time_plot, dpi = 300, units = 'mm', width = 240, height = 120)

#

# PCA complete vs imputed variables
# create data class subsets with complete variables only
x <- b_all_f
x1 <- b_all_fL


sumint <- x1 %>% 
  group_by(comp, class) %>%
  summarise(nonas = sum(is.na(peakArea) == TRUE))


complete_set <- function(CL){ 
  df <- x %>% filter(class == CL)
  rownames(df) <- df$Sample
  comp_vars <- sumint %>% filter(class == CL) %>% filter(nonas == 0)
  comp_df <- df[,comp_vars$comp]
}

std <- c('Acetone', 'Isoprene', 'Benzene', 'X3_Pentanone', 'X1.4_dioxane',
          'Pyridine', 'Toluene', 'Octane', 'P_Xylene', 'Nonane', 'Benzaldehyde',
          'X1_Heptanol', 'Decane', 'X3_Carene', 'Limonene', 'Undecane', 'Nonanal',
          'X1.2.3.4_tetrahydronaphthalene', 'Dodecane', 'X1_methylindole', 'Tridecane',
          'Pentadecane')

b_blank_f_c <- complete_set('Blank')
# for ES use spike-in data only
b_es_f_c <- b_all %>% filter(class == 'ES') %>% subset(select = colnames(.) %in% std) 
b_bg_f_c <- complete_set('BG')
b_s1_f_c <- complete_set('S1')
b_s2_f_c <- complete_set('S2')

# create data annotation subset for the experimental variables (group, class)
b_es_f <- b_all_f %>% filter(class == 'ES') 
b_blank_f <- b_all_f %>% filter(class == 'Blank') 
b_bg_f <- b_all_f %>% filter(class == 'BG') 
b_s1_f <- b_all_f %>% filter(class == 'S1') 
b_s2_f <- b_all_f %>% filter(class == 'S2') 

#
library(mixOmics)

plot_pca <- function(x, y, z){
  pca <- mixOmics::pca(log(x), scale = TRUE, center = TRUE)
  plot_pca <- plotIndiv(pca,
                        pch = 1,
                        cex = 1.5,
                        group = y$Batch,
                        #legend = TRUE,
                        title = z,
                        size.title = 10)#,
  #legend.title = 'Batch')
  
}


pc_blank_c <- plot_pca(b_blank_f_c, b_blank_f, 'Blank')
pc_blank_i <- plot_pca(assay_blank_imp, b_blank_f, 'Blank')

pc_es_c <- plot_pca(b_es_f_c, b_es_f, 'ES') # row 71 in B2 high number of NAs
pc_es_i <- plot_pca(assay_es_imp[-71,], b_es_f[-71,], 'ES')

pc_bg_c <- plot_pca(b_bg_f_c, b_bg_f, 'BG')
pc_bg_i <- plot_pca(assay_bg_imp, b_bg_f, 'BG')

pc_s1_c <- plot_pca(b_s1_f_c, b_s1_f, 'S1')
pc_s1_i <- plot_pca(assay_s1_imp, b_s1_f, 'S1')

pc_s2_c <- plot_pca(b_s2_f_c, b_s2_f, 'S2')
pc_s2_i <- plot_pca(assay_s2_imp, b_s2_f, 'S2')

pc_plots <- arrangeGrob(pc_blank_c$graph, pc_es_c$graph, pc_bg_c$graph, pc_s1_c$graph, pc_s2_c$graph,
                        pc_blank_i$graph, pc_es_i$graph, pc_bg_i$graph, pc_s1_i$graph, pc_s2_i$graph,
                        ncol = 5, nrow = 2)

dev.new()
plot(pc_plots)

ggsave('PCA_complete_vs_imputed_Train.tiff', pc_plots, dpi = 300, unit = 'mm', width = 310, height = 120)

#
#
#

# effect of missingness on technical replicates in breath samples
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
  
s_sum <- s_sum %>% left_join(s_sum %>% group_by(Batch) %>% summarise(BatchSum = sum(n))) %>%
  mutate(fraqBatch = n/BatchSum *100)

s_sum %>%
  ggplot(aes(x = as.factor(Batch), y = fraqBatch, fill = techNA)) + 
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_brewer(palette = 'Set2') +
  ylab('% of all observations') +
  xlab('Batch') +
  ggtitle('Missingness in breath sample replicates') +
  theme(plot.title = element_text(hjust = 0.5, size = 12))

ggsave('Missingness_S1S2_barchart_Train.tiff', dpi = 300, unit = 'mm', width = 100, height = 80)

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

pdf('Missing_pattern_sample_Train.pdf')
marrangeGrob(violin_plots$plots, nrow = 4, ncol = 4)
dev.off()

# following imputation
s_tech_imp <- b_imp_L %>% filter(class %in% c('S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Analysis_date, CoreVisit, RAD_ID),
            by = c('Sample' = 'Sample_ID')) %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  dplyr::select(!Time) %>%
  pivot_wider(id_cols = c(Sample, comp), names_from = class, values_from = peakArea) %>%
  left_join(s_tech %>% dplyr::select(Sample, comp, Batch, techNA))

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

pdf('Missing_pattern_sample_imp_B2.pdf')
marrangeGrob(violin_plots_imp$plots, nrow = 4, ncol = 4)
dev.off()

# compare distributions of breath data with and without imputation
s_tech1 <- s_tech %>% filter(techNA != 'Non-missing')
s_tech_imp1 <- s_tech_imp %>% filter(techNA != 'Non-missing')

View(s_tech1 %>% group_by(Sample) %>% summarise(n = n()))

s_tech_comp <- data.frame(comp = s_tech1$comp, 
                          Sample = s_tech1$Sample,
                          S1NA = s_tech1$S1,
                          S2NA = s_tech1$S2,
                          S1imp = s_tech_imp1$S1,
                          S2imp = s_tech_imp1$S2) %>%
  pivot_longer(cols = c(S1NA, S2NA, S1imp, S2imp), names_to = 'class', values_to = 'peakArea') %>%
  mutate(imputed = ifelse(grepl('imp', class), 'Yes', 'No')) %>%
  mutate(class = ifelse(grepl('S1', class), 'S1', 'S2'))

id <- '190524_RaDICA_RAD002'

comp_reg <- bind_rows(lapply(unique(s_tech_comp$Sample), function(id) {
  subset <- s_tech_comp %>%
    filter(Sample == id) %>%
    mutate(peakArea = peakArea/sd(peakArea, na.rm = TRUE)) %>%
    drop_na()
  
  model <- glm(peakArea ~ imputed,
                 data = subset,
                 family = 'Gamma'(link = 'identity'))
  
  coefs <- coef(summary(model)) %>% as.data.frame() %>%
    mutate(Sample = rep(id,2),
           cov = rownames(.)) 
}))

colnames(comp_reg)[4] <- 'p.value'

comp_reg_imp <- comp_reg %>% filter(cov == 'imputedYes') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

View(comp_reg_imp %>% filter(p.value < 0.05))

write.csv(comp_reg, 'Gamma_regression_breath_imputedVSobserved_SD_Train.csv')

## quantile regression
library(imputeLCMD)

assay2 <- b2[,-c(1:6)]
Lassay2 <- log(assay2)
assay2_imp <- impute.QRILC(Lassay2)
