## Missing value pattern evaluation and filtering

# author: Aggie Turlo
# project: RADicA
# date: 16/07/2024

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggfortify)
library(RColorBrewer)
library(scales)
library(gridExtra)

# load formatted study data
b_all <- read.csv('Train_dataset.csv')[,-1]

b_all <- b_all %>% filter(CoreVisit %ni% c('CV3', 'CV4'))


meta <- read.csv('RaDICA Breath sample metadata.csv')

b_all <- read.csv('RADicA_B2.csv')[,-1]

b_allL <- b_all %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                             names_to = 'comp', values_to = 'peakArea')

# EXPLORATION OF MISSING DATA PATTERNS
# plot frequency of NAs across features in each class
cnas <- function(d, cutoff){ # cutoff = % of NA per VOC
  fnas <- as.data.frame(colSums(is.na(d))/nrow(d)*100)
  colnames(fnas)[1] <- 'RatioNA'
  rownames(fnas) <- colnames(d)
  fnasF <- fnas %>% filter(RatioNA < cutoff)
  p2 <- fnas %>% ggplot(aes(x = RatioNA)) + 
    geom_histogram() +
    theme_bw(base_size = 12) +
    ylab('Sample number') +
    xlab('% Missing') +
    xlim(NA,100) +
    geom_vline(xintercept = cutoff, colour = 'red') +
    annotate('text', x = 12.5, y = 55, label = nrow(fnasF), colour = 'red',
             cex = 3)
  p2
  
}

cnas(b_all[,-c(1:5)], 20)

cnas1 <- function(cl) {
  cnas(b_all %>% filter(class == cl) %>% dplyr::select(!1:5), 20) +
    ggtitle(cl)}

cnas1('Blank')

cnas_plots <- lapply(unique(b_all$class), cnas1)
NAhists <- grid.arrange(grobs = cnas_plots, nrow = 2, ncol = 3)

ggsave('Histogram_NA_Train_class.tiff', NAhists, dpi = 300, unit = 'mm', width = 160, height = 100)

#
#
#

# plot frequency of NAs across samples (S1 and S2 only)
s1 <- cnas(t(b_all %>% filter(class %in% c('S1')) %>% dplyr::select(!1:5)),
           cutoff = 20) + 
  ggtitle('S1')

s2 <- cnas(t(b_all %>% filter(class %in% c('S2')) %>% dplyr::select(!1:5)),
           cutoff = 20) + 
  ggtitle('S2')

s <- grid.arrange(s1, s2, ncol = 2)

ggsave('Histogram_NA_Train_sample.tiff', s, dpi = 300, unit = 'mm', width = 105, height = 50)

table(b_all$class)


# identify outlying samples with > 50% NAs
obs_out <- function(CL) {
  d <- as.data.frame(t(b_all %>% filter(class %in% c(CL)) %>% dplyr::select(!1:5)))
  fnas <- as.data.frame(colSums(is.na(d))/nrow(d)*100)
  colnames(fnas)[1] <- 'RatioNA'
  out <- fnas %>% filter(RatioNA > 40) %>% mutate(obs = gsub('V', '', rownames(.))) 
  d1 <- b_all %>% filter(class %in% c(CL)) %>% dplyr::select(Sample)
  d2 <- d1[out$obs,]
  d2
  }

s1_out <- obs_out('S1')
s2_out <- obs_out('S2')

# remove the samples with high missing value ratio from the dataset
'%ni%' <- Negate('%in%')
b_all <- b_all %>% filter(Sample %ni% c(s1_out, s2_out))

#
#
#

# EFFECT OF EXPERIMENTAL FACTORS ON MISSINGNESS
# heatmap of NAs according to experimental condition and sample type
b_all <- b_all %>% dplyr::select(!c(RAD_ID, CoreVisit, Diagnosis))

data_heat <- b_all[,6:ncol(b_all)]
data_heat[is.na(data_heat)] <- 0
data_heat[data_heat > 0] <- 1
rownames(data_heat) <- b_all$Sample

cols1 <- (hue_pal()(5))
cols2 <- (hue_pal()(n_distinct(b_all$Batch))) 
names(cols1) <- unique(b_all$class)
names(cols2) <- unique(b_all$Batch)
annotCol <- list(class = cols1,
                 Batch = cols2)

b_all$Batch <- as.factor(b_all$Batch)
rownames(b_all) <- b_all$Sample

dev.new()
tiff('Heatmap_NA_Train_class_batch.tiff', res = 300, width = 250, height = 200, unit = 'mm')
heatmap <- pheatmap(as.matrix(t(data_heat)), 
                    annotation_col = b_all %>% 
                      dplyr::select(Batch, class), 
                    show_colnames = FALSE, cluster_cols =  T, 
                    show_rownames = FALSE, cluster_rows = FALSE,
                    fontsize = 10, main = 'Missing values across samples',
                    color = c('Tomato1', 'Black'),
                    treehight_col = 100,
                    annotation_colors = annotCol,
                    treeheight_col = 70) # use for analysis date mapping
heatmap
dev.off()

#

# in breath samples only, including diagnosis and treatment group
b_s <- b_all %>% filter(class %in% c('S1','S2')) %>% dplyr::select(!1:5)
data_heat1 <- b_s
data_heat1[is.na(data_heat1)] <- 0
data_heat1[data_heat1 > 0] <- 1
rownames(data_heat) <- b_s$Sample

b_all_an <- b_all %>% left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
                                  by = c('Sample' = 'Sample_ID')) %>% 
  mutate(Treatment = ifelse(CoreVisit %in% c('CV1', 'CV2'),
                                                                                                'Before', 'After'))
rownames(b_all_an) <- b_all_an$Sample 

cols3 <- (hue_pal()(n_distinct(b_all$Batch)))
cols <- (hue_pal()(5)) 
cols4 <- cols[1:2]
cols5 <- cols[4:5]
names(cols3) <- na.omit(unique(b_all_an$Batch))
names(cols4) <- na.omit(unique(b_all_an$Diagnosis))
names(cols5) <- unique(b_all_an$Treatment)
annotCol1 <- list(Batch = cols3,
                  Diagnosis = cols4,
                  Treatment = cols5)


dev.new()
tiff('Heatmap_NAs_Train_Breath.tiff', res = 300, width = 250, height = 200, unit = 'mm')
heatmap1 <- pheatmap(as.matrix(t(data_heat1)), 
                     annotation_col = b_all_an %>% 
                       filter(class %in% c('S1', 'S2')) %>%
                       dplyr::select(Batch, Diagnosis, Treatment), 
                     show_colnames = FALSE, cluster_cols =  T, 
                     show_rownames = FALSE, cluster_rows = FALSE,
                     fontsize = 10, main = 'Missing values across breath samples',
                     color = c('Tomato1', 'Black'),
                     annotation_colors = annotCol1,
                     treeheight_col = 70)

heatmap1
dev.off()

# count asthma cases between batches
b_all_an %>% group_by(Diagnosis, Batch) %>% summarise(n = n())

b_cases <- function(B){
  b <- b_all_an %>% filter(Batch == B) %>% filter(class == 'S2')
  b %>% group_by(Diagnosis, CoreVisit) %>% summarise(n = n())}

b_cases(6)

#
#

# FEATURE FILTERING
# filter features with > 20% NAs
class_na <- b_allL %>% 
  left_join(b_all_an %>% dplyr::select(Sample, Diagnosis, Treatment)) %>%
  mutate(class = ifelse(class == 'S1' & Diagnosis == 'Asthma', 'S1A', 
                        ifelse(class == 'S1' & Diagnosis == 'Not Asthma', 'S1NA',
                               ifelse(class == 'S2' & is.na(Diagnosis) == TRUE, 'S2',
                                      ifelse(class == 'S2' & Diagnosis == 'Not Asthma', 'S2NA', 
                                             ifelse(class == 'S2' & Diagnosis == 'Asthma', 'S2A', class)))))) %>%
  filter(class != 'S2') %>%  
  group_by(class, comp) %>% 
  summarise(nona = sum(is.na(peakArea) == TRUE)/length(peakArea))

class_na <- class_na %>% mutate(nona = ifelse(nona < 0.2, 1, 0))

# ALTERNATIVE FILTERING
# consider both breath measurement together depending on biological variable
class_na1 <- b_allL %>% 
  left_join(b1_all_an %>% dplyr::select(Sample, Diagnosis, Treatment)) %>%
  mutate(class = ifelse(class %in% c('S1', 'S2') & Diagnosis == 'Asthma', 'SAsth', 
                        ifelse(class %in% c('S1', 'S2') & Diagnosis == 'Not Asthma', 'SNAsth',
                               ifelse(class %in% c('S1','S2') & is.na(Diagnosis) == TRUE, 'S', class)))) %>%
  filter(class != 'S2') %>%  
  group_by(class, comp) %>% 
  summarise(nona = sum(is.na(peakArea) == TRUE)/length(peakArea))

class_na1 <- class_na1 %>% mutate(nona = ifelse(nona < 0.2, 1, 0))

# compare features retained across sample types
class_naF_wide <- class_na %>% pivot_wider(names_from = comp, values_from = nona)
assay_naF <- class_naF_wide[,-1]

# no of VOCs above threshold in all classes
sums <- as.data.frame(colSums(assay_naF))
colnames(sums) <- 'sum'
complete_var <- sums %>% filter(sum == 7) #7 if breath sub-classes included
reject_var <- sums %>% filter(sum == 0)

assay_naFf <- assay_naF %>% dplyr::select(!rownames(complete_var)) %>%
  dplyr::select(!rownames(reject_var))
rownames(assay_naFf) <- class_naF_wide$class

# heatmap of VOCs that passed 20/80 threshold in at least one class
tiff(filename = 'Heatmap_NA_Train_class_v2.tiff', width = 200, height = 170, unit = 'mm', res = 300)
dev.new()
heatmap <- pheatmap(t(as.matrix(assay_naFf)), 
                    show_rownames = TRUE, 
                    cluster_cols =  T, 
                    cluster_rows = T, 
                    show_colnames = TRUE, 
                    fontsize_row = 4,
                    main = 'Compounds with < 20% missing values across samples (black)',
                    legend  = F, color = c('Tomato1', 'Black'),
                    treeheight_col = 10)
heatmap
dev.off()

# no of VOCs above threshold in breath sample classes
assay_naFf_s <- assay_naFf %>% filter(rownames(.) %ni% c('BG', 'Blank', 'ES'))
sums_s <- as.data.frame(colSums(assay_naFf_s))
colnames(sums_s) <- 'sum'
retain_var <- sums_s %>% filter(sum > 0)

# retain features with < 20% NAs in S1 and/or S2 samples (in either diagnosis group)
b_all_f <- b_all %>% dplyr::select(c(1:5, rownames(retain_var), rownames(complete_var)))

write.csv(b_all_f, 'RADicA_B2_NAfiltered.csv') # NA filtered dataset
write.csv(b_all_f, 'RADicA_B1_NAfiltered_v2.csv')
write.csv(b_all_f, 'RADicA_Train_NAfiltered.csv')

b_all_fL <- b_all_f %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                       names_to = 'comp', values_to = 'peakArea')

# 267 features retained based on 20/80 threshold in breath samples

# compare retained VOCs between batches
b1 <- read.csv('RADicA_B1_NAfiltered_v2.csv')[,-1]
b2 <- read.csv('RADicA_B2_NAfiltered.csv')[,-1]

intersect(colnames(b1), colnames(b2))
setdiff(colnames(b1), colnames(b2))
setdiff(colnames(b2), colnames(b1))

library(ggvenn)

B1 <- colnames(b1) 
B1 <- B1[B1 %ni% c('class', 'Sample', 'Date', 'Time', 'Batch', 'Internal_Standard')]
B2 <- colnames(b2)
B2 <- B2[B2 %ni% c('class', 'Sample', 'Date', 'Time', 'Batch', 'Internal_Standard')]

ggvenn(list_venn <- list(B1 = B1,
                         B2 = B2),
set_name_size = 3,
       text_size = 4,
       auto_scale = TRUE,
       show_percentage = FALSE)

ggsave('Venn_NAfilter_B1B2.tiff', units = 'mm', dpi = 300, width = 50, height = 50)

#
#
#


# NAs according to peak intensity
sumint <- b_all_fL %>% 
  group_by(comp, class) %>%
  summarise(med = median(log(peakArea), na.rm = TRUE),
            nonas = sum(is.na(peakArea) == TRUE))

pdf('NAvsMedianIntensity_Train.pdf', width = 7, height = 4.5)

medVSna <- sumint %>% group_by(class) %>%
  arrange(med) %>%
  do(plots = ggplot(data = ., aes(x = med, y = nonas)) +
       geom_point(alpha = 0.6) + ggtitle(.$class) + theme_bw() +
       theme(axis.title.x = element_blank(),
             axis.title.y = element_blank()))

marrangeGrob(medVSna$plots, nrow = 2, ncol = 3)

dev.off()

# number of complete variables across groups
sumint %>% group_by(class) %>% summarise(complete_var = sum(nonas == 0))

#

#  NA patterns across technical replicates
techNA <- b_all_fL %>% filter(class %in% c('S1', 'S2', 'BG')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, RAD_ID, CoreVisit), 
            by = c('Sample' = 'Sample_ID')) %>%
  pivot_wider(id_cols = c(RAD_ID, comp, CoreVisit, Batch, Date), 
              names_from = class, values_from = peakArea,
              values_fn = mean)

techNA$techNA <- ifelse(is.na(techNA$S1) | is.na(techNA$S2), 'Missing', 'Non-missing')

techNAL <- techNA %>% dplyr::select(RAD_ID, comp, S1, S2, techNA, CoreVisit) %>%
  pivot_longer(cols = c(S1, S2), names_to = 'Sample_type', values_to = 'peakArea')

# per sample, across compounds
pdf('Missing_patterns_sample_Train.pdf')
pnas <- techNAL %>% group_by(RAD_ID) %>%
  do(plots = ggplot(data =.) + 
       (aes(x = techNA, y = log(peakArea))) + 
       geom_violin(aes(colour = Sample_type), position = position_dodge(0.9)) + 
       geom_boxplot(aes(colour = Sample_type), position = position_dodge(0.9), width = 0.1) + 
       theme_bw(base_size = 8) +
       theme(legend.position = 'none') +
       scale_color_brewer(palette = 'Dark2') +
       ggtitle(.$RAD_ID))

marrangeGrob(pnas$plots, nrow = 5, ncol = 4)
dev.off()

# per compound, across samples
pdf('Missing_patterns_compound.pdf')
pnac <- techNAL %>% group_by(comp) %>%
  do(plots = ggplot(data =.) + 
       aes(x = techNA, y = log(peakArea), colour = techNA) +
       geom_violin(position = position_dodge(0.9)) + 
       geom_boxplot(position = position_dodge(0.9), width = 0.1) +
       theme_bw(base_size = 8) +
       theme(legend.position = 'none',
             title = element_text(size = 6)) +
       scale_color_brewer(palette = 'Dark2') +
       ggtitle(.$comp))

marrangeGrob(pnac$plots, nrow = 5, ncol = 4)
dev.off()

#
#
#

# relationship between technical replicates
ps <- techNA %>% group_by(comp) %>%
  do(plots = ggplot(data = .) + (aes(y = log(S1), x = log(S2))) + geom_point() +
       ggtitle(.$comp) + coord_fixed() +
       theme_bw(base_size = 6) + geom_abline(intercept = 0, slope = 1, colour = 'red'))

pdf('Radica_S1vsS2_Train.pdf')
marrangeGrob(grobs = ps$plots, nrow = 4, ncol = 4)
dev.off()

