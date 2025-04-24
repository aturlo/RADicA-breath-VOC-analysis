## Missing value pattern evaluation and filtering

# author: Aggie Turlo
# project: RADicA
# date: 16/01/2022

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
library(ggvenn)
library(stringr)

## load data
# load formatted study data
b2_all <- read.csv('RADicA_B2.csv')[,-1]
b1_all <- read.csv('RADicA_B1.csv')[,-1]


# load patient metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

###############################

## EXPLORATION OF MISSING DATA PATTERNS

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
    annotate('text', x = 28, y = 33, label = nrow(fnasF), colour = 'red',
             cex = 2)
  p2
  
}

cnas(b1_all[,-c(1:5)], 20)
cnas(b2_all[,-c(1:5)], 20)

cnas1 <- function(d, cl, cutoff1) {
  cnas(d %>% filter(class == cl) %>% dplyr::select(!1:5), cutoff = cutoff1) + ggtitle(cl)}

cnas_plots_b1 <- lapply(unique(b1_all$class), cnas1, d = b1_all, cutoff1 = 20)
cnas_plots_b2 <- lapply(unique(b1_all$class), cnas1, d = b2_all, cutoff1 = 20)

NAhists_b1 <- grid.arrange(grobs = cnas_plots_b1, nrow = 2, ncol = 3)
NAhists_b2 <- grid.arrange(grobs = cnas_plots_b2, nrow = 2, ncol = 3)

ggsave('Histogram_NA_B1_class.tiff', NAhists_b1, dpi = 300, unit = 'mm', width = 160, height = 100)
ggsave('Histogram_NA_B2_class.tiff', NAhists_b1, dpi = 300, unit = 'mm', width = 160, height = 100)

#

# plot frequency of NAs across observations (breath samples only)
cnas_sample <- function(d, cl, cutoff2) {
  s1 <- cnas(t(d %>% filter(class %in% c(cl)) %>% dplyr::select(!1:5)),
           cutoff = cutoff2) + ggtitle(cl)
  s1
  }

s1_b1 <- cnas_sample(b1_all, 'S1', 35)
s1_b2 <- cnas_sample(b2_all, 'S1', 35)

s2_b1 <- cnas_sample(b1_all, 'S2', 35)
s2_b2 <- cnas_sample(b2_all, 'S2', 35)

b_b1 <- cnas_sample(b1_all, 'Blank', 35)
b_b2 <- cnas_sample(b2_all, 'Blank', 35)

bg_b1 <- cnas_sample(b1_all, 'BG', 35)
bg_b2 <-cnas_sample(b2_all, 'BG', 35)

es_b1 <- cnas_sample(b1_all, 'ES', 35)
es_b2 <- cnas_sample(b2_all, 'ES', 35)

s_b1 <- arrangeGrob(s1_b1, s2_b1, bg_b1, es_b1, b_b1, nrow = 2, ncol = 3)
plot(s_b1)

s_b2 <- arrangeGrob(s1_b2, s2_b2, bg_b2, es_b2, b_b2, nrow = 2, ncol = 3)
plot(s_b2)

ggsave('Histogram_NA_B1_sample.tiff', s_b1, dpi = 300, unit = 'mm', 
       width = 140, height = 100)
ggsave('Histogram_NA_B2_sample.tiff', s_b2, dpi = 300, 
       unit = 'mm', width = 140, height = 100)

#

# identify outlying samples with > 35% NAs
obs_out <- function(d, CL, cutoff3) {
  d1 <- as.data.frame(t(d %>% filter(class %in% c(CL)) %>% dplyr::select(!1:5)))
  fnas <- as.data.frame(colSums(is.na(d1))/nrow(d1)*100)
  colnames(fnas)[1] <- 'RatioNA'
  out <- fnas %>% filter(RatioNA >= cutoff3) %>% mutate(obs = gsub('V', '', rownames(.))) 
  d1 <- d %>% filter(class %in% c(CL)) %>% dplyr::select(Sample)
  d2 <- d1[out$obs,]
  d2
}

s1_out_b1 <- obs_out(b1_all, 'S1', 35)
s1_out_b2 <- obs_out(b2_all, 'S1', 35)

s2_out_b1 <- obs_out(b1_all, 'S2', 35)
s2_out_b2 <- obs_out(b2_all, 'S2', 35)

bg_out_b1 <- obs_out(b1_all, 'BG', 35)
bg_out_b2 <- obs_out(b2_all, 'BG', 35)

es_out_b1 <- obs_out(b1_all, 'ES', 35)
es_out_b2 <- obs_out(b2_all, 'ES', 35)

b_out_b1 <- obs_out(b1_all, 'Blank', 35)
b_out_b2 <- obs_out(b2_all, 'Blank', 35)

# remove matching breath sample replicates
s1_out_b2
s2_out_b2

s1_out_b1
s2_out_b1

s2_out_b1 <- b1_all %>% 
  filter(str_detect(b1_all$Sample, paste(str_sub(s1_out_b1, end = -11), '2', sep = ''))) %>%
  pull(Sample)

# remove the samples with high missing value ratio from the dataset
b1_all <- b1_all %>% filter(Sample %ni% c(s1_out_b1, s2_out_b1,
                                          bg_out_b1, es_out_b1, b_out_b1))

b2_all <- b2_all %>% filter(Sample %ni% c(s1_out_b2, s2_out_b2,
                                          bg_out_b2, es_out_b2, b_out_b2))



#
#
#


###############################

## FEATURE FILTERING
# specify dataset for analysis
data <- b1_all
dat <- 'B1'

#
#
#

# change dataset to long format
b_allL <- data %>% pivot_longer(cols =! c(class, Sample, Date, Time, Batch),
                                             names_to = 'comp', values_to = 'peakArea')

# filter features with > 20% NAs
class_na <- b_allL %>% 
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis), by = c('Sample' = 'Sample_ID')) %>%
  mutate(class = ifelse(class == 'S1' & Diagnosis == 'Asthma', 'S1A', 
                        ifelse(class == 'S1' & Diagnosis == 'Not Asthma', 'S1NA',
                               ifelse(class == 'S2' & is.na(Diagnosis) == TRUE, 'S2',
                                      ifelse(class == 'S2' & Diagnosis == 'Not Asthma', 'S2NA', 
                                             ifelse(class == 'S2' & Diagnosis == 'Asthma', 'S2A', class)))))) %>%
  filter(class != 'S2') %>%  
  group_by(class, comp) %>% 
  summarise(nona = sum(is.na(peakArea) == TRUE)/length(peakArea))

# dichotomise missingness based on 20% threshold
class_na <- class_na %>% mutate(nona = ifelse(nona < 0.2, 1, 0))

#

# compare misingness in features across sample types
class_naF_wide <- class_na %>% pivot_wider(names_from = comp, values_from = nona)
assay_naF <- class_naF_wide[,-1]

# no of classes where VOC met the missingness threshold
sums <- as.data.frame(colSums(assay_naF))
colnames(sums) <- 'sum'
complete_var <- sums %>% filter(sum == 7) # VOC missing in < 20% samples in all classes 
reject_var <- sums %>% filter(sum == 0) # VOC missing in > 20% samples in all classes

# filter VOCs with different missingness level in different sample classes
assay_naFf <- assay_naF %>% dplyr::select(!rownames(complete_var)) %>%
  dplyr::select(!rownames(reject_var))
rownames(assay_naFf) <- class_naF_wide$class

# presented in a heatmap
tiff(paste(filename = 'Heatmap_NA_', dat, '_class.tiff', sep =''), width = 200, height = 170, unit = 'mm', res = 300)
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

# no of VOCs above threshold in each of the breath sample classes
assay_naFf_s <- assay_naFf %>% filter(rownames(.) %ni% c('BG', 'Blank', 'ES'))
sums_s <- as.data.frame(colSums(assay_naFf_s))
colnames(sums_s) <- 'sum'
retain_var <- sums_s %>% filter(sum > 0)

# identify VOCs missing differentially between disease categories
retain_var1 <- retain_var %>% filter(sum == 2)
View(t(assay_naFf_s[,rownames(retain_var1)]))

# retain features with < 20% NAs in S1 and/or S2 samples (in either diagnosis group)
b1_all_f <- b1_all %>% dplyr::select(c(1:5, rownames(retain_var), rownames(complete_var)))

#
#
#

# change starting input!
b2_all_f <- b2_all %>% dplyr::select(c(1:5, rownames(retain_var), rownames(complete_var)))

#

# compare retained VOCs between datasets
retain_var2 <- intersect(colnames(b1_all_f), colnames(b2_all_f))

# Venn diagram
B1 <- colnames(b1_all_f) 
B1 <- B1[B1 %ni% c('class', 'Sample', 'Date', 'Time', 'Batch', 'Internal_Standard')]
B2 <- colnames(b2_all_f)
B2 <- B2[B2 %ni% c('class', 'Sample', 'Date', 'Time', 'Batch', 'Internal_Standard')]


ggvenn(list_venn <- list(Data_1 = B1,
                         Data_2 = B2),
set_name_size = 4,
       text_size = 4,
       auto_scale = TRUE,
       show_percentage = FALSE) +
  annotate("text", x = -1, y = 0, label = '3', size = 4)

ggsave('Venn_NAfilter_B1B2.tiff', units = 'mm', dpi = 300, width = 55, height = 65)

# retain only VOCs overlapping between the datasets
b1_all_f <- b1_all_f[,retain_var2]
b2_all_f <- b2_all_f[,retain_var2]

write.csv(b1_all_f, 'RADicA_B1_NAfiltered.csv') # NA filtered dataset
write.csv(b2_all_f, 'RADicA_B2_NAfiltered.csv')

#
#
#


