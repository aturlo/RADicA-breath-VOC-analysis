## RADicA study analysis v.02

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggfortify)

# Dataset formatting for all analyses

# load peak area file
data <- read.csv('240508_Radica_VOC_peak_areas_v5.1.csv', check.names = FALSE)
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta1 <- read.csv('RaDICA Breath sample metadata.csv')

# replace spaces in compound names with underscores
colnames(data)[192] <- 'Beta phellandrene'
colnames(data) <- gsub(' ', '_', colnames(data))
colnames(data) <- gsub('-', '_', colnames(data))

table(duplicated(data))

# filter Batch 1 data only
data <- data %>% separate(Acq_Date_Time, c('Date', 'Time'), sep = ' ')
data$Date <- as.Date(data$Date, format = "%d/%m/%Y")

b1_all <- data %>% filter(Date < '2020-10-01')
b2_all <- data %>% filter(Date > '2020-10-01')

# compare data distributions between batches
ggplot() + 
  geom_density(data = b1_pqn1 %>% 
                 filter(class %in% c('S1', 'S2', 'BG')) %>%
                          pivot_longer(cols =! c(1:4), 
                                       names_to = 'comp',
                                       values_to = 'peakArea'),
                       aes(x = log(peakArea), colour = class),
               size = 1) +
  geom_density(data = b_pqn1 %>% 
                 filter(class %in% c('S1', 'S2', 'BG')) %>%
                 pivot_longer(cols =! c(1:4), 
                              names_to = 'comp',
                              values_to = 'peakArea'),
               aes(x = log(peakArea), colour = class),
               linetype = 'dashed', size = 1) +
  theme_bw() +
  ggtitle('Density plot of raw peak areas in B1 (solid) and B2 (dashed)') +
  xlim(0,NA)

ggsave('Density_plot_B1vsB2_PQN.tiff', dpi = 300, unit = 'mm', width = 150, height = 100)


# assign sample classes
b_all <- data # change to b1_all if needed

b2_all <- b2_all %>% 
  mutate(class = ifelse(grepl('tracker|Tr', Sample),'Blank',
                              ifelse(grepl('QC', Sample), 'ES',
                                     ifelse(grepl('B1', Sample), 'BG',
                                            ifelse(grepl('S1', Sample), 'S1', 'S2'))))) %>%
  relocate(class)

b_all[b_all == 0] <- NA

write.csv(b_all, 'RADicA_B2.csv')

# split randomly into training and validation datasets
b_all <- b_all %>% left_join(meta %>% dplyr::select(Sample_ID, RAD_ID, Diagnosis, CoreVisit),
                           by = c('Sample' = 'Sample_ID')) 

ids <- b_all %>%
  filter(class %in% c('S1', 'S2')) %>%
  filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  dplyr::select(RAD_ID, Diagnosis, CoreVisit) %>%
  distinct()

View(ids %>% dplyr::group_by(Diagnosis, RAD_ID) %>% dplyr::summarise(n = n()) %>%
       filter(Diagnosis == 'Not Asthma'))

# keep the same ratio of asthma/not asthma patients
library(caret)

ids1 <- ids %>%  dplyr::select(!CoreVisit) %>% distinct()

set.seed(1410)
train_ids <- createDataPartition(ids1$Diagnosis,
                                 p = .8,
                                 list = FALSE,
                                 times = 1)

train_ids <- ids1[train_ids,]
train_all <- ids %>% filter(RAD_ID %in% train_ids$RAD_ID)

View(train_all %>% dplyr::group_by(Diagnosis, RAD_ID) %>% dplyr::summarise(n = n()) %>%
       filter(Diagnosis == 'Asthma'))

train <- b_all %>% filter(RAD_ID %in% train_all$RAD_ID) %>%
  rbind(b_all %>% filter(class %in% c('Blank', 'ES')))

test <- b_all %>% filter(RAD_ID %ni% train_all$RAD_ID) %>%
  rbind(b_all %>% filter(class %in% c('Blank', 'ES')))

write.csv(train, 'Train_dataset.csv')
write.csv(test, 'Test_dataset.csv')

#############################
#############################

#
#
#

# link to metadata
b_all1 <- b_all %>% left_join(meta %>% dplyr::select(Sample_ID, RAD_ID, Analysis_date),
                                  by = c('Sample' = 'Sample_ID')) %>%
  relocate(RAD_ID, Analysis_date)

b_all1$Analysis_date <- as.Date(b_all1$Analysis_date, format = "%d/%m/%Y")

# assign tracker samples (blanks) to study samples
meta1 <- meta1 %>% dplyr::select(Patient.ID,
                                 Sample.Tube.2...B1....Control...Glass.head, 
                                 Sample.Tube.3..S1....A2.Lower.airway,
                                 Sample.Tube.4..S2....B2.Lower.airway,
                                 Tracker.tube,
                                 Sample.analysis_date)

colnames(meta1)[1] <- 'RAD_ID'
colnames(meta1)[2] <- 'BG'
colnames(meta1)[3] <- 'S1'
colnames(meta1)[4] <- 'S2'
colnames(meta1)[5] <- 'Blank'
colnames(meta1)[6] <- 'Date'

# filter RAD_ID only in Batch 1/2
b_meta1 <- meta1 %>% filter(RAD_ID %in% b_all1$RAD_ID)

# match with peak area data
library(stringr)

b_meta1 <- b_meta1 %>% mutate(Date = as.Date(Date, format = '%d/%m/%Y')) %>%
  rename(Analysis_date = Date)

b_all1 <- b_all1 %>% mutate(Analysis_date = as.Date(Analysis_date, format = '%d/%m/%Y')) 

b_all1 <- b_all1 %>% left_join(b_meta1 %>% dplyr::select(RAD_ID, Analysis_date, Blank)) %>%
  relocate(Blank)

#
#
#

# long data layout
b_allL <- b_all %>% 
  pivot_longer(cols = c(9:ncol(b_all)), names_to = 'comp', values_to = 'peakArea') %>%
    mutate(logPeakArea = log(peakArea))
  
# density plots according to sample class
b_allL %>% ggplot(aes(x = logPeakArea, 
                       colour = factor(as.factor(class), levels = c('Blank', 'ES', 'BG', 'S1', 'S2')))) + 
                         geom_density() +
  theme_bw() + ggtitle('Density plot of log-transformed peak areas in post-covid dataset') +
  xlab('log(Peak area)') +
  theme(plot.title = element_text(hjust = 0.5, size = 11)) +
  labs(colour = 'Class')

ggsave('Density_class_B2.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

# descriptive statistics according to sample class
sumstat <- as.data.frame(
  b_allL %>% group_by(class) %>% 
  drop_na(logPeakArea) %>%
  summarise(
  min = min(logPeakArea),
  max = max(logPeakArea),
  median = median(logPeakArea),
  iqr = IQR(logPeakArea),
  mean = mean(logPeakArea),
  sd = sd(logPeakArea)))

write.csv(sumstat, 'Descriptive_stats_B2_class.csv')


#
#
#


# separate study samples from quality check runs
dataf <- data %>% filter(grepl('RAD', Sample)) %>%
  filter(!grepl('QC', Sample)) %>%
  filter(!grepl('Tr', Sample))

dataq <- data %>% filter(grepl('QC', Sample))


# check peak areas for duplication
table(duplicated(b1_all))

# identify samples with missing Internal_Standard observations
table(is.na(dataf$Internal_Standard))
View(dataf %>% filter(is.na(Internal_Standard) == TRUE))

# remove the samples from the dataset
dataf <- dataf %>% filter(is.na(Internal_Standard) == FALSE)
rownames(dataf) <- dataf$Sample

View(dataf %>% filter(Internal_Standard < 100000))
# filter samples with low Internal Standard value

data1 <- dataf %>% separate(Acq_Date_Time, c('Date', 'Time'), sep = ' ')
data1$Date = as.Date(data1$Date, format = "%d/%m/%Y")
data1 <- data1 %>% mutate(BatchB = ifelse(grepl('RaDICA', Sample), 'B1', 'B2'))

dataq <- dataq %>% separate(Acq_Date_Time, c('Date', 'Time'), sep = ' ')
dataq$Date = as.Date(dataq$Date, format = "%d/%m/%Y")
dataq <- dataq %>% mutate(BatchB = ifelse(Date > '2020-04-19', 'B2', 'B1'))

data1$Batch <- as.factor(data1$Batch)
dataq$Batch <- as.factor(dataq$Batch)

# flag poor quality QC samples
dataq$qc <- ifelse(rownames(dataq) %in% c(74,103,122,140), 'Bad', 'Good')
# remove poor quality QC samples
dataq <- dataq %>% filter(qc == 'Good') %>% dplyr::select(!qc)

# visualise internal and external standard peak areas across samples/QCs
dev.new()
p1 <- data1 %>% 
  filter(BatchB == 'B1') %>%
  #filter(Batch == 1) %>%
  ggplot(aes(x = Date, y = Internal_Standard, group = as.factor(Date))) +
  geom_boxplot(aes(fill = Batch), outliers = FALSE) + 
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('Batch 1') +
  geom_point(data = dataq %>% filter(BatchB == 'B1'), #%>%
             #filter(Batch == 1), 
             aes(x = Date, y = Internal_Standard, colour = qc),
             shape = 19) +
  geom_point(data = dataq %>% filter(BatchB == 'B1'), #%>%
              # filter(Batch == 1), 
             aes(x = Date, y = Octane, colour = qc),
             shape = 2) 

p1

# against behaviour of one ES

iso <- data1 %>%
  filter(BatchB == 'B1') %>%
  ggplot(aes(x = Date)) +
  #geom_boxplot(aes(x = Date, y = Isoprene, group = as.factor(Date)) 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  #geom_point(data = dataq %>% filter(BatchB == 'B1'), #%>%
             # filter(Batch == 1), 
             #aes(x = Date, y = Isoprene),
             #shape = 2, colour = 'black') +
  geom_point(data = dataq %>% filter(BatchB == 'B1') %>%
             pivot_longer(cols = c(Internal_Standard, Octane),
                          names_to = 'Standard',
                          values_to = 'peakArea'), #%>%
               #filter(Standard == 'Internal_Standard'),
             aes(x = Date, y = peakArea, 
                 colour = Batch, 
                 shape = Standard)) +
  scale_shape_manual(values = c(Internal_Standard = 19,
                                Octane = 2)) +
  ylab('Peak area') +
  ggtitle('Pre-covid dataset')

dev.new()
iso

ggsave('IS_ESvsDate_B1.tiff', dpi = 300, width = 180, height = 110, unit = 'mm')

# against behaviour of all ES
dataqL <- dataq %>% filter(BatchB == 'B1') %>%
  pivot_longer(cols = c(all_of(std)),
               names_to = 'Standard',
               values_to = 'peakArea') %>%
  mutate(RA = peakArea/Internal_Standard)

es <- dataqL %>% 
  filter(Standard != c('1_methylindole')) %>%
  ggplot(aes(x = Date, y = peakArea)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x = Date, y = log(RA), colour = Standard)) +
  #geom_line(data = dataqL %>% group_by(Date) %>% summarise(medRA = median(RA)),
            #aes(x = Date, y = medRA), colour = 'black', size = 1) +
  ylab('log (Relative abundance)') +
  ggtitle('Pre-covid dataset')

es

ggsave('RA_All_ESvsDate_B1.tiff', dpi = 300, width = 250, height = 110, unit = 'mm')

grid.arrange(p1, iso, nrow = 2)

p1a <- data1 %>% 
  filter(BatchB == 'B2') %>%
  ggplot(aes(x = Date, y = Internal_Standard, group = as.factor(Date))) +
  geom_boxplot(aes(fill = Batch), outliers = FALSE) + 
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('Post-covid') +
  geom_point(data = dataq %>% filter(BatchB == 'B2'), 
             aes(x = Date, y = Internal_Standard, colour = Batch),
             shape = 19) +
  geom_point(data = dataq %>% filter(BatchB == 'B2'), 
             aes(x = Date, y = Octane),
             shape = 2, colour = 'black')

p1a

p2 <- data1 %>% 
  filter(BatchB == 'B1') %>%
  ggplot(aes(x = Date, y = Internal_Standard, group = as.factor(Date))) +
  facet_wrap(~ Batch, scale = 'free') +
  geom_boxplot(outliers = FALSE) + 
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('Samples') +
  geom_point(data = dataq %>% filter(BatchB == 'B1'), 
             aes(x = Date, y = Internal_Standard, colour = Batch)) +
  geom_point(data = dataq %>% filter(BatchB == 'B1'), 
             aes(x = Date, y = Octane),
             shape = 2)

p2

p2a <- data1 %>% 
  filter(BatchB == 'B2') %>%
  ggplot(aes(x = Date, y = Internal_Standard, group = as.factor(Date))) +
  facet_wrap(~ Batch, scale = 'free') +
  geom_boxplot(outliers = FALSE) + 
  theme_bw() +
  theme(legend.position = 'none') +
  ggtitle('Samples') +
  geom_point(data = dataq %>% filter(BatchB == 'B2'), 
             aes(x = Date, y = Internal_Standard, colour = Batch)) +
  geom_point(data = dataq %>% filter(BatchB == 'B2'), 
             aes(x = Date, y = Octane),
             shape = 2)

p2a


dev.new()  
p2

grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)




# calculate peak area ratio  between internal and external standard (Octane)
# visualise across batches
dataq$StRatio <- dataq$Octane/dataq$Internal_Standard

dataq %>% filter(BatchB == 'B1') %>%
  ggplot(aes(x = Date, 
             #y = StRatio,
             y = StRatio, 
             colour = Batch)) +
  geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab('Relative abundance)') +
  ggtitle('Pre-covid dataset')


ggsave('RA_ESvsDate_B1.tiff', dpi = 300, width = 150, height = 110, unit = 'mm')


dataq %>% pivot_longer(cols = c(Internal_Standard, Octane),
                                names_to = 'Std', values_to = 'PeakArea') %>%
  ggplot(aes(x = Batch, y = PeakArea, fill = Std)) + geom_boxplot()

#
#
#

# PCA of relative abundance values in Batch1
# format QC dataset
data$Acq_Date_Time <- gsub(' .*', '', data$Acq_Date_Time)
data$Acq_Date_Time <- as.Date(data$Acq_Date_Time, 
                                   format = '%d/%m/%Y' )

dataq1 <- dataq %>% filter(BatchB == 'B1') 
assayq1 <- dataq1[,5:317]
rownames(assayq1) <- dataq1$Sample
assayq1RA <- assayq1/assayq1[,160]
std <- c('Acetone', 'Isoprene', 'Benzene', '3_Pentanone', '1,4_dioxane',
         'Pyridine', 'Toluene', 'Octane', 'P_Xylene', 'Nonane', 'Benzaldehyde',
         '1_Heptanol', 'Decane', '3_Carene', 'Limonene', 'Undecane', 'Nonanal',
         '1,2,3,4_tetrahydronaphthalene', 'Dodecane', '1_methylindole', 'Tridecane',
         'Pentadecane')
assayq1RA <- assayq1RA[,std]
nas <- nas %>% filter(column != 'Internal_Standard')
assayq1RA <- assayq1RA[,nas$column[-1]]
assayq1RAl <- log(assayq1RA)

dataq1RAl <- assayq1RAl %>% mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time))

pc_qc_ra <- pca(assayq1RAl, scale = TRUE, center = TRUE)
scores_qc_ra <- pc_qc_ra$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(dataq1RAl %>% dplyr::select(Sample, Batch, Acq_Date_Time))

plot_qc_ra <- scores_qc_ra %>% ggplot(aes(x = PC1, y = PC2, 
                            colour = as.factor(Batch))) +
                            #colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  #scale_colour_viridis_c(trans = 'date') +
  geom_point() +
  theme_bw() +
  xlab('PC1 (59.2%)') +
  ylab('PC2 (25%.3)') +
  #ggtitle('External standards') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(colour = 'Batch')

plot_qc_ra

# format sample dataset
datab1 <- data1 %>% filter(BatchB == 'B1') 
assayb1 <- datab1[,5:317]
assayb1RA <- assayb1/assayb1[,160]
assayb1RA <- assayb1RA[,-160]
assayb1RA[assayb1RA == 0] <- NA
assayb1RAl <- log(assayb1RA)

# remove variables with > 30% observations missing
assayAll <- rbind(assayb1, assayq1)

nas <- as.data.frame(colSums(is.na(assayb1RAl) == TRUE)/nrow(assayb1RAl)) %>%
  mutate(column = rownames(.))
colnames(nas)[1] <- 'naF'
nas <- nas %>% filter(naF < 0.3)

assayb1RAl <- assayb1RAl[,nas$column[-1]]
datab1RAl <- assayb1RAl %>% mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time)) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
            by = c('Sample' = 'Sample_ID')) %>%
  mutate(class = ifelse(grepl('B1', Sample), 'BG', 'Sample'))
rownames(datab1RAl) <- datab1RAl$Sample

pc_s_ra <- pca(datab1RAl[,-c(254:ncol(datab1RAl))], scale = TRUE, center = TRUE)

scores_s_ra <- pc_s_ra$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(datab1RAl %>% dplyr::select(Sample, Batch, Acq_Date_Time, class,
                                        Diagnosis, CoreVisit))

plot_s_ra <- scores_s_ra %>% ggplot(aes(x = PC1, y = PC2, shape = class,
                                          colour = as.factor(CoreVisit))) +
  #colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  #scale_colour_viridis_c(trans = 'date') +
  geom_point() +
  #geom_text(label = rownames(scores_s_ra)) +
  theme_bw() +
  xlab('PC1 (27.2%)') +
  ylab('PC2 (10.9%)') +
  #ggtitle('Study samples') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = 'Batch')

plot_s_ra

pc_ra <- arrangeGrob(plot_qc_ra, plot_s_ra, nrow = 1, widths = c(0.47, 0.53))

plot(pc_ra)

ggsave('PCA_RA.tiff', pc_ra, dpi = 300, units = 'mm', width = 280, height = 120)

# remove outlying observation and repeat PCA
o <- c(175, 176, 214, 215, 441)
pc_s_ra <- pca(datab1RAl[-o,-c(254:ncol(datab1RAl))], scale = TRUE, center = TRUE)

# look at Batch 1 only
b11 <- datab1RAl %>% filter(Batch == '1') 

pc_sb1_ra <- pca(b11[,-c(254:ncol(b11))], scale = TRUE, center = TRUE)

scores_sb1_ra <- pc_sb1_ra$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(datab1RAl %>% dplyr::select(Sample, Batch, Acq_Date_Time, class))

plot_sb1_ra <- scores_sb1_ra %>% ggplot(aes(x = PC1, y = PC2, shape = class,
                                        #colour = as.factor(Batch))) +
  colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  scale_colour_viridis_c(trans = 'date') +
  geom_point() +
  #geom_text(label = rownames(scores_s_ra)) +
  theme_bw() +
  xlab('PC1 (27.2%)') +
  ylab('PC2 (10.9%)') +
  #ggtitle('Study samples') +
  theme(plot.title = element_text(hjust = 0.5)) #+
  #labs(colour = 'Batch')

plot_sb1_ra


#
#
#

# PCA of raw peak areas in Batch 1
# QC samples
assayq1l <- log(assayq1[,std])

dataq1l <- assayq1l %>% mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time)) 

pc_qc_raw <- pca(assayq1l, scale = TRUE, center = TRUE)
scores_qc_raw <- pc_qc_raw$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(dataq1l %>% dplyr::select(Sample, Batch, Acq_Date_Time))

plot_qc_raw <- scores_qc_raw %>% ggplot(aes(x = PC1, y = PC2,
                                            #shape = class,
                                          colour = as.factor(Batch))) +
  #colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  #scale_colour_viridis_c(trans = 'date') +
  geom_point() +
  theme_bw() +
  xlab('PC1 (88.5%)') +
  ylab('PC2 (03.6%)') +
  #ggtitle('External standards') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  labs(colour = 'Batch')

plot_qc_raw  

# study samples
assayb1[assayb1 == 0] <- NA
assayb1l <- log(assayb1)
assayb1l <- assayb1l[,nas$column[-1]]

datab1l <- assayb1l %>% mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time)) %>%
  mutate(class = ifelse(grepl('B1', Sample), 'BG', 'Sample')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
            by = c('Sample' = 'Sample_ID'))
rownames(datab1l) <- datab1l$Sample

datab1l <- datab1l[-c(441, 215),]

pc_s_raw <- pca(datab1l[,-c(254:ncol(datab1l))], scale = TRUE, center = TRUE)
scores_s_raw <- pc_s_raw$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(datab1l %>% dplyr::select(Sample, Batch, Acq_Date_Time, class,
                                      Diagnosis, CoreVisit))

plot_s_raw <- scores_s_raw %>% ggplot(aes(x = PC1, y = PC2, 
                                          shape = class,
                                          colour = as.factor(CoreVisit))) +
  #colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  #scale_colour_viridis_c(trans = 'date') +
  #geom_text(label = rownames(scores_s_raw)) +
  geom_point() +
  theme_bw() +
  xlab('PC1 (41.7%)') +
  ylab('PC2 (9.5%)') +
  #ggtitle('Study samples') +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(colour = 'Batch')

plot_s_raw

pc_raw <- arrangeGrob(plot_qc_raw, plot_s_raw, nrow = 2, heights = c(0.4, 0.6))
plot(pc_raw)


ggsave('PCA_RAW.tiff', pc_raw, dpi = 300, units = 'mm', width = 120, height = 130)

# all PCA scores plots raw and ra data
plots_pc <- arrangeGrob(plot_qc_raw,  plot_s_raw, plot_qc_ra,  plot_s_ra,
                        widths = c(0.43, 0.57),
                        nrow = 2, ncol = 2)

plot(plots_pc)

ggsave('PCA_RAvsRAW.tiff', plots_pc, dpi = 300, units = 'mm', width = 185, height = 150)

# distribution of QC peak areas vs relative abundances

assayq <- dataq[,std]
rownames(assayq) <- dataq$Sample

assayqRA <- dataq[,5:317]/dataq[,164]
assayqRA <- assayqRA[,std]
rownames(assayqRA) <- dataq$Sample

assayqL <- assayq %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(cols =! Sample, names_to = 'comps', values_to = 'peakArea') %>%
  mutate(logPeakArea = log(peakArea)) %>%
  left_join(dataq %>% dplyr::select(Sample, Batch, BatchB))

assayqRAL <- assayqRA %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(cols =! Sample, names_to = 'comps', values_to = 'peakArea') %>%
  mutate(logPeakArea = log(peakArea)) %>%
  left_join(dataq %>% dplyr::select(Sample, Batch, BatchB))


mean(assayqL$logPeakArea)
sd(assayqL$logPeakArea)
quantile(assayqL$logPeakArea)

##
fn_cv <- function(sigma) {

  sqrt(
    exp( ((log(10))^2)*sigma^2) - 1
  ) * 100 
  
}


fn_mean <- function(mu,var) {
  exp(mu + var/2)
}
fn_var <- function(mu,var){
  #exp(var-1)*exp(2*mu + var)
  log(1 + (var)/(mu^2) )
}

fn_var(1.66,0.144^2)

##

mean(assayqRAL$logPeakArea, na.rm = TRUE)
sd(assayqRAL$logPeakArea, na.rm = TRUE)
quantile(assayqRAL$logPeakArea)

IQR(assayqlL$peakArea)
IQR(assayqRAlL$peakArea, na.rm = TRUE)

sd_assayql <- assayqL %>% 
  group_by(comps) %>%
  summarise(sd = sd(logPeakArea)) %>%
  mutate(cv = fn_cv(sd))


densq_raw <- assayqlL %>%
  ggplot(aes(x = peakArea)) + geom_density() + theme_bw() +
  xlab('log(PeakArea)')

densq_ra <- assayqRAlL %>%
  ggplot(aes(x = peakArea)) + geom_density() + theme_bw() +
  xlab('log(Relative abundance)')

densq_g <- arrangeGrob(densq_raw, densq_ra, ncol = 2, 
                      top = textGrob('Density plots of raw and relative abundance data'))
plot(densq_g)
ggsave('DensityPlotsQC.tiff', densq_g, dpi = 300, width = 125, height = 65, unit = 'mm')

# DRIFT CORRECTION

#
#
#

# evaluate noise features in QC samples
library(forcats)

assayq1 <- dataq1[,5:317]

nasq1 <- as.data.frame(colSums(is.na(assayq1) == TRUE)/nrow(assayq1)*100) %>%
  mutate(column = rownames(.)) 
colnames(nasq1)[1] <- 'naFq'
nasq1 %>% arrange(naFq) %>% ggplot(aes(x = naFq)) + geom_histogram() +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank()) +
  xlab('% NA') +
  geom_vline(xintercept = 30, colour = 'black') +
  geom_vline(xintercept = 50, colour = 'blue') +
  geom_vline(xintercept = 70, colour = 'red') +
  annotate("text", x=15, y=75, label= "< 30% NA \n p = 214", colour = 'black', cex = 3) +
  annotate("text", x=40.5, y=75, label= "< 50% NA \n p = 272", colour = 'blue', cex = 3) +
  annotate("text", x=60.5, y=75, label= "< 70% NA \n p = 298", colour = 'red', cex = 3) +
  ggtitle('Missing data in quality control dataset')
   
table(nasq1$naFq < 70)

ggsave('NA_ES.tiff', dpi = 300, unit = 'mm', width = 100, height = 70)


dataq1L <- dataq1 %>% pivot_longer(cols = c(5:317), names_to = 'comps', values_to = 'peakArea') %>%
  mutate(Std = ifelse(comps %in% std, 'Spike', 'Noise'),
         logPeakArea = log(peakArea))

dev.new()
dataq1L %>% ggplot(aes(x = Date, y = logPeakArea, group = comps, 
                       #colour = Std
                       colour = comps
                       )) + 
  geom_line() +
  theme_bw() +
  theme(legend.position = 'none') #+
  ylim(2, 8)

table(dataq1L$peakArea > 1000)

corES <- cor(assayq1, use = 'pairwise.complete.obs')
table(corES[corES > 0.5] == TRUE)

qnas <- dataq1L %>% group_by(comps) %>% summarise(nona = sum(is.na(peakArea) == TRUE)/69)
table(qnas$nona < 0.3)

# format Batch 1 database
b1_complete <- rbind(dataq1[,- ncol(dataq1)], datab1[,2:319]) %>%
  dplyr::select(!BatchB)
rownames(b1_complete) <- b1_complete$Sample
b1_assay <- b1_complete[,5:ncol(b1_complete)]

# remove variables with > 30% observations missing
# in sample dataset
assayb1 <- datab1[,6:318]

nasb1 <- as.data.frame(colSums(is.na(assayb1) == TRUE)/nrow(assayb1)) %>%
  mutate(column = rownames(.)) 

colnames(nasb1)[1] <- 'naFb'
nasb1f <- nasb1 %>% filter(naFb < 0.3)

setdiff(std, rownames(nas))

# in qc dataset
assayq1 <- dataq1[,5:317]

nasq1 <- as.data.frame(colSums(is.na(assayq1) == TRUE)/nrow(assayq1)) %>%
  mutate(column = rownames(.)) 

colnames(nasq1)[1] <- 'naFq'
nasq1f <- nasq1 %>% filter(naFq < 0.3)

# compare missing value ratio between qc and sample class dataset
nas <- nasb1 %>% left_join(nasq1)
View(nas %>% filter(naFb < 0.3))

b1_assay <- b1_assay[,nas$column]
b1_complete <- cbind(b1_complete[,1:4], b1_assay) 

# detection of outliers in QC samples based on robust PCA (Batch 1)
library(pcaPP)

dataq1 <- dataq %>% filter(BatchB == 'B1')

assayq1 <- dataq1[,5:317]

assayq1[assayq1 < 1000] <- NA

rownames(assayq1) <- dataq1$Sample

std <- c('Acetone', 'Isoprene', 'Benzene', '3_Pentanone', '1,4_dioxane',
'Pyridine', 'Toluene', 'Octane', 'P_Xylene', 'Nonane', 'Benzaldehyde',
'1_Heptanol', 'Decane', '3_Carene', 'Limonene', 'Undecane', 'Nonanal',
'1,2,3,4_tetrahydronaphthalene', 'Dodecane', '1_methylindole', 'Tridecane',
'Pentadecane')

assayq1 <- assayq1[,std]
assayq1 <- assayq1[,-160]
assayq1 <- assayq1[, nas$column]

assayq1l <- log(assayq1)

# identify QC samples of poor quality and remove from dataset
# all poor QCs from Batch 2
qna <- which(is.na(assayq1) == TRUE, arr.ind = TRUE) %>% as.data.frame()
unique(qna$row)
View(dataq[c(74,103,122,140),])
assayq <- assayq[-c(74,103,122,140),]

# impute missing values with NIPALS
assayq1l_imputed <- impute.nipals(assayq1l, ncomp = 10)

pc <- PCAgrid(assayq1l, scale = 'sd', center = 'mean')
pc1 <- prcomp(assayq1l, scale = TRUE, center = TRUE)
biplot(pc1)

scores <- as.data.frame(pc$scores) %>%
  mutate(Sample = rownames(.)) %>% 
  left_join(dataq %>% dplyr::select(Sample, Batch))

scores %>% ggplot(aes(x = Comp.1, y = Comp.2, colour = Batch)) + geom_point()

scores1 <- as.data.frame(pc1$x) %>%
  mutate(Sample = rownames(.)) %>% 
  left_join(dataq %>% dplyr::select(Sample, Batch))

scores1 %>% ggplot(aes(x = PC1, y = PC2, colour = Batch)) + geom_point()

sdod <- PCdiagplot(assayq1l, pc, plotbw = FALSE)

# identify outliers based on the critical OD and SD values (0.999 quantile)
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(assayq1l)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,3] | V2 > sdod$critSD[2,3])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(assayq1)
odOut <- od %>% filter(V1 > sdod$critOD[1,3] | V2 > sdod$critOD[2,3])

# remove outliers from the QC database
'%ni%' <- Negate('%in%')
dataq1 <- dataq1 %>% filter(Sample %ni% rownames(sdOut))

# potentially remove all Batch 4?

# merge QC and study samples from covid batch 1
datab1 <- data1 %>% filter(BatchB == 'B1') %>%
  mutate(class = 'Sample') %>% relocate(class)
datab1 <- rbind(datab1, dataq1[,-c(318:319)] %>% 
                  mutate(class = 'QC') %>% relocate(class))
rownames(datab1) <- datab1$Sample
datab1 <- datab1 %>% dplyr::select(!c(BatchB, Time))

# log-transform peak areas
datab1_t <- cbind(datab1[,1:4], log(datab1[,5:ncol(datab1)]))


# effect of class on raw data
library(mixOmics)

rownames(datab1_t) <- c(1:nrow(datab1_t))

pca_raw <- pca(datab1_t[,5:ncol(datab1_t)], scale = TRUE, center = TRUE)
plotIndiv(pca_raw,
          group = datab1_t$Batch,
          legend = TRUE)


# median fold change normalisation of log-transformed sample data
# as function
assayb1 <- datab1_t[,8:ncol(datab1_t)] 
d <- assayq1l

mfc <- function(d){
  mediansV <- sapply(1:ncol(d), function(x) {
    median(d[,x], na.rm = TRUE)
  }) # compute median of each variable
  assayb1n <- mapply('/', d, mediansV) %>% t()
  # divide each observation by variable median
  mediansS <- sapply(1:ncol(assayb1n), function(x) { 
    median(assayb1n[,x], na.rm = TRUE)
  }) # compute sample median of the rescaled dataset
  assayb1nn <- d/mediansS
  # divide each observation by sample median 
}

assayb1ln <- mfc(assayb1l)

pca_mf <- pca(assayb1ln, scale = TRUE, center = TRUE)
plotIndiv(pca_mf,
          group = datab1$Batch,
          legend = TRUE,
          pch = 1)


# CC normalisation
rownames(assayq1l) <- dataq1l$Sample

assayq1l <- assayq1l[,nas$column[-1]]

Yqc <- assayq1l

pca_yqc <- pca(Yqc, scale = TRUE, center = TRUE)
loadings1 <- pca_yqc$loadings[['X']] %>% as.data.frame() %>% dplyr::select(PC1) %>%
  as.matrix()

Y <- as.matrix(rbind(assayb1l, assayq1l))
Y <- as.matrix(assayb1l)

x0 <- Y
x0[is.na(Y)] <- 0 # managment of NAs in matrix multiplication

drift <- (x0 %*% loadings1) %*% t(loadings1)
drift[is.na(Y)] <- NA

Z <- Y - drift

#rownames(Z) <- 1:nrow(Z)

# overlay with Sample type category
pca_cc <- pca(Z, scale = TRUE, center = TRUE)

scores_cc <- pca_cc$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time)) %>%
  mutate(class = ifelse(grepl('QC', Sample), 'ES', 
                        ifelse(grepl('B1', Sample), 'BG', 'S'))) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
            by = c('Sample' = 'Sample_ID'))

plot_cc <- scores_cc %>% 
  filter(class != 'ES') %>%
  ggplot(aes(x = PC1, y = PC2,
             shape = class,
             colour = as.factor(Batch))) +
  #colour = Acq_Date_Time)) +
  #scale_colour_brewer(palette = 'Set1', name = 'Batch') +
  #scale_colour_viridis_c(trans = 'date') +
  geom_point() +
  theme_bw() +
  xlab('PC1 (20.5%)') +
  ylab('PC2 (10.2%)') +
  #ggtitle('External standards') +
  theme(plot.title = element_text(hjust = 0.5)) +
        #legend.position = 'none') +
  labs(colour = 'Batch')

plot_cc  

datanF$Analysis_date <- as.Date(datanF$Analysis_date, '%d/%m/%Y')

# compare with raw peaks
pca_raw <- pca(Y, scale = TRUE, center = TRUE)

scores_raw <- pca_raw$x %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Batch, Acq_Date_Time)) %>%
  mutate(class = ifelse(grepl('QC', Sample), 'ES', 
                        ifelse(grepl('B1', Sample), 'BG', 'S'))) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
            by = c('Sample' = 'Sample_ID'))

assayb1RAl <- assayb1RAl %>% dplyr::select(colnames(assayq1RAl))
Y1 <- rbind(assayb1RAl, assayq1RAl)

#
#
#


# normalise peak areas to the internal standard
dataf1 <- b1_all %>% dplyr::select(9:ncol(b1_all))

max(dataf1[is.na(dataf1) == FALSE])
min(dataf1[is.na(dataf1) == FALSE])

dataf1[dataf1 == 0] <- NA

datafn <- (dataf1/dataf1[,160])


dataf1L <- dataf1 %>% 
  mutate(Sample = rownames(.)) %>%
  left_join(b1_all %>% dplyr::select(Sample, class, Date)) %>%
  pivot_longer(cols = !c(Sample, class, Date), names_to = 'comp', values_to = 'peakArea') %>%
  mutate(logPeakArea = log(peakArea))

datafnL <- datafn %>%
  mutate(Sample = rownames(.)) %>%
  left_join(b1_all %>% dplyr::select(Sample, class, Date)) %>%
  pivot_longer(cols = !c(Sample, class, Date), names_to = 'comp', values_to = 'RA') %>%
  mutate(logRA = log(RA))

mean(dataf1L$logPeakArea, na.rm = TRUE)
sd(dataf1L$logPeakArea, na.rm = TRUE)
quantile(dataf1L$logPeakArea, na.rm = TRUE)

mean(datafnL$logPeakArea, na.rm = TRUE)
sd(datafnL$logPeakArea, na.rm = TRUE)
quantile(datafnL$logPeakArea, na.rm = TRUE)

# check normalisation and remove IS from dataset
table(datafn[160]) 
datafn <- datafn[,-160]


# relative standard deviation of raw and normalised data
rsd_f1 <- dataf1L %>%
  group_by(class, comp) %>%
  summarise(mean = mean(logPeakArea, na.rm = TRUE),
            sd = sd(logPeakArea, na.rm = TRUE),
            rsd = 100*sd/abs(mean))

rsd_f1 %>% group_by(class) %>%
  summarise(meanRSD = mean(rsd),
            minRSD = min(rsd),
            maxRSD = max(rsd),
            medianRSD = median(rsd, na.rm = TRUE))

rsd_fn <- datafnL %>%
  group_by(class, comp) %>%
  summarise(mean = mean(logRA, na.rm = TRUE),
            sd = sd(logRA, na.rm = TRUE),
            rsd = 100*sd/abs(mean))

rsd_fn %>% group_by(class) %>%
  summarise(meanRSD = mean(rsd, na.rm = TRUE),
            minRSD = min(rsd, na.rm = TRUE),
            maxRSD = max(rsd, na.rm = TRUE),
            medianRSD = median(rsd, na.rm = TRUE))

datafnL <- datafnL %>% mutate(RA100 = 1000*RA,
                              logRA100 = log(RA100))


# pca of relative abudnances
# remove variables with > 30% NAs
datafn1 <- datafn
datafn1 <- log(datafn1)
naf <- as.data.frame(colSums(is.na(datafn1) == TRUE)/nrow(datafn1))
colnames(naf) <- 'naf'
nafk <- naf %>% filter(naf < 0.3)
datafn1 <- datafn1[,rownames(nafk)]

library(mixOmics)
rownames(datafn) <- 1:nrow(datafn)
pc_ra <- mixOmics::pca(datafn, scale = TRUE, center = TRUE) 
plotIndiv(pc_ra)

datafn <- datafn[-c(215, 441, 214, 175, 176),]

# descriptive statistics
library(grid)
library(stats)

max(datafn[is.na(datafn) == FALSE])
min(datafn[is.na(datafn) == FALSE])

max(dataf1[is.na(dataf1) == FALSE])
min(dataf1[is.na(dataf1) == FALSE])

q <- quantile(dataf1, na.rm = TRUE)
q

q/q[3]

qn <- quantile(datafn, na.rm = TRUE)
qn
qn/qn[3]


quantile(datafn[is.na(datafn) == FALSE], probs = c(0.99))
sd(as.matrix(datafn), na.rm = TRUE)
mean(as.matrix(datafn), na.rm = TRUE)

colnames(data1)[1] <- 'Sample_ID'

dataf1L <- dataf1 %>% mutate(Sample_ID = rownames(.)) %>%
  pivot_longer(cols = colnames(dataf1),
                                   names_to = 'comp', values_to = 'peakArea') %>%
  left_join(meta) %>%
  left_join(data %>% dplyr::select(Sample, Batch), by = c('Sample_ID' = 'Sample')) %>%
  mutate(BatchB = ifelse(grepl('RaDICA', Sample_ID), 'B1', 'B2'))

datafnL <- datafn %>% mutate(Sample_ID = rownames(datafn)) %>% 
  pivot_longer(cols = colnames(datafn), names_to = 'comp', values_to = 'RA') %>%
  left_join(meta) %>% 
  left_join(data1 %>% dplyr::select(Sample, Batch, BatchB), by = c('Sample_ID' = 'Sample'))


datafnL <- datafnL %>% mutate(logRA = ifelse(RA == 0, NA, log(RA)))

datafnL$Analysis_date <- as.Date(datafnL$Analysis_date, format = '%d/%m/%Y')
datafnL %>% filter(comp == 'Isoprene') %>% 
filter(BatchB == 'B1') %>%  
#ggplot(aes(x = Analysis_date, y = logRA)) + geom_point()
ggplot(aes(x = comp, y = logRA)) + geom_boxplot()


max(datafnL$RA, na.rm = TRUE)
quantile(datafnL$RA, na.rm = TRUE, prob = 0.99)

datafnL1 <- datafnL %>% filter(BatchB == 'B1')
datafnL2 <- datafnL %>% filter(BatchB == 'B2')

quantile(datafnL1$RA, na.rm = TRUE, prob = 0.99)
quantile(datafnL2$RA, na.rm = TRUE, prob = 0.99)

# look at observation above 99% percentile
# in each covid batch separately
Hiobs <- datafnL1 %>%
  filter(RA > quantile(datafnL$RA, na.rm = TRUE, prob = 0.99))
View(table(Hiobs$comp))
n_distinct(Hiobs$CoreVisit)
View(table(Hiobs$Analysis_date))

quantile(datafnL$logRA[is.na(datafnL$logRA) == FALSE])
sd(datafnL$logRA, na.rm = TRUE)*2.5
mean(datafnL$logRA, na.rm = TRUE)

# look at outlier observations
View(datafnL %>% filter(logRA > 0.798))
View(datafnL %>% filter(logRA > 0.798))

dens_raw <- dataf1L  %>% filter(BatchB == 'B1') %>%
  ggplot(aes(x = log(peakArea))) + geom_density() + theme_bw() +
  xlab('log(Peak area)')
dens_raw

dens <- datafnL %>%
  ggplot(aes(x = RA)) + geom_density() + theme_bw() + xlab('relative abudnance')

log_dens <- datafnL %>% filter(BatchB == 'B1') %>%
  ggplot(aes(x = logRA), group = as.factor(Batch)) + geom_density() + theme_bw() + 
  xlab('log(Relative abundance)')

log_dens

dens_g <- arrangeGrob(dens_raw, log_dens, ncol = 2, 
                      top = textGrob('Density plots of raw and relative abundance data'))
plot(dens_g)
ggsave('DensityPlots.tiff', dens_g, dpi = 300, width = 125, height = 65, unit = 'mm')

which(datafn == 0, arr.ind = TRUE)

View(dataf %>% dplyr::select(Gamma_Butyrolactone))

#

# change in median relative abundance over time (across all compounds in study samples)
datafnL %>% filter(class %in% c('S1', 'S2', 'BG')) %>%
  group_by(Date, class) %>%
  summarise(median = median(logRA, na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = median)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~ class, scale = 'free', ncol = 1) +
  ylab('median log( Relative abundance)')
  
ggsave('Median_RAvsDate_B1.tiff', dpi = 300, units = 'mm', width = 90, height = 180)  

# and ES spike-ins
datafnL %>% filter(class %in% c('ES')) %>%
  filter(comp %in% std) %>%
  group_by(Date) %>%
  summarise(median = median(logRA, na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = median)) +
  geom_line() +
  theme_bw() +
  ylab('median log( Relative abundance)') +
  ggtitle('External Standard') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

ggsave('Median_RAvsDate_B1_ES.tiff', dpi = 300, units = 'mm', width = 90, height = 62)  

# change in peak area spike ins vs IS
p1 <- dataf1L %>%
  filter(class == 'ES') %>%
  filter(comp %in% c('Internal_Standard', std)) %>%
  ggplot(aes(x = Date, y = peakArea)) +
  geom_line(aes(colour = comp)) +
  ylim(0, 3.5*10^7) +
  scale_colour_manual(values = c('Internal_Standard' = 'red')) +
  theme_bw() +
  ylab('Peak area')

p2 <- datafnL %>%
  filter(class == 'ES') %>%
  filter(comp %in% std) %>%
  ggplot(aes(x = Date, y = logRA)) +
  geom_line(aes(colour = comp)) +
  ylim(-4, 3) +
  theme_bw() +
  ylab('log( Relative sbundance)')


ggsave('ES_peakArea_vs_Date.tiff', p1, dpi = 300, unit = 'mm', width = 175, height = 80)

ggsave('ES_logRA_vs_Date.tiff', p2, dpi = 300, unit = 'mm', width = 218, height = 80)

#

# heteroscedasticity
sumint1 <- datafnL %>% 
  group_by(comp) %>%
  drop_na() %>%
  summarise(meanInt = mean(RA),
            sdInt = sd(RA)) %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())

scatter <- sumint1 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point(alpha = 0.5) + 
  theme_bw() + ylim(0,4) + ylab('SD(relative abundance)') + xlab('rank(mean relative abundance)') +
  ggtitle('Non-transformed data') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

# log transformed
sumint2 <- datafnL %>% 
  group_by(comp) %>%
  drop_na() %>%
  summarise(meanInt = mean(log10(RA)),
            sdInt = sd(log10(RA))) %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())

log_scatter <- sumint2 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point(alpha = 0.5) + 
  theme_bw() + ylab('SD(relative abundance)') + xlab('rank(mean relative abundance)') +
  ggtitle('Log-transformed data') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))


scatter_g <- arrangeGrob(scatter, log_scatter, ncol = 2)
plot(scatter_g)
ggsave('ScatterPlots.tiff', scatter_g, dpi = 300, width = 125, height = 65, unit = 'mm')


# annotate with sample metadata
colnames(meta)[colnames(meta) == 'Sample'] <- 'Sample_type'
colnames(dataf)[colnames(dataf) == 'Sample'] <- 'Sample_ID'

datan <- cbind(datafn, dataf[,1:3]) %>% left_join(meta) 

setdiff(meta$Sample_ID, dataf$Sample_ID)
setdiff(dataf$Sample_ID, meta$Sample_ID)

# remove internal standard from variables
datan <- datan %>% dplyr::select(!Internal_Standard)

# divide into VOC and metadata datasets
rownames(datan) <- datan$Sample_ID

assay <- datan[,1:312]

meta1 <- datan[,313:ncol(datan)]

#
#
#

# MISSING DATA
# plot frequency of NAs across features
b1_all # Batch 1 dataset with all classes

cnas <- function(d, cutoff){
  fnas <- as.data.frame(colSums(is.na(d))/nrow(d)*100)
  colnames(fnas)[1] <- 'RatioNA'
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

cnas1 <- function(cl) {
  cnas(b1_all %>% filter(class == cl) %>% dplyr::select(!1:8), 20) +
  ggtitle(cl)}

cnas1('Blank')

cnas_plots <- lapply(unique(b1_all$class), cnas1)
NAhists <- grid.arrange(grobs = cnas_plots, nrow = 2, ncol = 3)

ggsave('Histogram_NA_B1-class.tiff', NAhists, dpi = 300, unit = 'mm', width = 160, height = 100)

# % missing per sample (BG, S1 and S2 only) - after summarising breath samples?
s1 <- cnas(t(b1_all %>% filter(class %in% c('S1')) %>% dplyr::select(!1:8)),
     cutoff = 40) + 
  ggtitle('S1')

s <- grid.arrange(s1, s2, ncol = 2)
ggsave('Histogram_NA_B1_sample.tiff', s, dpi = 300, unit = 'mm', width = 105, height = 50)

table(b1_all$class)

# effect of biological variable on missing values
b1_all_s <- b1_all %>% filter(class %in% c('S1', 'S2')) %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, Diagnosis),
            by = c('Sample' = 'Sample_ID')) %>%
  mutate(Treatment = ifelse(CoreVisit %in% c('CV1', 'CV2'), 'Before', 'After')) %>%
  relocate(CoreVisit, Diagnosis, Treatment)

cnas(b1_all_s %>%
#filter(class == 'S2') %>%
filter(Treatment == 'After') %>%
#filter(Treatment == 'After') %>%
dplyr::select(11:ncol(b1_all_s)),
cutoff = 20)


# plot NAs according to experimental condition and sample type
data_heat <- b1_all[,9:ncol(b1_all)]
data_heat[is.na(data_heat)] <- 0
data_heat[data_heat > 0] <- 1

rownames(data_heat) <- b1_all$Sample

library(RColorBrewer)
library(scales)

#newCols <- colorRampPalette(grDevices::rainbow(length(unique(meta_heat$Sampling_date))))
#annoCol <- newCols(length(unique(meta_heat$Sampling_date)))
#names(annoCol) <- unique(meta_heat$Sampling_date)
#annotCol <- list(Sampling_date = annoCol)

cols1 <- (hue_pal()(5))
cols2 <- (hue_pal()(5))
names(cols1) <- unique(b1_all$class)
names(cols2) <- unique(b1_all$Batch)
annotCol <- list(class = cols1,
                 Batch = cols2)

b1_all$Batch <- as.factor(b1_all$Batch)
rownames(b1_all) <- b1_all$Sample

b1_all_an <- b1_all %>% left_join(meta %>% dplyr::select(Sample_ID, Diagnosis, CoreVisit),
                                  by = c('Sample' = 'Sample_ID'))
rownames(b1_all_an) <- b1_all_an$Sample 

dev.new()
tiff('Heatmap_NA_B1_class_batch.tiff', res = 300, width = 250, height = 200, unit = 'mm')
heatmap <- pheatmap(as.matrix(t(data_heat)), 
                    annotation_col = b1_all %>% 
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

# in breath samples only
b1_s <- b1_all %>% filter(class %in% c('S1','S2')) %>% dplyr::select(!1:8)
data_heat1 <- b1_s
data_heat1[is.na(data_heat1)] <- 0
data_heat1[data_heat1 > 0] <- 1

rownames(data_heat) <- b1_s$Sample

b1_all_an <- b1_all_an %>% mutate(Treatment = ifelse(CoreVisit %in% c('CV1', 'CV2'),
                                                      'Before', 'After'))

cols3 <- (hue_pal()(5))
cols4 <- cols3[1:2]
cols5 <- cols3[4:5]
names(cols3) <- na.omit(unique(b1_all_an$Batch))
names(cols4) <- na.omit(unique(b1_all_an$Diagnosis))
names(cols5) <- unique(b1_all_an$Treatment)
annotCol1 <- list(Batch = cols3,
                 Diagnosis = cols4,
                 Treatment = cols5)


dev.new()
tiff('Heatmap_NAs_B1_Breath.tiff', res = 300, width = 250, height = 200, unit = 'mm')
heatmap1 <- pheatmap(as.matrix(t(data_heat1)), 
                    annotation_col = b1_all_an %>% 
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

# filter features with > 20% NAs
class_na <- b1_allL %>% 
  left_join(b1_all_an %>% dplyr::select(Sample, Diagnosis, Treatment)) %>%
  mutate(class = ifelse(class == 'S1' & Diagnosis == 'Asthma', 'S1A', 
                        ifelse(class == 'S1' & Diagnosis == 'Not Asthma', 'S1NA',
                               ifelse(class == 'S2' & is.na(Diagnosis) == TRUE, 'S2',
                                      ifelse(class == 'S2' & Diagnosis == 'Not Asthma', 'S2NA', 
                                             ifelse(class == 'S2' & Diagnosis == 'Asthma', 'S2A', class)))))) %>%
      filter(class != 'S2') %>%  
      group_by(class, comp) %>% 
      summarise(nona = sum(is.na(peakArea) == TRUE)/length(peakArea))

class_na <- class_na %>% mutate(nona = ifelse(nona < 0.2, 1, 0))

# consider both breath measurement together depending on biological variable
class_na1 <- b1_allL %>% 
  left_join(b1_all_an %>% dplyr::select(Sample, Diagnosis, Treatment)) %>%
  mutate(class = ifelse(class %in% c('S1', 'S2') & Diagnosis == 'Asthma', 'SAsth', 
                        ifelse(class %in% c('S1', 'S2') & Diagnosis == 'Not Asthma', 'SNAsth',
                               ifelse(class %in% c('S1','S2') & is.na(Diagnosis) == TRUE, 'S', class)))) %>%
  filter(class != 'S2') %>%  
  group_by(class, comp) %>% 
  summarise(nona = sum(is.na(peakArea) == TRUE)/length(peakArea))

class_na1 <- class_na1 %>% mutate(nona = ifelse(nona < 0.2, 1, 0))
# compare features retained across sample types
class_naF_wide <- class_na1 %>% pivot_wider(names_from = comp, values_from = nona)
assay_naF <- class_naF_wide[,-1]

# no of VOCs above threshold in all classes
sums <- as.data.frame(colSums(assay_naF))
colnames(sums) <- 'sum'
complete_var <- sums %>% filter(sum == 5) #7 if breath sub-classes included
reject_var <- sums %>% filter(sum == 0)

assay_naFf <- assay_naF %>% dplyr::select(!rownames(complete_var)) %>%
  dplyr::select(!rownames(reject_var))
rownames(assay_naFf) <- class_naF_wide$class

tiff(filename = 'Heatmap_NA_B1_class_v2.tiff', width = 200, height = 170, unit = 'mm', res = 300)
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
b1_all_f <- b1_all %>% dplyr::select(c(2:8, rownames(retain_var), rownames(complete_var))) #%>%
  dplyr::select(!'Internal_Standard')



# FILTERED DATASET
write.csv(b1_all_f, 'RADicA_B1_NAfiltered_v2.csv')
b1_all_f <- read.csv('RADicA_B1_NAfiltered_v2.csv')
b1_all_f <- b1_all_f[,-1]

cnas(b1_all_f[,8:ncol(b1_all_f)], cutoff = 20) + ggtitle('All samples')
table(colSums(is.na(b1_all_f[,8:ncol(b1_all_f)]))/nrow(b1_all_f)*100 < 20)

#
#
#

# see if peak intensities in QC sample follow the same trend as study samples
b1_all_fL <- b1_all_f %>% pivot_longer(cols = 8:ncol(b1_all_f),
                                       names_to = 'comp',
                                       values_to = 'peakArea') %>%
  mutate(logPeakArea = log(peakArea))

dev.new()

medClass <- b1_all_fL %>% 
  filter(class != 'ES') %>%
  group_by(class, Date) %>%
  summarise(medPeakArea = median(logPeakArea, na.rm = TRUE)) %>%
  left_join(b1_all_f %>% dplyr::select(Date, Batch)) %>%
  ggplot(aes(x = Date, y = medPeakArea)) +
  geom_line() +
  facet_wrap(~ factor(as.factor(class), levels = c('Blank', 'ES', 'BG', 'S1','S2')),
                      scale = 'free') +
  theme_bw() +
  ylab('median log( PeakArea)')

medClass

ggsave('Median_peakAVsDate.tiff', medClass, dpi = 300, unit = 'mm', 
       width = 200, height = 120)

medClassES <- b1_all_fL %>% 
  filter(class == 'ES') %>%
  mutate(spikeIn = ifelse(comp %in% std, 'Spike-in', 'Noise')) %>%
  group_by(spikeIn, Date) %>%
  summarise(medPeakArea = median(logPeakArea, na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = medPeakArea)) +
  facet_wrap(~spikeIn, scale = 'free', nrow = 2) +
  geom_line() +
  theme_bw() +
  ylab('median log( PeakArea)') 


medClassES

ggsave('Median_peakAVsDate_ES.tiff', medClassES, dpi = 300, unit = 'mm', 
       width = 100, height = 120)

# NAs according to peak intensity
sumint <- b1_all_fL %>% 
  group_by(comp, class) %>%
  summarise(med = median(logPeakArea, na.rm = TRUE),
            nonas = sum(is.na(logPeakArea) == TRUE))

pdf('NAvsMedianIntensity.pdf', width = 7, height = 4.5)

medVSna <- sumint %>% group_by(class) %>%
  arrange(med) %>%
  do(plots = ggplot(data = ., aes(x = med, y = nonas)) +
  geom_point(alpha = 0.6) + ggtitle(.$class) + theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()))

marrangeGrob(medVSna$plots, nrow = 2, ncol = 3)

dev.off()

sumint %>% group_by(class) %>% summarise(complete_var = sum(nonas == 0))

#  NA patterns across technical replicates
datanFl <- datanF %>% pivot_longer(cols = c(1:286), names_to = 'comp', values_to = 'RA') %>%
  pivot_wider(id_cols = c(RAD_ID, comp, CoreVisit, Batch, BatchB, Analysis_date), 
              names_from = Sample_type, values_from = RA,
              values_fn = mean)

datanFl$techNA <- ifelse(is.na(datanFl$S1) | is.na(datanFl$S2), 'Missing', 'Non-missing')

datanFl2 <- datanFl %>% select(RAD_ID, comp, S1, S2, techNA, CoreVisit) %>%
  pivot_longer(cols = c(S1, S2), names_to = 'Sample_type', values_to = 'RA')

library(gridExtra)

# per sample, across compounds
pdf('Missing_patterns_sample.pdf')
pnas <- datanFl2 %>% group_by(RAD_ID) %>%
  do(plots = ggplot(data =.) + 
       (aes(x = techNA, y = log2(RA))) + 
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
pnac <- datanFl2 %>% group_by(comp) %>%
  do(plots = ggplot(data =.) + 
       aes(x = techNA, y = log2(RA), colour = techNA) +
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
ps <- datanFl %>% group_by(comp) %>%
  do(plots = ggplot(data = .) + (aes(y = log2(S1), x = log2(S2))) + geom_point() +
       ggtitle(.$comp) + coord_fixed() +
       theme_bw(base_size = 6) + geom_abline(intercept = 0, slope = 1, colour = 'red'))

pdf('Radica_S1vsS2.pdf')
marrangeGrob(grobs = ps$plots, nrow = 4, ncol = 4)
dev.off()

corrs <- datanFl %>% dplyr::select(!MaskBG) %>% 
  drop_na() %>%
  group_by(comp) %>% summarise(cor = cor(S1, S2, method = 'pearson'))


#

# MISSING DATA IMPUTATION
# impute blank and ES datasets to be used for additive noise and drift correction
install.packages("C:/Users/q29634at/OneDrive - The University of Manchester/Documents/R/win-library/4.0/GMSimpute", 
                 repos = NULL, type = "source")

library(GMSimpute)

# proportion of missing values in Blank and ES datasets
# create assay dataset for each class
assay <- function(x){
  b1_all_f %>% filter(class == x) %>% dplyr::select(9:ncol(.))
}

rownames(b1_all_f) <- b1_all_f$Sample

assay_blank <- assay(x = 'Blank')
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

# evaluate effects of imputation 
imp_long <- function(x, y) {
  x %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(cols =! Sample, names_to = 'comp', values_to = 'peakAreaImp') %>%
    mutate(logPeakAreaImp = log(peakAreaImp)) %>%
    left_join(b1_all_fL %>% filter(class == y) %>% 
                dplyr::select(Sample, comp, logPeakArea, Date, class))
}

blank_imp <- imp_long(assay_blank_imp, 'Blank')
es_imp <- imp_long(assay_es_imp, 'ES')
bg_imp <- imp_long(assay_bg_imp, 'BG')
s1_imp <- imp_long(assay_s1_imp, 'S1')
s2_imp <- imp_long(assay_s2_imp, 'S2')

b1_imp_L <- rbind(es_imp, blank_imp, bg_imp, s1_imp, s2_imp) 

write.csv( b1_imp_L, 'RADicA_B1_NAfiltered_imputed_long.csv')
b1_imp_L <- read.csv('RADicA_B1_NAfiltered_imputed_long.csv')
b1_imp_L <- b1_imp_L[,-1]
b1_imp_L <- b1_imp_L %>% dplyr::select(!c(logPeakArea, logPeakAreaImp))

# density plot
dens_plot <-
  b1_imp_L %>% 
  filter(is.na(class) == FALSE) %>%
  ggplot() + 
  geom_density(aes(x = logPeakAreaImp, colour = 'Imputed values')) +
  geom_density(aes(x = logPeakArea, colour = 'Observed values')) +
  theme_bw() +
  #ggtitle('Density plot of log-transformed peak areas in blank samples') +
  xlab('log( Peak area)') +
  #theme(legend.position = 'none',
  facet_wrap(~ class, scale = 'free')

dens_plot

ggsave('Density_imputed.tiff', dens_plot, dpi = 300, units = 'mm', width = 240, height = 120)

# time-course change in median
time_plot <- 
  b1_imp_L %>% 
  drop_na(class) %>%
  group_by(class, Date) %>%
  summarise(medianObs = median(logPeakArea, na.rm = TRUE),
            medianImp = median(logPeakAreaImp)) %>%
  ggplot(aes(x = Date)) + 
  geom_line(aes(y = medianObs, colour = 'Observed values')) +
  geom_line(aes(y = medianImp, colour = 'Imputed values')) +
  theme_bw() +
  facet_wrap(~ class, scale = 'free') +
  ylab('median log( Peak area)') 

dev.new()
time_plot

ggsave('Time_imputed.tiff', time_plot, dpi = 300, units = 'mm', width = 240, height = 120)


# PCA complete vs imputed variables
library(ggplotify)

# create data class subsets with complete variables only
complete_set <- function(y){ 
  df <- b1_all_f %>% filter(class == y)
  rownames(df) <- df$Sample
  comp_vars <- sumint %>% filter(class == y) %>% filter(nonas == 0)
  comp_df <- df[,comp_vars$comp]
}

b1_blank_f_c <- complete_set('Blank')
b1_es_f_c <- b1_all %>% filter(class == 'ES') %>% subset(select = colnames(.) %in% std)
b1_bg_f_c <- complete_set('BG')
b1_s1_f_c <- complete_set('S1')
b1_s2_f_c <- complete_set('S2')

b1_es_f <- b1_all_f %>% filter(class == 'ES') 
b1_blank_f <- b1_all_f %>% filter(class == 'Blank') 
b1_bg_f <- b1_all_f %>% filter(class == 'BG') 
b1_s1_f <- b1_all_f %>% filter(class == 'S1') 
b1_s2_f <- b1_all_f %>% filter(class == 'S2') 

plot_pca <- function(x, y, z){
  pca <- pca(log(x), scale = TRUE, center = TRUE)
  plot_pca <- plotIndiv(pca,
          pch = 1,
          cex = 1.5,
          group = y$Batch,
          #legend = TRUE,
          title = z,
          size.title = 10)#,
          #legend.title = 'Batch')

}


pc_blank_c <- plot_pca(b1_blank_f_c, b1_blank_f, 'Blank')
pc_blank_i <- plot_pca(assay_blank_imp, b1_blank_f, 'Blank')

pc_es_c <- plot_pca(b1_es_f_c, b1_es_f, 'ES')
pc_es_i <- plot_pca(assay_es_imp, b1_es_f, 'ES')

pc_bg_c <- plot_pca(b1_bg_f_c, b1_bg_f, 'BG')
pc_bg_i <- plot_pca(assay_bg_imp, b1_bg_f, 'BG')

pc_s1_c <- plot_pca(b1_s1_f_c, b1_s1_f, 'S1')
pc_s1_i <- plot_pca(assay_s1_imp, b1_s1_f, 'S1')

pc_s2_c <- plot_pca(b1_s2_f_c, b1_s2_f, 'S2')
pc_s2_i <- plot_pca(assay_s2_imp, b1_s2_f, 'S2')

pc_plots <- arrangeGrob(pc_blank_c$graph, pc_es_c$graph, pc_bg_c$graph, pc_s1_c$graph, pc_s2_c$graph,
                        pc_blank_i$graph, pc_es_i$graph, pc_bg_i$graph, pc_s1_i$graph, pc_s2_i$graph,
                        ncol = 5, nrow = 2)

dev.new()
plot(pc_plots)

ggsave('PCA_complete_vs_imputed.tiff', pc_plots, dpi = 300, unit = 'mm', width = 310, height = 120)


#
#
#

# DRIFT CORRECTION THROUGH CPCA
# !! some variables missing from reference datasets!
# remove outliers from reference datasets (es and blank)
library(pcaPP)

pc <- PCAgrid(log(assay_es_imp), scale = 'sd', center = 'mean')
pc1 <- prcomp(log(assay_es_imp), scale = TRUE, center = TRUE)

summary(pc)
pc1$prop_expl_var

scores <- as.data.frame(pc$scores) %>%
  mutate(Sample = rownames(assay_es_imp)) %>% 
  left_join(b1_all_f %>% dplyr::select(Sample, Batch))

scores %>% ggplot(aes(x = Comp.1, y = Comp.2, colour = as.factor(Batch))) + geom_point()

scores1 <- as.data.frame(pc1$x) %>%
  mutate(Sample = rownames(.)) %>% 
  left_join(dataq %>% dplyr::select(Sample, Batch))

scores1 %>% ggplot(aes(x = PC1, y = PC2, colour = as.factor(Batch))) + geom_point()

sdod <- PCdiagplot(assay_es_imp, pc, plotbw = FALSE)

# identify outliers based on the critical OD and SD values (0.999 quantile)
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(assay_es_imp)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,3] | V2 > sdod$critSD[2,3])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(assay_es_imp)
odOut <- od %>% filter(V1 > sdod$critOD[1,3] | V2 > sdod$critOD[2,3])

# remove outliers from the QC database
'%ni%' <- Negate('%in%')
assay_es_imp1 <- assay_es_imp %>% filter(rownames(.) %ni% rownames(odOut))

#
#
#

# concatenate imputed database, with variable number matching the smallest var no
b1_imputed <- b1_imp_L %>% 
  #dplyr::select(!c(logPeakAreaImp, logPeakArea)) %>%
  pivot_wider(names_from = comp, values_from = peakArea)


diff_vars <- c(setdiff(colnames(assay_s1_imp), colnames(assay_es_imp),
               setdiff(colnames(assay_es_imp), colnames(assay_blank_imp))))

b1_imputed <- rbind(assay_es_imp[, colnames(assay_es_imp) %ni% diff_vars], 
                    assay_bg_imp[, colnames(assay_bg_imp) %ni% diff_vars],
                    assay_s1_imp[, colnames(assay_s1_imp) %ni% diff_vars],
                    assay_s2_imp[, colnames(assay_s2_imp) %ni% diff_vars],
                    assay_blank_imp[, colnames(assay_blank_imp) %ni% diff_vars])

# if loaded from csv file
b1_imputed <- b1_imputed[,1:230]

b1_all_no_blank <- b1_all_f %>% 
  filter(class != 'Blank') #%>%
  #filter(class %in% c('S1', 'S2', 'BG')) %>%
  #left_join(data %>% dplyr::select(Sample, Diagnosis, Treatment))
rownames(b1_imputed) <- b1_imputed$Sample

b1_imputed <- b1_imputed[b1_all_no_blank$Sample,]

b1_imputed <- b1_imputed[ order(match(rownames(b1_imputed), rownames(b1_all_no_blank))), ]

b1_imputed1 <- b1_imputed[,4:ncol(b1_imputed)]

b1_imputed_nIS <- b1_imputed1[, colnames(b1_imputed1) != 'Internal_Standard']


# run pca on concatenated log peak areas
pc_all <- pca(log2(b1_imputed_nIS), 
              scale = TRUE, center = TRUE, ncomp = 4)                    
Bpa_pc23 <- plotIndiv(pc_all,
          group = b1_all_no_blank$Batch,
          pch = 1,
          cex = 2,
          comp = c(1,2),
          legend = TRUE,
          title = 'Peak Area (non-normalised)',
          size.title = 10)


Bpa_pc <- arrangeGrob(Bpa_pc12$graph, Bpa_pc23$graph, nrow = 2, ncol = 1)
ggsave('PCA_raw_peak_area_imputed_Batch.tiff', Bpa_pc, dpi = 300, unit = 'mm', width = 100, height = 150)

colnames(b1_imp_L)[3] <- 'peakArea'
rsd_pa <- rsd(b1_imp_L)
imp <- rsd_pa %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()

# run PCA on IS normalised data
# compare with internal standard normalisation
b1_imp_IS <- b1_imputed1/b1_imputed1$Internal_Standard
b1_imp_IS <- b1_imp_IS[, colnames(b1_imp_IS) != 'Internal_Standard']

pca_cc2 <- pca(log(b1_imp_IS), scale = TRUE, center = TRUE, ncomp = 4)

Bis_pc23 <- plotIndiv(pca_cc2,
          group = b1_all_no_blank$Batch,
          pch = 1,
          cex = 2,
          legend = TRUE,
          comp = c(1,2),
          title = 'Relative abundance (IS normalised)',
          size.title = 10)

b1_imp_ISL <- cbind(b1_imputed[,1:3], b1_imp_IS) %>%
  pivot_longer(cols =! c(Sample, Date, class), names_to = 'comp', values_to = 'peakArea')
rsd_imp_IS <- rsd(b1_imp_ISL)
imp_IS <- rsd_imp_IS %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()

Bis_pc <- arrangeGrob(Bis_pc12$graph, Bis_pc23$graph, nrow = 2, ncol = 1)
ggsave('PCA_IS_normalised_imputed_Batch.tiff', Bis_pc, dpi = 300, unit = 'mm', width = 100, height = 150)

# run CC on spikes database
assay_es_imp1 <- b1_imputed %>% filter(class == 'ES') %>% dplyr::select(!c(Sample, Date, class))
pc_es <- pca(log(assay_es_imp1 %>% dplyr::select(!Internal_Standard)), scale = TRUE, center = TRUE, ncomp = 2) 

b1_es_f <- b1_all_f %>% filter(Sample %in% rownames(assay_es_imp1)) 

plotIndiv(pc_es,
          #group = b1_es_f$Batch,
          pch = 1)

loadings1 <- pc_es$loadings[['X']] %>% as.data.frame() %>% dplyr::select(PC1) %>%
  as.matrix()

loadings2 <- pc_es$loadings[['X']] %>% as.data.frame() %>% dplyr::select(PC2) %>%
  as.matrix()

rownames(b1_all_f) <- b1_all_f$Sample

Y <- b1_imputed %>% dplyr::select(!c(Sample, class, Date, Internal_Standard)) %>% 
  #filter(class != 'Blank') %>%
  #dplyr::select(!c(1:7)) %>% 
  log() %>% as.matrix()

# fill in missing loading values with min loading?
setdiff(colnames(Y), colnames(assay_es_imp))

mis_loads <- data.frame(comps = c(setdiff(colnames(Y), colnames(assay_es_imp))),
           PC1 = rep(-0.003, times = 8)) # change to PC1 when appropriate
rownames(mis_loads) <- mis_loads$comps
mis_loads <- mis_loads %>% dplyr::select(PC1) # change to PC1 when appropriate

loadings1m <- rbind(loadings1, mis_loads) %>% as.matrix()
loadings2m <- rbind(loadings2, mis_loads) %>% as.matrix()

#

x0 <- Y
x0[is.na(Y)] <- 0 # managment of NAs in matrix multiplication

drift <- (x0 %*% loadings1) %*% t(loadings1)
drift[is.na(Y)] <- NA

Z <- Y - drift


# overlay with Sample type category
Z1 <- Z[b1_all_no_blank$Sample,]

pca_cc <- pca(Z, scale = TRUE, center = TRUE, ncomp = 4)

groups <- b1_all_f %>% filter(class != 'Blank')

cc_pa_pc12 <- plotIndiv(pca_cc,
          #group = b1_all_no_blank$Batch,
          pch = 1,
          cex = 2,
          legend = TRUE,
          comp = c(1,2),
          title = 'Peak Area (CC normalised)',
          size.title = 10)

ggsave('PCA_CCnormalised_imputed_Batch.tiff', cc_pa_pc12$graph, dpi = 300, unit = 'mm', width = 100, height = 75)

rownames(Z) <- b1_imputed$Sample
ZL <- Z %>% as.data.frame() %>% 
  mutate(Sample = rownames(.)) %>% left_join(b1_imputed %>% dplyr::select(Sample, class, Date)) %>%
  pivot_longer(cols =! c(Sample, Date, class), names_to = 'comp', values_to = 'logPeakArea')

ZL <- ZL %>% mutate(peakArea = exp(logPeakArea))
rsd_imp_cc <- rsd(ZL)
imp_corr_cc <- rsd_imp_cc %>% ggplot(aes(x = class, y = rsd, fill = class)) + geom_boxplot()
imp_corr_cc
# 

# compare results of normalisation (density plots, RSD)
annotateL <- function(x) {
  x %>% mutate(Sample = rownames(.)) %>%
    left_join(b1_all_f %>% dplyr::select(Sample, class, Date)) %>%
    pivot_longer(cols = !c(Sample, class, Date), names_to = 'comp', values_to = 'value')%>%
    filter(comp != 'Internal_Standard') #%>%
    #filter(class == 'ES')
}

b1_impL <- annotateL(b1_imputed)
b1_imp_ISL <- annotateL(b1_imp_IS)



p1 <- b1_impL %>% group_by(Date) %>% summarise(median = median(logValue)) %>%
  ggplot(aes(x = Date, y = median)) + geom_line() + theme_bw() +
  ggtitle('Peak Area')

p <- arrangeGrob(p1, p2, p3, nrow = 3, ncol = 1)

ggsave('CC_median_time.tiff', p, dpi = 300, units = 'mm', width = 100, height = 180)


b1_impL <- b1_impL %>% mutate(logValue = log(value))
b1_imp_ISL <- b1_imp_ISL %>% mutate(logValue = log(value))
b1_imp_IS_logL <- b1_imp_IS_logL %>% mutate(logValue = value)
b1_imp_CCL <- b1_imp_CCL %>% mutate(logValue = value)



p11 <- b1_imp_CCL %>% ggplot(aes(x = logValue)) + geom_density() + theme_bw() +
  ggtitle('Peak Area + CC')
quantile(b1_imp_CCL$value)
mean(b1_imp_CCL$value)
sd(b1_imp_CCL$value)
mad(b1_imp_CCL$value)
IQR(b1_imp_CCL$value)

p12 <- b1_impL %>% ggplot(aes(x = logValue)) + geom_density() + theme_bw() +
  ggtitle('Peak Area')
quantile(b1_impL$logValue)
mean(b1_impL$logValue)
sd(b1_impL$logValue)
mad(b1_impL$logValue)
IQR(b1_impL$logValue)


b1_imp_ISL %>% ggplot(aes(x = log(value))) + geom_density() + theme_bw() +
  ggtitle('Relative abundance (IS normalised)')
quantile(b1_imp_ISL$logValue)
mean(b1_imp_ISL$logValue)
sd(b1_imp_ISL$logValue)
mad(b1_imp_ISL$logValue)
IQR(b1_imp_ISL$logValue)


p <- arrangeGrob(p12, p13, p11, nrow = 3, ncol = 1)

ggsave('CC_density.tiff', p, dpi = 300, units = 'mm', width = 100, height = 180)


#
#


# run CPC on two reference datasets (ES and blank)
install.packages("C:/Users/q29634at/OneDrive - The University of Manchester/Documents/R/win-library/4.0/cpca",  
                 repos = NULL, type = "source")

library(cpca)
# create covariance matrices 
# same number of variables required in both ref datasets!
diffnames <- c(setdiff(colnames(assay_blank_imp), colnames(assay_es_imp)),
               setdiff(colnames(assay_es_imp), colnames(assay_blank_imp)))

blank_cpc <- log(assay_blank_imp[, colnames(assay_blank_imp) %ni% diffnames])
es_cpc <- log(assay_es_imp[, colnames(assay_es_imp) %ni% diffnames])

# log transform, scale and center peak area datasets
blank_cpc1 <- scale(blank_cpc, scale = TRUE, center = TRUE)
es_cpc1 <- scale(es_cpc, scale = TRUE, center = TRUE)

cov_blank <- cov(blank_cpc1)
cov_es <- cov(es_cpc1)

ar <- array(c(cov_blank, cov_es), dim = c(225, 225, 2),
            dimnames = list(colnames(cov_blank), colnames(cov_blank), 
                         c('Blank', 'ES')))

ref_cpc <- cpc(ar, method = 'stepwise')

# eigenvalues
# variances associated with the components are allowed to vary across datasets
ev <- ref_cpc$D
ev

# % variance explained (if eigenvalues sum to the number of components)
var_exp <- ev/225*100 

# eigenvectors (scores) - same for both datasets 
library(stringi)

loadings_cpc <- as.data.frame(ref_cpc[['CPC']])
names <- stri_paste('PC', 1:225,sep = '')
colnames(loadings_cpc) <- names
rownames(loadings_cpc) <- colnames(blank_cpc)

# fill in missing loading values with min loading?
setdiff(colnames(Y), rownames(scores_cpc))

mis_loads <- data.frame(comps = c(setdiff(colnames(Y), rownames(scores_cpc))),
                        PC1 = rep(-0.016, times = 10))
rownames(mis_loads) <- mis_loads$comps
mis_loads <- mis_loads %>% dplyr::select(PC1)

loadings_cpcM <- rbind(loadings_cpc %>% dplyr::select(PC1), 
                       mis_loads) %>% as.matrix()

Y <- b1_all_f %>%
  dplyr::select(!c(1:7)) %>% log() %>% as.matrix()

x0 <- Y
x0[is.na(Y)] <- 0 # managment of NAs in matrix multiplication

drift <- (x0 %*% loadings_cpcM) %*% t(loadings_cpcM)
drift[is.na(Y)] <- NA

Z = Y - drift

# overlay with Sample type category
pca_cc1 <- pca(Z, scale = TRUE, center = TRUE)

plotIndiv(pca_cc1,
          group = b1_all_f$Batch,
          pch = 1,
          cex = 2,
          legend = TRUE)


#
#
#


# BATCH EFFECT
# data distribution across batches
bb_dens <- datanF %>%
  pivot_longer(cols = c(1:286), names_to = 'VOC', values_to = 'RA') %>%
  ggplot(aes(x = RA, colour = as.factor(BatchB))) + geom_density() +
  scale_x_log10() +
  labs(col = 'Batch', x = 'log10(relative abundance)') +
  theme_bw() +
  theme(axis.title.x = element_blank())

b_dens <- datanF %>%
  pivot_longer(cols = c(1:286), names_to = 'VOC', values_to = 'RA') %>%
  ggplot(aes(x = RA, colour = as.factor(Batch))) + geom_density() +
  scale_x_log10() +
  labs(col = 'Batch', x = 'log10(relative abundance)') +
  theme_bw() +
  theme(axis.title.x = element_blank())

# sample type effect
s_dens <- datanF %>%
  pivot_longer(cols = c(1:286), names_to = 'VOC', values_to = 'RA') %>%
  ggplot(aes(x = RA, colour = as.factor(Sample_type))) + geom_density() +
  scale_x_log10() +
  labs(col = 'Batch', x = 'log10(relative abundance)') +
  theme_bw() +
  theme(axis.title.x = element_blank())

s_dens

g_dens_b <- arrangeGrob(bb_dens, b_dens, s_dens, ncol = 3,
                        bottom = textGrob('log10(relative abundance)'))
plot(g_dens_b)

ggsave('Batch_effect_density.tiff', g_dens_b, dpi = 300, width = 300, height = 80, 
       unit = 'mm')



#
#
#


# principal component analysis
rownames(dataF) <- dataf$Sample_ID
dataFs <- dataF %>% filter((!grepl('B1', rownames(.)))) 

library(mixOmics)
dataFs <- dataFs[,-156]
dataFs <- dataFs %>% filter(rowSums(is.na(.)) < ncol(dataFs))

anF <- datan %>% filter(Sample_ID %in% rownames(dataFs))
rownames(anF) <- anF$Sample_ID

anF <- anF %>% mutate(ID = Sample_ID) %>% 
separate(ID, into = c('ID1', 'ID2', 'ID3', 'ID4', 'ID5', 'ID6'), sep = '_')


batch_pca <- pca(log10(dataFs), scale = TRUE, center = TRUE,
                 multilevel = anF$ID3)

plotIndiv(batch_pca,
          pch = as.factor(anF$BatchB),
          group = anF$CoreVisit,
          legend = TRUE)

plotIndiv(batch_pca,
          pch = as.factor(anF$BatchB),
          group = anF$Batch,
          legend = TRUE)

#
#
#

# Split dataset along covid batch for train-test datasets
train <- datanF %>% filter(BatchB == 'B1')
test <- datanF %>% filter(BatchB == 'B2')

n_distinct(test$RAD_ID)

trainL <- train %>% pivot_longer(cols = c(1:286), names_to = 'comp', values_to = 'RA') %>%
  pivot_wider(id_cols = c(RAD_ID, CoreVisit, comp), names_from = Sample_type, values_from = RA)

corrs <- trainL %>% dplyr::select(!MaskBG) %>% 
  drop_na() %>%
  group_by(comp) %>% summarise(cor = cor(S1, S2, method = 'pearson'))


#
#
#

# RELATIONSHIP WITH BACKGROUND VALUES
pdf('RADicA_BGvsS_test.pdf')

trainL1 <- trainL %>% pivot_longer(cols = c(S1, S2), 
                                   names_to = 'Breath_sample', 
                                   values_to = 'Breath_RA')

bg_plots <- trainL1 %>%
  group_by(comp) %>%
  do(plots = ggplot(data = .) + 
       (aes(x = log10(MaskBG), y = log10(Breath_RA), colour = Breath_sample)) + 
       geom_point(alpha = 0.7) + ggtitle(.$comp) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) +
       xlab('log10(background)') +
       ylab('log10(sample)'))

marrangeGrob(grobs = bg_plots$plots, nrow = 4, ncol = 4)

dev.off()

# only in Batch 3 
b3 <- datanF %>% filter(Batch == 3)

b3L <- b3 %>% pivot_longer(cols = c(1:286), names_to = 'comp', values_to = 'RA') %>%
  pivot_wider(id_cols = c(RAD_ID, CoreVisit, comp), names_from = Sample_type, values_from = RA)

b3L1 <- b3 %>% pivot_longer(cols = c(S1, S2), 
                                   names_to = 'Breath_sample', 
                                   values_to = 'Breath_RA')


# summarise technical replicates
datanFl$techNA <- ifelse(is.na(datanFl$S1) & is.na(datanFl$S2), 'Missing2', 
                    ifelse(is.na(datanFl$S1) | is.na(datanFl$S2), 'Missing1', 
                           'Non-missing'))

b3L <- datanFl %>% filter(Batch == 3)


b3L <- b3L %>% mutate(S1 = ifelse(is.na(S1) == TRUE & techNA == 'Missing1', 0, S1),
                            S2 = ifelse(is.na(S2) == TRUE & techNA == 'Missing1', 0, S2),
                            S = (S1+S2)/2,
                            log_S = log10(S),
                            MaskBG = ifelse(MaskBG == 0, NA, MaskBG),
                            log_MaskBG = log10(MaskBG))


nasums <- b3L %>% group_by(comp) %>% summarise(nona = sum(is.na(log_S) == TRUE)/60)
nasums1 <- nasums %>% filter(nona < 0.3)
b3L <- b3L %>% filter(comp %in% nasums1$comp) 

pdf('RADicA_BGvsS_B3_1.pdf')

b3_plots <- b3L %>%
  group_by(comp) %>%
  do(plots = ggplot(data = .) + 
       (aes(x = log_MaskBG, y = log_S)) + 
       geom_point(alpha = 0.7) + ggtitle(.$comp) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) +
       xlab('log10(background)') +
       ylab('log10(sample)'))


marrangeGrob(grobs = b3_plots$plots, nrow = 4, ncol = 4)
dev.off()

# plotting S-MaskBG relationship without the outlier observations
trainLout <- trainL %>% group_by(comp) %>% drop_na() %>%
  summarise(mean = mean(log_S),
  sd = sd(log_S))

trainL <- trainL %>% left_join(trainLout)

trainLoutF <- trainL %>%
  filter(log_S < mean + 2.5*sd & log_S > mean - 2.5*sd)

n_distinct(trainL$comp)
n_distinct(trainLoutF$comp)

pdf('RADicA_BGvsS_test1_Out.pdf')

bg_plots <- trainLoutF %>%
  group_by(comp) %>%
  do(plots = ggplot(data = .) + 
       (aes(x = log_MaskBG, y = log_S)) + 
       geom_point(alpha = 0.7) + ggtitle(.$comp) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) +
       xlab('log10(background)') +
       ylab('log10(sample)'))

marrangeGrob(grobs = bg_plots$plots, nrow = 4, ncol = 4)

dev.off()

# plot histograms for each compound
pdf('RADicA_S_hist.pdf')

bg_plots <- trainL %>%
  group_by(comp) %>%
  do(plots = ggplot(data = .) + 
       (aes(x = S)) + 
       geom_histogram() + 
       ggtitle(.$comp) + 
       theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) +
       ylab('log10(sample)'))

marrangeGrob(grobs = bg_plots$plots, nrow = 4, ncol = 4)

dev.off()


library(lme4)
library(lmerTest)
library(nlme)
library(MASS)
library(foreign)
library(sfsmisc)
library(robustlmm)
library(confintROB)

trainL <- b3L

# test effect of outliers on regression model fit
subset <- trainL %>% filter(comp == 'Acetophenone') %>%
  drop_na(log_S) %>% drop_na(log_MaskBG)

mod <- lmer(log_S ~ log_MaskBG + CoreVisit + (1 | RAD_ID), #+ Batch
          data = subset)

mod2 <- lm(log_S ~ log_MaskBG + CoreVisit,
           data = subset) 

summary(mod2)
opar <- par(mfrow = c(2,2), oma = c(0, 0, 1.1, 0))
plot(mod2, las = 1)



# Gamma regression
modG <- glm(S ~ MaskBG + CoreVisit, #+ Batch,
            data = subset,
            family = 'Gamma'(link = 'identity'))

summary(modG)

# view observations with highest Cook's distance
d1 <- cooks.distance(mod)
r <- stdres(mod)
a <- cbind(subset, d1, r)
a[d1 > 4/58, ]

# view observations with the highest absolute residual values
rabs <- abs(r)
a <- cbind(subset, d1, r, rabs)
asorted <- a[order(-rabs), ]
asorted[1:10, ]

# run robust regression
rmod <- rlm(log_S ~ log_MaskBG + CoreVisit, #+ Batch,
            data = subset)

summary(rmod)
f.robftest(rmod, var = "log_MaskBG")$p.value


plot(rmod)
# mixed-effect robust regression
rmod2 <- rlmer(log_S ~ log_MaskBG + CoreVisit  + (1 | RAD_ID) + Batch,
               data = subset)

summary(rmod2)
coef(summary(rmod2))
plot(rmod2)

# obtain confidence intervals through bootstrap
set.seed(123)
wild.rmod2 <- confintROB(object = rmod2,
                         boot.type = "wild",
                         nsim = 100)

wild.rmod2[1:5,]

#
#
#

# function applying robust regression to vocs (mixed-effect)
library(parallel)

reg_coefs <- function(voc){
  subset <- trainL %>%
    filter(comp == voc) %>%
    drop_na()
  
  model <- rlmer(log_S ~ log_MaskBG + CoreVisit + (1 | RAD_ID), #+Batch
              data = subset)
  
  coefs <- as.data.frame(coef(summary(model)))
  wild.rmod2 <- confintROB(object = model,
                           boot.type = "wild",
                           nsim = 100)
  coefs <- cbind(coefs, wild.rmod2[1:5,])
  
}

# robust regression fixed effects only
reg_coefs <- function(voc){
  subset <- trainL %>%
    filter(comp == voc) %>%
    drop_na()
  
  model <- rlm(log_S ~ log_MaskBG + CoreVisit, #+ Batch,
               data = subset,
               psi = psi.bisquare)
  
  coefs <- as.data.frame(coef(summary(model))) %>%
                           filter(rownames(.) == 'log_MaskBG')
  coefs <- cbind(coefs, f.robftest(model, var = "log_MaskBG")$p.value)
  
}

# Gamma regression
reg_coefs <- function(voc){
  subset <- trainL %>%
    filter(comp == voc) %>%
    drop_na()
  
  model <- glm(S ~ MaskBG + CoreVisit, #+ Batch,
               data = subset,
               family = 'Gamma'(link = 'inverse'))
  
  coefs <- as.data.frame(coef(summary(model))) %>%
    filter(rownames(.) == 'MaskBG')
  
}

#

phenol <- reg_coefs('Phenol')

# extract list of target protein names
vocs <- unique(trainL$comp)

# apply function to the list of target proteins
coefs_bg_G <- lapply(vocs, reg_coefs)

# transform the list of data frames to a data frame
coefs_bg1 <- bind_rows(coefs_bg1, .id = 'voc') %>%
  mutate(voc = rep(vocs, each = 1)) %>%
  remove_rownames() #%>%
  #rename(p.value = 'p-value')
  rename(p.value = 'Pr(>|t|)')

colnames(coefs_bg1)[5] <- 'p.value'

# apply Benjamini-Hochberg correction for multiple comparisons
# exclude p-values associated with intercept
coefs_bg2 <- coefs_bg1 %>% mutate(adjusted_p.value = p.adjust(p.value, method = 'BH'))
#coefs_bg3 <- coefs_bg2 %>% filter(predictor == 'log_MaskBG')

coefs_bg2 <- coefs_bg2 %>% 
  mutate(Adjust = ifelse(adjusted_p.value < 0.01, 'Yes', 'No'))

# adjust log_S values according to the regression coefficient
b3L <- b3L %>% left_join(coefs_bg2 %>% dplyr::select(voc, Value, Adjust),
                         by = c('comp' = 'voc'))

b3L <- b3L %>% mutate(adj_log_S = ifelse(Adjust == 'Yes', log_S - (Value * log_MaskBG),
                                         log_S))

write.csv(b3L, 'Batch3_BGadjusted.csv')


pdf('RADicA_BGvsS_B3_adj.pdf')

b3_plots <- b3L %>%
  group_by(comp) %>%
  do(plots = ggplot(data = .) + 
       (aes(x = log_MaskBG, y = adj_log_S, colour = Adjust)) + 
       geom_point(alpha = 0.7) + 
       scale_colour_manual(values = c('Yes' = 'blue', 
                                      'No' = 'black')) +
       ggtitle(.$comp) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) +
       xlab('log10(background)') +
       ylab('log10(sample)'))


marrangeGrob(grobs = b3_plots$plots, nrow = 4, ncol = 4)
dev.off()

dev.new()

b3L %>% filter(Adjust == 'Yes') %>%
  pivot_longer(cols = c(log_S, adj_log_S), names_to = 'measure', values_to = 'value') %>%
  ggplot(aes(x = measure, y = value, fill = measure)) + geom_boxplot(outliers = FALSE)


## MISSING VALUE IMPUTATION
