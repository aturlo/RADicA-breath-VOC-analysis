## Radica batch exploratory analysis
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(gridExtra)

batch1 <- read.csv('231208_Radica Batch 1_all samples_aligned S1 to Mask BG_updated 240125.csv')
batch2 <- read.csv('240105_Radica Batch 2_all samples_aligned S1 to MaskBG.csv')

list <- read.csv('Radica v4 VOC list.csv')
View(as.data.frame(unique(list$File.Name)))

# number of patients in each batch
n_distinct(na.omit(batch1$RAD_ID))
n_distinct(na.omit(batch2$RAD_ID))

# batch duiagnosis breakdown
batch1 %>% group_by(Diagnosis) %>% summarise(n = n_distinct(RAD_ID))
batch2 %>% group_by(Diagnosis) %>% summarise(n = n_distinct(RAD_ID))

# batch diagnosis and visit breakdown
n1 <- batch1 %>% group_by(Diagnosis, CoreVisit) %>% summarise(n = n_distinct(RAD_ID)) 
n2 <- batch2 %>% group_by(Diagnosis, CoreVisit) %>% summarise(n = n_distinct(RAD_ID))

cvsplit <- batch1 %>% drop_na(RAD_ID) %>% group_by(RAD_ID) %>% summarise(n = n()) %>% 
  group_by(n) %>% summarise(count = n())

cvsplit1 <- batch1 %>% drop_na(RAD_ID) %>% group_by(Diagnosis, RAD_ID) %>% summarise(n = n())

cvsplit2 <- batch1 %>% group_by(Diagnosis, RAD_ID, CoreVisit) %>% summarise(n = n())

View(cvsplit1 %>% filter(n == 9))


cv1 <- unique(batch2 %>% filter(CoreVisit == 'CV1' & Diagnosis == 'Asthma') %>% select(RAD_ID))
cv2 <- unique(batch2 %>% filter(CoreVisit == 'CV2' & Diagnosis == 'Asthma') %>% select(RAD_ID))
cv12 <- unique(batch2 %>% filter(CoreVisit %in% c('CV1', 'CV2') & Diagnosis == 'Asthma') %>% select(RAD_ID))

setdiff(cv12, cv1)
setdiff(cv12, cv2)

# extract intensity data
assay1 <- batch1[, 14:ncol(batch1)]
assay2 <- batch2[, 14:ncol(batch1)]

rownames(assay1) <- batch1$Sample_ID
rownames(assay2) <- batch2$Sample_ID

# duplicated sample_ID in batch 2
table(duplicated(batch2$Sample_ID))
batch2 %>% select(Sample_ID) %>% filter(duplicated(.))

# compare summary statistics between batches
assay1l <- assay1 %>% pivot_longer(cols = colnames(assay1), names_to = 'compound', values_to = 'intensity')
assay2l <- assay2 %>% pivot_longer(cols = colnames(assay2), names_to = 'compound', values_to = 'intensity')

assay1l <- assay1l %>% separate(compound, c('compound', 'RT'), sep = '_')
assay2l <- assay2l %>% separate(compound, c('compound', 'RT'), sep = '_')

med1 <- assay1l %>% group_by(compound) %>% summarise(median = median(na.omit(intensity))) %>% arrange(median)
med2 <- assay2l %>% group_by(compound) %>% summarise(median = median(na.omit(intensity))) %>% arrange(median)

quantile(na.omit(assay1l$intensity))
quantile(na.omit(assay2l$intensity))

quantile(med1$median)
quantile(med2$median)

med1$batch <- 'B1'
med2$batch  <- 'B2'

rbind(med1, med2) %>% ggplot(aes(x = batch, y = median, fill = batch)) + geom_boxplot() + 
  scale_y_continuous(trans = 'log2') + ylab('Log2_intensity') +
  theme_bw() + scale_fill_brewer(palette = 'Dark2') +
  ggtitle('Median compound intensity values between batches')

assay %>% ggplot(aes(x = batch, y = log2(intensity), fill = batch)) + geom_boxplot() +
  theme_bw(base_size =12) + scale_fill_brewer(palette = 'Dark2') + ylab('Log2_intensity') + xlab('Batch') +
  theme(legend.position = 'none')

assay %>% ggplot(aes(x = log2(intensity), colour = batch)) + geom_density() +
  theme_bw(base_size =12) + scale_colour_brewer(palette = 'Dark2') + xlab('Log2_intensity')

ggsave('Density_batches.tiff', width = 100, height = 67, units = 'mm', dpi = 300)

# group medians into bins of 100
bins <- sapply(1:20, function(x){
  rep(paste('b', x, sep =''), 51)
})

med1$bin <- bins[-c(1006:1020)]
med2$bin <- bins[-c(1006:1020)]

med <- rbind(med1, med2)

assay1l <- assay1l %>% full_join(med1)
assay2l <- assay2l %>% full_join(med2)

med %>% filter(bin == 'b7') %>% ggplot(aes(x = compound, y = log2(median), colour = batch)) + geom_point() +
  theme_bw(base_size = 12) + ylab('Log2_median_intensity') + xlab('compounds') + scale_colour_brewer(palette = 'Dark2') +
  theme(axis.text.x = element_blank())

ggsave('Medians_batches_b7.tiff', width = 100, height = 50, units = 'mm', dpi = 300)

pdf('Batch_medians.pdf')
box <- med %>% group_by(bin) %>% 
  do(plots = ggplot(data = .) + (aes(x = compound, y = log2(median), colour = batch)) + geom_point() +
       theme_bw() + theme(axis.text.x = element_blank()) + ggtitle(.$bin) + scale_colour_brewer(palette = 'Dark2') +
       ylab('Log2_intensity') + xlab('compounds'))

marrangeGrob(box$plots, nrow = 4, ncol =1)
dev.off()


# correlation plots for compounds close to each other in data frame
comps <- as.data.frame(colnames(batch2))
View(comps)

dm <- batch2 %>% ggplot(aes(x = X.2E.4E..3.7.Dimethylocta.2.4.diene_9.089, 
                      y = Cyclohexene..1.methyl.5..1.methylethenyl.....R.._9.084)) +
  geom_point() + theme_bw() +
  ylab('Monoterpene_9.084') + xlab('Diene_9.089') + coord_fixed()

dl <- batch2 %>% ggplot(aes(x = X.2E.4E..3.7.Dimethylocta.2.4.diene_9.089, 
                            y = D.Limonene_9.198)) +
  geom_point() + theme_bw() +
  ylab('Limonene_9.198') + xlab('Diene_9.089')
dl

dc <- batch2 %>% ggplot(aes(x = X.2E.4E..3.7.Dimethylocta.2.4.diene_9.089, 
                            y = X....3.Carene_8.643)) +
  geom_point() + theme_bw() +
  ylab('Carene_8.643') + xlab('Diene_9.089')
dc

db <- batch2 %>% ggplot(aes(x = Benzene..1.bromo.3.fluoro._6.615, 
                            y = X....3.Carene_8.643)) +
  geom_point() + theme_bw() +
  ylab('Bromofluorobenzene_6.615') + xlab('Diene_9.089')
db

d <- arrangeGrob(dm, dl, dc, db, nrow = 2, ncol = 2)
plot(d)
ggsave('VOCCorrRT.tiff', d, dpi = 300, width = 180, height = 120, units = 'mm')

# RI duplication plot
ri <- read.csv('RIs for Radica.csv')
View(as.data.frame(colnames(ri)))

ri1 <- ri %>% select(Sample, X.2E.4E..3.7.Dimethylocta.2.4.diene.Results, D.Limonene.Results)
ri1l <- ri1 %>% pivot_longer(cols = !Sample, names_to = 'compound', values_to = 'RI')
ri1l$RI <- as.numeric(ri1l$RI)

ri1l %>% ggplot(aes(x = Sample, y = RI, colour = compound)) +
  geom_point(alpha = 0.5) +  theme_classic() + ylim(1017, 1037) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_colour_brewer(palette = 'Dark2', labels = c('Limonene', 'Diene')) +
  xlab('Samples')

ggsave('RILimvsDiene.tiff', dpi = 300, units = 'mm', width = 200, height = 120) 
