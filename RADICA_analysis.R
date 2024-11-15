## RADICA study exploratory analysis

# title:
# author:
# output:

install.packages('dplyr', dependencies = TRUE)
install.packages('lifecycle', dependencies = TRUE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)

data <- read.csv('231208_Radica Batch 1_all samples_aligned S1 to Mask BG_updated 240125.csv')

# unfiltered dataset
assay <- data[,14:ncol(data)]
rownames(assay) <- data$Sample_ID


# filter first visit only (diagnosis)
data <- data %>% filter(CoreVisit == 'CV1')


# MF70 filtered dataset
mf_70w <- read.csv('MF70_intensities_Batch1.csv')
assay <- mf_70w
assay <- assay[,-1]
rownames(assay) <- mf_70w$Sample

# MF70 <30% NA filtered dataset
mf_70f <- read.csv('MF70_intensities_Batch1_NAfiltered.csv')
assay <- mf_70f
assay <- assay[,-1]
rownames(assay) <- mf_70f$X

# ! samples not included in individual sample MF score df
dif <- data %>% filter(Sample_ID %ni% mf2$Sample)

# add missing background values to RAD032 (copy of RAD031)
b32 <- assay[43,]
b32$Sample <- '190903_RaDICA_RAD032_B1_367147_1'
assay <- rbind(assay, b32)
assay <- assay[c(1:45, 171, 46:170), ]
assay <- select(assay, !X)

data$Sample_ID[52] <- '190903_RaDICA_RAD032_B1_367147_1'

# create metdata df
meta <- data[,1:13]

#####################################

## overview of filtered data (MF70 <30% NA)
assays <- assay %>% mutate(Sample_ID = rownames(assay)) %>%
  filter(Sample_ID %ni% bnames)

assaysan <- assays %>% left_join(meta %>% select(Sample_ID, Diagnosis, CoreVisit)) %>%
  relocate(Sample_ID, Diagnosis, CoreVisit) #%>%
  pivot_longer(cols = !c(Sample_ID, Diagnosis, CoreVisit), names_to = 'compound', values_to = 'intensity')

pdf('MF70filt_scatterplots.pdf')
scats <- assaysan %>% group_by(compound) %>%
  do(plots = ggplot(data =.) + 
       aes(x = CoreVisit, y = log2(intensity), colour = Diagnosis) + 
       geom_point(position = position_jitterdodge()) + theme_bw(base_size = 8) +
       theme(legend.position = 'none') +
       scale_color_brewer(palette = 'Dark2') +
       ggtitle(.$compound))

marrangeGrob(scats$plots, nrow = 3, ncol = 3)
dev.off()

## Missing data
assayt <- as.data.frame(t(assay))
assayt <- assayt[-1,]
assayt <- assayt %>% mutate_if(is.character, as.numeric)



# plot frequency of NAs across features

cnas <- function(x){
  fnas <- as.data.frame(table(rowSums(is.na(x))/ncol(x)*100))
  fnas$Var1 <- as.numeric(levels(fnas$Var1))[fnas$Var1]
  p2 <- fnas %>% ggplot(aes(y = Freq, x = Var1)) + 
  geom_bar(stat = 'identity') +
  theme_bw(base_size = 12) +
  ylab('Frequency') +
  xlab('% Missing') +
  geom_vline(xintercept = 30, colour = 'red')
  p2}

cnas(assayt) + ggtitle('Missing rate across compounds') + theme(plot.title = element_text(hjust = 0.5))
ggsave('Missing_rate_compounds.tiff', width = 100, height = 80, units = 'mm', res = 300)

# NAs across features in different samples types (BG, S1, S2)
bnames <- colnames(assayt)[grepl('B', colnames(assayt))]
s1names <- colnames(assayt)[grepl('S1', colnames(assayt))]
s2names <- colnames(assayt)[grepl('S2', colnames(assayt))]

bg <- assayt[,bnames]
s1 <- assayt[,s1names]
s2 <- assayt[,s2names]

cnas(bg)
cnas(s1)
cnas(s2)

# filter features with < 30% NAs
allf <- assayt %>% filter(rowSums(is.na(assayt))/ncol(assayt) < 0.3)

bgf <- bg %>% filter(rowSums(is.na(bg))/ncol(bg) < 0.3)
s1f <- s1 %>% filter(rowSums(is.na(s1))/ncol(s1) < 0.3)
s2f <- s2 %>% filter(rowSums(is.na(s1))/ncol(s1) < 0.3)

# compare features retained across sample types
allnames <- data.frame(comp = rownames(allf), All = rep(1, nrow(allf)))

bgnames <- data.frame(comp = rownames(bgf), BG = rep(1, nrow(bgf)))
s1fnames <- data.frame(comp = rownames(s1f), S1 = rep(1, nrow(s1f)))
s2fnames <- data.frame(comp = rownames(s2f), S2 = rep(1, nrow(s2f)))

af <- full_join(bgnames, s1fnames) %>% full_join(s2fnames) %>% full_join(allnames)
af[is.na(af)] <- 0
rownames(af) <- af$comp
af <- af[,-1]

tiff(filename = 'Heatmap_Missing_subsets.tiff', width = 180, height = 160, unit = 'mm', res = 300)
heatmap <- pheatmap(as.matrix(af), show_rownames = TRUE, cluster_cols =  F, 
                    cluster_rows = T, show_colnames = TRUE, fontsize_row = 6,
                    main = 'Compounds with >30% missing values across data subsets (blue)',
                    legend  = F)
heatmap
dev.off()

# retain features with < 30% NAs in S1 and S2 samples
assaytf <- assayt[rownames(assayt) %in% s1fnames$comp,]

assayf <- as.data.frame(t(assaytf))

write.csv(assayf, 'MF70_intensities_Batch1_NAfiltered.csv')

# NA across samples (only S1 and S2 values)
assaytf <- assayt
'%ni%' <- Negate('%in%')

assaytfs <- assaytf[colnames(assaytf) %ni% bnames]

assayf1 <- assayf %>% 
  mutate(Sample_ID = assay$Sample) %>% 
  left_join(meta[,c('Sample_ID', 'RAD_ID', 'Sample', 'Diagnosis', 'Analysis_date','Visit', 'Sampling_date')]) %>%
  relocate(Sample_ID, RAD_ID, Sample, Diagnosis, Analysis_date)

assayfs <- assayf1 %>% filter(Sample != 'MaskBG')
assayfs <- assayfs %>% mutate(nas = rowSums(is.na(assayfs[,6:73]))/nrow(assayfs)*100) %>%
  relocate(nas)

library(Polychrome)

p33 <- createPalette(33,  c("#ff0000", "#00ff00", "#0000ff"))
p33 <- sortByHue(p33)
p33 <- as.vector(t(matrix(p33, ncol=3)))
names(p33) <- NULL

assayfs$RAD_ID <- as.factor(assayfs$RAD_ID)

dev.new()
p0 <- assayfs %>% arrange(nas) %>% ggplot(aes(y = nas, x = RAD_ID, group = Visit)) + 
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6)) #+
  scale_fill_manual(values = p33)
  
p0

                  
snas <- as.data.frame(table(colSums(is.na(assaytfs))/nrow(assaytfs)*100))
snas$Var1 <- as.numeric(levels(snas$Var1))[snas$Var1]


p1 <- snas %>% ggplot(aes(y = Freq, x = Var1)) + 
  geom_bar(stat = 'identity') +
  ylab('Frequency') +
  xlab('% NAs') +
  ggtitle('Histogram of samples acording to %NAs')
p1


# Heatmap of NAs according to condition
rownames(meta) <- meta$Sample_ID
meta_heat <- meta %>% select(!c('RAD_ID', 'CoreVisit', 'Mask.BG', 'S1', 'S2', 'Tube.BG', 'Sample_ID'))

assayf <- as.data.frame(t(assaytf)) # filtered < 30% NAs

data_heat <- assayf
rownames(data_heat) <- meta$Sample_ID
rownames(data_heat) <- data_heat$Sample
data_heat <- data_heat[,-1]

data_heat <- assay[,-1] # non-filtered

'%ni%' <- Negate('%in%')

data_heat <- data_heat %>% filter(rownames(.) %ni% bnames) # exclude BG samples
meta_heat <- meta_heat %>% filter(Sample != 'MaskBG') 


data_heat[is.na(data_heat)] <- 0
data_heat[data_heat > 0] <- 1

library(RColorBrewer)

newCols <- colorRampPalette(grDevices::rainbow(length(unique(meta_heat$Sampling_date))))
annoCol <- newCols(length(unique(meta_heat$Sampling_date)))
names(annoCol) <- unique(meta_heat$Sampling_date)
annotCol <- list(Sampling_date = annoCol)

meta_heat$Visit <- as.factor(meta_heat$Visit)

tiff(filename = 'Heatmap_Missing_S.tiff', width = 180, height = 146, unit = 'mm', res = 300)

heatmap <- pheatmap(as.matrix(t(data_heat)), annotation_col = meta_heat %>% select(Diagnosis, age, Visit), 
                    show_rownames = FALSE, cluster_cols =  T, show_colnames = FALSE, cluster_rows = FALSE,
                    fontsize = 10, main = 'Missing values across compounds (breath only)')  #, annotation_colors = annotCol) # use for analysis date mapping
heatmap
dev.off()

# NAs according to retention time (only S1 and S2 values)
comps <- read.csv('Compound_names_RT.csv')

assayt1 <- assaytfs %>% mutate(compound = rownames(assaytfs))
assayt1$compound <- str_sub(assayt1$compound, end = -2)
assayt1 <- assayt1 %>% left_join(comps %>% select(compound, RT), by = 'compound')

assayt1$RT <- as.numeric(assayt1$RT)
assayt1$nonas <- rowSums(is.na(assayt1[,-(315:316)]))

dev.new()
assayt1 %>% ggplot(aes(x = RT, y = nonas/ncol(assayt1)*100)) + 
  geom_point() +
  ylab('% samples with NA')

# NAs according to peak intensity - breath samples
assaytna <- assayt1 %>%
  filter(nonas > 0)

assaytl <- assaytna %>% pivot_longer(cols = !c(compound, RT, nonas),
                                    names_to = 'sample',
                                    values_to = 'intensity') 

assaytl$intensity <- as.numeric(assaytl$intensity)

sumint <- assaytl %>% select(!RT) %>% 
  group_by(compound) %>%
  drop_na() %>%
  summarise(medint = median(intensity),
            nonas = mean(nonas)) %>%
  arrange(medint) %>%
  mutate(rankmed = row_number())

  
sumint %>% ggplot(aes(x = medint)) + geom_histogram(bins = 20)

sumint %>% 
  ggplot(aes(x = rankmed, y = nonas/314*100)) + geom_point() +
  ylab('% Missing values') + xlab('Rank of median intensity') + theme_bw(base_size = 13) +
  ggtitle('Missing values depending on median compound intensity (breath)') +
  theme(plot.title = element_text(hjust = 0.5, size = 13))

ggsave('NAvsMedianIntensity.tiff', width = 180, height = 100, unit = 'mm', dpi = 300)

# in background samples
assayb <- assay %>% filter(rownames(assay) %in% bnames)
assayb <- as.data.frame(t(assayb))
assayb <- assayb %>% mutate(nonas = rowSums(is.na(assayb))) %>% relocate(nonas)

assaybl <- assayb %>% mutate(compound = rownames(assayb)) %>% 
  pivot_longer(cols = !c(nonas, compound),
               names_to = 'sample',
               values_to = 'intensity') 

sumintb <- assaybl %>% 
  group_by(compound) %>%
  drop_na() %>%
  summarise(medint = median(intensity),
            nonas = mean(nonas)) %>%
  arrange(medint) %>%
  mutate(rankmed = row_number())

sumintb %>% 
  ggplot(aes(x = rankmed, y = nonas/156*100)) + geom_point() +
  ylab('% Missing values') + xlab('Rank of median intensity') + theme_bw(base_size = 13) +
  ggtitle('Missing values depending on median compound intensity (background)') +
  theme(plot.title = element_text(hjust = 0.5, size = 13)) + geom_hline(yintercept = 30, colour = 'red')

ggsave('NAvsMedianIntensityBG.tiff', width = 180, height = 100, unit = 'mm', dpi = 300)

#  NA patterns across technical replicates
assayfl1$techNA <- ifelse(is.na(assayfl1$S1) | is.na(assayfl1$S2), 'Missing', 'Non-missing')

assayfl2 <- assayfl1 %>% select(RAD_ID, compound, S1, S2, techNA, CoreVisit) %>%
  pivot_longer(cols = c(S1, S2), names_to = 'sample', values_to = 'intensity')

# per sample, across compounds
pdf('Missing_patterns_sample.pdf')
pnas <- assayfl2 %>% group_by(RAD_ID) %>%
  do(plots = ggplot(data =.) + 
       aes(x = techNA, y = log2(intensity), colour = sample) + 
  geom_boxplot(width = 0.1) + geom_violin() + theme_bw(base_size = 8) +
    theme(legend.position = 'none') +
    scale_color_brewer(palette = 'Dark2') +
  ggtitle(.$RAD_ID))

marrangeGrob(pnas$plots, nrow = 5, ncol = 4)
dev.off()

table(rowSums(is.na(assaytf)))

assayfl2 %>% filter(compound == 'Acetic.acid.') %>%
  ggplot(aes(x = techNA, y = log2(intensity))) + 
  geom_violin(aes(colour = sample), position = position_dodge(0.9)) + 
  geom_boxplot(aes(colour = sample), position = position_dodge(0.9), width = 0.1) + 
  theme_bw(base_size = 8) +
       theme(legend.position = 'none') +
       scale_color_brewer(palette = 'Dark2')

# per compound, across samples
pdf('Missing_patterns_compound.pdf')
pnac <- assayfl2 %>% group_by(compound) %>%
  do(plots = ggplot(data =.) + 
       aes(x = techNA, y = log2(intensity), colour = techNA) +
       geom_violin(position = position_dodge(0.9)) + 
       geom_boxplot(position = position_dodge(0.9), width = 0.1) +
       theme_bw(base_size = 8) +
       theme(legend.position = 'none',
             title = element_text(size = 6)) +
       scale_color_brewer(palette = 'Dark2') +
       ggtitle(.$compound))

marrangeGrob(pnac$plots, nrow = 5, ncol = 4)
dev.off()

# Complete missing breath measurements (S1 & S2 NA)
fullna <- assayfl1 %>% filter(is.na(S1) == TRUE & is.na(S2) == TRUE)
nasum <- fullna %>% group_by(compound) %>% summarise(nonafull = n())
nasum$compound <- str_sub(nasum$compound, end = -2)
fullnasum <- nasum %>% full_join(sumint) %>% pivot_longer(cols = c(nonas, nonafull), 
                                                          names_to = 'natype', values_to =  'nacount')

fullnasum %>% 
  ggplot(aes(x = rankmed, y = nacount/314*100, colour = natype)) + geom_point()

nasumcv <- fullna %>% group_by(compound, CoreVisit) %>% summarise(nasum = n())

# correlation matrix of compounds
library(ggcorrplot)

cormat <- cor(assays %>% select(!Sample_ID), use = 'pairwise.complete.obs', method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)

tiff(filename = 'CorrPlot_S.tiff', width = 100, height = 100, unit = 'mm', res = 300)
ggcorrplot(cormat, tl.cex = 4, hc.order = TRUE, legend.title = 'r', title = 'Pearsons correlation matrix of breath compounds')
dev.off()


################################



# Imputation of missing data


# NA imputation with GMSimpute
library(GMSimpute)

assaytfs_imp <- GMS.Lasso(assaytfs, alpha = 1, nfolds = 10, log.scale = TRUE, TS.Lasso = TRUE)

cormatI <- cor(t(assaytfs_imp1), method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)
ggcorrplot(cormatI, tl.cex = 6) #hc.order = TRUE)

# check assumption that log intensities follow multivariate normal distribution
x <- log10(assays_min)
cm <- colMeans(x)
S <- cov(x)
d <- apply(x, 1, function(x) t(x-cm) %*% solve(S) %*% (x - cm))

plot(qc <- qchisq((1:nrow(x) - 1/2) / nrow(x), df = 69), sort(d))
abline(a = 0, b =1)

# replace NAs with 0 (log transform first)
assaytfs_0 <- log10(assaytfs)
assaytfs_0[is.na(assaytfs_0)] <- 0

# replace NAs with compound minimal value
assays_min <- as.data.frame(t(assaytfs))

namin <- function(x) {
  assays_min %>% select(x) %>% mutate(x = ifelse(is.na(.), min(na.omit(.)), get(x)))}

sam <- colnames(assays_min)

assays_min <- as.data.frame(lapply(sam, namin))
  
assays_min <- assays_min[, seq(2, ncol(assays_min), 2)]
colnames(assays_min) <- sam

# Knn imputation (k = 10)
BiocManager::install("MSnbase")
BiocManager::install("MsCoreUtils")
library(MSnbase)
library(MsCoreUtils)

assaytfsK <- assaytfs
assaytfsK <- impute_matrix(as.matrix(assaytfsK), method = 'knn', k = 10)


# evaluate effect of different imputation methods
library(factoextra)
library(patchwork)

pcaL <- prcomp(t(log10(assaytfs_imp)), scale = TRUE)
pcamin <- prcomp(log10(assays_min), scale = TRUE)
pcaK <- prcomp(t(log10(assaytfsK)), scale = TRUE)

pcaL_scores <- as.data.frame(pcaL$x)
pcaL_load <- as.data.frame(pcaL$rotation)

pcamin_scores <- as.data.frame(pcamin$x)
pcamin_load <- as.data.frame(pcamin$rotation)

pcaK_scores <- as.data.frame(pcaK$x)
pcaK_load <- as.data.frame(pcaK$rotation)

metas <- meta %>% select(Diagnosis, Sample_ID, age, CoreVisit, Sampling_date, Analysis_date) %>%
  filter(Sample_ID %in% rownames(assays))


metas <- metas %>% separate(Analysis_date, c('Day', 'Month', 'Year'), sep = '/') %>%
  mutate(DateA = paste(Year, Month, Day)) %>%
  separate(Sampling_date, c('Day', 'Month', 'Year'), sep = '/') %>%
  mutate(DateS = paste(Year, Month, Day)) 

p2 <- ggplot(pcaL_scores, aes(x = PC1, y = PC2, colour = metas$DateS, shape = metas$Diagnosis)) + 
  geom_point() + 
  xlab(paste0('PC1: ', round(as.numeric(summary(pcaL)$importance[2,1]*100)), '% expl.var')) + 
  ylab(paste0('PC2: ', round(as.numeric(summary(pcaL)$importance[2,2]*100)), '% expl.var')) + 
  #scale_colour_brewer(palette = "Dark2") +
  theme_bw() + 
  ggtitle("Lasso-TS imputation")

p2

p3 <- ggplot(pcamin_scores, aes(x = PC1, y = PC2, colour = metas$DateA, shape = metas$Diagnosis)) + 
  geom_point() + 
  xlab(paste0('PC1: ', round(as.numeric(summary(pcamin)$importance[2,1]*100)), '% expl.var')) + 
  ylab(paste0('PC2: ', round(as.numeric(summary(pcamin)$importance[2,2]*100)), '% expl.var')) + 
  #scale_colour_brewer(palette = "Dark2") +
  theme_bw() + 
  ggtitle("Min compound value imputation")

p3

p4 <- ggplot(pcaK_scores, aes(x = PC1, y = PC2, colour = metas$CoreVisit, shape = metas$Diagnosis)) + 
  geom_point() + 
  xlab(paste0('PC1: ', round(as.numeric(summary(pcamin)$importance[2,1]*100)), '% expl.var')) + 
  ylab(paste0('PC2: ', round(as.numeric(summary(pcamin)$importance[2,2]*100)), '% expl.var')) + 
  #scale_colour_brewer(palette = "Dark2") +
  theme_bw() + 
  ggtitle("Knn imputation (k = 10)")

p4

grid.arrange(p2 + theme(legend.position = 'none'),
             p3 + theme(legend.position = 'none'), 
             p4 + theme(legend.position = 'none'), 
             nrow = 2)

eigval0 <- get_eigenvalue(pca0)
eigvalL <- get_eigenvalue(pcaL)
eigvalmin <- get_eigenvalue(pcamin)
eigvalK <- get_eigenvalue(pcaK)

p5 <- eigval0[1:10,] %>% ggplot(aes(x = fct_inorder(rownames(.)), y = variance.percent)) + 
  geom_bar(stat = 'identity') + ggtitle('Zero imputation')

p6 <- eigvalL[1:10,] %>% ggplot(aes(x = fct_inorder(rownames(.)), y = variance.percent)) + 
  geom_bar(stat = 'identity') + ggtitle('Lasso-TS imputation')

p7 <- eigvalmin[1:10,] %>% ggplot(aes(x = fct_inorder(rownames(.)), y = variance.percent)) + 
  geom_bar(stat = 'identity') + ggtitle('Min compound imputation')

p8 <- eigvalK[1:10,] %>% ggplot(aes(x = fct_inorder(rownames(.)), y = variance.percent)) + 
  geom_bar(stat = 'identity') + ggtitle('Knn imputation (k = 10)')


grid.arrange(p5, p6, p7, p8, nrow = 2)

################################



# Batch effects

# Multidimensional Scaling
library(vsn)
vsn <- vsn::vsnMatrix(as.matrix(assaytfs))
vsnout <- vsn@hx
vsnout <- as.data.frame(t(vsnout))

assays <- assays %>% select(!Sample_ID)
assaysL <- log10(assays_min)


Eudist <- dist(vsnout)
Mandist <- dist(assays, 'manhattan')

EuMds <- cmdscale(Eudist, k = 20, eig = TRUE)
ManMds <- cmdscale(Mandist, k = 20, eig = TRUE)

barplot(EuMds$eig[1:10])
barplot(ManMds$eig[1:10])

round(100*sum(EuMds$eig[1]/sum(EuMds$eig)))
round(100*sum(ManMds$eig[1]/sum(ManMds$eig)))

dfEu <- as_tibble(EuMds$points[,1:2]) %>% setNames(paste0('MDS', 1:2))
dfMan <- as_tibble(ManMds$points[,1:2]) %>% setNames(paste0('MDS', 1:2))

ggplot(dfMan, aes(x = MDS1, y = MDS2, colour = metas$DateS)) + geom_point()

# rank PCA
assayRank <- apply(assays, 2 , rank)
pcarank <- prcomp(assayRank)

pcaR_scores <- as.data.frame(pcarank$x)
pcaR_load <- as.data.frame(pcarank$rotation)

eigrank <- get_eigenvalue(pcarank)

eigrank[1:10,] %>% ggplot(aes(x = fct_inorder(rownames(.)), y = variance.percent)) + 
  geom_bar(stat = 'identity')

p2 <- ggplot(pcaR_scores, aes(x = PC1, y = PC2, colour = metas$DateA, shape = metas$Diagnosis)) + 
  geom_point() + 
  xlab(paste0('PC1: ', round(as.numeric(summary(pcaL)$importance[2,1]*100)), '% expl.var')) + 
  ylab(paste0('PC2: ', round(as.numeric(summary(pcaL)$importance[2,2]*100)), '% expl.var')) + 
  #scale_colour_brewer(palette = "Dark2") +
  theme_bw() + coord_fixed()

p2

###################################



# Relationship between S1 and S2 measurements
assayf <- assayf %>% mutate(Sample_ID = rownames(assayf)) %>% relocate(Sample_ID)
assayfl <- assayf %>% pivot_longer(cols = colnames(assayf[,2:ncol(assayf)]), names_to = 'compound', values_to = 'intensity')
assayfl <- assayfl %>% left_join(meta %>% select(RAD_ID, Sample_ID, Sample, CoreVisit))

assayfl1 <- assayfl %>% select(!Sample_ID) %>% pivot_wider(id_cols = c(RAD_ID, compound, CoreVisit), 
                                    names_from = Sample, values_from = intensity)

# plot S1 vs S2
ps <- assayfl1 %>% group_by(compound) %>%
  do(plots = ggplot(data = .) + (aes(y = log2(S1), x = log2(S2))) + geom_point() +
  ggtitle(.$compound) + coord_fixed() +
  theme_bw(base_size = 6) + geom_abline(intercept = 0, slope = 1, colour = 'red'))

pdf('Radica_S1vsS2.pdf')
marrangeGrob(grobs = ps$plots, nrow = 4, ncol = 4)
dev.off()

corrs <- assayfl1 %>% select(!MaskBG) %>% drop_na() %>%
  group_by(compound) %>% summarise(cor = cor(S1, S2, method = 'pearson'))

tiff('CorrHist_tech_reps.tiff', width = 100, height = 80, unit = 'mm', res = 300)
hist(corrs$cor, main = 'Correlation between technical replicates',
     xlab = 'Pearsons correlation coefficient (r)', cex.main = 0.9, cex.lab = 0.9, cex.axis = 0.9,
     breaks = 15)
dev.off()

# NAs in S1/S2 according to condition

assayfl1 <- assayfl1 %>% left_join(meta[, 1:2], by = 'RAD_ID') %>% unique()

################################



# Duplicate identifications

## Identify duplicated intensity values across all CV1 samples
library(tibble)

assayt <- assaytf
assayt <- as.data.frame(t(assay))

dupid <- function(y){
  dupval <- assayt %>% select(y) %>% 
    drop_na() %>% filter(duplicated(.))
  
  dup <- assayt %>% filter(get(y) %in% dupval[[1]]) %>% 
    select(y) %>%
    rename(quant = y) %>%
    mutate(sample = y)
  
  dup$quant <- as.numeric(dup$quant)
  
  dup <- dup %>% mutate(compound = rownames(dup)) %>%
    remove_rownames()
}

snos <- mf_70w$Sample

duplist <- lapply(snos, dupid)
dupframe <- bind_rows(duplist) %>%
  pivot_wider(names_from = sample, values_from = quant)


# dupframe - data frame with duplicated intensity values (cols = samples, rows = compounds)


## count number of duplications for each compound
rownames(dupframe) <- dupframe$compound
dupframet <- t(dupframe)
dupframet <- dupframet[-1,] %>% as.data.frame()

dupframe1 <- dupframe[,-1] 
dupframe1[is.na(dupframe1)] <- 0
dupframe1[dupframe1 > 0] <- 1
dupframe1 <- dupframe1 %>% mutate_if(is.character, as.numeric)

dupframe2 <- dupframe1 %>% 
  mutate(dupCount = rowSums(.),
         compound = dupframe$compound) %>% 
  arrange(desc(dupCount)) %>% relocate(dupCount, compound) # counts number of compound intensities that were also id'd as different compound

idCount <- assayt %>% transmute(idCount = ncol(.) - rowSums(is.na(.)),
                                compound = rownames(assayt))

fraqDup <- dupframe2 %>% select(dupCount, compound) %>% left_join(idCount) %>%
  mutate(fraqDup = dupCount/idCount*100)

# summary of duplication frequency across compounds
quantile(fraqDup$fraqDup)
hc <- hist(fraqDup$fraqDup, main = 'Duplication rate across compounds',
     xlab = '% values duplicated', breaks = 20, cex.main = 1)

## count number of duplications for each sample
dupframet1 <- as.data.frame(t(dupframe1)) %>% mutate(dupCount = rowSums(.)) %>% 
  relocate(dupCount)

SidCount <- assay %>% transmute(idCount = ncol(.) - rowSums(is.na(.)),
                                sample = rownames(assay))

SfraqDup <- dupframet1 %>% mutate(sample = rownames(dupframet1)) %>%
  select(sample, dupCount) %>% left_join(SidCount) %>%
  mutate(fraqDup = dupCount/idCount*100)
  
quantile(SfraqDup$fraqDup)
hs <- hist(SfraqDup$fraqDup, main = 'Duplication rate across samples',
     xlab = '% values duplicated', cex.lab = 1, cex.axis = 1, cex.main = 1, breaks = 20) 

# combine histograms
par(mfrow = c(2, 1))
hc <- hist(fraqDup$fraqDup, main = 'Duplication rate across compounds',
           xlab = '% values duplicated', breaks = 20, cex.main = 1)
hs <- hist(SfraqDup$fraqDup, main = 'Duplication rate across samples',
           xlab = '% values duplicated', cex.lab = 1, cex.axis = 1, cex.main = 1, breaks = 20) 

# plot duplicated values in the first sample
dupframe3 <- dupframe
dupframe3$compound <- str_sub(dupframe3$compound, end = -2)

comps <- read.csv('Compound_names_RT.csv')

dupframe3 <- dupframe3 %>% left_join(comps, by = 'compound')
dupframe3 <- dupframe3 %>% relocate(RT, compound_name) %>% select(!c(compound_name, X))
cnames <- colnames(dupframe3)
colnames(dupframe3)[3:467] <- sapply(1:465, function(x) {paste('S', x, sep = '')})

# link match factor scores with duplicated intensity values
s1 <- dupframe3 %>% select(RT, S1, compound) %>% drop_na() %>%
  left_join(mf2 %>% filter(Sample == '190524_RaDICA_RAD002_B1_353321_1') %>% 
              select(compound, MatchScore), by = 'compound')

duplot <- s1 %>%
  ggplot(aes(x = RT, y = S1, colour = as.factor(S1))) + 
  geom_point(position = position_dodge2(width = 0.3, preserve = 'single'), size = 1) +
  ylab('log10(intensity)') + xlab('Retention time [min]') +
  scale_y_log10() +
  ggtitle('Duplicated intensity values in one sample \n (13.2% of all sample values)') +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = 'none') +
  geom_text_repel(aes(x = RT, y = S1, label = MatchScore), max.overlaps = Inf, size = 4) +
  scale_color_manual(values = dark2)

duplot
ggsave('B1_duplicates.tiff', width = 160, height = 100, units = 'mm', dpi = 300)

write.csv(dupframe1, 'Radica_duplicates.csv')
dupframe1 <- read.csv('Radica_duplicates.csv')


# heatmap of duplicated values across compounds
library(pheatmap)
library(Polychrome)

dupframe1 <- dupframe1 %>% arrange(RT)

dupframet <- as.data.frame(t(dupframe1)) 
dupframet <- dupframet[-c(1:3), ]

dupframet <- dupframet %>% mutate_if(is.character, as.numeric)
dupframetL <- log10(dupframet)
dupframetL[is.na(dupframetL)] <- 0

tiff(filename = 'Heatmap_duplicates.tiff', width = 159, height = 129, unit = 'mm', res = 300)
tiff(filename = 'Heatmap_duplicates_all.tiff', width = 200, height = 140, unit = 'mm', res = 300)

dupHeat <- pheatmap(as.matrix(dupframetL), show_rownames = FALSE, cluster_cols =  FALSE, 
                    show_colnames = FALSE, cluster_rows = FALSE, 
                    main = 'Duplicated log10 intensity values across compounds',
                    border_color = 'black')
dupHeat
dev.off()

## compare match scores and delta RT of duplicated values
# try to choose a 'better match' or create compound groups matched to one peak intensity

# change column names back to sample IDs from 'data' df

# link full compound names with names from 'data' df
assayt$compound <- rownames(assayt)
datacomps <- as.data.frame(rownames(t(data[,14:1018])))
colnames(datacomps) <- 'compound'

comps <- read.csv('Compounds.csv')
comps <- cbind(comps, datacomps) %>% 
  separate(compound, c('compound', 'RT'), sep = '_',) %>% 
  mutate(RT = as.numeric(RT))

write.csv(comps, 'Compound_names_RT.csv')

# add full compound names to df with duplicated values
dupframe1 <- dupframe1 %>% left_join(comps) %>% 
  separate(compound_name, 'compound_name', sep = '_') # remove RT information




# format the match factor dataset (cols: Sample_ID, RT, Compounds, Match score)
mf <- read.csv('240206 Radica v4 VOC table with Library match scores per sample.csv')
mf <- mf[,-2]
compnames <- as.data.frame(colnames(mf))
colnames(compnames) <- 'Compound'

remove <- sapply(0:4020, function(i) {
  paste('X.', i, sep ='')
})

'%ni%' <- Negate('%in%')
compnames <- compnames %>% 
  filter(Compound %ni% remove) %>% 
  filter(Compound!= 'Sample')

compnames$Compound <- str_sub(compnames$Compound, end = -8)

compnames1 <- lapply(compnames$Compound, function(i) {
  c(paste(i, 'Area', sep = '_'), 
  paste(i, 'MatchScore', sep ='_'),
  paste(i, 'RT', sep = '_'),
  paste(i, 'RI', sep = '_'),
  paste(i, 'MZ', sep = '_'))
})

install.packages('purrr')
library(purrr)

compnames2 <- list_c(compnames1)

mf <- mf[-1,]
colnames(mf)[2:5026] <- compnames2
mf <- mf %>% pivot_longer(cols = 2:5026, names_to = 'compound', values_to = 'value') 
mf[,3] <- mf[,3] %>% mutate_if(is.character, as.numeric)
mf <- mf %>% separate(compound, c('compound', 'param'), sep ='_')
mf <- mf %>% pivot_wider(names_from = param, values_from = value)
mf$Sample <- str_sub(mf$Sample, end = -3)
mf$compound<- str_sub(mf$compound, end = -2)
mf2 <- mf %>% filter(Sample %in% data$Sample_ID)

mf_50 <- mf2 %>% filter(MatchScore >= 50)
mf_60 <- mf2 %>% filter(MatchScore >= 60)
mf_70 <- mf2 %>% filter(MatchScore >= 70)
mf_80 <- mf2 %>% filter(MatchScore >= 80)
mf_90 <- mf2 %>% filter(MatchScore >= 90)

fun <- function(x){
  mfw <- x %>% select(Sample, compound, Area) %>%
    pivot_wider(names_from = compound, values_from = Area)
  na <- as.data.frame(table(colSums(is.na(mfw[,-1]))/nrow(mfw)*100))
  na$Var1 <- as.numeric(levels(na$Var1))[na$Var1]
  nasum <- na %>% filter(Var1 < 30) %>% summarise(comps = sum(Freq))
  nasum[1,1]}

fun(mf_70)

mf_comp <- data.frame(MatchFactor = c(seq(from = 50, to = 90, by = 10)),
                  NoComps = c(n_distinct(mf_50$compound), 
                              n_distinct(mf_60$compound),
                              n_distinct(mf_70$compound),
                              n_distinct(mf_80$compound),
                              n_distinct(mf_90$compound)),
                  FiltComps = c(fun(mf_50), fun(mf_60), fun(mf_70), fun(mf_80), fun(mf_90))) %>%
  pivot_longer(cols = c(NoComps, FiltComps), names_to = 'Filter', values_to = 'Number')

mf_comp$Filter <- as.factor(mf_comp$Filter)
levels(mf_comp$Filter)

mf_comp$Filter <- factor(mf_comp$Filter, levels = c('NoComps', 'FiltComps'))

mfplot <- mf_comp %>% ggplot(aes(x = MatchFactor, y = Number, fill = Filter)) + 
  geom_bar(stat = 'identity', position = position_dodge2()) + 
  theme_bw(base_size = 10) + xlab('Match Factor cut-off') + ylab('Number of compounds') + ylim(0, 1005) +
  geom_text(aes(label = Number), vjust = -0.5, position = position_dodge(width = 8)) +
  scale_fill_manual(values = c('grey30', 'grey0'),
                    labels = c('No filter', 'Max 30% missing')) +
  ggtitle('Effect of match factor cut-off on identified compounds') +
  theme(plot.title = element_text(hjust = 0.5))

mfplot

ggsave('MFvsCompNo_plot.tiff', mfplot, width = 150, height = 84, units = 'mm', dpi = 300)


fun <- function(x){
mfw <- x %>% select(Sample, compound, Area) %>%
  pivot_wider(names_from = compound, values_from = Area)
na <- as.data.frame(table(colSums(is.na(mfw[,-1]))/nrow(mfw)*100))
na$Var1 <- as.numeric(levels(na$Var1))[na$Var1]
na %>% filter(Var1 < 30) %>% summarise(comps = sum(Freq))}

fun(mf_70)

write.csv(mf_70w, 'MF70_intensities_Batch1.csv')




# read match confidence scores and format to match the duplicate dataset
conf <- read.csv('Radica v4 VOC list.csv') %>% 
  rename(compound_name  = Compound.Name,
         sample = File.Name)
conf$sample <- str_sub(conf$sample, end = -3)

# create matching sample ID columns between match scores ('conf') and duplicate data ('dupframe1')

library(stringr)

dupnames <- as.data.frame(rownames(t(dupframe1))) %>%
  rename(ID = 'rownames(t(dupframe1))') %>%
  subset(!(ID %in% c('compound', 'RT', 'compound_name', 'X'))) %>% 
  cbind(data$Sample_ID) %>% 
  rename(Sample_ID = 'data$Sample_ID')

dupframel <- dupframe1 %>% pivot_longer(cols = 'B1':'S57.2', names_to = 'sample', 
                                        values_to =  'intensity') %>%
  mutate(sample = rep(dupnames$Sample_ID, 940)) %>% 
  select(!'compound') %>% drop_na()

# match duplicated values with match confidence scores
match <- dupframel %>% left_join(conf, by = c('compound_name', 'sample'))

# count missing match scores for duplicated intensities
nascores <- match %>% group_by(compound_name) %>% summarise(na = sum(is.na(Best.Hit)))
hist(nascores$na, breaks = 30)

# keep only duplicated intensities with associates match scores
match_complete <- match %>% drop_na()
match_dupv <- match_complete$intensity[duplicated(match_complete$intensity)]
match_dup <- match_complete %>% filter(intensity %in% match_dupv)

# plot duplicated intensities in one sample with match scores

library(ggrepel)

scoreplot <- match %>% filter(sample == '190524_RaDICA_RAD002_S1_353323_1') %>%
  ggplot(aes(x = as.factor(RT), y = intensity)) + 
  geom_point(size = 3, alpha = 0.5) +
  scale_y_continuous(trans = 'log2') +
  geom_text_repel(aes(label = round(Match.Factor, 2))) +
  theme_bw() +
  theme(axis.text.x = element_blank())

scoreplot

# all duplicated intensities with known match scores (across samples)
scoreplot1 <- match_dup %>%
  ggplot(aes(x = as.factor(round(RT, 2)), y = intensity, colour = sample)) + 
  geom_point(size = 3, alpha = 0.5) +
  scale_y_continuous(trans = 'log2') +
  geom_text_repel(aes(label = round(Delta.RI, 2))) +
  theme_bw() +
  theme(legend.position = 'none')

scoreplot1

# search for multiple matches to one comp/sample combination and choose one with better RT match
dup2 <- match_dup[duplicated(match_dup[,c('compound_name','sample')]),] %>% select(compound_name, sample)

dup3 <- match_dup %>% mutate(abs_RT_diff = abs(Component.RT - RT)) %>% 
  relocate(abs_RT_diff) %>% 
  filter(compound_name %in% dup2$compound_name & sample %in% dup2$sample) %>%
  filter(abs_RT_diff < 0.031)

# remove one repeat by hand
dup3 <- dup3[-2,]  

'%ni%' <- Negate('%in%')
dup4 <- match_dup %>% filter(compound_name %ni% dup2$compound_name & sample %ni% dup2$sample)

match_dup_corr <- rbind(dup4, dup3[,-1])
dupv2 <- match_dup_corr$intensity[duplicated(match_dup_corr$intensity)]
match_dup_corr <- match_dup_corr %>% filter(intensity %in% dupv2)
write.csv(match_dup_corr, 'Radica_duplicate_scores.csv')

# FINAL
# plot of true duplicates with associated match scores and delta RI

everysecond <- function(x){
  x <- sort(unique(x))
  x <- round(x, 2)
  x[seq(2, length(x), 2)] <- ""
  x
}

match_dup_corr <- match_dup_corr %>% 
  left_join(meta[,c('Sample_ID', 'Sample')], by = c('sample' = 'Sample_ID'))

scoreplot2 <- match_dup_corr %>%
  ggplot(aes(x = as.factor(RT), y = intensity, colour = Sample)) + 
  geom_point(size = 3) +
  scale_y_log10() +
  scale_x_discrete(labels = everysecond(match_dup_corr$RT)) +
  geom_text_repel(aes(label = round(Match.Factor, 1))) +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5, size = 17)) +
  ylab('Log10(intensity)') + xlab('Retention time') +
  ggtitle('Match factor of duplicated identifications') +
  scale_color_brewer(palette = 'Dark2')

scoreplot2  
ggsave('Match_Factor_duplicates.tiff', width = 282, height = 152, unit = 'mm', dpi = 300)


# score information about 70 duplicated intensity values out of 28945
# bigger difference between duplicates in Retention index than Match Score




#########################################

# Data transformation

# plot value distributions for each compound

library(gridExtra)
library(RColorBrewer)

BiocManager::install("vsn")
library(vsn)

assay <- assay %>% mutate(Sample_ID = rownames(assay))

assayl <- assay %>% pivot_longer(cols = !Sample_ID,
                                 names_to = 'compound',
                                 values_to = 'intensity')

assayl <- drop_na(assayl)
assayl <- assayl %>% mutate(logInt = log2(intensity))
assaylan <- assayl %>% left_join(meta[, c('Diagnosis', 'Sample', 'Sample_ID', 'RAD_ID')])  

pdf('Radica_hist_log_CV14.pdf')

p1 <- assayfl %>% group_by(compound) %>%
  do(plots = ggplot(data = .) + (aes(x = intensity)) + geom_histogram() +
       ggtitle(.$compound) + ylab('log_intensity') +
       theme_bw(base_size = 6) +
       scale_x_continuous(trans = 'log10'))

marrangeGrob(grobs = p1$plots, nrow = 5, ncol = 4)

dev.off()


assaylan %>% filter(compound == 'Cyclopropylacetylene_1.791') %>%
  ggplot(aes(x = logInt, fill = Sample)) + geom_histogram(position = 'identity', alpha = 0.6)


densint <- assaylan %>% ggplot(aes(x = intensity, colour = Sample)) + 
  geom_density(alpha = 0.5, linewidth = 0.7) + theme_bw(base_size = 12) + theme(legend.position = 'none')

densint
  
denslog <- assaylan %>% ggplot(aes(x = logInt, colour = Sample)) + 
  geom_density(alpha = 0.5, linewidth = 0.7) + theme_bw(base_size = 12) + 
  xlab('log2(intensity)')

denslog

g <- arrangeGrob(densint, denslog, widths = c(0.4, 0.6))
plot(g)

ggsave('Density_plots.tiff', g, width = 150, height = 67, units = 'mm', dpi = 300)

# plot mean vs sd for each compound
assayl$intensity <- as.numeric(assayl$intensity)

table(is.na(assayl$intensity))

sumint1 <- assayl %>% group_by(compound) %>%
  drop_na() %>%
  summarise(meanInt = mean(intensity),
            sdInt = sd(intensity))

colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

mraw <- sumint1 %>% ggplot(aes(x = meanInt, y = sdInt)) + 
  geom_point() + coord_fixed() + ggtitle('Raw')

sumint1 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_hex() + coord_fixed()

# wo outlier
mraw1 <- sumint1 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_point() + 
  coord_fixed() + xlim(0, 2.5*10^7) + ylim(0, 2.5*10^7) + ggtitle('Raw')

sumint1 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_hex() + coord_fixed() + xlim(0, 2.5*10^7) + ylim(0, 2.5*10^7)

sumint1 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_point() + coord_fixed() + xlim(0, 2.5*10^7) + ylim(0, 2.5*10^7)

# plot mean rank vs sd for each compound
sumint1 <- sumint1 %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())


mrank <- sumint1 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point() + ylim(0, 2.0*10^7)

sumint1 %>% ggplot(aes(x = meanRank, y = sdInt)) + stat_bin_hex(bins = 50) + ylim(0, 2.0*10^7)

rank <- sumint1 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point(alpha = 0.5) + ylim(0, 2.0*10^7) +
  theme_bw(base_size = 12) + ylab('sd(intensity)') + xlab('rank(mean_intensity)') 

# log2 transformation
assayl <- assayl %>% mutate(logInt = log2(intensity))

sumintLog <- assayl %>% group_by(compound) %>%
  drop_na() %>%
  summarise(meanInt = mean(logInt),
            sdInt = sd(logInt))

lraw <- sumintLog %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_point() + ggtitle('Log2')

sumintLog <- sumintLog %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())

lrank <- sumintLog %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point(alpha = 0.5) +
  theme_bw(base_size = 12) + ylab('sd(Log2_intensity)') + xlab('rank(mean_Log2_intensity)')

# raw vs log transformed plots

p <- arrangeGrob(rank, lrank, widths = c(0.53, 0.47))
ggsave('Sd_vs_mean_rank.tiff', p, width = 120, height = 67, units = 'mm', dpi = 300)

# size factor?

assay1 <- assay
assay1[is.na(assay1)] <- 0 
assay1[,2:1006] <- assay1[,2:1006] %>% mutate_if(is.character, as.numeric)

assay1 <- assay1 %>% mutate(sumInt = rowSums(assay1[,2:1006]))

barplot(sort(assay1$sumInt))

# vsn normalisation
# ?? NAs throwing transformation off?

library(limma)
library(gridExtra)

rownames(assay) <- assay$Sample_ID
assay2 <- assay[,-1]
assay2t <- t(assay2)
vsn <- vsn::vsnMatrix(assay2t)
vsnout <- vsn@hx
vsnout <- as.data.frame(vsnout)
vsnout$compound <- rownames(vsnout)

vsnl <- vsnout %>% pivot_longer(cols = !compound,
                                 names_to = 'Sample_ID',
                                 values_to = 'intensity')

vsnl <- drop_na(vsnl)

sumint2 <- vsnl %>% group_by(compound) %>%
  drop_na() %>%
  summarise(meanInt = mean(intensity),
            sdInt = sd(intensity))

vraw <- sumint2 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_point() + ggtitle('vsn')
vraw

sumint2 <- sumint2 %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())

vrank <- sumint2 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point() + ggtitle('Vsn')
vrank

# rlog normalisation
## sensitive to outliers!!
BiocManager::install("DESeq2")
library(DESeq2)

assay2 <- as.matrix(assay2)
assay2[is.na(assay2)] <- 0 # rlog does not accept NAs as input
rlog <- rlog(assay2)
write.csv(rlog, 'Radica_rlog.csv')

rlog <- as.data.frame(rlog)
rlog$SampleID <- rownames(rlog)
rlogl <- rlog %>% pivot_longer(cols = !c(SampleID),
                             names_to = 'compound',
                             values_to = 'intensity')

quantile(rlogl$intensity)

sumint3 <- rlogl %>% group_by(compound) %>%
  drop_na() %>%
  summarise(meanInt = mean(intensity),
            sdInt = sd(intensity))

rlograw <- sumint3 %>% ggplot(aes(x = meanInt, y = sdInt)) + geom_point() + ylim(0, 75) + xlim(0, 75)

sumint3 <- sumint3 %>%
  arrange(meanInt) %>%
  mutate(meanRank = row_number())

rlogrank <- sumint3 %>% ggplot(aes(x = meanRank, y = sdInt)) + geom_point()

# combine plots
grid.arrange(mraw1, lraw, vraw, ncol = 3, nrow = 1)
grid.arrange(mrank, lrank, vrank, ncol = 3, nrow = 1)

# plot transformed data against original data
vsnl <- vsnl %>% rename(intensity = 'vsnInt')
assayl1 <- assayl %>% left_join(vsnl)

assayl1 %>% pivot_longer(cols = c(logInt, vsnInt), names_to = 'transf', values_to = 'transInt') %>%
  filter(Sample_ID == '190524_RaDICA_RAD002_B1_353321_1') %>%
  ggplot(aes(x = intensity, y = transInt, colour = transf)) + geom_line()




#########################################

# Relationship between background and sample VOC levels
assayl <- assay %>% mutate(Sample_ID = rownames(assay)) %>%
  pivot_longer(!Sample_ID, names_to = 'compound', values_to = 'intensity')

assaylan <- assayl %>% left_join(meta[, c('Diagnosis', 'Sample', 'Sample_ID', 'RAD_ID', 'CoreVisit')])


sbg <- assaylan %>% 
  pivot_wider(id_cols = c(RAD_ID, compound, CoreVisit, Diagnosis), 
                         names_from = Sample, 
                         values_from = intensity) %>%
  pivot_longer(cols = c(S1,S2), 
               names_to = 'Sample', 
               values_to = 'intensity_sample') 

pdf('Radica_SvsBG_CV14.pdf')

p2 <- sbg %>% drop_na() %>% group_by(compound) %>%
  do(plots = ggplot(data = .) + (aes(x = log10(MaskBG), y = log10(intensity_sample), colour = Diagnosis)) + 
       geom_point(alpha = 0.7) + ggtitle(.$compound) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) + coord_fixed() +
       xlab('log10(background)') +
       ylab('log10(sample)'))

marrangeGrob(grobs = p2$plots, nrow = 4, ncol = 4)

dev.off()

sbg <- sbg %>% left_join(unique(meta[,c('RAD_ID', 'Diagnosis')]))

pdf('Radica_SvsBG_Diagnosis.pdf')

p2 <- sbg %>% drop_na() %>% group_by(compound) %>%
  do(plots = ggplot(data = .) + (aes(x = MaskBG, y = logInt_Sample, colour = Diagnosis)) + 
       geom_point(alpha = 0.7) + ggtitle(.$compound) + theme_bw(base_size = 8) + 
       theme(legend.position = 'none',
             plot.title = element_text(size = 6)) + coord_fixed())

marrangeGrob(grobs = p2$plots, nrow = 4, ncol = 4)

dev.off()


# Modelling predictive value of one compound
# Hexyl.octyl.ether_19.821 - S dep on BG, good predictor
# Spiro.2.4.hepta.4.6.diene. - S dep on BG diff across groups, good predictor
# Acetoin_2.722 - S non dep on BG, good predictor
# Butyrolactone_6.053 - S dep on BG, no diff across groups


library(lme4)
library(lmerTest)

# modelling compound level 

test3 <- assaylan %>% filter(compound == 'Hexyl.octyl.ether_19.821') %>% filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  drop_na(intensity) %>% mutate(intensityS = intensity/sd(intensity)) #%>% mutate(logInt  = log10(intensity)) 

hist(test3$intensityS)

test3 <- test3 %>% select(!Sample_ID) %>% pivot_wider(id_cols = c(RAD_ID, compound, CoreVisit), 
                                                      names_from = Sample, 
                                                      values_from = intensityS) %>% #logInt
  left_join(unique(meta[,c('RAD_ID', 'Diagnosis')])) %>%
  pivot_longer(cols = c(S1,S2), 
               names_to = 'Sample', 
               values_to = 'int_Sample') 

test3 <- test3 %>% drop_na(intensity_sample)
test3 <- test3 %>% mutate(logInt = log2(intensity_sample),
                          logBG = log2(MaskBG))

m0 <- lmer(logInt ~ logBG + Diagnosis + (1 | RAD_ID),
            data = test3)

summary(m0)

m1 <- glmer(int_Sample ~ Diagnosis + Diagnosis:CoreVisit + (1 | RAD_ID),
            data = test3,
            family = Gamma(link = 'log'))

summary(m1)

# with log transformation of intensity data and lmer

test3 <- assaylan %>% filter(compound == 'Aminoacetonitrile_3.946') %>% filter(CoreVisit %in% c('CV1', 'CV2')) %>%
  drop_na(intensity) %>% mutate(logInt  = log10(intensity)) 

hist(test3$logInt)

test3 <- test3 %>% select(!Sample_ID) %>% pivot_wider(id_cols = c(RAD_ID, compound, CoreVisit), 
                                                      names_from = Sample, 
                                                      values_from = logInt) %>%
  left_join(unique(meta[,c('RAD_ID', 'Diagnosis')])) %>%
  pivot_longer(cols = c(S1,S2), 
               names_to = 'Sample', 
               values_to = 'int_Sample') 

test3 <- test3 %>% drop_na(int_Sample)

m0 <- lmer(int_Sample ~ MaskBG + Diagnosis + Diagnosis:CoreVisit + (1 | RAD_ID),
            data = test3)

summary(m0)

m1 <- lmer(int_Sample ~ Diagnosis + Diagnosis:CoreVisit + (1 | RAD_ID),
            data = test3)

summary(m1)


# correcting for background based on linear regression
test1 %>% pivot_longer(cols = c(MaskBG, S1, S2), names_to = 'Sample', values_to = 'intensity') %>%
  ggplot(aes(x = Sample, y = intensity, group = RAD_ID)) + geom_point() + 
  geom_line(aes(colour = RAD_ID)) + scale_y_continuous(trans = 'log') + ylab('logInt')


test2 %>% ggplot(aes(x = MaskBG, y = int_Sample, colour = Diagnosis)) + geom_point() + 
  scale_x_continuous(trans = 'log2') + scale_y_continuous(trans = 'log2') + coord_fixed()


test3 <- sbg %>% filter(compound == 'Acetic.acid..TMS.derivative_2.772') #%>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

hist1 <- test3 %>% ggplot(aes(x = log2(intensity_sample))) +
  geom_histogram(aes(fill = Diagnosis), 
                 position = 'identity', alpha = 0.6,
                 bins = 15) +
  theme_bw() + scale_fill_brewer(palette = 'Dark2') +
  xlab('log2(intensity)') + ggtitle('No background correction') +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = 'none')

hist1

dot <- test3 %>% ggplot(aes(x = log2(MaskBG), y = log2(intensity_sample))) + 
  geom_point(aes(colour = Diagnosis)) + theme_bw() +
  scale_color_brewer(palette = 'Dark2') + 
  xlab('log2(intensity_background)') + ylab('log2(intensity_breath)') +
  geom_smooth(method = 'lm', se = FALSE) +
  annotate('text', label = 'slope = 0.4', x = 10, y = 20, colour = 'blue') +
  theme(legend.position = 'none')

dot

ggsave('BGvsBreath_TMS.tiff', width = 100, height = 90, unit = 'mm', dpi = 300)

m3a <- lm(log2(intensity_sample) ~ log2(MaskBG),
           data = test3)

summary(m3a)

c <- coef(m3a)[2]
c

test3 <- test3 %>% mutate(intensity_Corr = log2(intensity_sample) - (c * log2(MaskBG)))

hist2 <- test3 %>% ggplot(aes(x = log2(intensity_Corr))) +
  geom_histogram(aes(fill = Diagnosis), 
                 position = 'identity', alpha = 0.6,
                 bins = 15) +
  theme_bw() + scale_fill_brewer(palette = 'Dark2') +
  xlab('log2(intensity)') + ggtitle('With background correction') +
  theme(plot.title = element_text(hjust = 0.5))

hist2

hist1 <- hist1 + xlim(5,20)
hist2 <- hist2 + xlim(-2, 3)

g <- arrangeGrob(hist1, hist2, nrow = 1, ncol = 2, widths = c(0.43, 0.57))

plot(g)

ggsave('BackgroundCorrHistTMS.tiff', g, width = 230, height = 90, unit = 'mm', dpi = 300)

sd(log2(na.omit(test3$intensity_sample)))
sd(na.omit(test3$intensity_Corr))

# effect of correction on predicting diagnosis

test3 <- test3 %>% mutate(outcome = ifelse(Diagnosis == 'Asthma', 1, 0))

md <- glm(outcome ~ log2(intensity_sample),
            data = test3,
            family = binomial) 

summary(md)
confint(md)

md1 <- glm(outcome ~ intensity_Corr,
          data = test3,
          family = binomial) 
confint(md1)

summary(md1)

test3 %>% ggplot(aes(x = Diagnosis, y = log2(intensity_sample))) + geom_boxplot()
test3 %>% ggplot(aes(x = Diagnosis, y = intensity_Corr)) + geom_boxplot()

# modelling diagnosis category based on corrected intensity value

test3 %>% ggplot(aes(x = Diagnosis, y = int_SamCorr)) + geom_boxplot() + scale_y_continuous(trans = 'log2')

m4 <- glm(Diagnosis ~ int_SamCorr,
          data = test2,
          family = binomial)

summary(m4)

###

dates <- data %>% filter(Sampling_date %in% c('16/01/2020', '19/11/2019', '20/11/2019', '25/11/2019'))
dates <- dates %>% select(RAD_ID, Sampling_date, Diagnosis, Sample, 14:1018)
datesl <- dates %>% pivot_longer(cols = 5:1009, names_to = 'compound', values_to = 'intensity')
datesl %>% filter(Sample == 'MaskBG') %>% ggplot(aes(x = Sampling_date, y = intensity, colour = RAD_ID)) + 
  geom_boxplot() +  scale_y_continuous(trans = 'log2')

#############################

## Multilevel decomposition


assayst <- assays %>% select(!Sample_ID) %>% t() %>% as.data.frame()
assayst <- log2(assayst)
assayst1 <- assayst %>% filter(rowSums(is.na(assayst))/ncol(assayst) <= 0)

assaysc <- as.data.frame(t(assayst1)) %>% mutate(Sample_ID = rownames(.))
assayscan <- assaysc %>% left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, Diagnosis, RAD_ID, Sample)) %>%
  relocate(Sample_ID, CoreVisit, Diagnosis, RAD_ID, Sample)

complete <- assayscan %>% group_by(RAD_ID) %>% 
  summarise(n = n_distinct(CoreVisit)) %>%
  filter(n == 4)
  
a <- assaysc %>% dplyr::select(!Sample_ID) 

a2 <- assayscan %>% filter(Diagnosis == 'Not Asthma') %>% 
  filter(CoreVisit %in% c('CV1', 'CV4')) %>%
  filter(RAD_ID %in% complete$RAD_ID)

a2assay <- a2[,6:ncol(a2)]

library(mixOmics)

pca2 <- pca(a2assay, scale = TRUE)

p1 <- plotIndiv(pca2,
          group = a2$CoreVisit,
          ind.names = a2$RAD_ID,
          legend = TRUE,
          title = 'PCA on breath VOC data',
          size.title = rel(0.8),
          cex = 2.5)

ggsave('PCAexample.tiff', width = 140, height = 95, units = 'mm', dpi = 300)

# pca with multilevel argument
multipca2 <- pca(a2assay,
                 multilevel = a2$RAD_ID,
                 scale = TRUE)

p2 <- plotIndiv(multipca2,
          group = a2$CoreVisit,
          ind.names = a2$RAD_ID,
          title = 'Multilevel PCA on breath VOC data',
          legend = TRUE,
          size.title = rel(0.8),
          cex = 2.5)

ggsave('MultilevelPCAexample.tiff', width = 140, height = 95, units = 'mm', dpi = 300)
