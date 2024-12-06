## Normalisation of peak area data

# author: Aggie Turlo
# project: RADicA
# date: 17/07/2024

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(mixOmics)
library(pcaPP)
library(stringr)

'%ni%' <- Negate('%in%')

b_all$Date <- as.Date(b_all$Date)
b_imp$Analysis_date <- as.Date(b_imp$Analysis_date)

b1_praw <- b1_imp %>%  
  filter(class %in% c('BG', 'S', 'S2')) %>%
  ggplot(aes(y = log(Octane), x = Analysis_date, group = Analysis_date)) + 
  geom_boxplot(outliers = FALSE) +
  theme_bw(base_size = 8) +
  ggtitle('Raw peak area') +
  ylab('log(peak area)') +
  ylim(8, 13) +
  theme(plot.title = element_text(hjust = 0.5)) #+
  #geom_hline(yintercept = median(log(b1_imp_L$Octane)),
             #linetype = 'dashed', colour = 'blue')


b1_norm_L %>%  
  filter(class %in% c('ES')) %>%
  filter(comp == 'Octane') %>%
  filter(Batch != 5) %>%
  ggplot(aes(y = peakArea, x = as.Date(Analysis_date), group = as.Date(Analysis_date))) + 
  geom_point() +
  theme_bw(base_size = 12) +
  ggtitle('Raw peak area of external standard (Octane)') +
  ylab('peak area') +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('Date')

  

b1_praw
b1_pcc2
b1_pcc2_pqn

b1p <- grid.arrange(b1_praw, b1_pcc2, b1_pcc2_pqn, ncol = 1)

ggsave('ES_B1_time.tiff', b1p, dpi = 300, unit = 'mm', width = 130, height = 100)



# read datasets
b_all <- read.csv('RADicA_B2.csv')
b_all_L <- b_all %>% pivot_longer(cols =! c(1:6), names_to = 'comp', values_to = 'peakArea')

# imputed
b_imp <- read.csv('RADicA_B2_NAfiltered_imputed.csv')[,-1] # B2

b1_imp_L <- read.csv('RADicA_B1_NAfiltered_imputed_long.csv')[,-1] # B1

b_imp_L <- b_imp %>% pivot_longer(cols = c(6:ncol(b_imp)), names_to = 'comp', values_to = 'peakArea')

b1_imp_sum_corr_w1 <- read.csv('RADicA_BG_adjusted_outl_removed.csv')[,-1]

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

#
train_L <- read.csv('RADicA_Train_NAfiltered_blank_filtered_imputed_long.csv')[,-1]


# match features between B2 and finished B1 dataset
intersect(b1_imp_L$comp, b_imp_L$comp)
setdiff(b1_imp_L$comp, b_imp_L$comp)
setdiff(b_imp_L$comp, b1_imp_L$comp)
b2nab1 <- setdiff(b_imp_L$comp, b1_imp_L$comp)
b1nab2 <- setdiff(b1_imp_L$comp, b_imp_L$comp)

b_imp_L <- b_imp_L %>% filter(comp %ni% b2nab1)
n_distinct(b_imp_L$comp)

b1_imp_L <- b1_imp_L %>% filter(comp %ni% b1nab2)
n_distinct(b1_imp_L$comp)

train_L <- train_L %>% filter(comp %ni% b1nab2)

# move B5 samples from dataset 1 to dataset 2
b_imp_L <- b_imp_L %>% rbind(b1_imp_L %>% filter(Batch == 5)) 

b1_imp_L <- b1_imp_L %>% filter(Batch != 5) # in pre-covid dataset  

#

b_imp <- b_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea) %>%
  mutate(Analysis_date = paste(str_sub(Sample, start = 5, end = 6),
                               str_sub(Sample, start = 3, end = 4),
                               20, str_sub(Sample, start = 1, end = 2),
                               sep = '')) %>%
  relocate(Analysis_date) %>%
  dplyr::select(!c(Date, Time))

b_imp$Analysis_date <- as.Date(b_imp$Analysis_date, format = '%d%m%Y')

rownames(b_imp) <- b_imp$Sample


# replace imputed peakArea values < 1 with 1
b1_imp_L <- b1_imp_L %>% 
  mutate(peakArea = ifelse(peakArea < 1, 1, peakArea))

b_imp_L <- b_imp_L %>% 
  mutate(peakArea = ifelse(peakArea < 1, 1, peakArea))

b1_imp[b1_imp < 1] <- 1

b_imp[b_imp < 1] <- 1

train_L <- train_L %>% mutate(peakArea = ifelse(peakArea < 1, 1, peakArea))

train[train < 1] <-1
#
#
#

# train dataset - keep ES and Blank samples analysed on the same days
dates <- train %>% filter(class %ni% c('ES', 'Blank')) %>% 
  dplyr::select(Analysis_date) %>% distinct()

table(train$class)

train <- train %>% filter(Analysis_date %in% dates$Analysis_date)


# NORMALISATION OF IMPUTED DATA
# IS normalisation
b_IS <- cbind(b_imp[,1:4], b_imp[,5:ncol(b_imp)]/b_imp$Internal_Standard) %>%
  #filter(class != 'Blank') %>%
  dplyr::select(!Internal_Standard)

head(b_IS) # IS normalised dataset

#
#
#

# CC normalisation
# remove outliers from reference datasets (es and blank)
data <- train

data <- as.data.frame(data)
rownames(data) <- data$Sample

assay_es_imp <- data %>% filter(class %in% c('Blank')) %>% 
  dplyr::select(!c(Sample, Analysis_date, class, Batch)) # or Analysis_date

nacols <- which(colSums(is.na(assay_es_imp) == TRUE) > 0)

assay_es_imp <- assay_es_imp[,-nacols] # remove variables with NAs 


pc <- PCAgrid(log(assay_es_imp), scale = 'sd', center = 'mean')

as.data.frame(pc$scores) %>%
  mutate(Sample = rownames(assay_es_imp)) %>% 
  left_join(data %>% filter(class == 'Blank') %>% 
              dplyr::select(Batch) %>% mutate(Sample = rownames(.))) %>% 
  ggplot(aes(x = Comp.1, y = Comp.2, colour = as.factor(Batch))) + geom_point()


# identify outliers based on the critical OD and SD values (0.97 quantile)
sdod <- PCdiagplot(assay_es_imp, pc, plotbw = FALSE)
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(assay_es_imp)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,1] | V2 > sdod$critSD[2,1])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(assay_es_imp)
odOut <- od %>% filter(V1 > sdod$critOD[1,1] | V2 > sdod$critOD[2,1])

# remove outliers from the QC database
'%ni%' <- Negate('%in%')
assay_es_imp1 <- assay_es_imp %>% filter(rownames(.) %ni% rownames(odOut)) %>%
  filter(rownames(.) %ni% rownames(sdOut)) #%>%
  #filter(rownames(.) %ni% c('220207_15_RADicA_Tr_RAD198_362254')) # outlier blanks

#assay_es_imp1 <- assay_es_imp1[-51,] # outlier ES B2

# run CC on spikes database
pc_es <- mixOmics::pca(log(assay_es_imp1), scale = TRUE, center = TRUE, ncomp = 2) 

b_es <- data %>% filter(class == 'Blank') %>% filter(Sample %ni% rownames(odOut))  %>%
  filter(Sample %ni% rownames(sdOut)) #%>%
  #filter(rownames(.) %ni% c('220207_15_RADicA_Tr_RAD198_362254')) # outlier blanks B2


#b_es <- b_es[-51,] # outlier ES B2

plotIndiv(pc_es,
          group = b_es$Batch,
          pch = 1,
          legend = TRUE,
          )

loadings <- pc_es$loadings[['X']] %>% as.data.frame() %>%
  as.matrix()

loadings1 <- pc_es$loadings[['X']] %>% as.data.frame() %>% dplyr::select(PC1) %>%
  as.matrix()


Y <- data %>% 
  filter(class != 'Blank') %>%
  dplyr::select(colnames(assay_es_imp1)) %>% 
  log() %>% as.matrix() %>%
  scale(scale = TRUE, center = TRUE)

#

x0 <- Y
x0[is.na(Y)] <- 0 # managment of NAs in matrix multiplication

drift <- (x0 %*% loadings1) %*% t(loadings1)
drift[is.na(Y)] <- NA

# both components in one matrix

drift2 <- (x0 %*% loadings) %*% t(loadings)
drift2[is.na(Y)] <- NA

Z <- Y - drift
Z2 <- Y - drift2
#

Zres <- Z
for (i in c(1:dim(Y)[2])){
  Zres[,i] <- (Z[,i]*attributes(Y)$`scaled:scale`[i])+attributes(Y)$`scaled:center`[i]
}

#
Zres2 <- Z2
for (i in c(1:dim(Y)[2])){
  Zres2[,i] <- (Z2[,i]*attributes(Y)$`scaled:scale`[i])+attributes(Y)$`scaled:center`[i]
}

#

ann <- data %>% 
  filter(class != 'Blank')

train_cc <- cbind(ann[,1:4], exp(Zres)) 
train_cc2 <- cbind(ann[,1:4], exp(Zres2))

head(b_cc) # CC normalised dataset (1 component)
head(b_cc2) # CC normalised datase (2 components)

write.csv(b_cc, 'RADicA_B2_NAfiltered_imputed_CC.csv')
write.csv(b_cc2, 'RADicA_B2_NAfiltered_imputed_CC2.csv')

write.csv(b1_cc2, 'RADicA_B1_NAfiltered_imputed_CC2.csv')
write.csv(b1_cc, 'RADicA_B1_NAfiltered_imputed_CC2.csv')

b_cc2 <- read.csv('RADicA_B1_NAfiltered_imputed_CC2.csv')[,-1]
b_cc <- read.csv('RADicA_B2_NAfiltered_imputed_CC.csv')[,-1]

#

# compare density distributions of cc-normalised datasets
b_cc2 %>% mutate(data = 'B2') %>% 
  filter(class %in% c('S1', 'S2', 'BG')) %>%
  dplyr::select(!c(1:4)) %>%
  pivot_longer(cols =! data, names_to = 'comp', values_to = 'peakArea') %>%
  drop_na() %>%
  rbind(# compare density distributions of cc-normalised datasets
    b1_cc2 %>% mutate(data = 'B1') %>% 
      filter(class %in% c('S1', 'S2', 'BG')) %>%
      dplyr::select(!c(1:4)) %>%
      pivot_longer(cols =! data, names_to = 'comp', values_to = 'peakArea') %>%
      drop_na()) %>%
  ggplot(aes(x = log(peakArea), colour = data)) + geom_density() +
  theme_bw() +
  ggtitle('Normalised peak areas (CC2)')

ggsave('Density_B1_B2_CC2.tiff', unit = 'mm', dpi = 300, width = 100, height = 70)




#
#
#

# PQN normalisation
# using median peak intensity of BG,S1,S2 and ES datasets as reference 
# create reference spectrum
x <- b_imp

pqn <- function(x) { 
  xL <- x %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class),
                           names_to = 'comp', values_to = 'peakArea')
  ref <- xL %>% #filter(class %ni% c('Blank')) %>%
    group_by(comp) %>% summarise(median = median(peakArea, na.rm = TRUE))
# calculate quotients of each observation and reference spectrum
  quot <- xL %>% left_join(ref) %>%
  mutate(quotient = peakArea/median) %>% 
  group_by(Sample) %>%
  summarise(sizeEf = median(quotient, na.rm = TRUE)) # take median of quotients for each sample (size effect)
# divide each observation by the sample size effect
  b1_pqn <- xL %>% filter(class %ni% c('Blank')) %>%
  left_join(quot) %>% mutate(peakAreaN = peakArea/sizeEf) %>%
  dplyr::select(!c(peakArea, sizeEf)) %>%
  rename(peakArea = peakAreaN) %>%
  pivot_wider(names_from = comp, values_from = peakArea) #%>%
  #dplyr::select(!Internal_Standard)
  }

b_pqn <- pqn(b_imp) # PQN normalised dataset

write.csv(b1_pqn, 'RADicA_B1_NAfiltered_imputed_PQN.csv')

#

# use only ES/blanks as reference
pqn1 <- function(ref, data) { 
  xL <- ref %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class),
                           names_to = 'comp', values_to = 'peakArea')
  ref_spectrum <- xL %>% #filter(class %in% c('Blank', 
                          #          'ES'
                           #         )) %>%
    group_by(comp) %>% summarise(median = median(peakArea, na.rm = TRUE))
  # format dataset for normalisation
  yL <- data %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class), #Analysis_date
                             names_to = 'comp', values_to = 'peakArea')
  # calculate quotients of each observation and reference spectrum
  quot <- yL %>% left_join(ref_spectrum) %>%
    mutate(quotient = peakArea/median) %>% 
    group_by(Sample) %>%
    summarise(sizeEf = median(quotient, na.rm = TRUE)) # take median of quotients for each sample (size effect)
  # divide each observation by the sample size effect
  b1_pqn <- yL %>% 
    #filter(class %ni% c('Blank')) %>%
    left_join(quot) %>% mutate(peakAreaN = peakArea/sizeEf) %>%
    dplyr::select(!c(peakArea, sizeEf)) %>%
    rename(peakArea = peakAreaN) %>%
    pivot_wider(names_from = comp, values_from = peakArea) #%>%
  #dplyr::select(!Internal_Standard)
}

b1_pqn1 <- pqn1(b1_imp, b1_imp)
b_pqn1 <- pqn1(ref = b1_imp, data = b_imp)

write.csv(b1_pqn1, 'RADicA_B1_NAfiltered_imputed_PQN_noB5.csv')
write.csv(b_pqn1, 'RADicA_B2_NAfiltered_imputed_PQN.csv')



#
#
#

# CC + PQN
b1_cc2_pqn <- pqn1(b1_cc2, b1_cc2) 
b_cc2_pqn <- pqn1(b_cc2, b_cc2)
b_cc2_pqn_b1ref <- pqn1(b1_cc2, b_cc2) 

test_cc2_pqn <- pqn1(test_cc2, test_cc2)
train_cc2_pqn <- pqn1(train_cc2, train_cc2)

write.csv(b1_cc2_pqn, 'RADicA_B1_NAfiltered_imputed_CC2_PQN_noB5.csv')
write.csv(b_cc2_pqn_b1ref, 'RADicA_B2_NAfiltered_imputed_CC2_PQN.csv')

# ASSESSMENT OF NORMALISATION OUTCOMES
# remove IS column from raw data
b_imp1 <- b_imp %>% #dplyr::select(!Internal_Standard) %>%
  #filter(class != 'Blank') %>% 
  as.data.frame()

# Measures of relative dispersion
disp <- function(x, y){ # x is dataframe of normalised values, y is normalisation method
  x %>%
  pivot_longer(cols =! c(Sample, Analysis_date, Batch, class),
               names_to = 'comp', values_to = 'peakArea') %>%
  group_by(class, comp) %>%
  summarise(mean = mean(peakArea, na.rm = TRUE),
            sd = sd(peakArea, na.rm = TRUE),
            rsd = sd/abs(mean),
            median = median(peakArea, na.rm = TRUE),
            IQR = IQR(peakArea, na.rm = TRUE),
            rIQR = IQR/median,
            MAD = mad(peakArea, na.rm = TRUE),
            rMAD = MAD/median) %>%
    mutate(norm = y)
}


disp_raw <- disp(b_imp, 'Raw')
disp_IS <- disp(b_IS, 'IS')
disp_cc <- disp(b_cc, 'CC')
disp_cc2 <- disp(b_cc2, 'CC2')
#rsd_pqn <- RSD(b1_pqn, 'PQN')
disp_pqn1 <- disp(b_pqn1, 'PQN')
disp_cc_pqn <- disp(b_cc2_pqn_b1ref, 'CC PQN')
disp_cc_es <- disp(b_cc_es, 'CC ES')
disp_cc2_es <- disp(b_cc2_es, 'CC2 ES')

disp_conc <- rbind(disp_raw, #disp_IS, 
                   disp_cc, disp_cc2, #disp_cc_es, disp_cc2_es, disp_pqn1,
                   disp_cc_pqn)

# visualise distribution of relative dispersion values following different normalisation methods
plot_disp <- function(m) { # m = dispersion measure: RSD, rIQR, rMAD
  disp_conc %>% 
  filter(class %ni% c('Blank')) %>%
  ggplot(aes(x = .data[[m]], y = factor(norm, levels = c('CC2 ES', 'CC ES', 'CC2', 'CC', 'PQN', 'PQN QC', 'IS', 'Raw')), 
             fill = norm)) +
  geom_boxplot(outliers = FALSE) +
  facet_wrap(~ class, scale = 'free', 
             ncol = 1) +
  theme_bw() +
  ylab('Normalisation method') + xlab(m) +
  #ggtitle('Distribution of coefficients of variation in VOC abundance') +
  theme(plot.title = element_text(hjust = 0.4, size = 11))}


plot_disp('rsd')

ggsave('CV_boxplots_normalisation_v04_noB5_rMAD.tiff', dpi = 300, units = 'mm', width = 90, height = 170)

ggsave('CV_boxplots_normalisation_v03_noB5_S1S2.tiff', dpi = 300, units = 'mm', width = 90, height = 70)

ggsave('CV_boxplots_normalisation_B2_rIQR.tiff', dpi = 300, units = 'mm', 
       width = 90, height = 170)


#
#
#

# PRINCIPAL COMPONENT ANALYSIS
PCA_plot <- function(x, y) { # x is a dataframe of normalised values, y is the name of normalisation method
  #rownames(x) <- 1:nrow(x)
  plotIndiv(mixOmics::pca(log(x[,-c(1:4)]), scale = TRUE, center = TRUE, ncomp = 3),
            pch = 1,
            cex = 3,
            group = x$Batch,
            title = y,
            #legend = TRUE,
            legend.title = 'Batch',
            comp = c(1,2),
            
            )
  } 


dev.new()
raw_batch <- PCA_plot(test %>% filter(class %ni% c('Blank')), #%>%
                        #filter(Sample %ni% c('220207_35_RADicA_B1_RAD201_367133',
                                                 #'220207_11_RADicA_B1_RAD198_373824')),
                      'Raw')

is_batch <- PCA_plot(b_IS %>% filter(class != 'Blank'), # remove obs 144, 118, 295
         'IS')
cc_class <- PCA_plot(b_cc %>%
                       #filter(class == 'ES'),
                       filter(class %ni% c('Blank')),
         'CC')

cc2_batch <- PCA_plot(b_cc2 %>%
                    #filter(class == 'ES'),
                      filter(class %ni% c('Blank')), 
         'CC2')

pqn_batch <- PCA_plot(b_pqn %>% 
                         filter(class %ni% c('Blank')), 
                       'PQN')

pqn1_batch <- PCA_plot(b1_pqn1 %>% 
                         filter(class %ni% c('Blank')), 
                       'PQN1')

cc_pqn <- PCA_plot(b_cc2_pqn_b1ref %>% 
                     filter(class %ni% c('Blank')), 
                   'CC PQN')

cc2_test <- PCA_plot(test_cc2 %>% filter(class %ni% 'Blank'), 'CC')

# figure displaying class annotation
pca <- arrangeGrob(raw_class$graph, is_class$graph, cc_class$graph, pqn1_class$graph,
                   raw_batch$graph, is_batch$graph, cc_batch$graph, pqn1_batch$graph,
            ncol = 4, nrow = 2)

dev.new()
plot(pca)

ggsave('PCA_normalisation_comp_v03_noB5.tiff', pca, dpi = 300, unit = 'mm', width = 320, height = 150)
ggsave('PCA_normalisation_comp_B2.tiff', pca, dpi = 300, unit = 'mm', width = 320, height = 150)

#
#
#

# trajectory of median VOC level in time following normalisation
sum <- function(df) {
  df %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea') %>%
  group_by(class, Analysis_date) %>%
  summarise(median = median(log(peakArea), na.rm = TRUE))
}

imp1_sum <- sum(b_imp1) %>% mutate(norm = 'Raw')
is_sum <- sum(b_IS) %>% mutate(norm = 'IS')
cc_sum <- sum(b_cc) %>% mutate(norm = 'CC')
#pqn_sum <- sum(b1_pqn) %>% mutate(norm = 'PQN')
pqn1_sum <- sum(b_pqn1) %>% mutate(norm = 'PQN')

df_sum <- rbind(imp1_sum, is_sum, cc_sum, pqn1_sum)

plots <- df_sum %>% 
  filter(class == 'ES') %>%
  group_by(factor(norm, levels = c('Raw', 'IS', 'CC', 'PQN'))) %>%
  do(plots  = ggplot(data = ., aes(x = Analysis_date, y = median)) + geom_line() +
  theme_bw() + geom_line(aes(y = median(median)), colour = 'red', linetype = 'dashed') +
    ggtitle(.$norm) +
    theme(plot.title = element_text(hjust = 0.5)))

pdf('MedianvsTime_norm_comp_v03_noB5.pdf', width = 7, height = 5)
pdf('MedianvsTime_norm_comp_B2.pdf', width = 7, height = 5)
marrangeGrob(plots$plots, nrow = 2, ncol = 2)
dev.off()


# density plot VOC levels following normalisation
dens_plot <- function(df) {
  df %>% pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea') %>%
    filter(class == 'ES') %>% 
    ggplot(aes(x = log(peakArea))) + geom_density()
}

dens_plot(b_imp1)
dens_plot(b_IS)
dens_plot(b_cc2)
dens_plot(b_pqn)

#
#
#

# PQN using blank mean on the same day
blanks_sum <- b1_imp_corr %>%
  filter(class == 'Blank') %>%
  pivot_longer(cols =! c(Sample, Analysis_date, class, Batch),
               names_to = 'comp', values_to = 'Blank') %>%
  group_by(comp, Analysis_date) %>%
  summarise(Blank = mean(Blank, na.rm = TRUE)) %>%
  drop_na()


x <- b1_imp_corr 

pqn2 <- function(x) { 
  xL <- x %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class),
                           names_to = 'comp', values_to = 'peakArea')
  ref <- xL %>% filter(class %in% c('Blank')) %>% drop_na() %>%
    group_by(comp, Analysis_date) %>% summarise(Blank = mean(peakArea, na.rm = TRUE))
  ref2 <- ref %>% group_by(comp) %>% summarise(median = median(Blank))
  # calculate quotients of each observation and reference spectrum
  quot <- xL %>% filter(comp %in% unique(ref$comp)) %>% 
    filter(class %ni% 'Blank') %>%
    filter(Analysis_date %in% ref$Analysis_date) %>%
    left_join(ref) %>%
    #left_join(ref2) %>%
    mutate(quotient = peakArea/Blank) %>% 
    #group_by(Sample) %>%
    #summarise(sizeEf = median(quotient, na.rm = TRUE)) # take median of quotients for each sample (size effect)
    dplyr::select(!c(peakArea, Blank)) %>%
    pivot_wider(names_from = comp, values_from = quotient)
   
  
  
  # divide each observation by the sample size effect
  #b1_pqn <- xL %>%
    #filter(comp %in% unique(ref$comp)) %>% 
    #filter(class %ni% 'Blank') %>%
    #filter(Analysis_date %in% ref$Analysis_date) %>%
    #left_join(quot) %>% 
    #mutate(peakAreaN = peakArea/sizeEf) %>%
    #dplyr::select(!c(peakArea, sizeEf)) %>%
    #rename(peakArea = peakAreaN) %>%
    #pivot_wider(names_from = comp, values_from = peakArea) #%>%
  #dplyr::select(!Internal_Standard)
}

b1_pqn2 <- pqn2(b1_imp_corr)

# subtract background first
b1_imp_corr_L1 <- b1_imp_corr %>% filter(class %in% c('BG', 'S')) %>%
  pivot_longer(cols =! c(Sample, Analysis_date, class, Batch), 
               names_to = 'comp', values_to = 'peakArea') %>%
  pivot_wider(id_cols = c(Sample, Analysis_date, comp),
              names_from = class, values_from = peakArea) %>%
  mutate(Scorr = S - BG)

b1_imp_corr_L1 %>% ggplot() + 
  geom_histogram(aes(x = S), alpha = 0.6) + 
  geom_histogram(aes(x = Scorr), fill = 'blue', alpha = 0.6) +
  theme_bw()

sd(b1_imp_corr_L1$S)/mean(b1_imp_corr_L1$S)
sd(b1_imp_corr_L1$Scorr)/mean(b1_imp_corr_L1$Scorr)

# pqn norm
b1_imp_Scorr <- b1_imp_corr_L1 %>% dplyr::select(!S) %>%
  pivot_longer(cols = c(BG, Scorr), names_to = 'class', values_to = 'peakArea') %>%
  mutate(class = ifelse(class == 'Scorr', 'S', class)) %>%
  pivot_wider(names_from = comp, values_from = peakArea) %>%
  left_join(b1_imp_corr %>% dplyr::select(Sample, Analysis_date, class, Batch)) %>%
  relocate(Batch) %>%
  filter(Batch != 5) %>%
  rbind(b1_imp_corr %>% filter(class %in% c('Blank', 'ES')))

b1_pqn_corr <- pqn1(b1_imp_Scorr)

b1_pqn_corr_L <- b1_pqn_corr %>% 
  filter(class %in% c('BG','S')) %>%
  pivot_longer(cols = ! c(Sample, Batch, Analysis_date, class),
               names_to = 'comp', values_to = 'peakArea') %>%
  pivot_wider(names_from = class, values_from = peakArea)

sd(b1_pqn_corr_L$S)/mean(b1_pqn_corr_L$S)

#

corr_pqn_cv <- b1_pqn_corr_L %>% group_by(comp) %>%
  summarise(cvS = sd(S)/abs(mean(S)))

cv <- b1_imp_corr_L1 %>% group_by(comp) %>%
  summarise(cvS = sd(S)/abs(mean(S)),
            cvScorr = sd(Scorr)/abs(mean(Scorr)))

quantile(cv$cvS)
quantile(cv$cvScorr)
quantile(corr_pqn_cv$cvS)

#
#
#

# pca plot with analysis date
out <- c('220207_35_RADicA_B1_RAD201_367133',
         '220207_11_RADicA_B1_RAD198_373824')

out <- NULL

date_pca <- function(data, norm){
  data <- as.data.frame(data)
  rownames(data) <- data$Sample
  pc <- mixOmics::pca(data %>% 
                        filter(class %ni% c('ES','Blank')) %>%
                        filter(Sample %ni% out) %>%
                        dplyr::select(!c(1:4)) %>%
                        log(),
                      scale = TRUE, center = TRUE)
  scores <- pc$variates$X %>% as.data.frame() %>% 
    mutate(Sample = rownames(.)) %>%
    left_join(data %>% dplyr::select(Sample, Analysis_date, class)) %>%
    mutate(class = ifelse(class %in% c('S1', 'S2'), 'S', class))
  scores$Analysis_date <- as.Date(scores$Analysis_date)
  scores %>%
    ggplot(aes(x = PC1, y = PC2)) + 
    geom_point(aes(colour = as.numeric(Analysis_date),
                   shape = class)) +
    scale_colour_gradientn(colours = c('yellow', 'blue'), 
                           labels = as.Date, 
                           name = 'Analysis date') +
    theme_bw() +
    xlab(paste(round(pc$prop_expl_var$X[1], 3)*100, '%')) +
    ylab(paste(round(pc$prop_expl_var$X[2], 3)*100, '%')) +
    ggtitle(norm) +
    theme(plot.title = element_text(hjust = 0.5))
}

date_pca(b_imp, 'b2 raw') + xlim(-25, NA)

ggsave('PCA_date_B2_CC2_PQN.tiff', dpi = 300, width = 120, height = 90, unit = 'mm')

data <- b_cc2

b_cc2_loadings <- as.data.frame(pc$loadings$X)

write.csv(b_cc2_loadings, 'B2_CC2_Blank_loadings.csv')


