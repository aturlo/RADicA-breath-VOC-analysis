## Normalisation of VOC peak area data

# author: Aggie Turlo
# project: RADicA
# date: 20/01/2025

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(mixOmics)
library(pcaPP)
library(stringr)
library(grafify)
library(grid)
library(gridExtra)


## load data
# load study data with imputed values
b1_imp <- read.csv('RADicA_B1_NAfiltered_imputed.csv')[,-1]
b2_imp <- read.csv('RADicA_B2_NAfiltered_imputed.csv')[,-1]

rownames(b1_imp) <- b1_imp$Sample
rownames(b2_imp) <- b2_imp$Sample

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

################################

# format imputed dataset
# specify dataset for analysis
b_imp <- b2_imp

## replace analysis date recorded by GCMS with collection date
b_imp_L <- b_imp %>% pivot_longer(cols = c(6:ncol(b_imp)), names_to = 'comp', values_to = 'peakArea') %>%
  mutate(peakArea = ifelse(peakArea < 10, 10, peakArea)) # replace imputed peak area values < 10 with 10


b_imp <- b_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea) %>%
  mutate(Analysis_date = paste(str_sub(Sample, start = 5, end = 6),
                               str_sub(Sample, start = 3, end = 4),
                               20, str_sub(Sample, start = 1, end = 2),
                               sep = '')) %>%
  relocate(Analysis_date) %>%
  dplyr::select(!c(Date, Time))

b_imp$Analysis_date <- as.Date(b_imp$Analysis_date, format = '%d%m%Y')
rownames(b_imp) <- b_imp$Sample

#
#
#

b1_imp <- b_imp
b2_imp <- b_imp


#
#
#

###############################

# NORMALISATION OF IMPUTED DATA
# specify dataset for analysis
b_imp <- b2_imp
dat <- 'B2'

b_imp_L <- b_imp %>% pivot_longer(cols = c(6:ncol(b_imp)), names_to = 'comp', values_to = 'peakArea')

# Internal Standard normalisation (IS)
b_IS <- cbind(b_imp[,1:4], b_imp[,5:ncol(b_imp)]/b_imp$Internal_Standard) %>%
  dplyr::select(!Internal_Standard)

#
#
#

## Component correction (CC)
# remove VOCs not represented in reference sample class in either dataset
data <- b2_imp

data <- as.data.frame(data)
rownames(data) <- data$Sample

# subset reference dataset (here blank samples)
assay_es_imp <- data %>% filter(class %in% c('Blank')) %>% 
  dplyr::select(!c(Sample, Analysis_date, class, Batch))

# remove features missing in the reference dataset
nas <- which(is.na(assay_es_imp), arr.ind = TRUE)
table(nas[,2])

nacols <- which(colSums(is.na(assay_es_imp) == TRUE) > 0)
nas_voc <- colnames(assay_es_imp[nacols])

assay_es_imp <- assay_es_imp[,-nacols]


# remove outliers from reference datasets (es and blank)
# identify multivariate outliers using robust PCA
pc <- PCAgrid(log(assay_es_imp), scale = 'sd', center = 'mean')

#  flag outliers based on the critical OD and SD values (0.97 quantile)
sdod <- PCdiagplot(assay_es_imp, pc, plotbw = FALSE)
sd <- as.data.frame(sdod$SDist)
rownames(sd) <- rownames(assay_es_imp)
sdOut <- sd %>% filter(V1 > sdod$critSD[1,1] | V2 > sdod$critSD[2,1])

od <- as.data.frame(sdod$ODist)
rownames(od) <- rownames(assay_es_imp)
odOut <- od %>% filter(V1 > sdod$critOD[1,1] | V2 > sdod$critOD[2,1])

intersect(rownames(sdOut), rownames(odOut))

# remove outliers from the QC database
assay_es_imp1 <- assay_es_imp %>% filter(rownames(.) %ni% rownames(odOut)) %>%
  filter(rownames(.) %ni% rownames(sdOut))

# compare PCA scores plots before and after outlier removal
annot_es <- data %>% filter(class == 'Blank')
annot_es1 <- data %>% filter(class == 'Blank') %>%
  filter(rownames(.) %ni% rownames(odOut)) %>%
  filter(rownames(.) %ni% rownames(sdOut))

pc_es <- pca(log(assay_es_imp), scale = TRUE, center = TRUE, ncomp = 2)
pc_es1 <- pca(log(assay_es_imp1), scale = TRUE, center = TRUE, ncomp = 2)

#

plotIndiv(pc_es,
          group = annot_es$Batch,
          pch = 1,
          legend = TRUE)

plotIndiv(pc_es1,
          group = annot_es1$Batch, 
          pch = 1,
          legend = TRUE)

# extract loadings from PCA w/o outliers
loadings <- pc_es1$loadings[['X']] %>% as.data.frame() %>% 
  as.matrix()

# subset all observations apart from reference class
Y <- data %>% 
  filter(class %ni% c('Blank')) %>%
  dplyr::select(colnames(assay_es_imp1)) %>% 
  log() %>% as.matrix() %>%
  scale(scale = TRUE, center = TRUE)

#

x0 <- Y
x0[is.na(Y)] <- 0 # managment of NAs for matrix multiplication


# calculate drift using reference PCA loadings 
# 1 component
drift1 <- (x0 %*% loadings[,1]) %*% t(loadings[,1])
drift1[is.na(Y)] <- NA

# 2 components
drift <- (x0 %*% loadings) %*% t(loadings)
drift[is.na(Y)] <- NA

# subtract drift from the dataset
Z1 <- Y - drift1
Z <- Y - drift

# reverse scaling
Zres1 <- Z1
for (i in c(1:dim(Y)[2])){
  Zres1[,i] <- (Z1[,i]*attributes(Y)$`scaled:scale`[i]) + attributes(Y)$`scaled:center`[i]
}

Zres <- Z
for (i in c(1:dim(Y)[2])){
  Zres[,i] <- (Z[,i]*attributes(Y)$`scaled:scale`[i]) + attributes(Y)$`scaled:center`[i]
}


#
annot <- data %>% 
  filter(class %ni% c('Blank'))

b_cc <- cbind(annot[,1:4], exp(Zres1))
b_cc2 <- cbind(annot[,1:4], exp(Zres)) # reverse log transformation


#
#
#

# PQN normalisation
# using median peak intensity of Blank and ES datasets as reference 
# use only ES/blanks as reference
pqn <- function(ref, data) { 
  xL <- ref %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class),
                           names_to = 'comp', values_to = 'peakArea')
  ref_spectrum <- xL %>% filter(class %in% c('Blank', 'ES')) %>%
    group_by(comp) %>% summarise(median = median(peakArea, na.rm = TRUE))
  # format dataset for normalisation
  yL <- data %>% pivot_longer(cols =! c(Sample, Analysis_date, Batch, class), 
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
    pivot_wider(names_from = comp, values_from = peakArea) 
}


b_pqn <- pqn(b_imp, b_imp) # change to b2_imp for B1 - use the same reference for both

#
#
#

# CC + PQN
b2_cc2 <- b_cc2

b2_cc2_pqn <- pqn(b2_cc2, b2_cc2) # change to b2_cc2 for B1 - use the same reference for both
b1_cc2_pqn <- pqn(b2_cc2, b_cc2)

plotIndiv(pca(log(b2_cc2_pqn[,-c(1:4)]), scale = TRUE), pch = 1)

# save normalised dataset
write.csv(b1_cc2_pqn, 'RADicA_B1_NAfiltered_imputed_CC2_PQN.csv') # amend file name depending on dataset
write.csv(b2_cc2_pqn, 'RADicA_B2_NAfiltered_imputed_CC2_PQN.csv')
#
#
#

##################################

# ASSESSMENT OF NORMALISATION OUTCOMES
# Measures of relative dispersion
disp <- function(x1, y){ # x is dataframe of normalised values, y is normalisation method
  x1 %>%
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

#

disp_raw <- disp(b_imp, 'Raw')
disp_IS <- disp(b_IS, 'IS')
disp_cc <- disp(b_cc, 'CC')
disp_cc2 <- disp(b_cc2, 'CC2')
disp_pqn <- disp(b_pqn, 'PQN')
disp_cc2_pqn <- disp(b2_cc2_pqn, 'CC2 PQN')

#

disp_conc <- rbind(disp_raw, disp_IS, #disp_cc,
                   disp_cc2, disp_pqn, disp_cc2_pqn)

# visualise distribution of relative dispersion values following different normalisation methods
pals <- grafify::graf_palettes
my_pal = pals$fishy

plot_disp <- function(m) { # m = dispersion measure: RSD, rIQR, rMAD
  disp_conc %>% 
  filter(class %ni% c('Blank')) %>%
  ggplot(aes(x = round(.data[[m]], 1), y = factor(norm, levels = c('CC2 PQN', 'CC2','CC', 'PQN', 'IS', 'Raw')), 
             fill = norm)) +
  geom_boxplot(outliers = FALSE, lwd = 0.3) +
  facet_wrap(~ class, scale = 'free', 
             ncol = 1) +
  theme_bw(base_size = 10) +
  ylab('Normalisation method') + xlab(m) +
  theme(plot.title = element_text(hjust = 0.4, size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(0, 'mm'),
        axis.text = element_text(size = 6),
        panel.grid = element_blank()) +
    scale_fill_grafify(palette = 'fishy')}


rsd_p <- plot_disp('rsd') + theme(legend.position = 'none') #+ 
  ggtitle('Relative standard deviation')

rmad_p <- plot_disp('rMAD') + theme(legend.position = 'none') #+ 
                                      ggtitle('Relative median absolute deviation')

riqr_p <- plot_disp('rIQR') + theme(legend.position = 'none',
                                    axis.title.y = element_blank()) #+
  ggtitle('Relative interquartile range')

disp_plots <- arrangeGrob(rsd_p, rmad_p, riqr_p, nrow = 1, widths = c(0.29, 0.29, 0.42))

dev.new()
plot(disp_plots)

ggsave(paste('Boxplots_normalisation_', dat, '.tiff', sep =''), disp_plots, dpi = 300, units = 'mm', width = 240, height = 170)

#
#
#

###############################

# Figure S5B
s5b_b2 <- arrangeGrob(rmad_p, riqr_p, nrow = 1, top = textGrob('Training dataset', vjust = 1.2,
                                                               gp = gpar(fontsize = 10)),
                      widths = c(0.53, 0.47))
s5b_b1 <- arrangeGrob(rmad_p, riqr_p, nrow = 1, top = textGrob('Validation dataset', vjust = 1.2,
                                                               gp = gpar(fontsize = 10)),
                      widths = c(0.53, 0.47))

s5b_f <- ggplot() + 
  geom_rect(data = disp_conc, mapping = aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), 
            fill = 'white', alpha = 1) + theme(axis.text = element_blank(),
                                               axis.ticks = element_blank())
  

s5b <- arrangeGrob(s5b_b2, s5b_f, s5b_b1, nrow = 1, widths = c(0.475,  0.05, 0.475),
                   top = textGrob('Measures of dispersion in VOC datasets',
                                  gp = gpar(fontsize = 10, fontface = 'bold')))
plot(s5b)


#
#
#

# PRINCIPAL COMPONENT ANALYSIS
pals <- grafify::graf_palettes
my_pal = pals$fishy

PCA_plot <- function(x, y) { # x is a dataframe of normalised values, y is the name of normalisation method
  #rownames(x) <- 1:nrow(x)
  plot <- plotIndiv(mixOmics::pca(log(x[,-c(1:4)]), scale = TRUE, center = TRUE, ncomp = 3),
            pch = 1,
            cex = 0.4,
            group = x$class,
            title = y,
            #legend = TRUE,
            legend.title = 'class',
            comp = c(1,2),
            size.title = 8,
            size.xlabel = 8,
            size.ylabel = 8,
            size.legend.title = 6,
            size.legend = 6,
            size.axis = 6)
  
  plot_fin <- plot$graph + theme(strip.background = element_blank()) + 
    theme(strip.background = element_blank(),
          plot.title = element_text(face = 'plain'),
          plot.margin = unit(c(0, 0.1, 0, 0.1), 'cm'),
          panel.grid = element_blank()) + 
    scale_colour_manual(values = c('ES' = my_pal[[6]],
                                   'BG' = my_pal[[7]],
                                   'S1' = my_pal[[8]],
                                   'S2' = my_pal[[9]]))
  
  } 


raw_class <- PCA_plot(b_imp %>% filter(class %ni% c('Blank')), 
                      'Raw')

is_class <- PCA_plot(b_IS %>% filter(class != 'Blank'), 
         'IS')

cc_class <- PCA_plot(b_cc %>%
                       filter(class %ni% c('Blank')),
         'CC')

cc2_class <- PCA_plot(b_cc2 %>%
                      filter(class %ni% c('Blank')), 
         'CC2')

pqn_class <- PCA_plot(b_pqn %>% 
                         filter(class %ni% c('Blank')), 
                       'PQN')


cc_pqn_class <- PCA_plot(b2_cc2_pqn %>% 
                     filter(class %ni% c('Blank')), 
                   'CC2 PQN')



# figure displaying class annotation
pca <- arrangeGrob(raw_class, is_class, pqn_class, cc_class, 
                   cc2_class, cc_pqn_class,
            ncol = 3, nrow = 2)

dev.new()
plot(pca)

ggsave(paste('PCA_normalisation_comp_', dat, '.tiff', sep = ''), pca, dpi = 300, unit = 'mm', width = 250, height = 150)

#
#
#

################################

# Figure S5C
leg5c <- get_legend(cc_pqn_class)

s5c_b2 <- arrangeGrob(raw_class, is_class, pqn_class, 
                      cc2_class, cc_pqn_class, leg5c,
                      ncol = 2, nrow = 3,
                      top = textGrob('Training dataset', vjust = 1,
                                     gp = gpar(fontsize = 10)))

plot(s5c_b2)

s5c_b1 <- arrangeGrob(raw_class, is_class, pqn_class, 
                      cc2_class, cc_pqn_class, leg5c,
                      ncol = 2, nrow = 3,
                      top = textGrob('Validation dataset', vjust = 1,
                                     gp = gpar(fontsize = 10)))

plot(s5c_b1)

s5c <- arrangeGrob(s5c_b2, s5b_f, s5c_b1, nrow = 1, widths = c(0.475, 0.05, 0.475),
                   top = textGrob('PCA score plots of VOC datasets',
                                  gp = gpar(fontsize = 10, fontface = 'bold')))
plot(s5c)


#

###############################

# Assemble FigureS5

s5 <- plot_grid(s5c, s5b, nrow = 2, rel_heights = c(0.55, 0.45),
                labels = c('a)', 'b)'), label_size = 12)
plot(s5)

ggsave('figs5.tiff', s5, unit = 'mm', dpi = 300, width = 157, height = 220)

pdf('FigureS6.pdf', width = 6.18, height = 8.66)
plot(s5)
dev.off()


#
#
#

# PCA plot with analysis date
date_pca <- function(data, norm){
  data <- as.data.frame(data)
  rownames(data) <- data$Sample
  pc <- mixOmics::pca(data %>% 
                        filter(class %ni% c('Blank')) %>%
                        #filter(Sample %ni% out) %>%
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
    xlab(paste('PC1:', round(pc$prop_expl_var$X[1], 3)*100, '%')) +
    ylab(paste('PC2:', round(pc$prop_expl_var$X[2], 3)*100, '%')) +
    ggtitle(norm) +
    theme(plot.title = element_text(hjust = 0.5))
}

t_raw <- date_pca(b_imp, 'Raw') + theme(legend.position = 'none')
t_cc2_pqn <- date_pca(b_cc2_pqn, 'CC2 PQN')

t_plots <- arrangeGrob(t_raw, t_cc2_pqn, nrow = 1, widths = c(0.43, 0.57)) 
plot(t_plots)

ggsave(paste('PCA_date_norm_', dat, '.tiff', sep = ''), t_plots, dpi = 300, width = 200, height = 90, unit = 'mm')

#
#
#

# Effect of normalisation on standard peak area
std1 <- intersect(colnames(b_imp), std)

es_data <- b_imp[,c('Sample', 'Analysis_date', std1, 'class', 'Internal_Standard')] %>%
  filter(class == 'ES') %>% mutate(norm = 'Raw') %>%
  rbind(b_cc2_pqn[,c('Sample', 'Analysis_date', std1, 'class', 'Internal_Standard')] %>%
          filter(class == 'ES') %>% mutate(norm = 'CC2 PQN'))


# B2 only - remove poor quality outlier observations
es_data <- es_data %>% filter(Analysis_date %ni% c('2021-02-12', '2022-06-29', '2022-02-07', '2021-09-24'))

#

es_plots <- lapply(std1, function(voc) {
  es_voc <- es_data %>% rename(out = voc)
  
  es_voc %>%
  ggplot(aes(y = log(out), x = Analysis_date, group = Analysis_date)) + 
  geom_point() +
  facet_grid(~ factor(norm, levels = c('Raw', 'CC2 PQN'))) +
  theme_bw(base_size = 8) +
  ggtitle(voc) +
  ylab('log(peak area)') +
  theme(plot.title = element_text(hjust = 0.5))
  
})

pdf(paste('ES_normalisation_effect_time_', dat, '.pdf'))
marrangeGrob(es_plots, nrow = 3, ncol = 1)
dev.off()


#
#
#

# compare density distributions of CC2 and CC2 PQN normalised datasets
b1_cc2 <- b_cc2
b1_cc2_pqn <- b_cc2_pqn

b2_cc2 <- b_cc2
b2_cc2_pqn <- b_cc2_pqn

norm_dens <- function(data1, data2){
  data1 %>% mutate(data = 'B2') %>% 
  filter(class %in% c('S1', 'S2', 'BG')) %>%
  dplyr::select(!c(1:4)) %>%
  pivot_longer(cols =! data, names_to = 'comp', values_to = 'peakArea') %>%
  drop_na() %>%
  rbind(# compare density distributions of cc-normalised datasets
    data2 %>% mutate(data = 'B1') %>% 
      filter(class %in% c('S1', 'S2', 'BG')) %>%
      dplyr::select(!c(1:4)) %>%
      pivot_longer(cols =! data, names_to = 'comp', values_to = 'peakArea') %>%
      drop_na()) %>%
  ggplot(aes(x = log(peakArea), colour = data)) + geom_density() +
  theme_bw()
  }

cc2_dens <- norm_dens(b1_cc2, b2_cc2) + ggtitle('CC2') + theme(legend.position = 'none')
cc2_pqn_dens <- norm_dens(b1_cc2_pqn, b2_cc2_pqn) + ggtitle('CC2 PQN')

dens_plots <- arrangeGrob(cc2_dens, cc2_pqn_dens, nrow = 1, widths = c(0.43, 0.57))
plot(dens_plots)

ggsave('Density_B1_B2_CC2.tiff', dens_plots, unit = 'mm', dpi = 300, width = 140, height = 70)

#
#
#
