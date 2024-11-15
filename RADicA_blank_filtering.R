## Feature filtering based on comparison with blank samples

# author: Aggie Turlo
# project: RADicA
# date: 17/07/2024

#####################

library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(tibble)
library(factoextra)
library(ggfortify)
library(lme4)
library(gridExtra)

# load relevant datasets (in wide and long format)
b1_imp_L <- read.csv('RADicA_B1_NAfiltered_imputed_long.csv') # NA filtered, imputed

b1_imp_es <- b1_imp_L %>% filter(class == 'ES')
unique(b1_imp_es$comp)

b1_imp <- b1_imp_L %>% pivot_wider(names_from = comp, values_from = peakArea)

b1_all_fL # NA filtered
b1_all_f

# dataset with summarised breath values
b1_imp_L_sum <- read.csv('RADicA_B1_NAfiltered_imputed_PQN_summarised_long_ICC.csv')[,-1]

b1_imp_sum <- b1_imp_L_sum %>% pivot_wider(names_from = comp, values_from = peakArea) %>%
  dplyr::select(!Internal_Standard)

# metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

meta$Analysis_date <- as.Date(meta$Analysis_date, format = '%d/%m/%Y')

# VISUALISE CONDITIONAL DISTRIBUTIONS OF COMPOUNDS
x <- b1_imp_L_sum

pdf('Hist_PQN_sum_ICC_bgVSblank.pdf')
hists <- x %>%
  filter(comp %in% b_s_Ga_f$comp) %>%
  filter(class %in% c('Blank', 'S')) %>%
  group_by(comp) %>%
  do(plots = ggplot(data = ., aes(x = log(peakArea), fill = class, alpha = 0.5)) + 
       geom_histogram(position = 'identity') + 
       ggtitle(.$comp) + theme(#legend.position = 'none',
                               plot.title = element_text(size = 10)))

marrangeGrob(hists$plots, nrow = 3, ncol = 2)
dev.off()


#
#
#

# test if blank and sample observations come from the same Gamma distribution
# argument x is a dataframe in long format
# NOT CONVERGING
x <- b1_imp_L_sum

gamma_blank <- function(CL){
  subset <- x %>% filter(class %in% c('Blank', CL))
  blank_comps <- x %>% filter(class == 'Blank') %>% drop_na() %>% 
    dplyr::select(comp) %>% unique()
  results <- bind_rows(lapply(blank_comps$comp ,function(voc) {
    subset1 <- subset %>% filter(comp == voc) %>% drop_na() %>%
      mutate(peakArea = peakArea/sd(peakArea, na.rm = TRUE))
    test <- glmer(peakArea ~ class + (1 | Analysis_date),
                family = 'Gamma'(link = 'log'),
                data = subset1)
    df <- c(voc, coef(summary(test))[2,]) %>% as.data.frame() %>% t() %>% as.data.frame()
  }))
  colnames(results)[1] <- 'comp'
  colnames(results)[5] <- 'p.value'
  results <- results %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
}

#

b_bg_G <- gamma_blank('BG') %>% filter(adj.p.value > 0.05) 
b_s_G <- gamma_blank('S') %>% filter(adj.p.value > 0.05)

intersect(b_s2_G$comp, b_s1_G$comp)

remove <- intersect(b_s2_G$comp, b_s1_G$comp)
remove <- remove[remove != 'Internal_Standard']

# in summarised dataset
remove <- b_s_G$comp
remove <- remove[remove != 'Internal_Standard']

# remove variables with equal distribution in Blank/S
b1_imp_corr <- b1_imp %>% dplyr::select(!all_of(remove))
b1_all_f_corr <- b1_all_f %>% dplyr::select(!all_of(remove))

b1_imp_corr <- b1_imp_sum %>% dplyr::select(!all_of(remove))

write.csv(b1_imp_corr, 'RADicA_B1_NAfiltered_sum_blank_filtered.csv')

#
#
#

# 
# using linear regression and log-transformed values
# simple model
x <- b1_imp_L_sum

gauss_blank <- function(CL){
  subset <- x %>% filter(class %in% c('Blank', CL))
  blank_comps <- x %>% filter(class == 'Blank') %>% drop_na() %>% 
    dplyr::select(comp) %>% unique()
  results <- bind_rows(lapply(blank_comps$comp ,function(voc) {
    subset1 <- subset %>% filter(comp == voc) %>% drop_na() %>%
      mutate(logPeakArea = log(peakArea))
    test <- lmer(logPeakArea ~ class + (1 | Analysis_date),
                data = subset1)
    df <- coef(summary(test))[2,] %>% as.data.frame() %>% t() %>% as.data.frame() %>%
      mutate(comp = voc)
  }))
  colnames(results)[5] <- 'p.value'
  results 
}

b_s_Ga <- gauss_blank('S')

b_s_Ga_f <- b_s_Ga %>% 
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(adj.p.value > 0.05) 

# model using summarised paired blank values and patient information (ID, Diagnosis)

x <- wilcox_input %>% filter(class == 'S') %>%
  pivot_wider(names_from = class, values_from = peakArea) %>%
  mutate(RAD_ID = str_sub(Sample, start = 15)) %>%
  left_join(meta %>% dplyr::select(RAD_ID, Analysis_date, Diagnosis)) %>%
  pivot_longer(cols = c(S, Blank), names_to = 'class', values_to = 'peakArea')

gauss_blank1 <- function(CL){
  subset <- x %>% filter(class %in% c('Blank', CL))
  blank_comps <- x %>% filter(class == 'Blank') %>% drop_na() %>% 
    dplyr::select(comp) %>% unique()
  results <- bind_rows(lapply(blank_comps$comp ,function(voc) {
    subset1 <- subset %>% filter(comp == voc) %>% drop_na() %>%
      mutate(logPeakArea = log(peakArea)) %>%
      filter(logPeakArea > mean(logPeakArea) - 3*sd(logPeakArea)) %>%
      filter(logPeakArea < mean(logPeakArea) + 3*sd(logPeakArea))
    test <- lmer(logPeakArea ~ class*Diagnosis + (1 | Analysis_date) + (1 | RAD_ID),
                 data = subset1)
    df <- coef(summary(test))[2:4,] %>% as.data.frame() %>%
      mutate(comp = rep(voc, 3),
             coef = rownames(.))
  }))
  colnames(results)[5] <- 'p.value'
  results 
}

b_s_Ga1 <- gauss_blank1('S')

b_s_Ga1_f <- b_s_Ga1 %>% filter(coef == 'classS') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(adj.p.value > 0.05)

b_s_Ga1_f1 <- b_s_Ga1 %>% filter(coef == 'classS:DiagnosisNot Asthma') %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) %>%
  filter(comp %in% b_s_Ga1_f$comp) %>%
  filter(adj.p.value < 0.05)


results <- results %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

b1_imp_corr <- b1_imp_sum %>% dplyr::select(!all_of(b_s_Ga$comp))

write.csv(b1_imp_corr, 'RADicA_B1_NAfiltered_imputed_PQN_sum_blank_filtered.csv')

#

x %>%
  filter(comp %in% b_s_Ga_f1$voc) #%>%
  ggplot(aes(x = log(peakArea), fill = class)) + geom_histogram(position = 'identity') +
  facet_wrap(~comp)
  

#

# using Wilcoxon-Mann-Whitney test
# ASSUMES INDEPENDENCE BETWEEN OBSERVATION PAIRS

# input for paired test
blanks_sum <- b1_imp_sum %>%
  filter(class == 'Blank') %>%
  pivot_longer(cols =! c(Sample, Analysis_date, class, Batch),
               names_to = 'comp', values_to = 'Blank') %>%
  group_by(comp, Analysis_date) %>%
  summarise(Blank = mean(Blank, na.rm = TRUE)) %>%
  drop_na()

wilcox_input <- b1_imp_L_sum %>% 
  filter(class %in% c('S', 'BG')) %>%
  left_join(blanks_sum %>% dplyr::select(Analysis_date, comp, Blank),
            by = c('Analysis_date', 'comp')) #%>%
  #mutate(diff = peakArea - Blank) 


wilcox_blank <- function(CL){
  subset_cl <- x %>% filter(class == CL)
  subset_blank <- x %>% filter(class == 'Blank') %>%
    filter(CoreVisit == 'CV1')
  blank_comps <-  x %>% filter(class == 'Blank') %>% drop_na() %>% 
    dplyr::select(comp) %>% unique()
  results <- bind_rows(lapply(blank_comps$comp ,function(y) {
    subset_cl1 <- subset_cl %>% 
      filter(comp == y) %>%
      drop_na()
    subset_blank1 <- subset_blank %>% 
      filter(comp == y) %>%
      drop_na()
    test <- wilcox.test(subset_cl1$peakArea, subset_blank1$peakArea)
    df <- data.frame(comp = y, p.value = test[["p.value"]])
  }))
  results <- results %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
}

W_test <- wilcox_blank('S')
W_test_f <- W_test %>% filter(adj.p.value > 0.05)



# BH adjustment by hand
results <- results %>% arrange(p.value) %>% mutate(rank = c(1:nrow(.))) %>%
  mutate(BHp = p.value*nrow(.)/rank)
  





# plot relationship between breath, bg and blank distributions
hist_BGvsBlankvsS <- wilcox_input %>% 
  pivot_wider(names_from = class, values_from = peakArea) %>%
  pivot_longer(cols = c(Blank, BG, S), names_to = 'class', values_to = 'peakArea') %>%
  group_by(comp) %>%
  do(plots = ggplot(data =., aes(x =  log(peakArea), fill = class)) + 
       geom_histogram() +
       ggtitle(.$comp) +
       theme_bw(base_size = 8) +
       theme(legend.position = 'none'))

marrangeGrob(hist_BGvsBlankvsS$plots, ncol = 4, nrow = 4)

wilcox_input %>% 
  pivot_wider(names_from = class, values_from = peakArea) %>%
  pivot_longer(cols = c(Blank, BG, S), names_to = 'class', values_to = 'peakArea') %>%
  filter(comp == 'Azetidine') %>%
  ggplot(aes(x =  log(peakArea), fill = class)) + 
       geom_histogram() +
       theme_bw(base_size = 8)
