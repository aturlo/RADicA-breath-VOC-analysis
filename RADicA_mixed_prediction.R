## Predicting asthma based on clinical outcomes and breath VOCs

# author: Aggie Turlo
# project: RADicA
# date: 25/02/2025

#####################


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(ranger)
library(pROC)
library(caret)
library(matrixStats)
library(optRF)
library(ggbeeswarm)
library(scales)
library(irr)
library(forcats)
library(foreign)
library(stringr)

# load data
# load summarised BG corrected datasets w/o multivariate outliers
b1_corr_out <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b2_corr_out <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')[,-1]

# load clinical data
clin <- read.spss('RADicA Active 29052024.sav', to.data.frame = TRUE)

#
#
#

# custom functions
ber_fun <-  function(conf_mat) {
  ber <- 0.5*((conf_mat[1,2]/(conf_mat[1,2]+conf_mat[1,1])) + 
                (conf_mat[2,1]/(conf_mat[2,1] + conf_mat[2,2])))
  names(ber) <- 'BER'
  ber
}

er_fun <- function(conf_mat){
  er <- (sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
  names(er) <- 'ER'
  er
}

sens_fun <- function(conf_mat){
  sens <- conf_mat[1,1]/(conf_mat[1,1] + conf_mat[1,2])
  names(sens) <- 'Sensitivity'
  sens
}

spec_fun <- function(conf_mat){
  spec <- conf_mat[2,2]/(conf_mat[2,2] + conf_mat[2,1])
  names(spec) <- 'Specificity'
  spec
}

'%ni%' <- Negate('%in%')

#

# random partitioning of train data into M folds
fold_fun <- function(data, nfold){
  p <- n_distinct(data$RAD_ID)
  x <- sample(1:p)
  names(x) <- unique(data$RAD_ID)
  levels <- nfold
  folds <- split(x, x%%levels)
}

# train model on nfold - 1 and test on the left-out fold
train_test <- function(data, folds, m, nvar = 3, n = 1000){
  test_id <- folds[[m]]
  test_data <- data %>% filter(RAD_ID %in% names(test_id))
  train_data <- data %>% filter(RAD_ID %ni% names(test_id))
  
  model <-  forest <- ranger(Diagnosis ~ FeNOcat + BDRcat + BCTcat, #+ 
                               #Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                             data = train_data,
                             num.trees = n)
  
  model_pred <- predict(model,
                        data = test_data)
  
  test_conf <- table(factor(test_data$Diagnosis, levels = levels(as.factor(test_data$Diagnosis))), 
                     model_pred$predictions)
  
}

#
#
#

# keep clinical outcomes of interest
clin_core <- clin %>% dplyr::select(BBID, Ethnicity, MaleFemale, AgeInYears, SmokerPQ, 
                                    SmokerTypePQ, FeNOCV1, FeNOCV2, 
                                    FEV1PreCV1, FEV1PostCV1, FEV1PPrePredCV1, 
                                    FEV1PPostPredCV1, PD20CV2, ManPosChalOV1) 

# format clinical data 
clin_core$BBID <- str_sub(clin_core$BBID, end = -2)

# keep only patients with breath VOC data available
clin_voc <- clin_core %>% rename(RAD_ID = BBID) %>%
  filter(RAD_ID %in% c(b1_corr_out$RAD_ID, b2_corr_out$RAD_ID))

# examine availability of FeNO and BCT data
# FeNO recorded between core visits
View(clin_voc %>% mutate(FeNO = ifelse(is.na(FeNOCV1) == FALSE & is.na(FeNOCV2 == FALSE), 'both', 'one')))


# keep only patients with BCT data available
clin_voc <- clin_voc %>% filter(PD20CV2 != '         ')

# calculate bronchodilator reversibility
clin_voc <- clin_voc %>% mutate(BDR = ((FEV1PostCV1-FEV1PreCV1)/FEV1PreCV1)*100,
                                BDRPPre = FEV1PPostPredCV1-FEV1PPrePredCV1)

# create categorical outcomes based on NICE thersholds
clin_voc <- clin_voc %>% mutate(FeNOCV1cat = ifelse(FeNOCV1 > 50, 'Pos', 'Neg'),
                                FeNOCV2cat = ifelse(FeNOCV2 > 50, 'Pos', 'Neg'),
                                BDRcat = ifelse(BDR > 12, 'Pos', 'Neg'),
                                BDRPPrecat = ifelse(BDRPPre > 10, 'Pos', 'Neg'),
                                BCTcat = ifelse(PD20CV2 < 0.2, 'Pos', 'Neg'),
                                FeNOcat = ifelse((FeNOCV1 + FeNOCV2)/2 > 50, 'Pos', 'Neg'))


clin_voc$FeNOcat <- ifelse(is.na(clin_voc$FeNOcat) == TRUE & is.na(clin_voc$FeNOCV1cat) == TRUE, clin_voc$FeNOCV2cat, clin_voc$FeNOCV1cat)


# check agreement between FeNO and BDR categories
ag <- clin_voc %>% mutate(FeNOag = ifelse(FeNOCV1cat == FeNOCV2cat, 'Yes', 'No'),
                         BDRag = ifelse(BDRcat == BDRPPrecat, 'Yes', 'No'))
table(ag$FeNOag)
table(ag$BDRag)

#

# examine distribution of best VOC candidate markers
cols1 <- hue_pal()(8)

box_plot <- function(out, center = FALSE, sum = FALSE) {
  b2 <- b2_corr_out %>% rename(OUT = out) %>% mutate(dataset = 'Dataset 2')
  b1 <- b1_corr_out %>% rename(OUT = out) %>% mutate(dataset = 'Dataset 1')
  
  if (sum == TRUE) {
    b2a <- b2 %>% group_by(RAD_ID) %>% summarise(OUT = mean(OUT, na.rm = TRUE)) %>% as.data.frame() %>%
                                                   left_join(b2 %>% dplyr::select(Diagnosis, dataset, RAD_ID))
    b1a <- b1 %>% group_by(RAD_ID) %>% summarise(OUT = mean(OUT, na.rm = TRUE)) %>% as.data.frame() %>%
                                                   left_join(b1 %>% dplyr::select(Diagnosis, dataset, RAD_ID))
    
    ba <- rbind(b1a, b2a)
    
    ba %>% ggplot(aes(x = Diagnosis, y = OUT, fill = Diagnosis)) +
      geom_boxplot(outliers = FALSE) +
      geom_point(alpha = 0.4, size = 1, position = position_jitterdodge(jitter.width = 0.2)) +
      ggtitle(out) +
      facet_wrap(~dataset) +
      scale_fill_manual(values = c('Asthma' = cols1[2],
                                   'Not Asthma' = cols1[6])) +
      theme_bw(base_size = 8) +
      theme() +
      ylab('log (peak area)')
    } else {
    
  if (center == TRUE) {
    b2$OUT <- (b2$OUT-mean(b2$OUT))
    b1$OUT <- (b1$OUT-mean(b1$OUT))
  } else {
    b2 <- b2
    b1 <- b1
  }
  
  b <- rbind(b1, b2)
  
  b %>% ggplot(aes(x = CoreVisit, y = OUT, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_point(position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.2), alpha = 0.4, size = 1) +
    ggtitle(out) +
    facet_wrap(~dataset) +
    scale_fill_manual(values = c('Asthma' = cols1[2],
                                 'Not Asthma' = cols1[6])) +
    theme_bw(base_size = 8) +
    theme() +
    ylab('log (peak area)')
  }
  
}

eb <- box_plot('Ethyl_butanoate', sum = TRUE)

fm <- box_plot('Furan._2_methyl_', sum = TRUE)

mp <- box_plot('X3_methylpentane', center = TRUE)


bp <- grid.arrange(eb, mp, fm, nrow = 3)

ggsave('Top_VOC_marker_distribution_boxplots.tiff', bp, unit = 'mm', dpi = 300, width = 120, height = 155)

ggsave('Top_VOC_marker_distribution_boxplots_centered.tiff', bp, unit = 'mm', dpi = 300, width = 120, height = 155)

#

# add the best candidate VOC markers to the dataset
# center breath VOC levels within each dataset
b1_voc <- b1_corr_out %>% 
  transmute(RAD_ID = RAD_ID,
            CoreVisit = CoreVisit,
            Ethyl_butanoate = Ethyl_butanoate-mean(Ethyl_butanoate),
            X3_methylpentane = X3_methylpentane - mean(X3_methylpentane),
            Furan._2_methyl_ = Furan._2_methyl_ - mean(Furan._2_methyl_)) %>%
  pivot_wider(names_from = CoreVisit, 
              values_from = c(Ethyl_butanoate, X3_methylpentane, Furan._2_methyl_)) 

#

b2_voc <- b2_corr_out %>%
  transmute(RAD_ID = RAD_ID,
            CoreVisit = CoreVisit,
            Ethyl_butanoate = Ethyl_butanoate-mean(Ethyl_butanoate),
            X3_methylpentane = X3_methylpentane - mean(X3_methylpentane),
            Furan._2_methyl_ = Furan._2_methyl_ - mean(Furan._2_methyl_)) %>%
  pivot_wider(names_from = CoreVisit, 
              values_from = c(Ethyl_butanoate, X3_methylpentane, Furan._2_methyl_)) %>%
  filter(RAD_ID %ni% c('RAD076', 'RAD102', 'RAD103'))

b_voc <- rbind(b1_voc, b2_voc)

# summarise VOC levels between visits
b1_voc_sum <- b1_corr_out %>% group_by(RAD_ID) %>%
  summarise(Ethyl_butanoate = mean(Ethyl_butanoate, na.rm = TRUE),
            Furan._2_methyl_ = mean(Furan._2_methyl_, na.rm = TRUE),
            X3_methylpentane = mean(X3_methylpentane, na.rm = TRUE)) %>%
  transmute(RAD_ID = RAD_ID,
            Ethyl_butanoate = Ethyl_butanoate-mean(Ethyl_butanoate),
            X3_methylpentane = X3_methylpentane - mean(X3_methylpentane),
            Furan._2_methyl_ = Furan._2_methyl_ - mean(Furan._2_methyl_))
  
#

b2_voc_sum <- b2_corr_out %>% group_by(RAD_ID) %>%
  summarise(Ethyl_butanoate = mean(Ethyl_butanoate, na.rm = TRUE),
            Furan._2_methyl_ = mean(Furan._2_methyl_, na.rm = TRUE),
            X3_methylpentane = mean(X3_methylpentane, na.rm = TRUE)) %>%
  transmute(RAD_ID = RAD_ID,
            Ethyl_butanoate = Ethyl_butanoate-mean(Ethyl_butanoate),
            X3_methylpentane = X3_methylpentane - mean(X3_methylpentane),
            Furan._2_methyl_ = Furan._2_methyl_ - mean(Furan._2_methyl_))


#

b1_clin <- clin_voc %>% filter(RAD_ID %in% b1_voc$RAD_ID) %>%
  left_join(b1_voc) %>% 
  left_join(b1_corr_out %>% dplyr::select(RAD_ID, Diagnosis)) %>% 
  distinct() %>% mutate(Diagnosis = as.factor(Diagnosis))

b2_clin <- clin_voc %>% filter(RAD_ID %in% b2_voc$RAD_ID) %>%
  left_join(b2_voc) %>% 
  left_join(b2_corr_out %>% dplyr::select(RAD_ID, Diagnosis)) %>% 
  distinct() %>% mutate(Diagnosis = as.factor(Diagnosis))

# 

b1_clin_sum <- clin_voc %>% filter(RAD_ID %in% b1_voc_sum$RAD_ID) %>%
  left_join(b1_voc_sum) %>% 
  left_join(b1_corr_out %>% dplyr::select(RAD_ID, Diagnosis)) %>% 
  distinct() %>% mutate(Diagnosis = as.factor(Diagnosis))

b2_clin_sum <- clin_voc %>% filter(RAD_ID %in% b2_voc_sum$RAD_ID) %>%
  left_join(b2_voc_sum) %>% 
  left_join(b2_corr_out %>% dplyr::select(RAD_ID, Diagnosis)) %>% 
  distinct() %>% mutate(Diagnosis = as.factor(Diagnosis))


# replace missing category with result from CV2
b2_clin$FeNOCV1cat[is.na(b2_clin$FeNOCV1cat) == TRUE] <- 'Neg'

# remove patients with VOC observations missing in CV1
b1_clin1 <- b1_clin %>% filter(is.na(Ethyl_butanoate_CV1) == FALSE)
b2_clin1 <- b2_clin %>% filter(is.na(Ethyl_butanoate_CV1) == FALSE)

#
#
#

# Classification performance of individual variables
# clinical outcomes
pred_clin <- function(data, out) {
  pred_b1 <- data %>% mutate(Prediction = ifelse({{out}} == 'Pos', 'Asthma', 'Not_Asthma'))
  conf_mat <- table(factor(pred_b1$Diagnosis, levels = levels(as.factor(pred_b1$Diagnosis))), 
                  pred_b1$Prediction)
  perf <- c(round(er_fun(conf_mat),2), round(ber_fun(conf_mat),2), 
            round(sens_fun(conf_mat),2), round(spec_fun(conf_mat),1))
  print(perf)
  }

pred_clin(b1_clin1, FeNOCV1cat)
pred_clin(b2_clin1, FeNOCV1cat)

pred_clin(b1_clin1, BDRcat)
pred_clin(b2_clin1, BDRcat)

pred_clin(b1_clin1, BCTcat)
pred_clin(b2_clin1, BCTcat)


# ROC on clinical data subset
roc_fun <- function(voc, data){
  resp <- as.character(data$Diagnosis)
  resp[resp == 'Asthma'] <- 1
  resp[resp == 'Not Asthma'] <- 0
  pred <- data %>% dplyr::select(voc)
  colnames(pred) <- 'VOC'
  data1 <- cbind(pred, resp)
  
  roc_voc <- roc(response = data1$resp, predictor = data1$VOC)
  
  out_roc_voc <- list(as.data.frame(pROC::coords(roc_voc)),
                      auc_voc = auc(roc_voc)[[1]],
                      pROC::coords(roc_voc, 'best'))

  
}

feno_roc1 <- roc_fun('FeNOCV1', b1_clin1)
feno_roc2 <- roc_fun('FeNOCV1', b2_clin1)

eb_roc1 <- roc_fun(voc = 'Ethyl_butanoate_CV1', data = b1_clin1)
eb_roc2 <- roc_fun(voc = 'Ethyl_butanoate_CV1', data = b2_clin1)

#

mp_roc1 <- roc_fun(voc = 'X3_methylpentane_CV1', data = b1_clin1)
mp_roc2 <- roc_fun(voc = 'X3_methylpentane_CV1', data = b2_clin1)

#

fm_roc1 <- roc_fun(voc = 'Furan._2_methyl__CV1', data = b1_clin1)
fm_roc2 <- roc_fun(voc = 'Furan._2_methyl__CV1', data = b2_clin1)


#
#
#

# Classification performance of NICE guidelines
b1_clin_nice <- b1_clin1 %>% 
  mutate(Prediction = ifelse(FeNOCV1cat == 'Pos', 'Asthma', 
                             ifelse(BDRcat == 'Pos', 'Asthma',
                                    ifelse(BCTcat == 'Pos', 'Asthma', 'Not_Asthma'))))

b1_nice_conf <- table(factor(b1_clin_nice$Diagnosis, levels = levels(as.factor(b1_clin_nice$Diagnosis))),
                             b1_clin_nice$Prediction)

er_fun(b1_nice_conf)
ber_fun(b1_nice_conf)
sens_fun(b1_nice_conf)
spec_fun(b1_nice_conf)

b1_clin_nice %>% group_by(FeNOCV1cat, Diagnosis) %>% summarise(n = n())
b1_clin_nice %>% filter(FeNOCV1cat == 'Neg') %>% group_by(BDRcat, Diagnosis) %>% summarise(n = n())
b1_clin_nice %>% filter(FeNOCV1cat == 'Neg') %>% filter(BDRcat == 'Neg') %>%
  group_by(BCTcat, Diagnosis) %>% summarise(n = n())


#
b2_clin_nice <- b2_clin1 %>% 
  mutate(Prediction = ifelse(FeNOCV1cat == 'Pos', 'Asthma', 
                             ifelse(BDRcat == 'Pos', 'Asthma',
                                    ifelse(BCTcat == 'Pos', 'Asthma', 'Not_Asthma'))))

b2_nice_conf <- table(factor(b2_clin_nice$Diagnosis, levels = levels(as.factor(b2_clin_nice$Diagnosis))),
                      b2_clin_nice$Prediction)

er_fun(b2_nice_conf)
ber_fun(b2_nice_conf)
sens_fun(b2_nice_conf)
spec_fun(b2_nice_conf)

b2_clin_nice %>% group_by(FeNOCV1cat, Diagnosis) %>% summarise(n = n())
b2_clin_nice %>% filter(FeNOCV1cat == 'Neg') %>% group_by(BDRcat, Diagnosis) %>% summarise(n = n())
b2_clin_nice %>% filter(FeNOCV1cat == 'Neg') %>% filter(BDRcat == 'Neg') %>%
  group_by(BCTcat, Diagnosis) %>% summarise(n = n())

# draw dendrogram
library(ggraph)
library(igraph)
library(grafify)
library(cowplot)

pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

tree_df <- data.frame(main = rep('Entry', 4),
                      FeNO = c('Pos_FeNO', rep('Neg_FeNO', 3)),
                      BDR = c(NA, 'Pos_BDR', rep('Neg_BDR', 2)),
                      BCT = c(rep(NA, 2), 'Pos_BCT', 'Neg_BCT'),
                      Diag = c(rep('A', 3), 'NA'))


edges_level1_1 <- tree_df %>% select(main, FeNO) %>% unique %>% rename(from = main, to = FeNO)
edges_level1_2 <- tree_df %>% select(FeNO, BDR) %>% unique %>% rename(from = FeNO, to = BDR)
edges_level2_3 <- tree_df %>% select(BDR, BCT) %>% unique %>% rename(from = BDR, to = BCT)
edges_level3_4 <- tree_df %>% select(BCT, Diag) %>% unique %>% rename(from = BCT, to = Diag)

edge_list = rbind(edges_level1_1, edges_level1_2, edges_level2_3)

mygraph <- graph_from_data_frame(edge_list)


mygraph1 <- permute(mygraph, permutation = c(1:5, 7, 6))

nice_tree <- ggraph(mygraph, layout = 'dendrogram') + 
  geom_edge_elbow() +
  theme_void() +
  geom_node_text(label = c('FeNO > 50ppb', NA, 'BDR > 12% FEV1', NA, NA, 'BCT PD20 < 0.2mg', NA, NA, NA),
                 nudge_y = 0.2, nudge_x = -0.2, colour = my_pal[[7]], fontface = 'bold', size = 3) 

nice_tree

nice_tree1 <- plot_grid(nice_tree) +
  theme(plot.margin = unit(c(1, 0.3, 0.5, 0.1), 'cm')) +
  draw_label('Confirm diagnosis', x = 0.35, y = 0, size = 8, vjust = 0.8) +
  draw_label('Consider \n alternative', x = 0.95, y = 0, size = 8, vjust = 0.8) +
  draw_label('Decision tree for asthma diagnosis \n following NICE guideline',
             size = 8, fontface = 'bold', x = 0.5, y = 1, vjust = -0.5) +
  draw_line(x = c(0.05, 0.68), y = 0.02, linewidth = 0.3) +
  draw_label('Yes', x = 0.05, y = 0.8, size = 8) +
  draw_label('No', x = 0.65, y = 0.8, size = 8) +
  draw_label('Yes', x = 0.33, y = 0.5, size = 8) +
  draw_label('No', x = 0.86, y = 0.5, size = 8) +
  draw_label('Yes', x = 0.63, y = 0.25, size = 8) +
  draw_label('No', x = 0.99, y = 0.25, size = 8) 

plot(nice_tree1)

#

# NICE tree with FeNO replaced by Ethyl butanoate at best threshold specified by ROC in training dataset
b1_clin_nice1 <- b1_clin1 %>%
  mutate(Ethyl_butanoate_CV1cat = ifelse(Ethyl_butanoate_CV1 > eb_roc2[[3]][[1]], 'Neg', 'Pos')) %>%
  mutate(Prediction = ifelse(Ethyl_butanoate_CV1cat == 'Pos', 'Asthma', 
                             ifelse(BDRcat == 'Pos', 'Asthma',
                                    ifelse(BCTcat == 'Pos', 'Asthma', 'Not_Asthma'))))

b1_nice_conf <- table(factor(b1_clin_nice1$Diagnosis, levels = levels(as.factor(b1_clin_nice1$Diagnosis))),
                      b1_clin_nice1$Prediction)

er_fun(b1_nice_conf)
ber_fun(b1_nice_conf)
sens_fun(b1_nice_conf)
spec_fun(b1_nice_conf)



#
#
#

# Multivariate classification model
# random forest with clinical data only
set.seed(123)
base_forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat, 
                      data = b2_clin1,
                      mtry = 3,
                      num.trees = 1000)

base_forest
treeInfo(base_forest, tree = 12)


# random forest with VOC CV1 data only
set.seed(123)
base_forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                        Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                      data = b2_clin1,
                      mtry = 3,
                      num.trees = 1000)

base_forest
treeInfo(base_forest, tree = 12)

# with summarised data
set.seed(123)
base_forest <- ranger(Diagnosis ~ FeNOcat + BDRcat + BCTcat + 
                        Ethyl_butanoate + X3_methylpentane + Furan._2_methyl_, 
                      data = b2_clin_sum,
                      mtry = 3,
                      num.trees = 1000)


base_forest

# prediction error on training dataset
train_conf <- base_forest$confusion.matrix

er_fun(train_conf)
ber_fun(train_conf)
sens_fun(train_conf)
spec_fun(train_conf)

#

# prediction error in cross-validation
set.seed(123)
cv <- replicate(n = 50, expr = {
  cv_folds <- fold_fun(b2_clin1, nfold = 5)
  
  BER <- vector(mode = "numeric")
  ER <- vector(mode = 'numeric')
  SENS <- vector(mode = 'numeric')
  SPEC <- vector(mode = 'numeric')
  
  for (m in as.character(0:(length(cv_folds) - 1))) {
  cv_model_pred <- train_test(b2_clin1, cv_folds, nvar = 3, m, n = 1000)
  BER[m] <- if (ncol(cv_model_pred) == 2 & nrow(cv_model_pred) == 2) {
    ber_fun(cv_model_pred)
  } else {
      NA}
  }
  ER[m] <- if (ncol(cv_model_pred) == 2) {
    er_fun(cv_model_pred)
  } else {
    NA}
  SENS[m] <- if (ncol(cv_model_pred) == 2 & nrow(cv_model_pred) == 2) {
    sens_fun(cv_model_pred)
  } else {
    NA}
  SPEC[m] <- if (ncol(cv_model_pred) == 2 & nrow(cv_model_pred) == 2) {
    spec_fun(cv_model_pred)
  } else {
    NA}
  
  cv_ber <- mean(BER, na.rm = TRUE)
  cv_er <- mean(ER, na.rm = TRUE)
  cv_sens <- mean(SENS, na.rm = TRUE)
  cv_spec <- mean(SPEC, na.rm = TRUE)
  out <- c(cv_ber, cv_er, cv_sens, cv_spec)
  
  })

out_sum <- data.frame(mean = rowMeans(cv), sd = rowSds(cv))
rownames(out_sum) <- c('BER', 'ER', 'Sens', 'Spec')


# prediction error on test dataset
pred_forest <- predict(base_forest, data = b1_clin1)

test_conf <- table(factor(b1_clin1$Diagnosis, levels = levels(as.factor(b1_clin1$Diagnosis))), 
                  pred_forest$predictions)

er_fun(test_conf)
ber_fun(test_conf)
sens_fun(test_conf)
spec_fun(test_conf)

#

# prediction error using CV2 (in train and test data)
# remove patients with VOC observations missing in CV1
b1_clin2 <- b1_clin %>% filter(is.na(Ethyl_butanoate_CV2) == FALSE)
b2_clin2 <- b2_clin %>% filter(is.na(Ethyl_butanoate_CV2) == FALSE)

b1_clin2 <- b1_clin2 %>% 
  mutate(Ethyl_butanoate_CV1sc = (Ethyl_butanoate_CV2-mean(Ethyl_butanoate_CV2))/sd(Ethyl_butanoate_CV2),
         X3_methylpentane_CV1sc = (X3_methylpentane_CV2 - mean(X3_methylpentane_CV2))/sd(X3_methylpentane_CV2),
         Furan._2_methyl__CV1sc = (Furan._2_methyl__CV2 - mean(Furan._2_methyl__CV2))/sd(Furan._2_methyl__CV2))

b2_clin2 <- b2_clin2 %>% 
  mutate(Ethyl_butanoate_CV1sc = (Ethyl_butanoate_CV2-mean(Ethyl_butanoate_CV2))/sd(Ethyl_butanoate_CV2),
         X3_methylpentane_CV1sc = (X3_methylpentane_CV2 - mean(X3_methylpentane_CV2))/sd(X3_methylpentane_CV2),
         Furan._2_methyl__CV1sc = (Furan._2_methyl__CV2 - mean(Furan._2_methyl__CV2))/sd(Furan._2_methyl__CV2))

#
pred_forest2 <- predict(base_forest, data = b2_clin2)

test_conf2 <- table(factor(b2_clin2$Diagnosis, levels = levels(as.factor(b2_clin2$Diagnosis))), 
                   pred_forest2$predictions)

er_fun(test_conf2)
ber_fun(test_conf2)
sens_fun(test_conf2)
spec_fun(test_conf2)

#
pred_forest3 <- predict(base_forest, data = b1_clin2)

test_conf3 <- table(factor(b1_clin2$Diagnosis, levels = levels(as.factor(b1_clin2$Diagnosis))), 
                    pred_forest3$predictions)

er_fun(test_conf3)
ber_fun(test_conf3)
sens_fun(test_conf3)
spec_fun(test_conf3)


#
#
#


#
#
#

# MODEL OPTIMISATION
# Optimise number of trees, number of variables available at each split and tree size (number of terminal nodes)
# prepare argument combinations
args <- expand.grid(nvar = c(2, 3, 6),
                    n = c(seq(50, 2000, 50)))

# training set
fun <- function(nvar, n) {
  forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
              data = b2_clin1,
              mtry = nvar,
              num.trees = n)
  train_conf <- forest$confusion.matrix
  ber <- ber_fun(train_conf)
}

train_opt <- replicate(n = 50, expr = {mapply(FUN = fun,
                    nvar = args$nvar,
                    n = args$n)
  })

train_op <- args %>%
  mutate(ber_mean = rowMeans(train_opt),
         ber_sd = rowSds(train_opt))

# visualise results
train_op %>% 
  ggplot(aes(x = n, y = ber_mean, group = nvar, colour = as.factor(nvar))) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = ber_mean - ber_sd, ymax = ber_mean + ber_sd)) +
  theme_bw()

ggsave('RF_opt_train_treeNum.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

#

# Cross-validation
# test model parameter grid with CV
args1 <- expand.grid(n = c(seq(50, 2000, 50)),
                    nvar = c(2,3))

cv_opt <- replicate(n = 50, expr = {
  cv_folds <- fold_fun(b2_clin1, nfold = 8)
  BER <- vector(mode = "numeric")
  
  cv_opt <- mapply(n = args1$n,
                   nvar = args1$nvar,
                   FUN = function(n, nvar) {
           for (m in as.character(0:(length(cv_folds) - 1))) {
             cv_model_pred <- train_test(b2_clin1, cv_folds, nvar = nvar, m, n = n)
             BER[m] <- if (ncol(cv_model_pred) == 2 & nrow(cv_model_pred) == 2) {
               ber_fun(cv_model_pred)
             } else {
                 NA}
             print(BER)
             }
           ber <- mean(BER, na.rm = TRUE)
           })
  })
  


#

cv_op <- args1 %>% mutate(ber_mean = rowMeans(cv_opt),
                          ber_sd = rowSds(cv_opt))

write.csv(cv_op, 'RF_CV_parameter_grid.csv')

cv_op %>% 
  #filter(nvar == 3) %>%
  ggplot(aes(x = n, y = ber_mean, colour = as.factor(nvar))) + 
  geom_point() + geom_errorbar(aes(ymin = ber_mean - ber_sd, ymax = ber_mean + ber_sd)) +
  geom_line() +
  theme_bw()

ggsave('RF_opt_CV_treeNum1.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)


#

# Test dataset
args2 <- expand.grid(nvar = c(2, 3),
                     n = c(seq(50, 2000, 50)))

fun2 <- function(nvar, n, train = b2_clin1, test = b1_clin1) {
  forest1 <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                     Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                   data = train,
                   mtry = nvar,
                   num.trees = n)
  
  pred_forest1 <- predict(forest1, data = test)
  
  test_conf1 <- table(factor(test$Diagnosis, levels = levels(as.factor(test$Diagnosis))), 
                     pred_forest1$predictions)
  
  ber <- ber_fun(test_conf1)
}


test_opt <- replicate(n = 50, expr = {
  mapply(FUN = fun2,
         nvar = args2$nvar,
         n = args2$n)
})

test_op <- args2 %>% 
  mutate(ber_mean = rowMeans(test_opt),
         ber_sd = rowSds(test_opt))

#

test_op %>% 
  #filter(nvar == 3) %>%
  ggplot(aes(x = n, y = ber_mean, group = nvar, colour = as.factor(nvar))) + 
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin = ber_mean - ber_sd, ymax = ber_mean + ber_sd)) +
  theme_bw()

ggsave('RF_opt_test_treeNum1.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

#
#
#

################################

# Variable importance
# through permutation
set.seed(123)
base_forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                        Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                      data = b2_clin1,
                      mtry = 3,
                      num.trees = 1000,
                      importance = 'permutation')

imp <- importance(base_forest) %>% as.data.frame()

importance_pvalues(base_forest, method = 'altmann',
                   formula = Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                     Ethyl_butanoate + X3_methylpentane + Furan._2_methyl_, 
                   data = b2_clin_sum)


# variable importance stability
# using package RFopt
var_stab <- opt_importance(y = b2_clin1$Diagnosis, X = b2_clin1 %>%
                             dplyr::select(FeNOCV1cat, BDRcat, BCTcat, 
                                           Ethyl_butanoate_CV1, X3_methylpentane_CV1,
                                           Furan._2_methyl__CV1))

var_stab1 <- opt_importance(y = b2_clin1$Diagnosis, X = b2_clin1 %>%
                              dplyr::select(FeNOCV1cat, BDRcat, BCTcat, 
                                            Ethyl_butanoate_CV1, X3_methylpentane_CV1,
                                            Furan._2_methyl__CV1),
                            recommendation = 'selection',
                            alpha = 3)

# manually
fun_stab <- function(n) {
  forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                     Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                   data = b2_clin1,
                   mtry = 3,
                   num.trees = n,
                   importance = 'permutation')
  importance(forest) %>% as.data.frame() %>%
    mutate(pred = rownames(.))
}


set.seed(1410)
train_stab <- replicate(n = 50, expr = {
  lapply(c(seq(100, 2000, 100)), FUN = fun_stab)
  }, simplify = FALSE)


stab_res <- data.frame()
for (m in 1:length(train_stab)){
  rep <- as.data.frame(train_stab[[m]])
  colnames(rep) <- c(seq(100, 2000, 100))
  rep <- rep %>% mutate(pred = rownames(.),
                        rep = m) %>%
    pivot_longer(cols = !c(pred, rep), names_to = 'n', values_to = 'imp')
  stab_res <- stab_res %>% rbind(rep)
}


stab_res_sum <- stab_res %>% group_by(pred, n) %>% 
  summarise(mean_imp = mean(imp),
            sd_imp = sd(imp),
            cv_imp = sd_imp/mean_imp) 

stab_res_sum %>%
  ggplot(aes(x = as.numeric(n), y = mean_imp, colour = pred)) +
  geom_point() + geom_errorbar(aes(ymin = mean_imp - sd_imp, ymax = mean_imp + sd_imp)) +
  geom_line(aes(group = pred)) +
  theme_bw(base_size = 10) +
  xlab('Tree number') +
  ylab('Variable importance')

ggsave('RF_opt_VIM_treeNum_var.tiff', dpi = 300, unit = 'mm', width = 140, height = 90)

#

stab_res_sum %>% group_by(n) %>%
  summarise(mean_sd = mean(sd_imp)) %>%
  ggplot(aes(x = as.numeric(n), y = mean_sd)) +
  geom_line() +
  theme_bw() +
  ylab('Standard deviation') +
  xlab('Tree number')

ggsave('RF_opt_VIM_treeNum_sd.tiff', dpi = 300, unit = 'mm', width = 80, height = 60)

#

stab_res_w <- stab_res %>% pivot_wider(names_from = rep, values_from = imp) %>%
  mutate(n = as.numeric(n))

stab_icc <- sapply(c(seq(100, 2000, 100)), function(ntree){
input <- stab_res_w %>% filter(n == ntree)
icc_stab <- icc(input[,3:50], model = 'twoway', type = 'agreement') #columns and rows random
output <- data.frame(ntree = ntree, 
                     ICC = icc_stab$value, 
                     p.value = icc_stab$p.value,
                     CI_lwr = icc_stab$lbound, CI_upr = icc_stab$ubound)
})

stab_icc_t <- as.data.frame(t(stab_icc))

stab_icc_t %>% ggplot(aes(x = as.numeric(ntree), y = as.numeric(ICC))) + 
  geom_point() + geom_line() +
  theme_bw(base_size = 8) +
  xlab('Tree number') +
  ylab('ICC') +
  ggtitle('Intraclass correlation coefficient of variable importance (50 repetitions)')

ggsave('RF_opt_VIM_treeNum.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

#
#
#
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

vip_plot <- stab_res_sum %>% filter(n == 1000) %>%
  arrange(desc(mean_imp)) %>%
  ggplot(aes(x = mean_imp, y = fct_inorder(as.factor(pred)), fill = pred)) +
  geom_col(colour = 'black', lwd = 0.3) +
  geom_errorbar(aes(xmin = mean_imp - sd_imp, xmax = mean_imp + sd_imp),
                width = 0.4) +
  theme_bw(base_size = 8) +
  theme(legend.position = 'none') +
  ylab('') +
  xlab('Variable importance') +
  ggtitle('Variable importance in random forest model') +
  theme(plot.title = element_text(hjust = 0.9, face = 'bold', size = 8),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), 'cm')) +
  scale_fill_manual(values = c('BCTcat' = my_pal[[6]],
                               'BDRcat' = my_pal[[5]],
                               'Ethyl_butanoate_CV1' = my_pal[[7]],
                               'Furan._2_methyl__CV1' = my_pal[[3]],
                               'X3_methylpentane_CV1' = my_pal[[8]],
                               'FeNOCV1cat' = my_pal[[9]])) +
  scale_y_discrete(labels = c('BCT',
                               'BDR',
                               'Ethyl butanoate',
                               'Furan-2-methyl',
                               '3_methylpentane',
                               'FeNO')) +
  annotate(geom = 'text', label = '*', x = 0.2, y = 1)  +
  annotate(geom = 'text', label = '*', x = 0.1, y = 2)

vip_plot

ggsave('RF_VIM_1000tr_50rep_sum.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

##################################

# Figure 4
fig4 <- plot_grid(nice_tree1, vip_plot, nrow = 1, labels = c('a)', 'b)'), label_size = 12)
plot(fig4)

ggsave('fig4.tiff', unit = 'mm', dpi = 300, width = 157, height = 60)

pdf('Figure4.pdf', width = 6.1, height = 2.36)
plot(fig4)
dev.off()

# using repeated 2-fold cross validation
set.seed(123)
cv_vim_rep <- replicate(50, expr = {
  cv_vim <- holdoutRF(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                        Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                      data = b2_clin1,
                      mtry = 3,
                      num.trees = 1000)
  cv_vim_out <- cv_vim[['variable.importance']]
  
  })

cv_vim_sum <- data_frame(pred = rownames(cv_vim_rep),
                         mean_vim = rowMeans(cv_vim_rep),
                         sds_vim = rowSds(cv_vim_rep))


cv_vim_sum %>%
  arrange(desc(mean_vim)) %>%
  ggplot(aes(x = mean_vim, y = fct_inorder(as.factor(pred)), fill = pred)) +
  geom_col(colour = 'black') +
  geom_errorbar(aes(xmin = mean_vim - sds_vim, xmax = mean_vim + sds_vim),
                width = 0.4) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none') +
  ylab('') +
  xlab('Variable importance') +
  ggtitle('Variable importance in random forest model (50 repetitions)') +
  theme(plot.title = element_text(hjust = 1)) 


ggsave('RF_VIM_1000tr_50rep_holdout.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)

#
#
#

# Partial Dependance plots
library(pdp)

set.seed(123)
base_forest <- ranger(Diagnosis ~ FeNOCV1cat + BDRcat + BCTcat + 
                        Ethyl_butanoate_CV1 + X3_methylpentane_CV1 + Furan._2_methyl__CV1, 
                      data = b2_clin1,
                      mtry = 3,
                      num.trees = 1000,
                      probability = TRUE)

base_pred <- base_forest$predictions %>% as.data.frame() %>%
  mutate(prediction = ifelse(Asthma > `Not Asthma`, 'Asthma', 'Not Asthma')) %>%
  cbind(b2_clin1 %>% dplyr::select(Diagnosis))

train_conf <- table(factor(base_pred$prediction, levels(as.factor(base_pred$Diagnosis))),
                           base_pred$Diagnosis)

er_fun(train_conf)
ber_fun(train_conf)
sens_fun(train_conf)
spec_fun(train_conf)


par_eb <- partial(base_forest, pred.var = c("Ethyl_butanoate_CV1"), chull = TRUE)
autoplot(par_eb, contour = TRUE)

par_mp <- partial(base_forest, pred.var = c("X3_methylpentane_CV1"), chull = TRUE,
                  grid.resolution = 20)
autoplot(par_mp, contour = TRUE)

par_fm <- partial(base_forest, pred.var = c("Furan._2_methyl__CV1"), chull = TRUE,
                  grid.resolution = 20)
autoplot(par_fm, contour = TRUE)

#
dev.new()
b2_clin1 %>% 
  pivot_longer(cols = c(Ethyl_butanoate_CV1, Furan._2_methyl__CV1, X3_methylpentane_CV1),
               names_to = 'VOC', values_to = 'level') %>%
  ggplot(aes(x = VOC, y = level)) +
  geom_boxplot(aes(fill = Diagnosis)) +
  geom_beeswarm(aes(group = Diagnosis), dodge.width = 0.7) +
  theme_bw() +
  ggtitle('Dataset 2')

# pull the best thresholds from ROC
eb_th <- eb_roc2[[3]][,1]
mp_th <- mp_roc2[[3]][,1]
fm_th <- fm_roc2[[3]][,1]

b2_clin2 <- b2_clin1 %>% mutate(Ethyl_butanoate_cat = ifelse(Ethyl_butanoate_CV1 > eb_th, 'Neg', 'Pos'),
                                Furan._2_methyl_cat = ifelse(Furan._2_methyl__CV1 > fm_th, 'Neg', 'Pos'),
                                X3_methylpentane_cat = ifelse(X3_methylpentane_CV1 > mp_th, ' Neg', 'Pos'))

b1_clin2 <- b1_clin1 %>% mutate(Ethyl_butanoate_cat = ifelse(Ethyl_butanoate_CV1 > eb_th, 'Neg', 'Pos'),
                                Furan._2_methyl_cat = ifelse(Furan._2_methyl__CV1 > fm_th, 'Neg', 'Pos'),
                                X3_methylpentane_cat = ifelse(X3_methylpentane_CV1 > mp_th, ' Neg', 'Pos'))

eb_roc2[[3]]
mp_roc2[[3]]
fm_roc2[[3]]

View(b2_clin2 %>% group_by(Ethyl_butanoate_cat, Furan._2_methyl_cat, FeNOcat,
                           Diagnosis) %>% summarise(n = n()))

b2_clin %>% group_by(FeNOcat, Diagnosis) %>% summarise(n = n())
b2_clin2 %>% #filter(FeNOcat == 'Neg') %>% 
  group_by(Ethyl_butanoate_cat, Diagnosis) %>% summarise(n = n())
b2_clin2 %>% filter(Ethyl_butanoate_cat == 'Neg') %>% 
  group_by(Furan._2_methyl_cat, Diagnosis) %>% summarise(n = n())
b2_clin2 %>% #filter(FeNOcat == 'Neg') %>% 
  filter(Ethyl_butanoate_cat == 'Neg') %>% filter(Furan._2_methyl_cat == 'Neg') %>%
  group_by(BDRcat, Diagnosis) %>% summarise(n = n())
b2_clin2 %>% #filter(FeNOcat == 'Neg') %>% 
  filter(Ethyl_butanoate_cat == 'Neg') %>% 
  filter(Furan._2_methyl_cat == 'Neg') %>%
  filter(BDRcat == 'Neg') %>% 
  group_by(BCTcat, Diagnosis) %>% summarise(n = n())

