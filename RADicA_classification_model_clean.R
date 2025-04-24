## Predicting asthma based on breath VOC samples

# author: Aggie Turlo
# project: RADicA
# date: 10/02/2025

#####################


library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(vsn)
library(mixOmics)
library(forcats)
library(irr)
library(ggcorrplot)
library(pROC)
library(scales)
library(mitoticFigureCounts)
library(shades)
library(cowplot)
library(purrr)

# load data
# load summarised BG corrected datasets 
b1_corr_w <- read.csv('RADicA_B1_BG_adjusted.csv')[,-1]
b2_corr_w <- read.csv('RADicA_B2_BG_adjusted.csv')[,-1]

# load summarised BG corrected datasets w/o multivariate outliers
b1_corr_out <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b2_corr_out <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]

# load metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')[,-1]

# load VOC origin annotation based on logFC analysis
endo_exo <- endo_exo <- read.csv('Endo_Exo_filters.csv')[,-1] 

#

# data formatting
# set columns in the same order
b1_corr_out <- b1_corr_out[names(b2_corr_out)]

# remove patient IDs repeated between datasets
b2_corr_out <- b2_corr_out %>% filter(RAD_ID %ni% intersect(b1_corr_out$RAD_ID, b2_corr_out$RAD_ID))

#
#
#

# custom functions
# balanced error rate
ber_fun <-  function(conf_mat) {
  ber <- 0.5*((conf_mat[1,2]/(conf_mat[1,2]+conf_mat[2,2])) + 
                (conf_mat[2,1]/(conf_mat[2,1] + conf_mat[1,1])))
  print(ber)
}


# error rate
er_fun <- function(conf_mat){
  er <- (sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
  print(er)
}

# test sensitivity
sens_fun <- function(conf_mat){
  sens <- conf_mat[1,1]/(conf_mat[1,1] + conf_mat[1,2])
  print(sens)
}

# test specificity
spec_fun <- function(conf_mat){
  spec <- conf_mat[2,2]/(conf_mat[2,2] + conf_mat[2,1])
  print(spec)
}

#
'%ni%' <- Negate('%in%')

#
#
#


# examine reproducibility of VOC levels across time points
# calculate intraclass correlation coefficient between Core Visits
# specify dataset for analysis
data <- b2_corr_out
vocs <- colnames(data)[-c(1:4)]

#

icc_b2 <- bind_rows(lapply(vocs, function(voc){
  input <-  data %>%
    pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea') %>%
    filter(comp == voc) %>%
    dplyr::select(!Sample) %>%
    pivot_wider(names_from = CoreVisit, values_from = peakArea) %>%
    as.data.frame()
  
  rownames(input) <- input$RAD_ID
  input <- input[,c(4:5)]
  
  icc_res <- icc(input, model = 'twoway', type = 'consistency', unit = 'single') #columns and rows random
  output <- data.frame(comp = voc, 
                       ICC = icc_res$value, 
                       p.value = icc_res$p.value,
                       CI_lwr = icc_res$lbound, CI_upr = icc_res$ubound)
  
})) 

# save to adequately names objects
icc_b1 <- icc_b1 %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
hist(icc_b1$ICC)
quantile(icc_b1$ICC)

icc_b2 <- icc_b2 %>% mutate(adj.p.value = p.adjust(p.value, method = 'BH'))
hist(icc_b2$ICC)
quantile(icc_b2$ICC)

# save results to csv file
write.csv(icc_b1, 'ICC_CoreVisits_B1.csv')
write.csv(icc_b2 , 'ICC_CoreVisits_B2.csv')

#
#
#

###############################

# clustered AUROC for individual VOCs
# center and scale vocs 
b1_corr_out_sc <- cbind(b1_corr_out[,1:4], scale(b1_corr_out[, - c(1:4)], center = TRUE))
b2_corr_out_sc <- cbind(b2_corr_out[,1:4], scale(b2_corr_out[, - c(1:4)], center = TRUE))

b1_corr_out_sc$Diagnosis <- as.factor(b1_corr_out_sc$Diagnosis)
b2_corr_out_sc$Diagnosis <- as.factor(b2_corr_out_sc$Diagnosis)

roc_fun <- function(voc, data){
  resp <- data$Diagnosis
  pred <- data %>% dplyr::select(voc)
  colnames(pred) <- 'VOC'
  clusID <- data$RAD_ID
  data1 <- cbind(pred, resp, clusID)
  
  out_roc_voc <- as.data.frame(doAUCcluster(predictor1 = data1$VOC,
                                            response = data1$resp,
                                            clusterID = data1$clusID)) %>%
    mutate(comp = voc)
  
}

# calculate AUROC for each dataset
voc_roc_b2 <- bind_rows(lapply(colnames(b2_corr_out_sc[, -c(1:4)]), roc_fun, data = b2_corr_out_sc))
voc_roc_b1 <- bind_rows(lapply(colnames(b1_corr_out_sc[, -c(1:4)]), roc_fun, data = b1_corr_out_sc))

b2_top <- voc_roc_b2 %>% filter(auc > 0.60)
b1_top <- voc_roc_b1 %>% filter(auc > 0.60)

intersect(b2_top$comp, b1_top$comp)

#

voc_roc_b1_cv1 <- bind_rows(lapply(colnames(b1_corr_out_cv1[, -c(1:4)]), roc_fun, data = b1_corr_out_cv1))
voc_roc_b1_cv2 <- bind_rows(lapply(colnames(b1_corr_out_cv2[, -c(1:4)]), roc_fun, data = b1_corr_out_cv2))

# concatenate results from two datasets at each time point
comp_auc <- voc_roc_b2 %>% dplyr::select(comp, auc) %>% distinct() %>%
  rename(auc_b2 = auc) %>%
  left_join(voc_roc_b1 %>% dplyr::select(comp, auc) %>% distinct()) %>%
  rename(auc_b1 = auc) %>%
  mutate(lab = ifelse(auc_b2 > 0.6 & auc_b1 > 0.6, 'yes', 'no'),
         name = ifelse(lab == 'yes', comp, NA)) %>%
  left_join(endo_exo)

#

# plot relationship between AUROC from two datasets
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

plot2b <- comp_auc %>%
  ggplot(aes(x = auc_b2, y = auc_b1, fill = lab)) +
  geom_point(size = 0.9, shape = 21) +
  geom_vline(xintercept = 0.6, linetype = 'dashed', colour = 'darkgrey', lwd = 0.3) +
  geom_hline(yintercept = 0.6, linetype = 'dashed', colour = 'darkgrey', lwd = 0.3) +
  theme_bw(base_size = 8) +
  scale_fill_manual(values = c('yes' = my_pal[[7]], 'no' = 'white')) +
  geom_text_repel(aes(label = name), colour = 'purple', size = 2,
                  nudge_x = 0.1, direction = 'y', nudge_y = 0.015, fontface = 'bold') +
  theme(plot.title = element_text(hjust = 0.5, size = 8, face = 'bold'),
        legend.position = 'none',
        panel.grid = element_blank()) +
  xlab('Training dataset') + ylab('Validation dataset') +
  ggtitle('Area under the Receiver Operating Characteristic \n curve for breath VOCs')

ggsave('AUROC_B1_B2_indi_VOCs.tiff', unit = 'mm', width = 90, height = 80)

##################################

# Figure 2b
plot2b <- plot2b +
  theme(plot.margin = unit(c(0.2, 2, 0, 2), 'cm'))

plot2b

# figure 2 top panel
plot2ab <- plot_grid(plot2a, plot2b, nrow = 2, rel_heights = c(0.6, 0.4), labels = c('a)', 'b)'), label_size = 12)
dev.new()
plot(plot2ab)



# Figure 2c
cols <- c('Diagnosis', 'RAD_ID', 'CoreVisit', 'Ethyl_butanoate', 'Furan._2_methyl_', 'X3_methylpentane')

comm_out <- b1_corr_out[, colnames(b1_corr_out) %in% cols] %>%
  mutate(dataset = 'Validation dataset') %>%
  rbind(b2_corr_out[, colnames(b2_corr_out) %in% cols] %>%
          mutate(dataset = 'Training dataset'))
  
plot2c_fun <- function(voc) {
    comm_out %>% dplyr::select(voc, Diagnosis, CoreVisit, dataset) %>%
    rename(VOC = voc) %>%
    ggplot(aes(x = CoreVisit, y = VOC, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE, lwd = 0.3) + 
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, size = 0.5) +
    theme_bw(base_size = 8) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 8, hjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = unit(c(0.5, 0.15, 0.1, 0.1), 'cm')) +
    facet_wrap(~dataset, nrow = 1) +
    scale_fill_manual(values = c('Asthma' = my_pal[[2]],
                               'Not Asthma' = my_pal[[1]])) +
    xlab('Collection visit')
  }
  

plot2c1 <- plot2c_fun('Ethyl_butanoate') + ggtitle('Ethyl butanoate')
plot2c2 <- plot2c_fun('Furan._2_methyl_') + ggtitle('Furan-2-methyl') +
  theme(plot.margin = unit(c(0.2, 0.15, 0.1, 0.1), 'cm'))
plot2c3 <- plot2c_fun('X3_methylpentane') + ggtitle('3-methylpentane') +
  theme(axis.title.x = element_text(),
        legend.position = 'bottom',
        plot.margin = unit(c(0.2, 0.15, 0.1, 0.1), 'cm'))


plot2c <- plot_grid(plot2c1, plot2c2, plot2c3, ncol = 1, rel_heights = c(0.3, 0.3, 0.4)) +
  theme(plot.margin = unit(c(0.45, 0, 0.1, 0.5), 'cm')) +
  draw_label(label = 'Log peak area', x = 0, y = 0.5, angle = 90, size = 8, vjust = 0) +
  draw_label(label = 'VOC abundance in breath', 
             x = 0.5, y = 1, size = 8, vjust = -0.5, fontface = 'bold')  +
  draw_label(label = 'Training dataset', x = 0.25, y = 0.97, vjust = -0.5, size = 8) +
  draw_label(label = 'Validation dataset', x = 0.75, y = 0.97, vjust = -0.5, size = 8)

dev.new()
plot(plot2c)

# Figure 2
fig2 <- plot_grid(plot2ab, plot2c, ncol = 2, rel_widths = c(0.6, 0.4), labels = c(NA, 'c)'), label_size = 12)
plot(fig2)

ggsave('fig2.tiff', fig2, units = 'mm', dpi = 300, width = 157, height = 130)

pdf('Figure2.pdf', width = 6.18, height = 5.15)
plot(fig2)
dev.off()

#
#
#

################################

# MULTIVARIATE CLASSIFICATION MODEL
# MINT PLS-DA
# specify data for analysis
train <- b1_corr_out

#
rownames(train) <- train$Sample

train <- train %>% mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>%
  relocate(study)

#

base <- mint.plsda(X = train[,-c(1:5)],  
                   Y = train$Diagnosis, study = as.factor(train$study), ncomp = 1, scale = TRUE)

base$prop_expl_var$X

# trained on dataset 1
base1 <- mint.plsda(X = train[,-c(1:5)], Y = train$Diagnosis, study = as.factor(train$study), ncomp = 1, scale = TRUE)

base1$prop_expl_var$X

#
#
#

# MODEL PERFORMANCE
# prediction on training or test dataset
test <- b1_corr_out
mod <- base

rownames(test) <- test$Sample

#
test <- test %>% mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>%
  relocate(study)

train_pred <- predict(mod,
                      newdata = test[,-c(1:5)],
                      dist = 'max.dist',
                      study.test = as.factor(test$study))

# sample-level prediction
dummy_train <- train_pred[["predict"]][,,1] %>% as.data.frame()
out_train <- train_pred[["class"]][["max.dist"]] %>% as.data.frame() %>% dplyr::select(1)
conf_train <- table(factor(out_train$comp1, levels = levels(as.factor(test$Diagnosis))), 
                  test$Diagnosis)

er_fun(conf_train)
ber_fun(conf_train)
sens_fun(conf_train)
spec_fun(conf_train)

resp <- test$Diagnosis
resp[resp == 'Asthma'] <- 1
resp[resp == 'Not Asthma'] <- 0

roc_train <- roc(response = resp, predictor = dummy_train$Asthma)
auc(roc_train)
coords(roc_train, "best", ret = "threshold")

# save ROC curve plots
tiff('ROC_MINT.tiff', unit = 'mm', res = 300, width = 90, height = 45)
par(mfrow = c(1,2))
plot(roc_train, legacy.axes = TRUE, main = 'Training dataset', cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7,
     mar = c(2,2,1,0.1), mgp = c(1,0.2,0), tck = -0.02)
plot(roc_train, legacy.axes = TRUE, main = 'Test dataset', cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7,
     mar = c(2,2,1,0.1), mgp = c(1,0.2,0), tck = -0.02)
dev.off()


# agreement between predictions for the same patient
out_train1 <- out_train %>% 
  mutate(Sample = rownames(.),
         RAD_ID = str_sub(Sample, end = 6),
         CoreVisit = str_sub(Sample, start = 8)) %>%
  dplyr::select(!Sample) %>%
  pivot_wider(names_from = CoreVisit, values_from = comp1) %>%
  drop_na() %>% mutate(agree = ifelse(CV1 == CV2, 'yes', 'no'))

table(out_train1$agree)

# patient-level prediction
dummy_train_sum <- dummy_train %>% 
  mutate(Sample = rownames(.),
         RAD_ID = str_sub(Sample, end = 6)) %>%
  group_by(RAD_ID) %>% summarise(Asthma = mean(Asthma),
                                 Not_Asthma = mean(`Not Asthma`)) %>%
  mutate(comp1 = ifelse(Asthma > Not_Asthma, 'Asthma', 'Not Asthma'))
  
out_train_pat <- dummy_train_sum %>% dplyr::select(comp1) %>% as.data.frame()
rownames(out_train_pat) <- dummy_train_sum$RAD_ID

input_pat <- test %>% dplyr::select(RAD_ID, Diagnosis) %>% distinct()
input_pat <- input_pat[match(rownames(out_train_pat), input_pat$RAD_ID), ]     

conf_train_pat <- table(factor(out_train_pat$comp1, levels = levels(as.factor(input_pat$Diagnosis))), 
                    input_pat$Diagnosis)

er_fun(conf_train_pat)
ber_fun(conf_train_pat)
sens_fun(conf_train_pat)
spec_fun(conf_train_pat)

resp <- input_pat$Diagnosis
resp[resp == 'Asthma'] <- 1
resp[resp == 'Not Asthma'] <- 0

roc_train <- roc(response = resp, predictor = dummy_train_sum$Asthma)
plot(roc_train, legacy.axes = TRUE)
auc(roc_train)

#
#
#

# prediction performance in cross-validation

# create M-fold cross-validation method for final MINT model
# preserve all observations from one patient in either train or held-out fold
fold_fun <- function(data, nfold){
  p <- n_distinct(data$RAD_ID)
  x <- sample(1:p)
  names(x) <- unique(data$RAD_ID)
  levels <- nfold
  folds <- split(x, x%%levels)
  }

# train model on nfold-1 and test on the left-out fold
train_test <- function(data, folds, m){
  test_id <- folds[[m]]
  test_data <- data %>% filter(RAD_ID %in% names(test_id)) %>%
    mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
  train_data <- data %>% filter(RAD_ID %ni% names(test_id)) %>%
    mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
  model <- mint.plsda(X = train_data[,-c(1:5)], Y = train_data$Diagnosis, 
                       study = as.factor(train_data$study), ncomp = 1)
  model_pred <- predict(model,
                        newdata = test_data[,-c(1:5)],
                        dist = 'max.dist',
                        study.test = as.factor(test_data$study))
  model_pred1 <- append(model_pred, list(test_data$Diagnosis)) %>% append(list(test_data$RAD_ID))
  names(model_pred1)[[9]] <- 'truth'
  names(model_pred1)[[10]] <- 'RAD_ID'
  model_pred1
  }


# use cross-validation to calculate prediction errors in train dataset
set.seed(1410)

cv_perf <- function(nrep, nfold, input_cv, pred_level){
  replicate(n = nrep, expr = {
    cv_folds <- fold_fun(input_cv, nfold)
    
    BER <- vector(mode = "numeric")
    ER <- vector(mode = 'numeric')
    SENS <- vector(mode = 'numeric')
    SPEC <- vector(mode = 'numeric')
    
    for (m in as.character(0:(length(cv_folds) - 1))) {
      cv_model_pred <- train_test(input_cv, cv_folds, m)
      
      if (pred_level == 'patient') {
      # for subject-level prediction use code below
      cv_pred <- cbind(cv_model_pred[['RAD_ID']], as.data.frame(cv_model_pred[["predict"]]))
      colnames(cv_pred) <- c('RAD_ID', 'Asthma_pred', 'Not_Asthma_pred')
      sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% 
        summarise(Asthma_pred = mean(Asthma_pred),
                  Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
        left_join(input_cv %>% dplyr::select(RAD_ID, Diagnosis)) %>%
        mutate(pred = ifelse(Asthma_pred > Not_Asthma_pred, 'Asthma', 'Not Asthma')) %>%
        distinct()
      cv_out_pred <- sum_cv_pred$pred
      cv_conf_mat <- table(factor(cv_out_pred, levels = levels(as.factor(sum_cv_pred$Diagnosis))), sum_cv_pred$Diagnosis)
      
      } else {
      
      # for sample-level prediction use code below
      cv_out_pred <- cv_model_pred[["class"]][["max.dist"]] %>% as.data.frame() %>% dplyr::select(1)
      cv_truth <- cv_model_pred[['truth']]
      cv_conf_mat <- table(factor(cv_out_pred$comp1, levels = levels(as.factor(cv_truth))), cv_truth)
      }
      
      BER[m] <- if (ncol(cv_conf_mat) == 2 & nrow(cv_conf_mat) == 2) {
        ber_fun(cv_conf_mat)
      } else {
        NA}
      ER[m] <- if (ncol(cv_conf_mat) == 2) {
        er_fun(cv_conf_mat)
      } else {
        NA}
      SENS[m] <- if (ncol(cv_conf_mat) == 2 & nrow(cv_conf_mat) == 2) {
        sens_fun(cv_conf_mat)
      } else {
        NA}
      SPEC[m] <- if (ncol(cv_conf_mat) == 2 & nrow(cv_conf_mat) == 2) {
        spec_fun(cv_conf_mat)
      } else {
        NA}
    }
    
    cv_ber <- mean(BER, na.rm = TRUE)
    cv_er <- mean(ER, na.rm = TRUE)
    cv_sens <- mean(SENS, na.rm = TRUE)
    cv_spec <- mean(SPEC, na.rm = TRUE)
    out <- c(cv_ber, cv_er, cv_sens, cv_spec)
  })
  }

# sample-level prediction
cv_perf_sam <- cv_perf(nrep = 100, nfold = 6, input_cv = b2_corr_out, pred_level = 'sample')

ber_sam <- mean(cv_perf_sam[1,], na.rm = TRUE)
er_sam <- mean(cv_perf_sam[2,], na.rm = TRUE)
sens_sam <- mean(cv_perf_sam[3,], na.rm = TRUE)
spec_sam <- mean(cv_perf_sam[4,], na.rm = TRUE)

ber_sam
er_sam
sens_sam
spec_sam

sd(cv_perf_sam[1,], na.rm = TRUE)
sd(cv_perf_sam[2,], na.rm = TRUE)
sd(cv_perf_sam[3,], na.rm = TRUE)
sd(cv_perf_sam[4,], na.rm = TRUE)

# patient-level prediction
cv_perf_pat <- cv_perf(nrep = 100, nfold = 6, input_cv = b2_corr_out, pred_level = 'patient')

ber_pat <- mean(cv_perf_pat[1,], na.rm = TRUE)
er_pat <- mean(cv_perf_pat[2,], na.rm = TRUE)
sens_pat <- mean(cv_perf_pat[3,], na.rm = TRUE)
spec_pat <- mean(cv_perf_pat[4,], na.rm = TRUE)

ber_pat
er_pat
sens_pat
spec_pat

sd(cv_perf_pat[1,], na.rm = TRUE)
sd(cv_perf_pat[2,], na.rm = TRUE)
sd(cv_perf_pat[3,], na.rm = TRUE)
sd(cv_perf_pat[4,], na.rm = TRUE)

#

# use cross-validation to calculate AUROC in train dataset
cv_auc <- function(nrep, nfold, input_cv, pred_level){
  replicate(n = nrep, expr = {
    cv_folds <- fold_fun(input_cv, nfold)
    
    SEN <- list()
    SPE <- list()
    AUC <- vector(mode = 'numeric')
    
    for (m in as.character(0:(length(cv_folds) - 1))) {
      cv_model_pred <- train_test(input_cv, cv_folds, m)
      
      if (pred_level == 'patient') {
        # for subject-level prediction use code below
        cv_pred <- cbind(cv_model_pred[['RAD_ID']], as.data.frame(cv_model_pred[["predict"]]))
        colnames(cv_pred) <- c('RAD_ID', 'Asthma_pred', 'Not_Asthma_pred')
        sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% 
          summarise(Asthma_pred = mean(Asthma_pred),
                    Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
          left_join(input_cv %>% dplyr::select(RAD_ID, Diagnosis)) %>%
          mutate(pred = ifelse(Asthma_pred > Not_Asthma_pred, 'Asthma', 'Not Asthma')) %>%
          distinct()
        dummy_pred_out <- sum_cv_pred %>% dplyr::select(Asthma_pred, Diagnosis) %>%
          mutate(Diagnosis = ifelse(Diagnosis == 'Asthma', 1, 0))
        
      } else {
        
        #for sample-level prediction use code below
        dummy_pred_out <- cbind(as.data.frame(cv_model_pred[['predict']][,,1]), cv_model_pred['truth']) %>%
          dplyr::select(Asthma, truth) %>% rename(Asthma_pred = Asthma, Diagnosis = truth) %>%
          mutate(Diagnosis = ifelse(Diagnosis == 'Asthma', 1, 0))
      }
      
      if (nlevels(as.factor(sum_cv_pred$Diagnosis)) == 2) {
        roc_cv <- roc(response = dummy_pred_out$Diagnosis, predictor = dummy_pred_out$Asthma_pred)
        sen <- roc_cv[['sensitivities']] %>% as.data.frame()
        spe <- roc_cv[['specificities']] %>% as.data.frame()
        auc <- auc(roc_cv)
      } else{
        sen <- c(rep(NA, nrow(test_data))) %>% as.data.frame()  
        spe <- c(rep(NA, nrow(test_data))) %>% as.data.frame()
      }
      
      SEN[[m]] <- sen # DEAL WITH UNEVEN TEST FOLD SIZES (different length of sens/spec vectors)
      SPE[[m]] <- spe
      AUC[m] <- auc
    }
    
    names(SEN) <- c(1:length(SEN))
    fold_names <- list_c(lapply(names(SEN), function(m) {
      rep(m, nrow(SEN[[m]]))
    }))
    names(AUC) <- c(1:length(AUC))
    auc <- list_c(lapply(names(AUC), function(x){
      rep(AUC[x], nrow(SEN[[x]]))
    }))
    
    SenSpec <- list_c(SEN) %>% cbind(list_c(SPE), fold_names, auc) %>% drop_na()
    
  }, simplify = FALSE)
}

# sample-level prediction
set.seed(1410)

cv_auc_sam <- cv_auc(nrep = 100, nfold = 2, input_cv = b2_corr_out, pred_level = 'patient')
      
names(cv_auc_sam) <- c(1:length(cv_auc_sam))
rep_names <- list_c(lapply(names(cv_auc_sam), function(m) {
  rep(m, nrow(cv_auc_sam[[m]]))
}))

rep_SenSpec <- list_c(cv_auc_sam) %>% cbind(rep_names) 

# patient-level prediction
set.seed(1410)

cv_auc_pat <- cv_auc(nrep = 100, nfold = 2, input_cv = b2_corr_out, pred_level = 'patient')

names(cv_auc_pat) <- c(1:length(cv_auc_pat))
rep_names <- list_c(lapply(names(cv_auc_pat), function(m) {
  rep(m, nrow(cv_auc_pat[[m]]))
}))

rep_SenSpec <- list_c(cv_auc_pat) %>% cbind(rep_names)

# results interpretation and visualisation for any level prediction
colnames(rep_SenSpec)[1:2] <- c('Sen', 'Spec')
rep_SenSpec$rep_fold = paste(rep_SenSpec$rep, rep_SenSpec$fold, sep = '_')
rep_SenSpec$fpr = 1-rep_SenSpec$Spec

auc_unique <- rep_SenSpec %>% dplyr::select(rep_fold, auc) %>% distinct()
mean(auc_unique$auc)
sd(auc_unique$auc)
median(auc_unique$auc)
IQR(auc_unique$auc)

# plot ROC from each cv repetition
rep_SenSpec %>%
  ggplot(aes(x = fpr, y = Sen, group = rep_fold)) +
  scale_x_reverse() +
  geom_step(colour = 'gray40', size = 0.1) + 
  geom_abline(aes(intercept = 0, slope = -1), colour = 'black', alpha = 1, linewidth = 1) +
  theme_bw() +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) +
  ylab('Sensitivity') +
  xlab('1 - Specificity') +
  ggtitle('Train data (2-fold cross-validation)') +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave('AUROC_CV_B2.tiff', unit = 'mm', dpi = 300, width = 120, height = 100)

#
#
#

# Model loadings
# extract loadings from the model and annotate with VOC origin categories
# specify model for analysis
mod <- base1

#
loads_glob <- mod$loadings$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>% 
  left_join(endo_exo %>% dplyr::select(comp, FC_filter)) %>%
  distinct()

# exctract loadings vip
vip <- vip(mod) %>% as.data.frame() %>% mutate(comp = rownames(.)) %>%
  rename(vip = comp1)

loads_glob <- loads_glob %>% left_join(vip) %>%
  mutate(FC_filter = ifelse(FC_filter == 'Endo_both_datasets', 'Breath-enriched (both datasets)',
                            ifelse(FC_filter == 'Exo', 'Ambiguous', 'Breath-enriched (one dataset)')))

# plot most important loadings 
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

Lglob1 <- loads_glob  %>%
  arrange(desc(abs(comp1))) %>% 
  slice(1:20) %>%
  mutate(comp = str_sub(comp, end = 32)) %>%
  ggplot(aes(x = comp1, y = fct_inorder(as.factor(comp)), 
             fill = factor(FC_filter, levels = c('Ambiguous',
                                                 'Breath-enriched (one dataset)',
                                                 'Breath-enriched (both datasets)')))) +
  geom_col(colour = 'grey30', lwd = 0.3) +
  theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.2), 'cm'),
        legend.key.size = unit(0.5, 'cm')) +
  scale_fill_manual(name = 'VOC classification',
                    values = c('Ambiguous' = my_pal[[4]],
                               'Breath-enriched (both datasets)' = my_pal[[3]],
                               'Breath-enriched (one dataset)' = my_pal[[7]]))

Lglob1

l1_mint <- Lglob1 + ggtitle('Validation dataset') +
  scale_y_discrete(labels = c('X1_Propene._1_.methylthio._._.E.' = '1_Propene_1_(methylthio)_(E)_',
                              'X1_Propene._1_.methylthio._._.Z.' = '1_Propene_1_(methylthio)_(Z)_',
                              'Bicyclo.3.1.0.hexane._4_methylen' = 'Bicyclo.3.1.0.hexane._4_',
                              'X1_Propanol._2_ethoxy_' = '1_Propanol_2_ethoxy',
                              'X3_methylbutan_1_ol' = '3_methylbutan_1_ol',
                              'X2_Butanone' = '2_Butanone',
                              'X3_Octene._.Z._' = '3_Octene_Z_',
                              'X1_pentanol' = '1_Pentanol'))

l2_mint <- Lglob1 + ggtitle('Training dataset') +
  scale_y_discrete(labels = c('X.S._..._6_Methyl_1_octanol' = '(S)_(+)_6_Methyl_1_octanol',
                              'X4_ethyl_2.2_dimethylhexane' = '4_ethyl_2,2_dimethylhexane',
                              'Acetic_acid._1.7.7_trimethyl_bic' = 'Acetic_acid._1.7.7_trimethyl_',
                              'Ethane._1.1.2_trichloro_1.2.2_tr' = 'Ethane._1.1.2_trichloro_'))

leg3 <- get_legend(Lglob1)

plot(leg3)

# assemble multipanel figure
l_plots_mint <- plot_grid(l2_mint, l1_mint, nrow = 1) +
  theme(plot.margin = unit(c(0.4, 0.1, 0.3, 0.1), 'cm')) +
  draw_label(label = 'MINT PLS-DA loadings', x = 0.6, fontface = 'bold', size = 8, y = 1, vjust = -0.5) +
  draw_label(label = 'Component 1', x = 0.6, y = 0, size = 8)

plot(l_plots_mint)

ggsave('Loadings_MINT.tiff', 
       l_plots_mint, dpi = 300, unit = 'mm', width = 110, height = 75)

# ratio of endo to exogenous VOCs among top loadings
loads_glob1 <- loads_glob %>% arrange(desc(abs(comp1))) %>% filter(vip > 1) %>% 
  mutate(comp1_class = ifelse(comp1 > 0, 'pos', 'neg')) 

loads_glob_sum <- loads_glob1 %>% group_by(comp1_class, FC_filter) %>% summarise(n = n()) %>%
  mutate(ratio = ifelse(comp1_class == 'pos', n/table(loads_glob1$comp1_class)[2], 
                        n/table(loads_glob1$comp1_class)[1]))

#

Lglob2 <- loads_glob_sum %>%
  ggplot(aes(x = comp1_class, y = ratio*100, fill = factor(FC_filter, levels = c('Ambiguous',
                                                 'Breath-enriched (one dataset)',
                                                 'Breath-enriched (both datasets)')))) +
  geom_col(colour = 'grey30', lwd = 0.3) +
  scale_x_discrete(labels = c('neg' = '< 0',
                              'pos' = '> 0')) +
  ylab('% VOCs') +
  theme_bw(base_size = 8) +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        legend.position = 'none',
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))  +
  scale_fill_manual(name = 'VOC classification',
                    values = c('Ambiguous' = my_pal[[4]],
                               'Breath-enriched (both datasets)' = my_pal[[3]],
                               'Breath-enriched (one dataset)' = my_pal[[7]]))

Lglob2

l_sum_b1 <- Lglob2 + ggtitle('Validation \n dataset') + theme(axis.title.y = element_blank())
l_sum_b2 <- Lglob2 + ggtitle('Training \n dataset') 

l_sum_plots <- plot_grid(l_sum_b2, l_sum_b1, nrow = 1, rel_widths = c(0.52, 0.48)) +
  theme(plot.margin = unit(c(1, 0.1, 0.3, 0.1), 'cm')) +
  draw_label(label = 'Loadings', size = 8, x = 0.5, y = 0) +
  draw_label(label = 'Classification of the \n most important loadings \n (VIP > 1) in MINT PLS-DA', size = 8,
             fontface = 'bold', x = 0.6, y = 1, vjust = -0.1)

plot(l_sum_plots)

ggsave('Bar_plot_class_loadings.tiff', l_sum_plots, unit = 'mm', dpi = 300, width = 60, height = 50)

# Figure 3d
plot3d <- plot_grid(l_sum_plots, leg3, ncol = 1, rel_heights = c(0.65, 0.35))
plot(plot3d)

plot3cd <- plot_grid(l_plots_mint, plot3d, labels = c('c)', 'd)'), label_size = 12,
                  rel_widths = c(0.7, 0.3))
plot(plot3cd)

################################


#
#
#

# Model variates
# extract and plot global variates
pals <- grafify::graf_palettes
my_pal = pals$fishy
swatch(my_pal)

plot_variates <- function(variates, data) {
  smint_scores <- variates %>% as.data.frame() %>% mutate(Sample = rownames(.)) %>%
    left_join(data %>% dplyr::select(Sample, Diagnosis, CoreVisit))
  
  colnames(smint_scores)[1] <- 'comp1'
  
  sp1 <- smint_scores %>% ggplot(aes(x = Diagnosis, y = comp1, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE, lwd = 0.3) +
    geom_jitter(size = 0.5, alpha = 0.5) +
    theme_bw(base_size = 8) +
    scale_fill_manual(values = c('Asthma' = my_pal[[2]],
                                 'Not Asthma' = my_pal[[1]])) +
    theme(plot.title =  element_text(hjust = 0.5, size = 8),
          panel.grid = element_blank(),
          plot.margin = unit(c(0.5, 0.1, 0.2, 0.1), 'cm'),
          axis.title.x = element_blank(),
          legend.key.width = unit(0.5, 'cm'))

}

sp1 <- plot_variates(base$variates$X, b2_corr_out) +
  ylab(paste('Component 1 (', round(as.numeric(base$prop_expl_var$X[3])*100, 2), '% )')) +
  theme(legend.position = 'none')  +
  ggtitle('Training dataset')

sp2 <- plot_variates(train_pred$variates, b1_corr_out) +
  theme(axis.title.y = element_blank()) +
  ggtitle('Validation dataset') +
  ylab('')

sp <- plot_grid(sp1, sp2, nrow = 1, rel_widths = c(0.39, 0.61)) +
  theme(plot.margin = unit(c(0.2, 0, 0.2, 0), 'cm')) +
  draw_label(label = 'MINT PLS-DA score plots', y = 1, x = 0.5, fontface = 'bold', size = 8) +
  draw_label(label = 'Diagnosis', y = 0, x = 0.5, size = 8)

plot(sp)

ggsave('Score_plots_MINT.tiff', sp, unit = 'mm', dpi = 300, width = 120, height = 60)

##################################

# Figure 3a
plot3a <- sp

# Figure 3b
plot_roc <- function(roc) {
  roc1 <- roc
  roc_df <- data.frame(sens = roc1$sensitivities,
                       spec = roc1$specificities) %>%
    mutate(fpr = 1-spec) %>%
    arrange(sens)
  
  roc_df %>%
  ggplot(aes(x = fpr, y = sens)) +
  geom_step() + 
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  geom_abline(aes(intercept = 0, slope = 1), colour = 'black', alpha = 1, linewidth = 0.3) +
  theme_bw(base_size = 8) +
  ylab('Sensitivity') +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.4, 0.1, 0.2, 0.1), 'cm')) +
  annotate(geom = 'text', label = paste('AUROC =', round(auc(roc), 2)), x = 0.7, y = 0.05, size = 2)
}

plot3b1 <- plot_roc(roc_train) + ggtitle('Training dataset') 
plot3b2 <- plot_roc(roc_test) + ggtitle('Validation dataset') + ylab('')

plot3b <- plot_grid(plot3b1, plot3b2, nrow = 1) +
  theme(plot.margin = unit(c(0.3, 0, 0.2, 0.1), 'cm')) +
  draw_label(label = 'Receiver Operating Characteristics curve for \n MINT-PLSDA prediction scores', x = 0.55, y = 1, size = 8, fontface = 'bold') +
  draw_label(label = '1 - Specificity', x = 0.5, y = 0, size = 8)

plot(plot3b)

# assemble top figure panel
plot3ab <- plot_grid(plot3a, plot3b, nrow = 1, labels = c('a)', 'b)'), label_size = 12, rel_widths = c(0.55, 0.45)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.4, 0.1), 'cm'))
plot(plot3ab)

###################################

# Figure 3
fig3 <- plot_grid(plot3ab, plot3cd, ncol = 1, rel_heights = c(0.4, 0.6))

fig3

ggsave('fig3.tiff', unit = 'mm', dpi = 300, width = 157, height = 120)

pdf('Figure3.pdf', width = 6.18, height = 4.72)
plot(fig3)
dev.off()
