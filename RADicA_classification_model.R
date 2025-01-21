## Predicting asthma based on breath VOC samples

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

# load summarised BG corrected dataset 
b1_imp_sum_corr_w <- read.csv('RADicA_BG_adjusted.csv')[,-1]

# w/o multivariate outliers
b1_corr_w1 <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b_corr_w1 <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]

# set columns in the same order
b1_corr_w1 <- b1_corr_w1[names(b_corr_w1)]

#
b1_corr_w1 <- b1_corr_w1 %>% dplyr::select(!'X2_butanone')
write.csv(b1_corr_w1, 'RADicA_BG_adjusted_B1_outl_removed.csv')
b_corr_w1 <- b_corr_w1 %>% dplyr::select(!'X2_butanone')
write.csv(b_corr_w1, 'RADicA_BG_adjusted_B2_outl_removed.csv')
#

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')[,-1]

reg_valid_join <- read.csv('RegValid_results.csv')[,-1]

endo <- reg_valid_join %>% filter(FC_filter != 'Exo') %>%
  filter(comp != 'Phenylethyne')

b1_endo <- b1_corr_w1 %>% dplyr::select(1:4, endo$comp) 
b2_endo <- b_corr_w1 %>% dplyr::select(1:4, endo$comp)

n_distinct(b_corr_w1$Sample)
n_distinct(b1_corr_w1$Sample)

#
#
#

# custom functions
ber_fun <-  function(conf_mat) {
  ber <- 0.5*((conf_mat[1,2]/(conf_mat[1,2]+conf_mat[1,1])) + 
                (conf_mat[2,1]/(conf_mat[2,1] + conf_mat[2,2])))
  print(ber)
}

er_fun <- function(conf_mat){
  er <- (sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
  print(er)
}

sens_fun <- function(conf_mat){
  sens <- conf_mat[1,1]/(conf_mat[1,1] + conf_mat[1,2])
  print(sens)
}

spec_fun <- function(conf_mat){
  spec <- conf_mat[2,2]/(conf_mat[2,2] + conf_mat[2,1])
  print(spec)
}

'%ni%' <- Negate('%in%')

#
#
#

# correlation structure of the datasets
vocs <- colnames(b1_corr_w1[,-c(1:4)])

b1_cm <- b1_corr_w1[,-c(1:4)]
colnames(b1_cm) <- sapply(1:145, FUN = function(x){
  paste('V', x, sep = '')})

b_cm <- b_corr_w1[,-c(1:4)]
colnames(b_cm) <- sapply(1:145, FUN = function(x){
  paste('V', x, sep = '')})

cor_b1 <- cor(b1_cm)
p_mat_b1 <- cor_pmat(b1_cm)

cor_b2 <- cor(b_cm)
p_mat_b2 <- cor_pmat(b_cm)

dev.new()
cp1 <- ggcorrplot(cor_b1, lab = FALSE, hc.method = 'ward.D2', hc.order = TRUE,
           p.mat = p_mat_b1, insig = 'blank', tl.cex = 3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Dataset 1 breath VOC correlation matrix')
ggsave('B1_correlation_plot.tiff', cp1, dpi = 300, unit = 'mm', width = 150, height = 150)

dev.new()
cp2 <- ggcorrplot(cor_b2, lab = FALSE, hc.method = 'ward.D2', hc.order = TRUE,
           p.mat = p_mat_b2, insig = 'blank', tl.cex = 3) +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Dataset 2 breath VOC correlation matrix')
ggsave('B2_correlation_plot.tiff', cp2, dpi = 300, unit = 'mm', width = 150, height = 150)

# retrieve hierarchical clusters
hc1 <- hclust(dist(cor_b1), method = 'ward.D2')
dev.new()
plot(hc1, cex = 0.5)
hc1_cut <- cutree(hc1, h = 5)
n_distinct(hc1_cut)

write.csv(hc1_cut, 'VOC_corr_hclust_B1.csv')

hc2 <- hclust(dist(cor_b2), method = 'ward.D2')
dev.new()
plot(hc2, cex = 0.5)
hc2_cut <- cutree(hc2, h = 5)
n_distinct(hc2_cut)

write.csv(hc2_cut, 'VOC_corr_hclust_B2.csv')

names(hc1_cut) <- vocs
names(hc2_cut) <- vocs

# number of variables overlapping between clusters
hc_ag <- matrix(nrow = 11, ncol = 10)

for(n in 1:10){
  for(i in 1:11) {
    agreement <- intersect(names(hc1_cut[hc1_cut == i]), names(hc2_cut[hc2_cut == n]))
    print(n_distinct(agreement))
    hc_ag[i,n] <- n_distinct(agreement)
    }
  }

reg_valid_join %>% filter(model == 'Model1') %>%
  left_join(as.data.frame(vocs) %>% mutate(VN = rownames(.)) %>% rename(comp = vocs)) %>%
  filter(comp %in% intersect(names(hc1_cut[hc1_cut == 7]), names(hc2_cut[hc2_cut == 4]))) %>%
  dplyr::select(c(sign, wilcox_filter, FC_filter, comp, VN))


#
#
#

# MINT
rownames(b1_corr_w1) <- b1_corr_w1$Sample

X <- b1_corr_w1[,-c(1:4)]
Y <- b1_corr_w1$Diagnosis
study <- b1_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2

study <- as.factor(study)

set.seed(1410)
base <- mint.plsda(X, Y, study = study, ncomp = 1, scale = TRUE)
base$prop_expl_var$X

# 1 component
# global variates
plot_variates <- function(model, variates, data) {
  smint_scores <- variates %>% as.data.frame() %>% mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Diagnosis, CoreVisit))
  
  colnames(smint_scores)[1] <- 'comp1'
  
  sp1 <- smint_scores %>% ggplot(aes(x = Diagnosis, y = comp1, fill = Diagnosis)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 1, alpha = 0.7) +
  theme_bw(base_size = 10) +
  scale_fill_manual(values = c('Not Asthma' = 'darkorange1',
                               'Asthma' = 'dodgerblue2')) +
  ylab(paste('Component 1 (', round(as.numeric(model$prop_expl_var$X[3])*100, 2), '% )')) +
  ggtitle('Global scores plot') +
  theme(plot.title =  element_text(hjust = 0.5))
  
  print(sp1)
}

sp1 <- plot_variates(base$variates$X, b_corr_w1)
sp1_b1 <- plot_variates(base1$variates$X, b1_corr_w1)

# partial variates    
plot_variates_part <- function(variates, data) {
  smint_scores_part <- variates[1] %>% as.data.frame() %>%
  rbind(variates[2] %>% as.data.frame()) %>%
  mutate(Sample = rownames(.)) %>%
  left_join(data %>% dplyr::select(Sample, Diagnosis, CoreVisit))
  
  smint_scores_part %>% 
    ggplot(aes(x = Diagnosis, y = comp1, fill = Diagnosis)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(size = 1, alpha = 0.7) +
    theme_bw(base_size = 10) +
    facet_wrap(~ CoreVisit) +
    scale_fill_manual(values = c('Not Asthma' = 'darkorange1',
                                 'Asthma' = 'dodgerblue2')) +
    ggtitle('Partial scores plots') +
    theme(plot.title =  element_text(hjust = 0.5)) +
    ylab('Component 1')
}

sp2 <- plot_variates_part(base$variates.partial$X, b_corr_w1)
sp2_b1 <- plot_variates_part(base1$variates.partial$X, b1_corr_w1)
  
sp <-  plot_grid(sp1, sp2, nrow = 2)
sp

ggsave('Scores_plots_sMINT.tiff', sp, unit = 'mm', dpi = 300, width = 100, height = 130)

sp_b1 <- plot_grid(sp1_b1, sp2_b1, nrow = 2)
sp_b1

ggsave('Scores_plots_sMINT_B1.tiff', sp_b1, unit = 'mm', dpi = 300, width = 100, height = 130)

#

plotVar(base, cutoff = 0.6)

plotLoadings(base, ndisplay = 15, comp = 1, 
             study = 'global', 
             method = 'median',
             contrib = 'max')


# tuning
base_opt <- perf(base, dist = 'all')
plot(base_opt)

base_opt$global.error$BER
base_opt$global.error$overall
base_opt$global.error$error.rate.class
base_opt$study.specific.error[['2']]

opt_pred <- base_opt[["class"]][["centroids.dist"]] %>% as.data.frame() 

conf_mat <- table(factor(b_corr_w1$Diagnosis, levels = levels(as.factor(b_corr_w1$Diagnosis))), 
                  opt_pred$comp1)

ber_fun(conf_mat)
er_fun(conf_mat)

#

# create M-fold cross-validation method for final MINT model
# preserve all observations from one patient in either train or held-out fold
data <- b1_corr_w1

set.seed(1410)

cv_perf <- 
  replicate(n = 10, expr = {
    p <- n_distinct(data$RAD_ID)
    x <- sample(1:p)
    names(x) <- unique(data$RAD_ID)
    levels <- 8
    folds <- split(x, x%%levels)
    
    BER <- vector(mode = "numeric")
    ER <- vector(mode = 'numeric')
    SENS <- vector(mode = 'numeric')
    SPEC <- vector(mode = 'numeric')
    
    for (m in as.character(0:(levels - 1))) {
      test_id <- folds[[m]]
      test_data <- data %>% filter(RAD_ID %in% names(test_id)) %>%
        mutate(study = '3')
      train_data <- data %>% filter(RAD_ID %ni% names(test_id)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2'))
      conc_data <- rbind(test_data, train_data) %>%
        relocate(study)
      study <- as.factor(conc_data$study)
      test <- which(study == '3')
      X <- conc_data[,-c(1:5)]
      Y <- conc_data$Diagnosis
      model <- mint.splsda(X = X[-c(test),], Y = Y[-c(test)], 
                          study = droplevels(study[-c(test)]), ncomp = 1)
      model_pred <- predict(model,
                            newdata = X[test,],
                            dist = 'max.dist',
                            study.test = factor(study[test]))
      # for subject-level prediction use code below
      #cv_pred <- cbind(test_data$RAD_ID, as.data.frame(model_pred[["predict"]])) %>%
      #rename(RAD_ID = 'test_data$RAD_ID',
       #      Asthma_pred = "Asthma.dim1",
        #    Not_Asthma_pred = "Not Asthma.dim1")
      #sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% 
      #summarise(Asthma_pred = mean(Asthma_pred),
       #         Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
      #left_join(data %>% dplyr::select(RAD_ID, Diagnosis)) %>%
      #mutate(pred = ifelse(abs(Asthma_pred - mean(Asthma_pred)) > # centroids dist?
      #                           abs(Not_Asthma_pred - mean(Not_Asthma_pred)), 
       #                              'Asthma', 'Not Asthma')) %>%
      #distinct()
      #out_pred <- sum_cv_pred$pred
      #conf_mat <- table(factor(out_pred, levels = levels(as.factor(out_pred))), 
      #sum_cv_pred$Diagnosis)
      # for sample-level prediction use code below
      out_pred <- model_pred[["class"]][["max.dist"]] %>% as.data.frame() %>% dplyr::select(1)
      conf_mat <- table(factor(out_pred$comp1, levels = levels(as.factor(test_data$Diagnosis))), 
                        test_data$Diagnosis)
      BER[m] <- if (ncol(conf_mat) == 2 & nrow(conf_mat) == 2) {
        ber_fun(conf_mat)
      } else {
        NA}
      ER[m] <- if (ncol(conf_mat) == 2) {
        er_fun(conf_mat)
      } else {
        NA}
      SENS[m] <- if (ncol(conf_mat) == 2 & nrow(conf_mat) == 2) {
        sens_fun(conf_mat)
      } else {
        NA}
      SPEC[m] <- if (ncol(conf_mat) == 2 & nrow(conf_mat) == 2) {
        spec_fun(conf_mat)
      } else {
        NA}
    }
    cv_ber <- mean(BER, na.rm = TRUE)
    cv_er <- mean(ER, na.rm = TRUE)
    cv_sens <- mean(SENS, na.rm = TRUE)
    cv_spec <- mean(SPEC, na.rm = TRUE)
    out <- c(cv_ber, cv_er, cv_sens, cv_spec)
  })


ber <- mean(cv_perf[1,], na.rm = TRUE)
er <- mean(cv_perf[2,], na.rm = TRUE)
sens <- mean(cv_perf[3,], na.rm = TRUE)
spec <- mean(cv_perf[4,], na.rm = TRUE)

ber
er
sens
spec

# M-fold CV for variable number selection
data <- b_corr_w1
lv <- 5
nvar <- 10

set.seed(123)

cv_perf_vs <-
  replicate(n = 1, expr = {
    p <- n_distinct(data$RAD_ID)
    x <- sample(1:p)
    names(x) <- unique(data$RAD_ID)
    levels <- lv
    folds <- split(x, x%%levels)
    
    vs_BER <- data.frame()
    vs_vocs <- list()
  
    
    for (m in as.character(0:(levels - 1))) {
      test_id <- folds[[m]]
      test_data <- data %>% filter(RAD_ID %in% names(test_id)) %>%
        mutate(study = '3')
      train_data <- data %>% filter(RAD_ID %ni% names(test_id)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2'))
      conc_data <- rbind(test_data, train_data) %>%
        relocate(study)
      study <- as.factor(conc_data$study)
      test <- which(study == '3')
      X <- conc_data[,-c(1:5)]
      Y <- conc_data$Diagnosis
      
      vs_ber <- function(vn) {
        BER <- vector(mode = "numeric", vn)
        VOCS <- list()
        for(n in seq_len(vn)) {
        model <- mint.splsda(X = X[-c(test),], Y = Y[-c(test)], 
                        study = droplevels(study[-c(test)]), ncomp = 1,
                        keepX = c(n)) 
        model_pred <- predict(model,
                            newdata = X[test,],
                            dist = 'centroids.dist',
                            study.test = factor(study[test]))
      # for subject-level prediction use code below
      #cv_pred <- cbind(test_data$RAD_ID, as.data.frame(model_pred[["predict"]])) %>%
      # rename(RAD_ID = 'test_data$RAD_ID',
      #       Asthma_pred = "Asthma.dim1",
      #      Not_Asthma_pred = "Not Asthma.dim1")
      #sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% 
      # summarise(Asthma_pred = mean(Asthma_pred),
      #          Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
      #left_join(data %>% dplyr::select(RAD_ID, Diagnosis)) %>%
      #mutate(pred = ifelse(Asthma_pred > Not_Asthma_pred, 'Asthma', 'Not Asthma')) %>%
      #distinct()
      #out_pred <- sum_cv_pred$pred
      #conf_mat <- table(factor(out_pred, levels = levels(sum_cv_pred$Diagnosis)), 
      #sum_cv_pred$Diagnosis)
      # for sample-level prediction use code below
        
      loads <- model[['loadings']][['X']] %>% as.data.frame() %>% filter(comp1 != 0) 
      vocs <- rownames(loads)

      out_pred <- model_pred[["class"]][["centroids.dist"]] %>% as.data.frame() %>% dplyr::select(1)
      conf_mat <- table(factor(out_pred$comp1, levels = levels(as.factor(test_data$Diagnosis))), 
                        test_data$Diagnosis)
      VOCS <- append(VOCS, list(vocs))
      BER[n] <- if (ncol(conf_mat) == 2) {
        ber_fun(conf_mat)
        } else {
        NA}
      
      }
        out <- list(BER, VOCS)
      }
      vs_out <- vs_ber(nvar)
      vs_BER = rbind(vs_out[[1]], vs_BER)
      colnames(vs_BER) <- sapply(1:nvar, function(g){paste('VN', g, sep = '')})
      vs_vocs <- append(vs_vocs, list(unlist(vs_out[[2]])))
    }
    out <- list(vs_BER, vs_vocs)
  }, simplify = F)


cv_perf_vs_df <- bind_rows(cv_perf_vs[[1]][[1]])
cv_perf_vs_df_sum <- cv_perf_vs_df %>% pivot_longer(cols = everything(), 
                        names_to = 'VN', values_to = 'BER') %>%
  mutate(VN = as.factor(VN)) %>%
  group_by(VN) %>%
  summarise(mean = mean(BER, na.rm = TRUE), sd = sd(BER, na.rm = TRUE))

cv_perf_vs_df_sum$VN <- factor(cv_perf_vs_df_sum$VN, levels = c(str_sort(cv_perf_vs_df_sum$VN, numeric = T)))

dev.new()
cv_perf_vs_df_sum %>%
  ggplot(aes(x = fct_inorder(VN), y = mean, group = VN)) + geom_point(size = 1) +
  geom_errorbar(aes(min = mean-sd, max = mean+sd), linewidth = 0.3) +
  scale_x_discrete(limits  = levels(cv_perf_vs_df_sum$VN),
                   labels = c(1, rep('', 8), 10, rep('', 9), 20, rep('', 9), 30,
                              rep('', 9), 40, rep('', 9), 50, rep('', 9), 60, 
                              rep('', 9), 70, rep('', 9), 80, rep('', 9), 90, 
                              rep('', 9), 100)) +
  ggtitle('Tuning of variable number in sparse multi-group PLS-DA model') +
  xlab('Number of variables') +
  theme_bw() +
  ylab('Balanced Error Rate') 

ggsave('sMINT_tuning_var_number_B2_50vars.tiff', unit = 'mm', dpi = 300, width = 180, height = 80)

#

vocs_freq <- table(unlist(cv_perf_vs[[1]][[2]])) %>% as.data.frame() %>%
  mutate(stab = Freq/(lv*nvar))

View(vocs_freq %>% filter(stab > 0.5))

write.csv(vocs_freq, 'sMINT_variable_stability.csv')

#
#
#

# tune variable number using validation dataset
vs_ber <- function(vn) {
  BER <- vector(mode = "numeric", vn)
  for(n in seq_len(vn)) {
    train <- b_corr_w1
    test_b1_df <- b1_corr_w1 %>% #filter(CoreVisit == 'CV2') %>%
      mutate(CoreVisit = 'B1_CV')
    conc_corr_w1 <- rbind(train, test_b1_df)
    rownames(conc_corr_w1) <- conc_corr_w1$Sample
    X <- conc_corr_w1[,-c(1:4)]
    Y = conc_corr_w1$Diagnosis
    study <- conc_corr_w1$CoreVisit
    study[study == 'CV1'] <- 1
    study[study == 'CV2'] <- 2
    study[study == 'B1_CV'] <- 3
    study <- as.factor(study)
    test_b1 <- which(study == '3')
    base <- mint.splsda(X = X[-c(test_b1),],
                       Y = Y[-c(test_b1)],
                       study = droplevels(study[-c(test_b1)]),
                       ncomp = 1,
                       keepX = c(n))
    model_pred <- predict(base,
                          newdata = X[test_b1,],
                          dist = 'centroids.dist',
                          study.test = factor(study[test_b1]))
    # for sample-level prediction use code below
    out_pred <- model_pred[["class"]][["centroids.dist"]] %>% as.data.frame() %>% dplyr::select(1)
    conf_mat <- table(factor(out_pred$comp1, levels = levels(as.factor(test_b1_df$Diagnosis))), 
                      test_b1_df$Diagnosis)
    BER[n] <- if (ncol(conf_mat) == 2) {
      ber_fun(conf_mat)
    } else {
      NA}
    
  }
  BER
}

set.seed(123)
vs_test <- vs_ber(50)
names(vs_test) <- c(1:50)
plot(names(vs_test), vs_test)
vs_test[vs_test == min(vs_test)]

vs_test %>% as.data.frame() %>% rename(BER = '.') %>%
  ggplot(aes(x = as.integer(rownames(.)), y = BER)) + geom_point() +
  xlab('Number of variables') +
  theme_bw() +
  ylab('Balanced Error Rate') 
ggsave('sMINT_tuning_var_number_test_B1.tiff', unit = 'mm', dpi = 300, width = 180, height = 80)


# CV for AUROC calculation
set.seet(1410)

data <- b_corr_w1

auc_perf <- 
  replicate(n = 1000, expr = {
    p <- n_distinct(data$RAD_ID)
    x <- sample(1:p)
    names(x) <- unique(data$RAD_ID)
    levels <- 2
    folds <- split(x, x%%levels)
    
    SEN <- list()
    SPE <- list()
    AUC <- vector(mode = 'numeric')
    
    for (m in as.character(0:(levels - 1))) {
      test_id <- folds[[m]]
      test_data <- data %>% filter(RAD_ID %in% names(test_id)) %>%
        mutate(study = '3')
      train_data <- data %>% filter(RAD_ID %ni% names(test_id)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2'))
      conc_data <- rbind(test_data, train_data) %>%
        relocate(study)
      study <- as.factor(conc_data$study)
      test <- which(study == '3')
      X <- conc_data[,-c(1:5)]
      Y <- conc_data$Diagnosis
      model <- mint.splsda(X = X[-c(test),], Y = Y[-c(test)], 
                           study = droplevels(study[-c(test)]), ncomp = 1)
      model_pred <- predict(model,
                            newdata = X[test,],
                            dist = 'max.dist',
                            study.test = factor(study[test]))
      # for subject-level prediction use code below
      #cv_pred <- cbind(test_data$RAD_ID, as.data.frame(model_pred[["predict"]])) %>%
      #rename(RAD_ID = 'test_data$RAD_ID',
      #      Asthma_pred = "Asthma.dim1",
      #    Not_Asthma_pred = "Not Asthma.dim1")
      #sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% 
      #summarise(Asthma_pred = mean(Asthma_pred),
      #         Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
      #left_join(data %>% dplyr::select(RAD_ID, Diagnosis)) %>%
      #mutate(pred = ifelse(abs(Asthma_pred - mean(Asthma_pred)) > # centroids dist?
      #                           abs(Not_Asthma_pred - mean(Not_Asthma_pred)), 
      #                              'Asthma', 'Not Asthma')) %>%
      #distinct()
      #out_pred <- sum_cv_pred$pred
      #conf_mat <- table(factor(out_pred, levels = levels(as.factor(out_pred))), 
      #sum_cv_pred$Diagnosis)
      # for sample-level prediction use code below
      dummy_pred <- model_pred[['predict']][,,1] 
      if (nlevels(as.factor(test_data$Diagnosis)) == 2) {
        roc_cv <- roc(response = test_data$Diagnosis, predictor = dummy_pred[,1])
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

names(auc_perf) <- c(1:length(auc_perf))
rep_names <- list_c(lapply(names(auc_perf), function(m) {
  rep(m, nrow(auc_perf[[m]]))
}))

rep_SenSpec <- list_c(auc_perf) %>% cbind(rep_names) 
colnames(rep_SenSpec)[1:2] <- c('Sen', 'Spec')
rep_SenSpec$rep_fold = paste(rep_SenSpec$rep, rep_SenSpec$fold, sep = '_')
rep_SenSpec$fpr = 1-rep_SenSpec$Spec

rep_SenSpec %>%
ggplot(aes(x = fpr, y = Sen, group = rep_fold)) +
  scale_x_reverse() +
  geom_step(colour = 'gray40', size = 0.1) + 
  geom_abline(aes(intercept = 0, slope = -1), colour = 'grey90', alpha = 0.6, linewidth = 0.5) +
  theme_bw() +
  coord_fixed(xlim=c(0,1), ylim=c(0,1)) +
  ylab('Sensitivity') +
  xlab('1 - Specificity') +
  ggtitle('Train data (2-fold cross-validation)') +
  theme(plot.title = element_text(hjust = 0.5)) 

ggsave('AUROC_CV_B2.tiff', unit = 'mm', dpi = 300, width = 120, height = 100)


auc_unique <- rep_SenSpec %>% dplyr::select(rep_fold, auc) %>% distinct()
mean(auc_unique$auc)
sd(auc_unique$auc)
auc_unique$auc_null <- rep(0.5, 300)
t.test(auc_unique$auc, auc_unique$auc_null)

spec95 <- rep_SenSpec %>% filter(Spec > 0.895 & Spec < 0.905) %>%
  group_by(rep_fold) %>% summarise(Sen = mean(Sen))
mean(spec95$Sen)
sd(spec95$Sen)

sen95 <- rep_SenSpec %>% filter(Sen > 0.895 & Sen < 0.905) %>%
  group_by(rep_fold) %>% summarise(Spec = mean(Spec))
mean(sen95$Spec)
sd(sen95$Spec)

#
#
#

# compare LOGOCV error rates to plsda error rates
cv1 <- b_corr_w1 %>% filter(CoreVisit == 'CV1')
cv2 <- b_corr_w1 %>% filter(CoreVisit == 'CV2')

pls.da <- plsda(X = cv1[,-c(1:4)], Y = cv1$Diagnosis, ncomp = 1, scale = TRUE)
pred_pls.da <- predict(pls.da, newdata = cv2[,-c(1:4)])

pred_class <- pred_pls.da[["class"]][["max.dist"]] %>% as.data.frame()

conf_mat <- table(factor(cv2$Diagnosis, levels = levels(as.factor(cv2$Diagnosis))), 
      pred_class$comp1)

er_fun(conf_mat)
ber_fun(conf_mat)

#
rownames(b_corr_w1) <- b_corr_w1$Sample

X <- b_corr_w1[,-c(1:4)]
Y <- b_corr_w1$Diagnosis
study <- b_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2

study <- as.factor(study)

base <- mint.plsda(X = X, Y = Y, study = study, ncomp = 1)

base_opt <- perf(base)
base_opt$study.specific.error[['2']]

#
ber_fun(conf_mat)
er_fun(conf_mat)
sens_fun(conf_mat)
spec_fun(conf_mat)

#
#
#

# compare MINT to PLSDA on stacked data
cv1 <- b1_corr_w1 %>% filter(CoreVisit == 'CV1')
cv2 <- b1_corr_w1 %>% filter(CoreVisit == 'CV2')

set.seed(1410)
pls.da <- plsda(X = b_corr_w1[,-c(1:4)], Y = b_corr_w1$Diagnosis, ncomp = 1, scale = TRUE)
pls.da_perf <- perf(pls.da, folds = 8, nrepeat = 100)
pls.da_perf$error.rate

pred_pls.da <- predict(pls.da, newdata = b_corr_w1[,-c(1:4)])

pred_class <- pred_pls.da[["class"]][["centroids.dist"]] %>% as.data.frame() %>% dplyr::select(1)

conf_mat <- table(factor(b1_corr_w1$Diagnosis, levels = levels(as.factor(b1_corr_w1$Diagnosis))), 
                  pred_class$comp1)

ber_fun(conf_mat)
er_fun(conf_mat)
sens_fun(conf_mat)
spec_fun(conf_mat)

plotLoadings(pls.da, ndisplay = 10,
             title = 'Stacked',
             contrib = 'max',
             legend = FALSE)

# PLS-DA + mixed effect model
plsda_pred_vals <- data.frame(Y.pred = pred_pls.da$predict[,,1][,1],
                              out = b_corr_w1$Diagnosis,
                              RAD_ID = b_corr_w1$RAD_ID,
                              CoreVisit = b_corr_w1$CoreVisit) %>%
  mutate(out = ifelse(out == 'Asthma', 1 ,0))

glmMod <- glm(out ~ Y.pred,
              data = plsda_pred_vals,
              family = 'binomial'(link = 'logit'))

summary(glmMod)

mixMod <- glmer(out ~ Y.pred + (1 | CoreVisit),
              data = plsda_pred_vals,
              family = 'binomial'(link = 'logit'))


summary(mixMod)

dev.new()
plsda_pred_vals %>% 
  filter(out == 0) %>%
  ggplot(aes(x = as.factor(CoreVisit), y = Y.pred, group = RAD_ID)) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_line() + 
  geom_point() 

#
#
#

# test MINT on Dataset 1
train <- b1_corr_w1
test_b1_df <- b_corr_w1 %>% filter(CoreVisit == 'CV2') %>%
  mutate(CoreVisit = 'B1_CV')

conc_corr_w1 <- rbind(train, test_b1_df)

rownames(conc_corr_w1) <- conc_corr_w1$Sample

X <- conc_corr_w1[,-c(1:4)]
Y = conc_corr_w1$Diagnosis
study <- conc_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2
study[study == 'B1_CV'] <- 3

study <- as.factor(study)

test_b1 <- which(study == '3')

base <- mint.plsda(X = X[-c(test_b1),],
                   Y = Y[-c(test_b1)],
                   study = droplevels(study[-c(test_b1)]),
                   ncomp = 1)

mint_pred_b1 <- predict(base,
                        newdata = X[test_b1,],
                        dist = 'centroids.dist',
                        study.test = factor(study[test_b1]))

#
dummy_pred_b1_int <- mint_pred_b1[['predict']][,,1]
pred_class_b1_int <- mint_pred_b1[["class"]][["centroids.dist"]] %>% as.data.frame() 

conf_mat <- table(factor(test_b1_df$Diagnosis, levels = levels(as.factor(test_b1_df$Diagnosis))), 
                  pred_class_b1_int$comp1)

conf_mat

er_fun(conf_mat)
ber_fun(conf_mat)
sens_fun(conf_mat)
spec_fun(conf_mat)

roc_b1_mint <- roc(response = test_b1_df$Diagnosis, predictor = dummy_pred_b1_int[,2])
plot(roc_b1_mint, legacy.axes = TRUE)
out_roc_b1_mint <- roc_b1_mint[['sensitivities']] %>% cbind(roc_b1_mint[['specificities']],
                                              roc_b1_mint[['thresholds']]) %>%
  as.data.frame()
colnames(out_roc_b1_mint) <- c('Sens', 'Spec', 'Thresh')
auc(roc_b1_mint)
out_roc_b1_mint %>% filter(Spec < 0.92 & Spec > 0.89) %>%
  summarise(mean = mean(Sens, na.rm = TRUE))
out_roc_b1_mint %>% filter(Sens < 0.92 & Sens > 0.88) %>%
  summarise(mean = mean(Spec))

plot(roc_b1_mint, legacy.axes = TRUE, main = 'Test data (Core Visit 1)') +
  text(label = round(auc(roc_b1_mint_int),2), x = 0.25, y = 0.2, col = 'blue')

View(roc_b1_mint_int[['sensitivities']] %>% cbind(roc_b1_mint_int[['specificities']],
                                        roc_b1_mint_int[['thresholds']]))

#

sp_t2 <- plot_variates(mint_pred_b1$variates, b1_corr_w1) + ggtitle('CV2') + 
  ylab('Component 1') + ylab('') + xlab('') + theme(plot.title = element_text(size = 10))

sp_t1 <- plot_variates(mint_pred_b1$variates, b1_corr_w1) + ggtitle('CV1') + 
  ylab('Component 1') + theme(legend.position = 'none') + xlab('') + theme(plot.title = element_text(size = 10))

sp_t1
sp_t2

sp_t <- arrangeGrob(sp_t1, sp_t2, widths = c(0.375, 0.625),
                    top = textGrob('Global scores plot'))
plot(sp_t)

ggsave('Scores_plots_MINT_test.tiff', sp_t, dpi = 300, unit = 'mm',
       width = 110, height = 65)

# AUROC
# mixOmics function using LOGOCV to calculate ROC
auc <- auroc(base, roc.comp = 1, roc.study = '2')
View(auc)
View(auc[["graph.Comp1"]][["data"]][["Sensitivity"]] %>% as.data.frame())

b_corr_w1 %>% ggplot(aes(x = Diagnosis, y = D_Menthone)) + geom_boxplot()

eb_roc <- roc(b_corr_w1$Diagnosis, b_corr_w1$D_Menthone)
plot.roc(eb_roc, legacy.axes = TRUE)
View(eb_roc)
coords(eb_roc, "best", ret = "threshold")

View(as.data.frame(eb_roc[['thresholds']]))

View(as.data.frame(eb_roc[['sensitivities']]))
View(as.data.frame(eb_roc[['specificities']]))

# using predictions for the training dataset
pred_t_mint <- predict(base, newdata = X, study.test = study, dist = 'max.dist')
dummy_pred <- pred_t_mint[["predict"]][,,1]
rowSums(dummy_pred)

#
roc_t_mint <- roc(response = Y, predictor = dummy_pred[,1])
plot(roc_t_mint, legacy.axes = TRUE)
roc_t_mint[['sensitivities']] %>% cbind(roc_t_mint[['specificities']],
                                        roc_t_mint[['thresholds']])

auc(roc_t_mint)
#
pred_t_class <- pred_t_mint$class$max.dist
conf_mat <- table(factor(Y, levels = levels(as.factor(Y))), 
      pred_t_class[,1])

sens_fun(conf_mat)
spec_fun(conf_mat)

# using predictions for the test dataset
cv1_b1 <- b1_corr_w1 %>% filter(CoreVisit == 'CV1')

study_b1 <- as.factor(rep('1', 48))

pred_b1_mint <- predict(base, newdata = cv1_b1[,-c(1:4)], study.test = study_b1, 
                       dist = 'max.dist')
dummy_pred_b1 <- pred_b1_mint[["predict"]][,,1]
class_pred_b1 <- pred_b1_mint[['class']][['max.dist']]

conf_mat_b1 <- table(factor(cv1_b1$Diagnosis, levels = levels(as.factor(cv1_b1$Diagnosis))),
                     class_pred_b1)

er_fun(conf_mat_b1)
ber_fun(conf_mat_b1)
sens_fun(conf_mat_b1)
spec_fun(conf_mat_b1)

roc_b1_mint <- roc(response = cv1_b1$Diagnosis, predictor = dummy_pred_b1[,2])
plot(roc_b1_mint, legacy.axes = TRUE)
View(roc_b1_mint[['sensitivities']] %>% cbind(roc_b1_mint[['specificities']],
                                        roc_b1_mint[['thresholds']]))
auc(roc_b1_mint)

#
#
#

# INTEGRATE EVERYTHING
conc_corr_w1 <- rbind(b_corr_w1, 
                      b1_corr_w1 %>%
                        mutate(CoreVisit = ifelse(CoreVisit == 'CV1', 'B1_CV1', 'B1_CV2')))

X <- conc_corr_w1[,-c(1:4)]
Y <- conc_corr_w1$Diagnosis
study <- conc_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2
study[study == 'B1_CV1'] <- 3
study[study == 'B1_CV2'] <- 4

study <- as.factor(study)

full <- mint.plsda(X = X, Y = Y, study = study, ncomp = 1, scale = TRUE)
plotIndiv(full)
dev.new()
plotLoadings(full, ndisplay = 15, study = 'global')

# tuning number of variables with LOGOCV
base_tune <- tune(X = X, Y = Y, study = as.factor(study),
                  ncomp = 1,
                  test.keepX = seq(1,50,1),
                  method = 'mint.splsda',
                  measure = 'BER',
                  dist = 'centroids.dist')

plot(base_tune)
base_tune$choice.keepX

#
loads_glob <- base$loadings$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>% 
  left_join(reg_valid_join %>% dplyr::select(comp, sign)) %>%
  left_join(endo_exo %>% dplyr::select(comp, FC_filter)) %>%
  distinct()

loads_glob1 <- loads_glob %>% arrange(desc(abs(comp1))) %>%
  slice(1:46) %>% mutate(comp1_class = ifelse(comp1 > 0, 'pos', 'neg')) 

loads_glob1_sum <- loads_glob1 %>% group_by(comp1_class, FC_filter) %>% summarise(n = n()) %>%
  mutate(ratio = ifelse(comp1_class == 'pos', n/table(loads_glob1$comp1_class)[2], 
                                                      n/table(loads_glob1$comp1_class)[1]),
         model = 'global')

#
loads_glob_b1 <- base1$loadings$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>% 
  left_join(reg_valid_join %>% dplyr::select(comp, sign)) %>%
  left_join(endo_exo %>% dplyr::select(comp, FC_filter)) %>%
  distinct()


loads_glob_b11 <- loads_glob_b1 %>% arrange(desc(abs(comp1))) %>%
  slice(1:42) %>% mutate(comp1_class = ifelse(comp1 > 0, 'pos', 'neg')) 

loads_glob_b11_sum <- loads_glob_b11 %>% group_by(comp1_class, FC_filter) %>% summarise(n = n()) %>%
  mutate(ratio = ifelse(comp1_class == 'pos', n/table(loads_glob_b11$comp1_class)[2], 
                        n/table(loads_glob_b11$comp1_class)[1]),
         model = 'global')

#
#

loads_part <- base$loadings.partial$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  rename(CV1_comp1 = comp1,
         CV2_comp1 = comp1.1) 

loads_part1_cv1 <- loads_part %>% 
  left_join(reg_valid_join %>% dplyr::select(comp, sign)) %>%
  left_join(endo_exo %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>% arrange(desc(abs(CV1_comp1))) %>%
  slice(1:46) %>% 
  mutate(comp1_class = ifelse(CV1_comp1 > 0, 'pos', 'neg')) 
  
loads_part1_cv1_sum <- loads_part_b11_cv1 %>% group_by(comp1_class, FC_filter) %>% summarise(n = n()) %>%
  mutate(ratio = ifelse(comp1_class == 'pos', n/table(loads_part_b11_cv1$comp1_class)[2], 
                        n/table(loads_part_b11_cv1$comp1_class)[1]),
         model = 'partial_cv1') 
  
# visualise VOC classes among most important loadings
loads_glob1_sum %>% rbind(loads_part_b11_cv1_sum, loads_part_b11_cv2_sum) %>%
  mutate(FC_filter = ifelse(is.na(FC_filter), 'Unknown', FC_filter)) %>%
  ggplot(aes(x = comp1_class, y = ratio, fill = FC_filter)) +
  facet_wrap(~model) +
  #geom_col(position = position_dodge(preserve = 'single')) +
  geom_col() +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  xlab('') +
  scale_x_discrete(labels = c('neg' = 'comp1 < 0',
                              'pos' = 'comp1 > 0')) +
  ggtitle('Classification of component 1 loadings (Dataset 1)') +
  ylab('')

ggsave('Bar_plot_class_loadings_B1.tiff', unit = 'mm', dpi = 300, width = 170, height = 80)

# variable importance
loads_vip_b <- vip(base) %>% as.data.frame()

nrow(loads_vip %>% filter(comp1 > 1))
hist(loads_vip$comp1)

colnames(vocs_freq)[1] <- 'comp'
vocs_perf_base <- loads_vip_b %>% mutate(comp = rownames(.)) %>%
  left_join(vocs_freq) %>%
  rename(vip = comp1)

plot(vocs_perf_base$vip, vocs_perf_base$stab)

top_stab <- vocs_perf_base %>% filter(stab > 0.7)

#
#

Lglob1 <- loads_glob %>% 
  mutate(sign = ifelse(is.na(sign), 'no', sign)) %>%
  #left_join(icc_b12w %>% dplyr::select(comp, ICC_B2, CI_lwr_B2)) %>%
  #mutate(ICC_B2 = ifelse(ICC_B2 < 0, 0, ICC_B2),
   #      sign_icc = ifelse(CI_lwr_B2 < 0, 'no', 'yes'),
    #     comp = str_trunc(comp, 33)) %>%
  rename(BG_corr = sign) %>%
  arrange(desc(abs(comp1))) %>% 
  slice(1:50) %>% 
  ggplot(aes(x = comp1, y = fct_inorder(as.factor(comp)), 
             #fill = ICC_B1,
             fill = FC_filter,
             colour = BG_corr
             #colour = sign_icc
             )) + 
  geom_col() +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  #scale_fill_viridis_b(limits = c(0, 0.7)) +
  theme(axis.title.y = element_blank())

Lglob1

Lglob1 + ggtitle('Global loadings')

ggsave('MINT_global_loadings_VS.tiff', dpi = 300, unit = 'mm', width = 135, height = 155)

Lglob2

Lglob <- arrangeGrob(Lglob1, Lglob2, ncol = 2, 
                     top = textGrob('Global loadings of multi-group PLS-DA'),
                     widths = c(0.55, 0.45))
plot(Lglob)
ggsave('MINT_global_loadings_B1.tiff', Lglob, unit = 'mm', dpi = 300, width = 240, height = 80)

#
dev.new()

Lpart11 <- loads_part %>% 
  left_join(reg_valid_join %>% dplyr::select(comp, sign)) %>%
  left_join(endo_exo_2 %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>%
  left_join(icc_b12w %>% dplyr::select(comp, ICC_B2, CI_lwr_B2)) %>%
  mutate(ICC_B2 = ifelse(ICC_B2 < 0, 0, ICC_B2),
         sign_icc = ifelse(CI_lwr_B2 < 0, 'no', 'yes'),
         comp = str_trunc(comp, 33)) %>%
  rename(BG_corr = sign) %>%
  arrange(desc(abs(CV1_comp1))) %>% 
  slice(1:46) %>% 
  ggplot(aes(x = CV1_comp1, y = fct_inorder(as.factor(comp)), 
             #fill = ICC_B1,
             fill = FC_filter,
             colour = BG_corr
             #colour = sign_icc
  )) + 
  geom_col() +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  #scale_fill_viridis_b(limits = c(0, 0.7)) +
  theme(axis.title.y = element_blank(),
        legend.position = 'none') +
  xlab('comp1')

Lpart11 + ggtitle('Partial loadings (CV1)') + theme(plot.title = element_text(hjust = 0.8))

ggsave('MINT_partial_loadings_CV1_VS.tiff', dpi = 300, unit = 'mm', width = 90, height = 155)

Lpart2 <- arrangeGrob(Lpart21, Lpart22, ncol = 2, 
                     top = textGrob('Partial loadings (Core Visit 2)'))
plot(Lpart2)
ggsave('MINT_partial_CV1_loadings_B1.tiff', Lpart2, unit = 'mm', dpi = 300, 
       width = 155, height = 80)


#
cv1_topL <- loads_part %>% arrange(desc(abs(CV1_comp1))) %>% 
  slice(1:20)
cv2_topL <- loads_part %>% arrange(desc(abs(CV2_comp1))) %>% 
  slice(1:20)

intersect(cv1_topL$comp, cv2_topL$comp)


dev.new()
cim(base, comp = 1,
    row.sideColors = color.mixo(as.numeric(as.factor(Y))), legend = TRUE)

auroc(base, roc.comp = 1, print = FALSE)

#
#
#

# PLS-DA
cv1 <- b_corr_w1 %>% filter(CoreVisit %in% c('CV1'))
cv2 <- b_corr_w1 %>% filter(CoreVisit %in% c('CV2'))

intersect(cv1$RAD_ID, cv2$RAD_ID)

# APPROACH 1 - two models
assay_c1 <- cv1[,-c(1:4)]
cv1

assay_c2 <- cv2[,-c(1:4)]
cv2

Y <- as.factor(cv1$Diagnosis)
Y2 <- as.factor(cv2$Diagnosis)

X <- assay_c1 
X2 <- assay_c2

sub_all <- plsda(X, Y, scale = TRUE, ncomp = 5)
sub_all2 <- plsda(X2, Y2, scale = TRUE, ncomp = 5)


plotIndiv(sub_all, legend = TRUE,
          pch = 1, ellipse = TRUE,
          title = 'CV1')

plotIndiv(sub_all2, legend = TRUE,
          pch = 1, ellipse = TRUE,
          title = 'CV2')


vars <- plotVar(sub_all, cutoff = 0.6)
vars2 <- plotVar(sub_all2, cutoff = 0.6)

plotLoadings(sub_all, ndisplay = 10, method = 'median', contrib = 'max',
             title = 'CV1')
plotLoadings(sub_all2, ndisplay = 10, method = 'median', contrib = 'max',
             title = 'CV2')

# tuning number of components
tune_comps <- perf(sub_all, validation = "Mfold", 
                   folds = 8, nrepeat = 100)

plot(tune_comps, sd = TRUE,
     legend.position = "horizontal")

tune_comps$choice.ncomp
# 1 comp, max dist

tune_comps2 <- perf(sub_all2, validation = "Mfold", 
                   folds = 6, nrepeat = 100)

plot(tune_comps2, sd = TRUE,
     legend.position = "horizontal")

tune_comps$choice.ncomp
# 1 comp, max dist

#

# tuning number of variables 
keepX <- c(1:20, seq(20, 30, 5))

set.seed(35)
tune_sub_all <- tune.splsda(X, Y, 
                            ncomp = 2,
                            validation = 'Mfold',
                            folds = 8, 
                            dist = 'max.dist',
                            test.keepX = keepX,
                            nrepeat = 100,
                            measure = "BER")

head(tune_sub_all$error.rate)
plot(tune_sub_all, sd = TRUE)

tune_sub_all$choice.ncomp$ncomp
tune_sub_all$choice.keepX
# 1 comp, 14 vars

#
tune_sub_all2 <- tune.splsda(X2, Y2, 
                            ncomp = 2,
                            validation = 'Mfold',
                            folds = 6, 
                            dist = 'max.dist',
                            test.keepX = keepX,
                            nrepeat = 100,
                            measure = "BER")

head(tune_sub_all2$error.rate)
plot(tune_sub_all2, sd = TRUE)

tune_sub_all2$choice.ncomp$ncomp
tune_sub_all2$choice.keepX
# 1 comp, 11 vars

# sparse model
sub_all_s <- splsda(X, Y, ncomp = 1,
                    keepX = 14, 
                    scale = TRUE)

sub_all_s2 <- splsda(X2, Y2, ncomp = 1,
                     keepX = 11,
                     scale = TRUE)

#

scores <- data.frame(C1 = sub_all_s$variates$X,
                     Y = Y)

write.csv(scores, 'Splsda_scores.csv')

c1_scores <- scores %>% ggplot(aes(x = Y, y = comp1, group = Y, fill = Y)) + 
  geom_boxplot(outliers = FALSE) + geom_jitter() +
  theme_bw() +
  ggtitle('CV1') +
  ylab('Comp1 (10.5%)') +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none') +
  scale_fill_manual(name = 'Diagnosis',
                    values = c('Asthma' = 'goldenrod1',
                               'Not Asthma' = 'olivedrab3')) +
  geom_hline(yintercept = 0, linetype = 'dashed')

c1_scores


scores2 <- data.frame(C1 = sub_all_s2$variates$X,
                    Y = Y2)

c2_scores <- scores2 %>% ggplot(aes(x = Y, y = comp1, group = Y, fill = Y)) + 
  geom_boxplot(outliers = FALSE) + geom_jitter() +
  theme_bw() +
  ggtitle('CV2') +
  ylab('Comp1 (9.3%)') +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Diagnosis',
                    values = c('Asthma' = 'goldenrod1',
                               'Not Asthma' = 'olivedrab3')) +
  geom_hline(yintercept = 0, linetype = 'dashed')

c2_scores

c_scores <- arrangeGrob(c1_scores, c2_scores, nrow = 1,
                        widths = c(0.41, 0.59),
                        top = textGrob('sPLS-DA scores plots'))
plot(c_scores)
ggsave('SPLSDA_scores_plots_CV1_2.tiff', c_scores, unit = 'mm', dpi = 300, 
       width = 150, height = 70)


#

sub_all_s$prop_expl_var
# 0.0714
sub_all_s2$prop_expl_var
# 0.0574

#
loadings <- as.data.frame(sub_all_s$loadings$X) %>%
  mutate(Var1 = rownames(.)) %>%
  filter(comp1 > 0 | comp1 < 0)

loadings2 <- as.data.frame(sub_all_s2$loadings$X) %>%
  mutate(Var1 = rownames(.)) %>%
  filter(comp1 > 0 | comp1 < 0)


# 5 VOCs selected in both models
plotLoadings(sub_all_s, ndisplay = 15,
             title = 'CV1')

plotLoadings(sub_all_s2, ndisplay = 15,
             title = 'CV2')

#



# stability
stab_sub <- perf(sub_all_s, 
                 ncomp = 1,
                 validation = 'Mfold',
                 folds = 8, 
                 dist = 'max.dist',
                 nrepeat = 100)

loadings_stab <- loadings %>% left_join(as.data.frame(stab_sub$features$stable[[1]]))

#
stab_sub2 <- perf(sub_all_s2, 
                 ncomp = 1,
                 validation = 'Mfold',
                 folds = 6, 
                 dist = 'max.dist',
                 nrepeat = 100)

loadings_stab2 <- loadings2 %>% left_join(as.data.frame(stab_sub2$features$stable[[1]]))

write.csv(loadings_stab, 'Splsda_loadings_stab.csv')

# plot loadings and stability
library(forcats)

cv1_L <- loadings_stab %>% 
  arrange(abs(comp1)) %>%
  ggplot(aes(x = fct_inorder(as.factor(Var1)), y = comp1, fill = Freq)) + 
  geom_col() +
  coord_flip() +
  scale_fill_continuous(type = 'viridis', name = 'Stability',
                        limits = c(0.4,1)) +
  ylab('Comp1') +
  ggtitle('CV1') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = 'bottom')

cv1_L

cv2_L <- loadings_stab2 %>% 
  arrange(abs(comp1)) %>%
  ggplot(aes(x = fct_inorder(as.factor(Var1)), y = comp1, fill = Freq)) + 
  geom_col() +
  coord_flip() +
  scale_fill_continuous(type = 'viridis', name = 'Stability',
                        limits = c(0.4,1)) +
  ylab('Comp1') +
  ggtitle('CV2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = 'bottom',
        plot.margin = margin(l = 155,
                             r = 5))

cv2_L

dev.new()
c_loadings <- arrangeGrob(cv1_L, cv2_L, ncol = 1,
                        top = textGrob('sPLS-DA loadings plots'),
                        heights = c(0.6, 0.4))
plot(c_loadings)
ggsave('SPLSDA_loadings_plots_CV1_2.tiff', c_loadings, unit = 'mm', dpi = 300, 
       height = 150, width = 150)


# loadings agreement
cons_voc <- intersect(loadings$Var1, loadings2$Var1)

View(loadings_stab %>% filter(Var1 %in% cons_voc))
View(loadings_stab2 %>% filter(Var1 %in% cons_voc))

v0 <- cbind(b1_imp_sum_corr_w1[,(1:4)], assay) %>% 
  ggplot(aes(x = CoreVisit,
             y = Pyridine, 
             fill = Diagnosis)) +
  geom_boxplot(outliers = FALSE) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Diagnosis',
                    values = c('Asthma' = 'goldenrod1',
                               'Not Asthma' = 'olivedrab3')) +
  ylab('Peak area') +
  geom_point(position = position_jitterdodge()) +
  ggtitle('Pyridine')
  

v0

selected_vars <- arrangeGrob(v1, v2, v3, v4, v5, v0, nrow = 2, ncol = 3)
plot(selected_vars)

ggsave('sPLSDA_selected_VOCs_agreement.tiff', selected_vars, units = 'mm',
       dpi = 300, width = 340, height = 150)

#

# compare performance of pls-da dn spls-da models for CV1
perf_sub_all <- perf(sub_all, folds = 8, validation = 'Mfold',
                     dist = 'max.dist', nrepeat = 100)

perf_sub_all$error.rate
perf_sub_all$error.rate.class

perf_sub_all_s <- perf(sub_all_s, folds = 8, validation = 'Mfold',
                      dist = 'max.dist', nrepeat = 100)

perf_sub_all_s$error.rate
perf_sub_all_s$error.rate.class

# compare performance of spls-da from CV1 and CV2
perf_sub_all_s2 <- perf(sub_all_s2, folds = 6, validation = 'Mfold',
                       dist = 'max.dist', nrepeat = 100)

perf_sub_all_s2$error.rate
perf_sub_all_s2$error.rate.class


#
#
#

# APPROACH 2
# predict CV2
predict_cv2 <- predict(sub_all_s, X2, dist = 'max.dist')

predict_c1 <- predict_cv2$class$max.dist
table(factor(predict_c1, levels = levels(Y)), Y2)

sens <- 21/(21+4)
spec <- 8/(9+8)

auc_cv1 <- auroc(sub_all_s)



#
#
#

# APPROACH 3
# remove patient effect

# mixOmics multilevel decomposition - IF PERFORMED ON CV1, CV2 SUBSET, PLSDA DOES NOT CONVERGE
# load summarised BG corrected dataset 
b1_sub_corr_w <- read.csv('RADiCA_breath_BGcorrected.csv')[,-1]

# keep only sample with min 2 observations present
c_reps <- b1_imp_sum_corr_w1 %>% group_by(RAD_ID) %>% summarise(n = n()) %>%
  filter(n >= 2)

b1_multi <- b1_imp_sum_corr_w1 %>% filter(RAD_ID %in% c_reps$RAD_ID)

#
# variance stabilise
assay_multi <- b1_sub_corr_w_multi[,-c(1:6)]
rownames(assay_multi) <- b1_sub_corr_w_multi$Sample

vsn <- vsn::vsnMatrix(as.matrix(assay_multi))
vsnout <- vsn@hx
assay_multi <- as.data.frame(vsnout)

# remove patient effect through multilevel decomposition
design <- data.frame(cv = b1_sub_corr_w$CoreVisit) %>%
  mutate(cv = ifelse(cv == 'CV1',1,2))

Xw <- withinVariation(X = assay, design = design)

# APPROACH 4 - one model
Xsum <- b1_multi[,-c(1:4)]
Ysum <- b1_multi$Diagnosis

mod_sum <- plsda(X = Xsum, Y = Ysum, ncomp = 2, scale = TRUE)
plotIndiv(mod_sum, pch = 1)
plotLoadings(mod_sum, ndisplay = 15)

keepX <- c(1:20, seq(20, 50, 5))

mod_sum_s <- splsda(X = Xsum, Y = Ysum, ncomp = 1, scale = TRUE,
                    keepX = 10)

plotIndiv(mod_sum_s, pch = 1)

as.data.frame(mod_sum_s$variates) %>% mutate(Sample = rownames(.)) %>%
  left_join(b1_multi %>% dplyr::select(Sample, Diagnosis)) %>%
  ggplot(aes(x = Diagnosis, y = comp1)) + geom_boxplot() +
  geom_jitter()

plotLoadings(mod_sum_s)

stab_sum <- perf(mod_sum_s, 
                  ncomp = 1,
                  validation = 'Mfold',
                  folds = 10, 
                  dist = 'max.dist',
                  nrepeat = 100)

loadings_sum <- as.data.frame(mod_sum_s$loadings) %>%
  mutate(comp = rownames(.)) %>%
  left_join(as.data.frame(stab_sum$features$stable[[1]]) %>%
              rename(comp = Var1))

loadings_sum %>% 
  arrange(abs(comp1)) %>%
  filter(comp1 != 0) %>%
  ggplot(aes(x = fct_inorder(as.factor(comp)), y = comp1, fill = Freq)) + 
  geom_col() +
  coord_flip() +
  scale_fill_continuous(type = 'viridis', name = 'Stability') +
  ylab('Comp1') +
  ggtitle('CV1 + CV2 (complete obs)') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = 'bottom')

# test interactionof diagnosis and cv on selected features
markers <- loadings_sum %>% filter(comp1 != 0) %>% filter(Freq > 0.5)

unitest <- function(df, voc){
  subset <- df %>%
    dplyr::select(c(voc, Diagnosis, CoreVisit, RAD_ID)) %>%
    rename(comp = voc)
  model <- lmer(comp ~ Diagnosis*CoreVisit + (1 | RAD_ID),
                data = subset)
  coefs <- cbind(as.data.frame(coef(summary(model))),
                 as.data.frame(confint(model)[3:6,]))
  print(coefs)}


unitest(b1_multi, 'X1_Propanol._2_ethoxy_')

box_comp <- function(voc) {
  b1_multi %>% rename(comp = {{voc}}) %>%
  ggplot(aes(x = CoreVisit, y = comp, fill = Diagnosis)) +
  geom_boxplot(outliers = FALSE) +
  ggtitle(voc) +
  theme_bw() +
  scale_fill_brewer(palette = 'Set2')
}

box_comp('Ethyl_Acetate')

marker_boxplots <- lapply(markers$comp, box_comp)
marrangeGrob(marker_boxplots, nrow = 3, ncol = 3)

b1_multi %>% rename(comp = voc) %>%
  ggplot(aes(x = CoreVisit, y = comp, colour = Diagnosis)) +
  geom_point() + geom_line(aes(group = RAD_ID))

voc_sum <- bind_rows(lapply(colnames(b1_multi)[5:148], sum_voc <- function(voc) {
  summary <- b1_multi %>% rename(comp = voc) %>%
  group_by(Diagnosis, CoreVisit) %>% 
  summarise(mean = mean(comp),
            sd = sd(comp)) %>%
  as.data.frame() %>%
  mutate(voc = rep(voc, 4))}))

voc_sum %>% ggplot(aes(x = CoreVisit, y = sd, colour = Diagnosis)) + 
  geom_point() + geom_line(aes(group = interaction(Diagnosis, voc)))

voc_sum %>% ggplot(aes(x = CoreVisit, y = sd, fill = Diagnosis)) + 
  geom_boxplot()
  
 
#
#
#

# vegan RDA 
library(vegan)

data <- b1_corr_w1

rownames(data) <- data$Sample
assay <- data[,-c(1:4)]

rda <- rda(assay ~ Diagnosis + CoreVisit,
           data = data,
           scale = TRUE)


rda

vif.cca(rda)
RsquareAdj(rda)
envfit(rda, data[,c(2,3)])

inertia <- summary(rda)[[7]]
rda_var_expl <- round(eigenvals(rda)[1:2]/inertia*100, 2)

rda_scores <- scores(rda, choices = 1:3, display = 'sites') %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(data[,c(1:4)])

rda_samples_plot1 <- rda_scores %>% 
  ggplot(aes(x = RDA1, y = RDA2, colour = Diagnosis, shape = CoreVisit)) + 
  geom_point(size = 2) + theme_bw(base_size = 10) +
  xlab(paste('RDA1:', rda_var_expl[1], '%')) +
  ylab(paste('RDA2:', rda_var_expl[2], '%')) +       
  scale_colour_manual(values = c('Asthma' = 'dodgerblue3',
                                 'Not Asthma' = 'darkorange1')) +
  scale_shape_manual(values = c('CV1' = 1,
                                'CV2' = 2)) +
  ggtitle('Dataset 1') +
  theme(plot.title = element_text(hjust = 0.5))

rdap <- arrangeGrob(rda_samples_plot1, rda_samples_plot, ncol = 1,
            top = textGrob('Redundancy Analysis sample scores plot'))

plot(rdap)

ggsave('RDA_sample_scores_plots.tiff', rdap, units = 'mm', dpi = 300, width = 110, height = 155)

#
rda_loadings_voc <- scores(rda, choices = 1:3, display = 'species') %>%
  as.data.frame() %>%
  mutate(comp = rownames(.))

b1rda1 <- rda_loadings_voc %>%
  left_join(reg_valid_join %>% dplyr::select(comp, FC_filter, sign)) %>%
  distinct() %>%
  mutate(FC_filter = ifelse(is.na(FC_filter), 'Unknown', FC_filter)) %>%
  rename(BGcorr = sign) %>%
  arrange(desc(abs(RDA1))) %>%
  slice(1:20) %>%
  mutate(comp = str_trunc(comp, 34)) %>%
  ggplot(aes(x = RDA1, y = fct_inorder(as.factor(comp)),
             fill = FC_filter, colour = BGcorr)) + 
  geom_col() +
  theme_bw(base_size = 9) +
  #ggtitle('Diagnosis') +
  ggtitle('Core Visit') +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank()#,
        #legend.position = 'none' 
        ) +
  scale_fill_brewer(palette = 'Set2')

b1rda1

b1rdap <- arrangeGrob(b1rda1, b1rda2, ncol = 2, widths = c(0.4, 0.6),
                      top = textGrob('RDA loadings (Dataset 1)')) 
plot(b1rdap)

ggsave('RDA_loadings_B1.tiff', b1rdap, unit = 'mm', dpi = 300, width = 200, 
       height = 80)

#
b2rdap <- arrangeGrob(b2rda1, b2rda2, ncol = 2, widths = c(0.4, 0.6),
                      top = textGrob('RDA loadings (Dataset 2)')) 
plot(b2rdap)

ggsave('RDA_loadings_B2.tiff', b2rdap, unit = 'mm', dpi = 300, width = 200, 
       height = 80)


#

rda1_scores <- scores(rda1, choices = 1, display = 'sites') %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(b1_sub_corr_w[,c(1:6)])

rda1_plot <- rda1_scores %>% ggplot(aes(x = Diagnosis, y = RDA1, fill = Diagnosis)) + 
  geom_boxplot() +
  geom_jitter() + theme_bw() +
  ggtitle('RDA Diagnosis + Conditional(CoreVisit)')

#
rda0_scores <- scores(rda0, choices = 1, display = 'sites') %>%
  as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(b1_sub_corr_w[,c(1:6)])

rda0_plot <- rda0_scores %>% ggplot(aes(x = Diagnosis, y = RDA1, fill = Diagnosis)) + 
  geom_boxplot() +
  geom_jitter() + theme_bw() +
  ggtitle('RDA Diagnosis')

rda_plots <- arrangeGrob(rda0_plot, rda_plot, rda1_plot, nrow = 1, ncol = 3)
dev.new()
plot(rda_plots)

#

plot(rda)
plot(rda, display = 'wa', type = 'points')
ordispider(rda, col = 'red')
text(rda, display = 'cn', col = 'blue')

# example from package authors
data(dune, dune.env, package = "vegan")
ord <- cca(dune ~ Moisture, data = dune.env)
ord
plot(ord, display = "wa", type = "points")
ordispider(ord, col = "red")
text(ord, display = "cn", col = "blue")

#

# PERMANOVA
perm <- adonis2(assay ~ Diagnosis + RAD_ID,
                permutations = 99,
                data = b1_imp_sum_corr_w1,
                method = 'euclidean')

perm



#
#
#

# Evaluating outcomes
# PCA
pca_multi <- pca(Xw, scale = TRUE, center = TRUE, ncomp = 4)
plotIndiv(pca_multi,
          pch = 1,
          group = b1_sub_corr_w_multi1$Diagnosis)



#

# PLS-DA
Y_multi <- as.factor(b1_sub_corr_w$Diagnosis)

X_multi <- Xw

sub_multi <- plsda(X_multi, Y_multi, scale = TRUE, ncomp = 5)

plotIndiv(sub_multi, legend = TRUE,
          pch = 1, ellipse = TRUE,
          title = 'CV1 + CV2')

plotVar(sub_multi, cutoff = 0.7)

plotLoadings(sub_multi, ndisplay = 10, method = 'median', contrib = 'max',
             title = 'CV1 + CV2')

# tuning number of components
tune_comps_multi <- perf(sub_multi, validation = "Mfold", 
                   folds = 10, nrepeat = 100)

plot(tune_comps_multi, sd = TRUE,
     legend.position = "horizontal")

tune_comps$choice.ncomp
# 1 comp, centroid distance?

# tuning number of variables 
keepX <- c(1:15, seq(20, 30, 5))

set.seed(507)
tune_sub_multi <- tune.splsda(X_multi, Y_multi, 
                            ncomp = 3,
                            validation = 'Mfold',
                            folds = 10, 
                            dist = 'max.dist',
                            test.keepX = keepX,
                            nrepeat = 100,
                            measure = "BER")

head(tune_sub_multi$error.rate)
plot(tune_sub_multi, sd = TRUE)

tune_sub_all$choice.ncomp$ncomp
tune_sub_all$choice.keepX

# 1 comp, 30 vars

sub_multi_s <- splsda(X_multi, Y_multi, scale = TRUE,
                      keepX = 30, ncomp = 1)

sub_multi$prop_expl_var

scores_multi <- data.frame(C1 = sub_multi_s$variates$X,
                      Y = Y_multi)

multi_scores <- scores_multi %>% ggplot(aes(x = Y, y = comp1, group = Y, fill = Y)) + 
  geom_boxplot(outliers = FALSE) + geom_jitter() +
  theme_bw() +
  ggtitle('CV1 + CV2') +
  ylab('Comp1 (11.5%)') +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Diagnosis',
                    values = c('Asthma' = 'goldenrod1',
                               'Not Asthma' = 'olivedrab3')) +
  geom_hline(yintercept = 0, linetype = 'dashed')

multi_scores

ggsave('SPLSDA_scores_plots_CV1&2.tiff', unit = 'mm', dpi = 300, 
       width = 100, height = 70)

#
plotLoadings(sub_multi_s, ndisplay = 10,
             method = 'median', contrib = 'max',
             title = 'CV1 + CV2')

#
loadings_multi <- as.data.frame(sub_multi_s$loadings$X) %>%
  mutate(Var1 = rownames(.)) %>%
  filter(comp1 > 0 | comp1 < 0)

# stability
stab_multi <- perf(sub_multi_s, 
                 ncomp = 1,
                 validation = 'Mfold',
                 folds = 10, 
                 dist = 'max.dist',
                 nrepeat = 100)

loadings_multi_stab <- loadings_multi %>% left_join(as.data.frame(stab_multi$features$stable[[1]]))

multi_L <- loadings_multi_stab %>% 
  arrange(abs(comp1)) %>%
  filter(Freq > 0.7) %>%
  ggplot(aes(x = fct_inorder(as.factor(Var1)), y = comp1, fill = Freq)) + 
  geom_col() +
  coord_flip() +
  scale_fill_continuous(type = 'viridis', name = 'Stability',
                        limits = c(0.4,1)) +
  ylab('Comp1') +
  ggtitle('CV1 + CV2') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        legend.position = 'bottom')

multi_L

ggsave('SPLSDA_loadings_plots_CV1&2.tiff', unit = 'mm', dpi = 300, 
       height = 120, width = 150)

# compare performance of spls-da from CV1 and CV2
perf_sub_multi <- perf(sub_multi_s, folds = 10, validation = 'Mfold',
                        dist = 'max.dist', nrepeat = 100)

perf_sub_multi$error.rate
perf_sub_multi$error.rate.class

#
#
#

# RANDOM FORESTS EXTENSION
library(devtools)
install_github("bcjaeger/bimm")
library(bimm)

input <- b_corr_w1 %>% dplyr::select(!c(Sample, CoreVisit))

input1 <- input %>% mutate(Diagnosis = ifelse(Diagnosis  == 'Asthma', 1, 0))

#input1 <- input1 %>% dplyr::select(!RAD_ID)

#dup_id <- as.data.frame(table(input$RAD_ID)) %>% filter(Freq > 1)

#input1 <- input %>% filter(RAD_ID %in% dup_id$Var1)

forest <- bimm_fit(data = input1,
                   formula = Diagnosis ~ . + (1 | RAD_ID),
                   n_iteration = 2,
                   verbose = TRUE)
summary(forest)

View(as.data.frame(forest[["model_ml"]][["variable.importance"]]))

input2 <- b1_corr_w1 %>% dplyr::select(!c(Sample, CoreVisit))
input2 <- input2 %>% mutate(Diagnosis = ifelse(Diagnosis  == 'Asthma', 1, 0))
input2$Diagnosis <- as.factor(input2$Diagnosis)

preds_train <- bimm_predict(forest, new_data = input2, type = 'new_sub') %>% as.data.frame()
summary(preds_train)

colnames(preds_train) <- 'Pred'
preds_train <- preds_train %>% mutate(Pred = ifelse(Pred > 0.5, 1 , 0)) 
preds_train$Pred <- as.factor(preds_train$Pred)

conf_mat <- table(factor(input2$Diagnosis, levels = levels(input2$Diagnosis)), preds_train$Pred)
conf_mat

ber_fun(conf_mat)

data <- input1

cv_perf <- 
  replicate(n = 10, expr = {
    p <- n_distinct(data$RAD_ID)
    x <- sample(1:p)
    names(x) <- unique(data$RAD_ID)
    levels <- 5
    folds <- split(x, x%%levels)
    
    BER <- vector(mode = "numeric")
    
    for (m in as.character(0:(levels - 1))) {
      test_id <- folds[[m]]
      test_data <- data %>% filter(RAD_ID %in% names(test_id)) 
      train_data <- data %>% filter(RAD_ID %ni% names(test_id))
      forest <- bimm_fit(data = train_data,
                         formula = Diagnosis ~ . + (1 | RAD_ID),
                         n_iteration = 1,
                         verbose = TRUE)
      preds_train <- bimm_predict(forest, new_data = test_data, type = 'new_sub') %>% as.data.frame()
      colnames(preds_train) <- 'Pred'
      preds_train <- preds_train %>% mutate(Pred = ifelse(Pred > 0.5, 1 , 0)) 
      preds_train$Pred <- as.factor(preds_train$Pred)
      test_data$Diagnosis <- as.factor(test_data$Diagnosis)
      conf_mat <- table(factor(test_data$Diagnosis, levels = levels(test_data$Diagnosis)), preds_train$Pred)
      BER[m] <- if (ncol(conf_mat) == 2 & nrow(conf_mat) == 2) {
        ber_fun(conf_mat)
      } else {
        NA
      }
      print(BER)
    }
    cv_ber <- mean(BER, na.rm = TRUE)
  })


ber <- mean(cv_perf, na.rm = TRUE)
ber

# RF with ranger without random effect
library(ranger)
base_forest <- ranger(Diagnosis ~., 
                      data = input1 %>% dplyr::select(!RAD_ID),
                      probability = TRUE)

pred_test <- predict(base_forest, data = input2 %>% dplyr::select(!c(RAD_ID)))
pred_test_df <- pred_test$predictions %>% as.data.frame() 
colnames(pred_test_df) <- c('Not_Asthma', 'Asthma')
pred_test_df <- pred_test_df %>% mutate(Diagnosis = ifelse(Asthma > 0.5, 1, 0))

conf_mat <- table(factor(input2$Diagnosis, levels = levels(input2$Diagnosis)), pred_test_df$Diagnosis)
conf_mat

ber_fun(conf_mat)

#
#
#

# examine reproducibility of VOC levels across time points
data <- b_corr_w1
vocs <- colnames(b_corr_w1)[-c(1:4)]

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

write.csv(icc_b1, 'ICC_CoreVisits_B1.csv')
write.csv(icc_b2 , 'ICC_CoreVisits_B2.csv')

icc_b1 <- read.csv('ICC_CoreVisits_B1.csv')[,-1]
icc_b2 <- read.csv('ICC_CoreVisits_B2.csv')[,-1]

icc_b1 <- icc_b1 %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH')) 

icc_b2 <- icc_b2 %>%
  mutate(adj.p.value = p.adjust(p.value, method = 'BH'))

hist(icc_b2$ICC, breaks = 8)

icc_b12 <- icc_b1 %>% mutate(dataset = 'B1') %>%
  rbind(icc_b2 %>% mutate(dataset = 'B2'))

icc_b12w <- icc_b12 %>% 
  pivot_wider(names_from = dataset, 
              values_from = c(ICC, p.value, adj.p.value, CI_lwr, CI_upr)) %>%
  mutate(p.value = ifelse(p.value_B1 < 0.05 & p.value_B2 < 0.05, 'both', 
                          ifelse(p.value_B1 > 0.05 & p.value_B2 > 0.05, 'neither' ,' one'))) %>%
  mutate(adj.p.value = ifelse(adj.p.value_B1 < 0.05 & adj.p.value_B2 < 0.05, 'both', 
                          ifelse(adj.p.value_B1 > 0.05 & adj.p.value_B2 > 0.05, 'neither' ,' one'))) %>%
  mutate(CI_lwr = ifelse(CI_lwr_B1 > 0 & CI_lwr_B2 > 0, 'both', 
                              ifelse(CI_lwr_B1 < 0 & CI_lwr_B2 < 0, 'neither' ,' one')))


icc_p <- icc_b12w %>% 
  dplyr::select(ICC_B1, ICC_B2, p.value, adj.p.value, CI_lwr) %>%
  pivot_longer(cols = c(p.value, adj.p.value, CI_lwr), 
               names_to = 'para', values_to = 'outcome') 

icc_p$para <- as.factor(icc_p$para)
levels(icc_p$para) <- c('adj.p.value < 0.05',
                        '95_CI_lwr > 0',
                        'p.value < 0.05')
icc_p$para <- factor(icc_p$para, levels = c('95_CI_lwr > 0',
                                            'p.value < 0.05',
                                            'adj.p.value < 0.05'))

icc_p %>%
  ggplot(aes(x = ICC_B1, y = ICC_B2)) +
  facet_wrap(~para, ncol = 3) +
  geom_point(aes(colour = outcome)) +
  scale_colour_brewer(palette = 'Set2') +
  geom_abline(intercept = 0, slope = 1, colour = 'gray33') +
  theme_bw() +
  ggtitle('Intraclass Correlation Coefficient of breath VOCs between Core Visits')

icc_p

ggsave('ICC_CV_B1B2.tiff', unit = 'mm', dpi = 300, width = 205, height = 80)

# Ratio of endo- to exogenous VOCs in breath samples
endo_exo <- read.csv('Endo_Exo_filters.csv')[,-1]
endo_exo_2 <- read.csv('Endo_Exo_filters_2FC.csv')[,-1]

data <- b_corr_w1

endo_exo_fc <- function(data, filters){
  data_L <- data %>% 
    pivot_longer(cols =! c(1:4), names_to = 'comp', values_to = 'peakArea') %>%
    left_join(loads_glob %>% dplyr::select(comp1, comp, sign)) %>%
    drop_na() %>%
    left_join(filters %>% dplyr::select(comp, FC_filter, wilcox_filter)) %>%
    distinct()
  
  endo_exo_fc <- 
    data_L %>%
    #filter(comp %in% rownames(loads_vip_top)) %>%
    mutate(FC_filter = ifelse(FC_filter == 'Exo', 'Exo', 'Endo')) %>%
    group_by(Sample, FC_filter) %>% summarise(median = median(peakArea)) %>%
    pivot_wider(names_from = FC_filter, values_from = median) %>%
    left_join(data %>% dplyr::select(Sample, Diagnosis, CoreVisit)) %>%
    mutate(FC = Endo - Exo
           ) 
  }

data_fc <- endo_exo_fc(b_corr_w1, endo_exo)

hist(data_fc$FC_yes)
hist(data_fc$FC_no)

data_fc %>% ggplot(aes(x = Diagnosis, y = FC)) + geom_boxplot(outliers = FALSE)

data_fc %>% ggplot(aes(x = Diagnosis, y = Endo_yes)) + geom_boxplot(outliers = FALSE)
data_fc %>% ggplot(aes(x = Diagnosis, y = Exo_no)) + geom_boxplot(outliers = FALSE)

data_fc %>% ggplot(aes(x = Exo_yes, Endo_yes, colour = Diagnosis)) + geom_point()
data_fc %>% ggplot(aes(x = Exo_no, y = Endo_no, colour = Diagnosis)) + geom_point()


# best central measure for endo / exo?
dev.new()
data_L %>% 
  #filter(RAD_ID == 'RAD002') %>%
  #filter(Sample == 'RAD002_CV1') %>%
  ggplot(aes(x = peakArea)) + 
  geom_histogram(aes(fill = Diagnosis), position = 'identity', alpha = 0.6, bins = 10) + 
  facet_wrap(~FC_filter, scale = 'free')

dev.new()
data_L %>%
  filter(Sample %in% sample(Sample, 30)) %>%
  mutate(FC_filter = ifelse(FC_filter == 'Exo', 'Exo', 'Endo')) %>%
  #filter(FC_filter == 'Exo') %>%
  #filter(comp %in% rownames(loads_vip_top)) %>%
  ggplot(aes(x = peakArea)) +
  geom_histogram(aes(fill = sign), alpha = 0.5, position = 'identity') +
  theme_bw() +
  facet_wrap(~FC_filter, ncol = 1)
  #ggplot(aes(x = sign, y = peakArea, colour = comp)) + geom_jitter(size = 1) +
  #theme_bw()


# ROC
# fold change
library(pROC)
roc_fc_both <- roc(response = data_fc$Diagnosis, predictor = data_fc$FC_2)
plot(roc_fc_both, legacy.axes = TRUE)
View(as.data.frame(roc_fc_both['specificities']) %>%
  cbind(as.data.frame(roc_fc_both['sensitivities']), as.data.frame(roc_fc_both['thresholds'])))
coords(roc_fc_both, "best", ret = "threshold")

data_fc <- data_fc %>% mutate(pred = ifelse(FC_both < 2.25, 'Asthma', 'Not Asthma'))

conf_mat <- table(factor(data_fc$Diagnosis, levels = levels(as.factor(data_fc$Diagnosis))), 
               data_fc$pred)
sens_fun(conf_mat)
spec_fun(conf_mat)


#
#
#

# AUROC for most important VOCs from MINT trained on dataset 2
loads_vip_top <- loads_vip %>% arrange(desc(comp1)) %>% filter(comp1 > 1.5)
loads_vip_b1_top <- loads_vip_b1 %>% arrange(desc(comp1)) %>% filter(comp1 > 1.5)

voc_perf_b1 <- as.data.frame(sapply(rownames(loads_vip), function(voc, data = b1_corr_w1){
  resp <- data$Diagnosis
  pred <- data %>% dplyr::select(voc)
  colnames(pred) <- 'VOC'
  roc_voc <- roc(response = resp, predictor = pred$VOC)
  out_roc_voc <- as.data.frame(roc_voc['specificities']) %>%
         cbind(as.data.frame(roc_voc['sensitivities']), 
               as.data.frame(roc_voc['thresholds']))
  auc_voc <- auc(roc_voc)[[1]]
  sen_95spec <- out_roc_voc %>% 
    filter(specificities > 0.94 & specificities < 0.96) %>%
    summarise(mean(sensitivities))
  sen_90spec <- out_roc_voc %>% 
    filter(specificities > 0.88 & specificities < 0.91) %>%
    summarise(mean(sensitivities))
  spec_95sens <- out_roc_voc %>% 
    filter(sensitivities > 0.94 & sensitivities < 0.96) %>%
    summarise(mean(specificities))
  spec_90sen <- out_roc_voc %>% 
    filter(sensitivities > 0.88 & sensitivities < 0.91) %>%
    summarise(mean(specificities))
  out_voc <- c(auc_voc, sen_95spec, sen_90spec, spec_95sens, spec_90sen)
}))

rownames(voc_perf_b1) <- c('AUROC', 'sen_95spec', 'sen_90spec', 'spec_95sen', 'spec_90sen')
voc_perf_bt <- t(voc_perf_b) %>% as.data.frame() %>% 
  mutate(comp = rownames(.),
         Top = ifelse(comp %in% rownames(loads_vip_top) & comp %in% rownames(loads_vip_b1_top), 'Both', 
                      ifelse(comp %in% rownames(loads_vip_top), 'Dataset 2', 
                             ifelse(comp %in% rownames(loads_vip_b1_top), 'Dataset 1', 'None'
                        )
                      )))

voc_perft <- voc_perf_bt %>% left_join(voc_perf_b1t, by = 'comp') %>%
  mutate(name = ifelse(AUROC.x > 0.601 & AUROC.y > 0.601, comp, NA))

voc_perft %>% ggplot(aes(x = as.numeric(AUROC.x), y = as.numeric(AUROC.y), colour = Top.x)) + 
  geom_point()  + theme_bw() + coord_fixed() + 
  geom_hline(yintercept = 0.6, linetype = 'dashed') +
  geom_vline(xintercept = 0.6, linetype = 'dashed') +
  geom_text_repel(aes(label = name), colour = 'black') +
  scale_colour_discrete(name = 'Top PLS-DA loadings') +
  xlab('AUROC in Dataset 2') +
  ylab('AUROC in Dataset 1') +
  ggtitle('Area under ROC curve for individual VOCs')

ggsave('AUROC_comparison.tiff', unit = 'mm', width = 150, height = 120, dpi = 300)

voc_perft %>% ggplot(aes(x = as.numeric(sen_90spec.x), y = as.numeric(sen_90spec.y), colour = Top.x)) + 
  geom_point() + coord_fixed() + theme_bw() +
  ggtitle('Sensitivity at 90% specificity for individual VOC') +
  xlab('Dataset 2') + ylab('Dataset 1') +
  theme(legend.position = 'none')

ggsave('Senat90spec_comparison.tiff', unit = 'mm', width = 130, height = 100, dpi = 300)

voc_perft %>% 
  #mutate(names = ifelse(spec_95sen.x > 18 & spec_95sen.y > 18, comp, NA)) %>%
  ggplot(aes(x = as.numeric(spec_90sen.x), y = as.numeric(spec_90sen.y), colour = Top.x)) + 
  geom_point() + coord_fixed() + theme_bw() +
  ggtitle('Specificity at 90% sensitivity for individual VOCs') +
  #geom_text_repel(aes(label = name), colour = 'black') +
  xlab('Dataset 2') + ylab('Dataset 1') +
  theme(legend.position = 'none')

ggsave('Specat90sen_comparison.tiff', unit = 'mm', width = 130, height = 100, dpi = 300)

hist(as.numeric(voc_perft$sen_90spec.y), breaks = 17)
hist(as.numeric(voc_perft$spec_90sen.y), breaks = 15)

plot(voc_perf1_t$AUROC, voc_perf_t$AUROC) + abline(0,1)
plot(voc_perf1_t$sen_90spec, voc_perf_t$sen_90spec) + abline(0,1)

write.csv(voc_perf_t, 'MINT_B2_top_VOC_perf_B1_S90.csv')

View(endo_exo %>% filter(comp %in% rownames(loads_vip_top)))

#
#

model <- glm(as.factor(Diagnosis) ~ X3_methylpentane +  Ethyl_butanoate + CoreVisit,
              data = b_corr_w1,
              family = 'binomial')

summary(model)

box2 <- b_corr_w1 %>% ggplot(aes(x = CoreVisit, y = X3_methylpentane, fill = Diagnosis)) +
  geom_violin(alpha = 0.7) + 
  geom_boxplot(width = 0.3, position = position_dodge(width = 0.9), alpha = 0) + 
  theme_bw() +
  #theme(legend.position = 'none') +
  ylab('log(peak area)') +
  scale_fill_manual(values = c('Not Asthma' = 'orange', 'Asthma' = 'blue'))

box1 <- box1 + ggtitle('Dataset 1')
box1

box2 <- box2 + ggtitle('Dataset 2')
box2

boxes <- arrangeGrob(box1, box2, widths = c(0.42, 0.58))
plot(boxes)

ggsave('3methylpentane_boxplots_diagnosis.tiff', boxes, dpi = 300, unit = 'mm', width = 160, height = 60)

b2_uni <- b_imp_sum_c %>% filter(comp == 'X3_methylpentane') %>%
  left_join(b_corr_w1 %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  left_join(clin_dyn_imp_b2 %>% dplyr::select(Sample, FVCPre))

outl <- infl_obs %>% filter(comp == 'X3_methylpentane')

b2_uni1 <- b2_uni %>% filter(Sample %ni% outl$Sample)

#

b1_uni <- b1_imp_sum_c %>% filter(comp == 'X3_methylpentane') %>%
  left_join(b1_corr_w1 %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  left_join(clin_dyn_imp_b1 %>% dplyr::select(Sample, FVCPre))


uniMod <- lmer(log(S) ~ log(BG) + Diagnosis + CoreVisit + (1|RAD_ID),
               data = b2_uni1)

summary(uniMod)

uniMod1 <- lmer(log(S) ~ log(BG) + Diagnosis + CoreVisit + FVCPre + (1|RAD_ID),
                data = b2_uni1)
summary(uniMod1)
