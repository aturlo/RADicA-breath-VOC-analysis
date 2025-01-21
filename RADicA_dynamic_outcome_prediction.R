## Predicting asthma based on breath VOC samples

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(mixOmics)
library(forcats)
library(foreign)
library(stringr)
library(ggcorrplot)
library(ggfortify)
library(irr)

# custom functions
'%ni%' <- Negate('%in%')

# load normalised bg-corrected data w/o multivariate outliers
b1_corr_w1 <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b_corr_w1 <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]
endo_exo <- read.csv('Endo_Exo_filters.csv')[,-1]

# set columns in the same order
b1_corr_w1 <- b1_corr_w1[names(b_corr_w1)]

rownames(b1_corr_w1) <- b1_corr_w1$Sample
rownames(b_corr_w1) <- b_corr_w1$Sample

# load clinical data
clin <- read.spss('RADicA Active 29052024.sav', to.data.frame = TRUE)
clin_core <- clin %>% dplyr::select(BBID, Ethnicity, MaleFemale, AgeInYears, SmokerPQ, SmokerTypePQ,
                                    ACQPointsCV1, ACQPointsCV2, FeNOCV1, FeNOCV2, FEV1PreCV1, FEV1PreCV2,
                                    FEV1PPrePredCV1, FEV1PPrePredCV2, FVCPreCV1, FVCPreCV2,
                                    FVCPPrePredCV1, FVCPPrePredCV2, FEV1PMPreCV1, FEV1PMPreCV2,
                                    EosinophilsCV1, EosinophilsCV2) 

# format clinical data 
clin_core$BBID <- str_sub(clin_core$BBID, end = -2)
clin_core <- clin_core %>% rename(RAD_ID = BBID) %>%
  filter(RAD_ID != '      ')

# static data
clin_stat <- clin_core %>% dplyr::select(RAD_ID, Ethnicity, MaleFemale, AgeInYears, SmokerPQ, 
                                         SmokerTypePQ, EosinophilsCV1, EosinophilsCV2) 

clin_dyn <- clin_core %>% dplyr::select(!colnames(clin_stat)) %>% cbind(clin_core$RAD_ID) %>%
  rename(RAD_ID = 'clin_core$RAD_ID') %>% relocate(RAD_ID)
  
clin_dyn_cv1 <- clin_dyn %>% dplyr::select(c(RAD_ID, contains('CV1'))) %>%
  mutate(CoreVisit = 'CV1') %>% relocate(CoreVisit)
colnames(clin_dyn_cv1)[-c(1:2)] <- gsub('CV1', '', colnames(clin_dyn_cv1[-c(1:2)]))

clin_dyn_cv2 <- clin_dyn %>% dplyr::select(c(RAD_ID, contains('CV2'))) %>%
  mutate(CoreVisit = 'CV2') %>% relocate(CoreVisit)
colnames(clin_dyn_cv2)[-c(1:2)] <- gsub('CV2', '', colnames(clin_dyn_cv2[-c(1:2)]))

clin_dyn <- rbind(clin_dyn_cv1, clin_dyn_cv2)
  
clin_stat <- clin_stat %>%
  mutate(Eosinophils = ifelse(is.na(EosinophilsCV2) == TRUE, EosinophilsCV1, EosinophilsCV2)) %>%
  dplyr::select(!c(EosinophilsCV2,EosinophilsCV1))

# keep only patients with VOC data available
clin_stat <- clin_stat %>% filter(RAD_ID %in% c(b_corr_w1$RAD_ID, b1_corr_w1$RAD_ID))
clin_dyn <- clin_dyn %>% filter(RAD_ID %in% c(b_corr_w1$RAD_ID, b1_corr_w1$RAD_ID))

clin_dyn_voc <- clin_dyn %>% mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  filter(Sample %in% c(b_corr_w1$Sample, b1_corr_w1$Sample))

#
#
#

# explore missing data pattern
table(is.na(clin_stat))
which(is.na(clin_stat) == TRUE, arr.ind = TRUE)
colnames(clin_stat)[2]

table(is.na(clin_dyn))

# exclude rows with all 7 dynamic data types missing
excl <- table(which(is.na(clin_dyn) == TRUE, arr.ind = TRUE)[,1]) %>% as.data.frame() %>% 
  filter(Freq == 7) %>% dplyr::select(Var1) %>% mutate(Var1  = as.numeric(levels(Var1))[Var1])

clin_dyn <- clin_dyn[-excl$Var1,]
table(is.na(clin_dyn))

#
#
#

###############################

# EXPLORATORY ANALYSIS OF DYNAMIC CLINICAL OUTCOMES
# explore correlation structure of dynamic outcomes
clin_cor <- cor(clin_dyn[,-c(1:2)], use = 'complete.obs')
clin_p <- cor_pmat(clin_dyn[,-c(1:2)], use = 'complete.obs')
ggcorrplot(clin_cor, lab = TRUE, hc.method = 'ward.D2', hc.order = TRUE,
           p.mat = clin_p, insig = 'blank', tl.cex = 10, lab_size = 3) +
  ggtitle('Pearsons correlation coefficients')

ggsave('Clin_dyn_cor.tiff', unit = 'mm', dpi = 300, width = 120, height = 100)

# unsupervised analysis of dynamic outcomes
clin_dyn1 <- clin_dyn %>% drop_na() %>%
  left_join(rbind(b1_corr_w1, b_corr_w1) %>% 
              dplyr::select(RAD_ID, Diagnosis)) %>%
  distinct() %>%
  mutate(Dataset = ifelse(RAD_ID %in% b1_corr_w1$RAD_ID, 'B1', 'B2')) %>%
  relocate(Dataset, Diagnosis) #%>%
  #dplyr::select(!c(FEV1Pre, FVCPre))

clin_pca <- prcomp(clin_dyn1[,-c(1:4)], scale = TRUE, center = TRUE)
summary(clin_pca)

autoplot(clin_pca, data = clin_dyn1, 
         #label.label = 'RAD_ID', label.size = 2, shape = FALSE,
         loadings = TRUE, loadings.label = TRUE, loadings.colour = 'black', loadings.label.colour = 'black',
         loadings.label.repel = TRUE,
         colour = 'Diagnosis', shape = 'CoreVisit', size = 2) +
  scale_colour_manual(values = c('Asthma' = 'dodgerblue2',
                                 'Not Asthma' = 'darkorange1')) +
  theme_bw() + #guides(col="none") 
  ggtitle('Principal Component Analysis of dynamic clinical outcomes')

ggsave('PCA_dynamic_outcomes.tiff', dpi = 300, unit = 'mm', width = 150, height = 100)


# agreement between core visits
icc_clin_dyn <- bind_rows(lapply(colnames(clin_dyn)[-c(1:2)], function(out) {
  input <- clin_dyn %>% dplyr::select(RAD_ID, CoreVisit, out) %>%
  pivot_wider(names_from = CoreVisit, values_from = out) %>%
  dplyr::select(!RAD_ID)
  icc_res <- icc(input, model = 'twoway', type = 'agreement', unit = 'single')
  output <- data.frame(var = out, 
                     ICC = icc_res$value, 
                     p.value = icc_res$p.value,
                     CI_lwr = icc_res$lbound, CI_upr = icc_res$ubound)
  }))


write.csv(icc_clin_dyn, 'ICC_dynamic_clinical_outcomes.csv')

#
#
#

##############################

# INTEGRATION OF CLINICAL AND BREATH VOC DATA
# impute missing outcome values
which(is.na(clin_dyn_voc) == TRUE, arr.ind = TRUE)
clin_dyn_voc1 <- clin_dyn_voc[-105,] %>%
  relocate(Sample)

clin_dyn_imp <- impute.nipals(clin_dyn_voc1[,-c(1:3)], ncomp = 7) %>%
  cbind(clin_dyn_voc1[,1:3]) %>%
  relocate(Sample, RAD_ID, CoreVisit)

clin_dyn_imp_b1 <- clin_dyn_imp %>% filter(Sample %in% b1_corr_w1$Sample) %>%
  relocate(Sample) %>%
  dplyr::select(Sample, FeNO, FEV1PMPre, FEV1PPrePred, FVCPre) 

clin_dyn_imp_b1 <- clin_dyn_imp_b1[match(b1_corr_w1$Sample, clin_dyn_imp_b1$Sample),] %>%
  drop_na()

clin_dyn_imp_b2 <- clin_dyn_imp %>% filter(Sample %in% b_corr_w1$Sample) %>%
  relocate(Sample) %>%
  dplyr::select(Sample, FeNO, FEV1PMPre, FEV1PPrePred, FVCPre)

clin_dyn_imp_b2 <- clin_dyn_imp_b2[match(b_corr_w1$Sample, clin_dyn_imp_b2$Sample),]

rownames(clin_dyn_imp_b2) <- clin_dyn_imp_b2$Sample
rownames(clin_dyn_imp_b1) <- clin_dyn_imp_b1$Sample

#
#
#

# Mapping new outcomes to unsupervised analysis
data_b <- b1_corr_w1
meta_b <- clin_dyn_imp_b1
var1 <- 'FVCPre'

pc_b <- mixOmics::pca(data_b[,-c(1:4)], scale = TRUE, center = TRUE, ncomp = 5)
barplot(pc_b$prop_expl_var$X)

b_scores <- pc_b$variates$X %>% as.data.frame() %>% 
  mutate(Sample = rownames(.),
         CoreVisit = str_sub(Sample, start = 8),
         RAD_ID = str_sub(Sample, end = -5)) %>%
  left_join(
    meta_b)
    #clin_stat %>% dplyr::select(RAD_ID, Eosinophils))
b_scores <- b_scores %>%
  arrange(desc(b_scores[,var1])) %>%
  drop_na()

# create quantile colour scale 
out <- b_scores[,var1]

colvalquant <- base::seq(from = 0, to = 1, length.out = 5)
quants <- quantile(out, colvalquant)
b_scores$ptile_var <- ecdf(out)(out)

b1_fvc <- 
  b_scores %>% ggplot(aes(x = PC1, y = PC2, colour = ptile_var, shape = CoreVisit)) + 
  geom_point(size = 2) +
  theme_bw() +  scale_colour_gradient2(low = '#0000FFFF',
                                     mid = '#009E61FF',
                                     high = '#FFFF00FF',
                                     midpoint = 0.5,
                                     labels = round(quants,1),
                                     name = 'FVC') +
  xlab(paste('PC1 (', round(pc_b$prop_expl_var$X[1], 3)*100, '%)')) +
  ylab(paste('PC2 (', round(pc_b$prop_expl_var$X[2], 3)*100, '%)'))


b1_feno
b1_fev1pm
b1_fev1
b2_fev1
b2_fev1pm
b2_feno

b2_fvc
b1_fvc

pca_dynout <- arrangeGrob(b1_feno, b1_fev1pm, b1_fev1, 
            b2_feno, b2_fev1pm, b2_fev1, 
            ncol = 3, nrow = 2)

plot(pca_dynout)

ggsave('PCA_VOC_dynamic_outcomes.tiff', pca_dynout, unit = 'mm', dpi = 300, 
       width = 310, height = 140)


#
#
#

clin_dyn_imp_b2 <- clin_dyn_imp_b2 %>% mutate(FeNO2 = log(FeNO))
clin_dyn_imp_b1 <- clin_dyn_imp_b1 %>% mutate(FeNO2 = log(FeNO))

#
#
#

# MINT with continuous response (FEV1PPred, FeNO, FEV1PMPre)

b1_corr_w1 <- b1_corr_w1 %>% filter(Sample != 'RAD010_CV2')

rownames(b_corr_w1) <- b_corr_w1$Sample
rownames(b1_corr_w1) <- b1_corr_w1$Sample

voc_train <- b_corr_w1 #%>% filter(Sample != 'RAD010_CV2')
data_train <- clin_dyn_imp_b2 
out <- 'FVCPre'

X <- voc_train[,-c(1:4)] 
Y <- data_train %>% dplyr::select(out) %>% rename(outcome = out) 

study <- voc_train %>% dplyr::select(CoreVisit)

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2

study <- as.factor(study$CoreVisit)

base_reg2 <- mint.pls(X, Y$outcome, study = study, ncomp = 1, scale = TRUE)

base_mod <- base_reg2

base_mod$prop_expl_var$X

variates <- base_mod$variates$X
  #base_reg$variates.partial$X[2] %>% as.data.frame()

mint_scores1 <- variates %>% as.data.frame() %>% 
  mutate(Sample = rownames(.),
  RAD_ID = str_sub(Sample, end = -5)) %>%
  left_join(data_train) %>%
  left_join(voc_train %>% dplyr::select(Sample, Diagnosis))

b_scoresp <- 
  mint_scores2 %>% 
  ggplot(aes(x = comp1, y = mint_scores2[,out], colour = Diagnosis)) + geom_point() +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        #legend.position = 'none'
        ) +
  ylab('FVC (L)') +
  xlab(paste('Component 1 (', round(base_mod$prop_expl_var$X[[3]]*100, 2), '% )')) +
  ggtitle('Dataset 2') +
  scale_colour_manual(values = c('Asthma' = 'dodgerblue2',
                                 'Not Asthma' = 'darkorange1')) #+
  #geom_hline(yintercept = 75, linetype = 'dashed', colour = 'grey20') #+
  ylim(NA,301)

cor(mint_scores$comp1, log(mint_scores[,out]), use = 'complete.obs')

b_scoresp
b1_scoresp

# feno_scoresp

fvcpred_scoresp <- arrangeGrob(b1_scoresp, b_scoresp, nrow = 1, widths = c(0.42, 0.58))
plot(fvcpred_scoresp)

ggsave('MINT_PLS_FVCPre_scores_plots.tiff', fvcpred_scoresp, unit = 'mm', dpi = 300, 
       width = 180, height = 70)

#
#
#

# loadings
vip_b1 <- vip(base_reg) %>% as.data.frame() %>% mutate(comp = rownames(.))
vip_b2 <- vip(base_reg2) %>% as.data.frame() %>% mutate(comp = rownames(.))

View(vip_b1 %>% filter(comp1 > 1))
nrow(vip_b2 %>% filter(comp1 > 1))

shared_voc <- intersect(rownames(vip_b1 %>% filter(comp1 > 1)),
          rownames(vip_b2 %>% filter(comp1 > 1)))

View(loads_glob_b1 %>% filter(comp %in% shared_voc) %>% left_join(loads_glob_b2, by = 'comp'))

loads_glob_b1 <- base_reg$loadings$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>% 
  left_join(endo_exo %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>%
  mutate(FC_filter = ifelse(FC_filter == 'Exo', FC_filter, 
                            ifelse(FC_filter == 'Endo_both_datasets', FC_filter,
                                   'Endo_one_dataset')))

rownames(loads_glob_b1) <- loads_glob_b1$comp
rownames(loads_glob_b2) <- loads_glob_b2$comp


#
#
#

Lglob1 <- loads_glob_b1 %>%
  arrange(desc(abs(comp1))) %>% 
  mutate(comp = str_trunc(comp, 33)) %>%
  slice(1:40) %>% 
  ggplot(aes(x = comp1, y = fct_inorder(as.factor(comp)), 
             fill = FC_filter
  )) + 
  geom_col() +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw() +
  theme(axis.title.y = element_blank())
dev.new()

Lglob1 <- Lglob1 + ggtitle('Dataset 1') + theme(plot.title = element_text(hjust = 0.5), 
                                                legend.position = 'none')
Lglob2 <- Lglob2 + ggtitle('Dataset 2') + theme(plot.title = element_text(hjust = 0.5),
                                                legend.position = 'none')
Lglob1
Lglob2


load_plots_fvcpred <- arrangeGrob(Lglob1, Lglob2, nrow = 1, widths = c(0.5, 0.5)) 
plot(load_plots_fvcpred)

ggsave('MINT_PLS_FVCPre_loading_plots.tiff', load_plots_fvcpred, dpi = 300, unit = 'mm', 
       width = 155, height = 80)

# correlation of outcome with individual VOCs
voc <- 'Unknown_C12H24'

cor(b1_corr_w1[,voc], clin_dyn_imp_b1$FVCPre)
plot(b1_corr_w1[,voc], clin_dyn_imp_b1$FVCPre)

cor(b1_corr_w1[, voc], clin_dyn_imp_b1$FVCPre)
plot(b1_corr_w1[,voc], clin_dyn_imp_b1$FVCPre)

test <- b1_corr_w1 %>% left_join(clin_dyn_imp_b1 %>% dplyr::select(Sample, FVCPre)) %>%
  mutate(Diagnosis = ifelse(Diagnosis == 'Asthma', 1, 0))

test %>% ggplot(aes(x = Diagnosis, y = FVCPre, group = Diagnosis)) + geom_boxplot()

finMod <- glm(Diagnosis ~ X3_methylpentane + Ethyl_butanoate,
                data = test,
                family = 'binomial'(link = 'logit'))

summary(finMod)

cor(b_corr_w1$Ethyl_butanoate, b_corr_w1$X3_methylpentane)
cor(b1_corr_w1$Furan._2_methyl_, b1_corr_w1$X3_methylpentane)

#
#
#

# predictive performance
# train data
base_reg_pred <- predict(base_mod, newdata = X, study = study)

prediction_fun <- function(pred_data, out_data, out){
  pred_data[["predict"]] %>% as.data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  left_join(out_data %>% dplyr::select(out, Sample)) %>%
  rename(outcome = out) %>%
  mutate(error = Y.dim1 - outcome)
}

prediction <- prediction_fun(base_reg_pred, data_train, out)

cor(prediction$Y.dim1, prediction$outcome)
  
mse_fun <- function(pred_data) {mean((pred_data$error)^2, na.rm = TRUE)}
rmse <- sqrt(mse)
mae_fun <- function(pred_data) {mean(abs(pred_data$error),  na.rm = TRUE)}
mape_fun <- function(pred_data) {mean(abs(pred_data$error/pred_data$outcome)*100, na.rm = TRUE)}
me_fun <- function(pred_data) {mean(pred_data$error, na.rm = TRUE)}

mse <- mse_fun(prediction)
mae <- mae_fun(prediction)
me <- me_fun(prediction)
mape <- mape_fun(prediction)

mse
mae
me
mape

#
#
#

# cross-validation on train data
out <- 'FVCPre'

tmp_data <- voc_train %>%
  left_join(data_train %>% 
              dplyr::select(Sample, out)) %>% relocate(out)

rownames(tmp_data) <- tmp_data$Sample

set.seed(1410)

cv_perf <- 
  replicate(n = 100, expr = {
  
    p <- n_distinct(tmp_data$RAD_ID)
    x <- sample(1:p)
    names(x) <- unique(tmp_data$RAD_ID)
    levels <- 3
    folds <- split(x, x%%levels)
    
    MAE <- vector(mode = "numeric")
    ME <- vector(mode = 'numeric')
    MAPE <- vector(mode = 'numeric')
    COR <- vector(mode = 'numeric')
    
    for (m in as.character(0:(levels - 1))) {
      test_id <- folds[[m]]
      test_data_cv <- tmp_data %>% filter(RAD_ID %in% names(test_id)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
      train_data_cv <- tmp_data %>% filter(RAD_ID %ni% names(test_id)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
      
      model_cv <- mint.pls(X = train_data_cv[,-c(1:6)], 
                        Y = train_data_cv[,out], 
                        study = as.factor(train_data_cv$study), 
                        ncomp = 1, scale = TRUE)
      
      model_cv_pred <- predict(model_cv,
                            newdata = test_data_cv[,-c(1:6)],
                            study.test = as.factor(test_data_cv$study),
                            scale = TRUE)
      
      prediction_cv <- prediction_fun(model_cv_pred, test_data_cv, out)
      
    MAE[m] <- mae_fun(prediction_cv)
    MAPE[m] <- mape_fun(prediction_cv)
    ME[m] <- me_fun(prediction_cv)
    COR[m] <- cor(prediction_cv$outcome, prediction_cv$Y.dim1, use = 'complete.obs')
    
    }
    
    cv_mae <- mean(MAE, na.rm = TRUE)
    cv_mape <- mean(MAPE, na.rm = TRUE)
    cv_me <- mean(ME, na.rm = TRUE)
    cv_cor <- mean(COR, na.rm = TRUE)
    
    out <- c(cv_mae, cv_mape, cv_me, cv_cor) 
    
  })


mae <- c(mean(cv_perf[1,], na.rm = TRUE), sd(cv_perf[1,], na.rm = TRUE))
mape <- c(mean(cv_perf[2,], na.rm = TRUE), sd(cv_perf[2,], na.rm = TRUE))
me <- c(mean(cv_perf[3,], na.rm = TRUE), sd(cv_perf[3,], na.rm = TRUE))
cor <- c(mean(cv_perf[4,], na.rm = TRUE), sd(cv_perf[4,], na.rm = TRUE))

mae
mape
me
cor

#
#
#

# test data
voc_train <- b_corr_w1
voc_test <- b1_corr_w1

data_train <- clin_dyn_imp_b2
data_test <- clin_dyn_imp_b1
out <- 'FVCPre'

train <- voc_train %>% 
  left_join(data_train %>% dplyr::select(Sample, out)) %>%
  mutate(study = ifelse(CoreVisit == 'CV1', 1 ,2)) %>%
  relocate(out, study) %>%
  rename(outcome = out) 

test_b1_df <- voc_test %>% #filter(CoreVisit == 'CV2') %>%
  left_join(data_test %>% dplyr::select(Sample, out)) %>%
  mutate(study = ifelse(CoreVisit == 'CV1', 1, 2)) %>%
  relocate(out, study) %>%
  rename(outcome = out) %>%
  drop_na()

rownames(test_b1_df) <- test_b1_df$Sample
rownames(train) <- train$Sample

# method from the textbook
conc_corr_w1 <- rbind(train, test_b1_df)
rownames(conc_corr_w1) <- conc_corr_w1$Sample

X <- conc_corr_w1[,-c(1:6)]
Y <- conc_corr_w1$outcome
study <- conc_corr_w1$study

study <- as.factor(study)

test_b1 <- which(study == '3')

base_reg1 <- mint.pls(X = X[-c(test_b1),],
                   Y = Y[-c(test_b1)],
                   study = droplevels(study[-c(test_b1)]),
                   ncomp = 1,
                   scale = TRUE)

base_reg_pred_test <- predict(base_reg1,
                              newdata = X[c(test_b1),],
                              study.test = study[c(test_b1)],
                              scale = TRUE)

# entering train and test data separately
base_reg_2 <- mint.pls(X = train[,-c(1:6)],
                      Y = train$outcome,
                      study = as.factor(train$study),
                      ncomp = 1,
                      scale = TRUE)

base_reg_pred_test2 <- predict(base_reg_2,
                               newdata = test_b1_df[,-c(1:6)],
                               study.test = as.factor(test_b1_df$study))

prediction_test <- prediction_fun(base_reg_pred_test2, data_test, out)
plot(prediction_test$Y.dim1, prediction_test$outcome) + abline(0,1,)
cor(prediction_test$Y.dim1, prediction_test$outcome, use = 'complete.obs')

mae_fun(prediction_test)
mape_fun(prediction_test)
me_fun(prediction_test)

#
#
#

# tuning number of variables using cross-validation
# M-fold CV for variable number selection
out <- 'FVCPre'

data_vs <- b1_corr_w1 %>%
  left_join(clin_dyn_imp_b1 %>% 
              dplyr::select(Sample, out)) %>% relocate(out) %>%
  rename(outcome = out)

rownames(data_vs) <- data_vs$Sample
lv <- 3
nvar <- 50

set.seed(1410)

cv_perf_vs <-
  replicate(n = 50, expr = {
    p_vs <- n_distinct(data_vs$RAD_ID)
    x_vs <- sample(1:p)
    names(x_vs) <- unique(data_vs$RAD_ID)
    levels_vs <- lv
    folds_vs <- split(x_vs, x_vs%%levels_vs)
    
    vs_MAE <- data.frame()
    vs_COR <- data.frame()
    vs_vocs <- list()
    
    for (m in as.character(0:(levels_vs - 1))) {
      test_id_vs <- folds_vs[[m]]
      test_data_vs <- data_vs %>% filter(RAD_ID %in% names(test_id_vs)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
      train_data_vs <- data_vs %>% filter(RAD_ID %ni% names(test_id_vs)) %>%
        mutate(study = ifelse(CoreVisit == 'CV1', '1', '2')) %>% relocate(study)
      
      vs_ber <- function(vn) {
        
        MAE1 <- vector(mode = "numeric", vn)
        COR1 <- vector(mode = "numeric", vn)
        VOCS <- list()
        
        for(n in seq_len(vn)) {
          model_vs <- mint.spls(X = train_data_vs[,-c(1:6)], 
                               Y = train_data_vs$outcome, 
                               study = as.factor(train_data_vs$study), 
                               ncomp = 1, 
                               scale = TRUE,
                               keepX = c(n))
          
          model_vs_pred <- predict(model_vs,
                                   newdata = test_data_vs[,-c(1:6)],
                                   study.test = as.factor(test_data_vs$study),
                                   scale = TRUE)
          
          prediction_vs <- prediction_fun(model_vs_pred, test_data_vs, 'outcome')
          
          loads <- model_vs[['loadings']][['X']] %>% as.data.frame() %>% filter(comp1 != 0) 
          vocs <- rownames(loads)
        
          VOCS <- append(VOCS, list(vocs))
          MAE1[n] <- mae_fun(prediction_vs)
          COR1[n] <- cor(prediction_vs$outcome, prediction_vs$Y.dim1, use = 'complete.obs')
        }
        out_vs <- list(MAE1, COR1, VOCS)
      }
      
      vs_out <- vs_ber(nvar)
      
      vs_MAE <- rbind(vs_out[[1]], vs_MAE)
      colnames(vs_MAE) <- sapply(1:nvar, function(g){paste('VN', g, sep = '')})
      
      vs_COR <- rbind(vs_out[[2]], vs_COR)
      colnames(vs_COR) <- sapply(1:nvar, function(g){paste('VN', g, sep = '')})
      
      vs_vocs <- append(vs_vocs, list(unlist(vs_out[[3]])))
    }
    out_all_vs <- list(vs_MAE, vs_COR, vs_vocs)
  }, simplify = F)
      



MAE2 <- sapply(1:length(cv_perf_vs), function(l){
  t(cv_perf_vs[[l]][[1]]) %>% as.data.frame() %>% rowMeans()
  })

COR2 <- sapply(1:length(cv_perf_vs), function(l){
  t(cv_perf_vs[[l]][[2]]) %>% as.data.frame() %>% rowMeans()
})

dev.new()
MAE2l <- t(MAE2) %>% as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'VN', values_to = 'MAE') 

MAE2l %>%
  ggplot(aes(x = fct_inorder(as.factor(VN)), y = MAE)) + geom_boxplot(outliers = FALSE) +
  theme_bw() +
  scale_x_discrete(labels = c(1, rep('', 8), 10, rep('', 9), 20, rep('', 9), 30, rep('', 9), 40, rep('',9), 50)) +
  xlab('VOC number') +
  ggtitle('Mean absolute error') 

ggsave('MINT_FeNO_VS_B2_CV_MAE.tiff', dpi = 300, unit = 'mm', width = 110, height = 70)

View(MAE2l %>% group_by(VN) %>% summarise(mean = mean(MAE),
                                     sd = sd(MAE),
                                     rsd = sd/mean))


# 1-4 V lowest error

COR2l <- t(COR2) %>% as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'VN', values_to = 'COR') 

COR2l %>%
  ggplot(aes(x = fct_inorder(as.factor(VN)), y = COR)) + geom_boxplot(outliers = FALSE) +
  theme_bw() +
  scale_x_discrete(labels = c(1, rep('', 8), 10, rep('', 9), 20, rep('', 9), 30, rep('', 9), 40, rep('',9), 50)) +
  xlab('VOC number') +
  ylab('r') +
  ggtitle('Correlation (predicted vs true outcome)') 

ggsave('MINT_FeNO_VS_B2_CV_COR.tiff', dpi = 300, unit = 'mm', width = 110, height = 70)

View(COR2l %>% group_by(VN) %>% summarise(mean = mean(COR),
                                                      sd = sd(COR),
                                                      rsd = sd/mean))

# 1-7 V highest correlation
FeNO_cv_vs <- cbind(MAE2l, COR2l)

write.csv(FeNO_cv_vs, 'VarSelection_B2_FeNO.csv')

# pull only models with 4 variables from each CV repeat
vim <- lapply(1:length(cv_perf_vs), function(l){
  sapply(1:lv, function(p){
    table(cv_perf_vs[[l]][[3]][[p]] %>% as.data.frame() %>% slice(7:10)/1) 
    })
  })

# summarise how often within CV repeat VOC was selected
vim_sum <- bind_rows(lapply(1:length(vim), function(r) {
  vim[[r]] %>% as.data.frame() %>% rowMeans() %>% as.data.frame() %>%
    rename(Freq = '.') %>% mutate(comp = rownames(.))
  })) 

# transform into frequency
vim_freq <- vim_sum %>% group_by(comp) %>% summarise(Freq = sum(Freq)/length(vim))

write.csv(vim_sum, 'VarSelection_VIM_B2_FeNO.csv')

#
level_order <- c("X.S._..._6_Methyl_1_octanol" , "X2.6_difluorobenzaldehyde",
                 "Unknown_C12H24", "X2_Butenal",
                   "Alpha_terpinene", "Toluene","Octane._2.4.6_trimethyl_")

vim_freq %>%
  arrange(desc(Freq)) %>%
  ggplot(aes(x = fct_inorder(as.factor(comp)), y = Freq)) + geom_col() +
  theme_bw(base_size = 9) +
  ggtitle('Variables most frequently selected in cross-validation') +
  scale_x_discrete(labels = c('Unknown_C12H24', 'Toluene', "X.S._..._6_Methyl_1_octanol",
  'Alpha_terpinene', 'Styrene', rep('', 25), 31)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        plot.margin = margin(l = 15))


ggsave('MINT_FeNO_VS_B2_VIM_barplot.tiff', dpi = 300, unit ='mm', width = 110, height = 80)

#
#
#

# tuning number of variable using test dataset
vs_mae <- function(vn) {
  MAE <- vector(mode = 'numeric', vn)
  for(n in seq_len(vn)) {
    vs_model <- mint.spls(X = train[,-c(1:6)],
                          Y = train$outcome,
                          study = as.factor(train$study),
                          ncomp = 1,
                          scale = TRUE,
                          keepX = n)
    
    base_reg_pred_test1 <- predict(vs_model,
                                   newdata = test_b1_df[,-c(1:6)],
                                   study.test = as.factor(test_b1_df$study))
    vn_prediction <- prediction_fun(base_reg_pred_test1, data_test, out)
    MAE[n] <- mae_fun(vn_prediction)
    #MAE[n] <- cor(vn_prediction$Y.dim1, vn_prediction$outcome)
    }
  MAE
}

  
  
set.seed(1410)
vs_test <- vs_mae(100)
names(vs_test) <- c(1:100)
dev.new()
plot(names(vs_test), vs_test, xlab = 'VOC number', ylab = 'MAE')
vs_test[vs_test == min(vs_test)]

#
#
#

# MINT with continuous response binned according to diagnostic test thresholds
clin_dyn_imp_b2cat <- clin_dyn_imp_b2 %>% 
  mutate(FeNOc = ifelse(FeNO > 50, 'Hi', 'Lo'),
         FEV1PMPrec = ifelse(FEV1PMPre > 75, 'Hi', 'Lo'))

clin_dyn_imp_b1cat <- clin_dyn_imp_b1 %>% 
  mutate(FeNOc = ifelse(FeNO > 50, 'Hi', 'Lo'),
         FEV1PMPrec = ifelse(FEV1PMPre > 75, 'Hi', 'Lo'))

# check how often the same ID is assigned different label between visits
cats_2_wide <- clin_dyn_imp_b2cat %>% left_join(b_corr_w1 %>% dplyr::select(Sample, RAD_ID, CoreVisit)) %>%
  dplyr::select(RAD_ID, CoreVisit, FeNOc, FEV1PMPrec) %>%
  pivot_wider(names_from = CoreVisit, values_from = c(FeNOc, FEV1PMPrec)) %>%
  mutate(FeNOc_agree = ifelse(FeNOc_CV1 == FeNOc_CV2, 'Yes', 'No'),
         FEV1PMPrec_agree = ifelse(FEV1PMPrec_CV1 == FEV1PMPrec_CV2, 'Yes', 'No'))

table(cats_2_wide$FeNOc_agree)
table(cats_2_wide$FEV1PMPrec_agree)

X <- b_corr_w1[,-c(1:4)]
Y <- clin_dyn_imp_b2cat$FeNOc
study <- b_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2

study <- as.factor(study)

set.seed(1410)
base_cat <- mint.plsda(X, Y, study = study, ncomp = 1, scale = TRUE)

variates <- base_reg$variates$X
data <- clin_dyn_imp_b2cat
out <- 'FeNOc'

mint_scores <- variates %>% as.data.frame() %>% 
  mutate(Sample = rownames(.),
         RAD_ID = str_sub(Sample, end = -5)) %>%
  left_join(data) 

dev.new()
mint_scores %>% ggplot(aes(x = FeNOc, y = comp1)) + geom_boxplot() +
  theme_bw() +
  theme(legend.position = 'none')

# prediction on train dataset (no CV)
base_cat_pred_train <- predict(base_cat, newdata = X, study = study, scale = TRUE)

out_pred <- base_cat_pred_train[["class"]][["max.dist"]] %>% as.data.frame() %>% dplyr::select(1)
conf_mat <- table(factor(out_pred$comp1, levels = levels(as.factor(clin_dyn_imp_b2cat$FeNOc))), 
                  clin_dyn_imp_b2cat$FeNOc)

ber_fun(conf_mat)
sens_fun(conf_mat)
spec_fun(conf_mat)

fenoc_roc <- roc(response = clin_dyn_imp_b2cat$FeNOc, predictor = base_cat_pred_train[["predict"]][,,1][,1])
plot(fenoc_roc)

View(fenoc_roc[['sensitivities']] %>% cbind(fenoc_roc[['specificities']],
                                              fenoc_roc[['thresholds']]))


# prediction on test dataset
train <- b_corr_w1 %>% left_join(clin_dyn_imp_b2cat %>% dplyr::select(Sample, FeNOc)) %>% relocate(FeNOc)
test_b1_df <- b1_corr_w1 %>% #filter(CoreVisit == 'CV2') %>%
  mutate(CoreVisit = 'B1_CV') %>% left_join(clin_dyn_imp_b1cat %>% dplyr::select(Sample, FeNOc)) %>% relocate(FeNOc)

conc_corr_w1 <- rbind(train, test_b1_df)

rownames(conc_corr_w1) <- conc_corr_w1$Sample

X <- conc_corr_w1[,-c(1:5)]
Y = conc_corr_w1$FeNOc
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

pred_class_b1 <- mint_pred_b1[["class"]][["centroids.dist"]] %>% as.data.frame() 

conf_mat <- table(factor(test_b1_df$FeNOc, levels = levels(as.factor(test_b1_df$FeNOc))), 
                  pred_class_b1$comp1)

conf_mat

ber_fun(conf_mat)
sens_fun(conf_mat)
spec_fun(conf_mat)

#
#

# Mixed effect random forest
library(LongituRF)
data <- DataLongGenerator(n=20)

b1_corr_w1a <- b1_corr_w1 %>% filter(Sample != "RAD010_CV2")

train_merf <- b_corr_w1
test_merf <- b1_corr_w1a

train_merf_data <- clin_dyn_imp_b2
test_merf_data <- clin_dyn_imp_b1

res1 <- MERF(X = as.matrix(train_merf[,-c(1:4)]),
             Y = train_merf_data$FEV1PPrePred,
             id = as.numeric(as.factor(train_merf$RAD_ID)),
             time = as.numeric(as.factor(train_merf$CoreVisit)),
             Z = as.matrix(c(rep(1, nrow(train_merf)))),
             sto = 'none')

pred_merf <- predict(res1,
                X = as.matrix(test_merf[,-c(1:4)]),
                id = as.numeric(as.factor(test_merf$RAD_ID)),
                time = as.numeric(as.factor(test_merf$CoreVisit)),
                Z = as.matrix(c(rep(1, nrow(test_merf)))))

plot(pred_merf, test_merf_data$FEV1PPrePred)
cor(pred_merf, test_merf_data$FEV1PPrePred, use = 'complete.obs')
mean(abs(pred_merf - test_merf_data$FEV1PPrePred), na.rm = TRUE)
mean(pred_merf - test_merf_data$FEV1PPrePred, na.rm = TRUE)
mean((abs(pred - test_merf_data$FEV1PPrePred)/clin_dyn_imp_b1$FEV1PPrePred)*100, na.rm = TRUE)

#

res <- MERF(X = as.matrix(b_corr_w1[,-c(1:4)]),
            Y = clin_dyn_imp_b2$FeNO,
            id = c(rep(1,101)),
            time = c(rep(1, 101)),
            Z = as.matrix(c(rep(1, 101))),
            sto = 'none')

plot(res$predicted, clin_dyn_imp_b2$FeNO)

#
#
#

# PLS stacked data + GLM
library(lme4)
library(lmerTest)

rownames(b_corr_w1) <- b_corr_w1$Sample
rownames(b1_corr_w1) <- b1_corr_w1$Sample

voc_train <- b_corr_w1 #%>% filter(Sample != 'RAD010_CV2')
data_train <- clin_dyn_imp_b2 
out <- 'FVCPre'

X <- voc_train[,-c(1:4)] 
Y <- data_train %>% dplyr::select(out) %>% rename(outcome = out) 

pls <- pls(X = X, Y = Y, ncomp = 1, scale = TRUE)
pred_train_pls <- predict(pls, newdata = X) 
prediction_pls <- cbind(pred_train_pls$predict, Y, voc_train$RAD_ID)
colnames(prediction_pls) <- c('PLSpred', 'out', 'RAD_ID')
cor(prediction_pls$PLSpred, prediction_pls$out)
plot(prediction_pls$PLSpred, prediction_pls$out)

linModel <- lm(out ~ PLSpred,
               data = prediction_pls)

mixModel <- lmer(out ~ PLSpred + (1 | RAD_ID),
                 data = prediction_pls)

summary(linModel)
summary(mixModel)

linModelpred <- predict(linModel, newdata = prediction_pls)
cor(linModelpred, prediction_pls$out)

mixModelpred <- predict(mixModel, newdata = prediction_pls)
cor(mixModelpred, prediction_pls$out)

# prediction of test data
voc_test <- b1_corr_w1 %>% filter(Sample != 'RAD010_CV2')
data_test <- clin_dyn_imp_b1 

pred_test_pls <- predict(pls, newdata = voc_test[,-c(1:4)])
prediction_pls_test <- cbind(pred_test_pls$predict, data_test[, out], voc_test$RAD_ID) %>%
  as.data.frame()
colnames(prediction_pls_test) <- c('PLSpred', 'outcome', 'RAD_ID')
prediction_pls_test <- prediction_pls_test %>%
  mutate(PLSpred = as.numeric(PLSpred),
         outcome = as.numeric(outcome),
         error = PLSpred - outcome)

cor(prediction_pls_test$PLSpred, prediction_pls_test$outcome)
plot(prediction_pls_test$PLSpred, prediction_pls_test$outcome)

mae_fun(prediction_pls_test)
me_fun(prediction_pls_test)
mape_fun(prediction_pls_test)

#

mixModelpred_test <- predict(mixModel, newdata = prediction_pls_test,
                             allow.new.levels = TRUE)

prediction_lmer_test <- cbind(mixModelpred_test, prediction_pls_test$outcome) %>%
  as.data.frame()
colnames(prediction_lmer_test) <- c('LMERpred', 'outcome')
prediction_lmer_test <- prediction_lmer_test %>%
  mutate(error = LMERpred - outcome)

cor(prediction_lmer_test$LMERpred, prediction_lmer_test$outcome)
mae_fun(prediction_lmer_test)
me_fun(prediction_lmer_test)
mape_fun(prediction_lmer_test)

#
#
#

# 3-methylpentane
b2_uni1 <- b_corr_w1 %>% left_join(clin_dyn_imp_b2 %>% dplyr::select(Sample, FVCPre))

b1_uni <- b1_imp_sum_c %>% filter(comp == 'X3_methylpentane') %>%
  left_join(b1_corr_w1 %>% dplyr::select(RAD_ID, Diagnosis)) %>% distinct() %>%
  mutate(Sample = paste(RAD_ID, CoreVisit, sep = '_')) %>%
  left_join(clin_dyn_imp_b1 %>% dplyr::select(Sample, FVCPre))

scat2 <- b2_uni1 %>% drop_na() %>%
  ggplot(aes(x = FVCPre, y = log(S), colour = Diagnosis)) + geom_point() +
  theme_bw() +
  scale_colour_manual(values = c('Asthma' = 'blue', 'Not Asthma' = 'orange')) +
  ylab('log(peak area)') #+
  theme(legend.position = 'none')

scat2 <- scat2 + ggtitle('Dataset 2')
scat2

scat1 <- scat1 + ggtitle('Dataset 1')
scat1

cor(log(b1_uni$S), b1_uni$FVCPre, use = 'complete.obs')
cor(log(b2_uni1$S), b2_uni1$FVCPre, use = 'complete.obs')

scats <- arrangeGrob(scat1, scat2, widths = c(0.42, 0.58))
plot(scats)

ggsave('3-methylpentane_scatterplots_FVCPre.tiff', scats, dpi = 300, units = 'mm',
       width = 180, height = 70)

#

boxFvc1 <- b1_uni %>% drop_na() %>% ggplot(aes(x = Diagnosis, y = FVCPre, fill = Diagnosis)) + 
  geom_violin(alpha = 0.7) + geom_boxplot(width = 0.3, alpha = 0) + theme_bw() +
  scale_fill_manual(values = c('Asthma' = 'blue', 'Not Asthma' = 'orange')) +
  theme(legend.position = 'none')
  
boxFvc2 <- boxFvc2 + ggtitle('Dataset 2')
boxFvc2

boxFvc1 <- boxFvc1 + ggtitle('Dataset 1')
boxFvc1

boxFvc <- arrangeGrob(boxFvc1, boxFvc2, widths = c(0.4, 0.6))
plot(boxFvc)

ggsave('FVCPre_vs_Diagnosis_boxplots.tiff', boxFvc, dpi = 300, units = 'mm', width = 125, height = 60)
