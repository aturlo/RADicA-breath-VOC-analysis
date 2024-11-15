# RADicA classification method comparison

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(bimm)
library(mixOmics)

# load summarised BG corrected dataset 
b1_imp_sum_corr_w <- read.csv('RADicA_BG_adjusted.csv')[,-1]

# w/o multivariate outliers
b1_imp_sum_corr_w1 <- read.csv('RADicA_BG_adjusted_outl_removed.csv')[,-1]

b1_imp_sum_corr_w1 <- b1_imp_sum_corr_w1 %>% dplyr::select(!c('Ethyl_Acetate'))

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')
#
#
#

# BASELINE MODEL
# observations stacked vertically
data <- b1_imp_sum_corr_w1[,-1]
data$Diagnosis <- as.factor(data$Diagnosis)

out <- as.factor(data$Diagnosis)
assay <- data %>% dplyr::select(!c(1:3))

# without clustering
# PLS-DA
base_plsda <- plsda(X = assay, Y = out, scale = TRUE, ncomp = 1)

# error for the train set 
# without cross-validation
base_plsda_pred <- predict(base_plsda, assay, dist = 'max.dist')
out_pred <- base_plsda_pred$class$max.dist
conf_mat <- table(factor(out_pred, levels = levels(out)), out)
print(conf_mat)
ber <- 0.5*(conf_mat[1,2]/(conf_mat[1,2]+conf_mat[1,1]) + 
              conf_mat[2,1]/(conf_mat[2,1]+conf_mat[2,2]))
ber

# agreement between predictions at each time point
ev_pred <- cbind(data[,1:3], out_pred) %>%
  rename(Diagnosis_pred = comp1) %>%
  pivot_wider(names_from = CoreVisit, values_from = Diagnosis_pred) %>%
  mutate(Diagnosis = as.factor(Diagnosis),
         CV1 = as.factor(CV1),
         CV2 = as.factor(CV2),
         agreement = ifelse(CV1 == CV2, 'Yes', 'No')) 

table(ev_pred$agreement)

# prediction at specific time point
confusion_tp <- function(timePoint) {
  ev_pred_tp <- ev_pred %>% dplyr::select(Diagnosis, timePoint) %>% 
    drop_na() %>% rename(CV = timePoint)
  conf_mat <- table(factor(ev_pred_tp$Diagnosis, levels = levels(ev_pred_tp$Diagnosis)), 
        ev_pred_tp$CV)
  print(conf_mat)
  ber <- 0.5*(conf_mat[1,2]/(conf_mat[1,2]+conf_mat[1,1]) + 
                conf_mat[2,1]/(conf_mat[2,1]+conf_mat[2,2]))
  print(paste('BER', ber))
  }

confusion_tp('CV1')
confusion_tp('CV2')

# with cross-validation
perf_base_plsda <- perf(base_plsda, 
                        validation = 'Mfold',
                        folds = 8,
                        nrepet = 100)

perf_base_plsda[["error.rate"]]
perf_base_plsda[["error.rate.class"]]

# calculate subject-level prediction
cv_pred <- cbind(data$RAD_ID, as.data.frame(perf_base_plsda[["predict"]][["comp1"]])) %>%
  rename(RAD_ID = 'data$RAD_ID',
         Asthma_pred = "Asthma.nrep1",
         Not_Asthma_pred = "Not Asthma.nrep1")

sum_cv_pred <- cv_pred %>% group_by(RAD_ID) %>% summarise(Asthma_pred = mean(Asthma_pred),
                                           Not_Asthma_pred = mean(Not_Asthma_pred)) %>%
  left_join(data %>% dplyr::select(RAD_ID, Diagnosis)) %>%
  mutate(pred = ifelse(Asthma_pred > Not_Asthma_pred, 'Asthma', 'Not Asthma'),
         correct = ifelse(pred == Diagnosis, ' yes', 'no')) %>%
  distinct()

# classification error 
sub_conf_mat <- 
  sum_cv_pred %>%
  group_by(Diagnosis, correct) %>% summarise(n = n()) %>%
  as.data.frame()

(sub_conf_mat[2,3] + sub_conf_mat[4,3])/ sum(sub_conf_mat$n) 

# balanced error
0.5 * ((sub_conf_mat[2,3]/(sub_conf_mat[1,3] + sub_conf_mat[2,3])) +
         (sub_conf_mat[4,3]/(sub_conf_mat[3,3] + sub_conf_mat[4,3])))

# with clustering
# create M-fold cross-validation method 
# preserve all observations from one patient in either train or held-out fold
cv_BER <- 
  replicate(n = 100, expr = {
  p <- n_distinct(data$RAD_ID)
  x <- sample(1:p)
  names(x) <- unique(data$RAD_ID)
  levels <- 8
  folds <- split(x, x%%levels)
  
  #BER <- vector(mode = "numeric")
  ER <- vector(mode = 'numeric')
  
  for (m in as.character(0:(levels - 1))) {
    test_id <- folds[[m]]
    test_data <- data %>% filter(RAD_ID %in% names(test_id))
    train_data <- data %>% filter(RAD_ID %ni% names(test_id))
    model <- plsda(train_data[,-c(1:3)], Y = train_data$Diagnosis, ncomp = 1, scale = TRUE)
    model_pred <- predict(model, test_data[,-c(1:3)], dist = 'max.dist')
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
    out_pred <- model_pred$class$max.dist
    conf_mat <- table(factor(out_pred, levels = levels(train_data$Diagnosis)), test_data$Diagnosis)
    #ber <- 0.5*(conf_mat[1,2]/(conf_mat[1,2]+conf_mat[1,1]) + 
                #conf_mat[2,1]/(conf_mat[2,1]+conf_mat[2,2]))
    er <- (conf_mat[1,2] + conf_mat[2,1])/sum(conf_mat)
    #BER[m] <- ber
    ER[m] <- er
  }
  #cv_ber <- mean(BER, na.rm = TRUE)
  cv_er <- mean(ER, na.rm = TRUE)
  })


mean(cv_BER)

# Random forest ++
# data transformation for RF++ software
rf_input <- b1_imp_sum_corr_w1 %>% 
  dplyr::select(!c(Sample, CoreVisit)) %>%
  relocate(Diagnosis, .after = Dl_Menthol)

rf_input$RAD_ID <- as.factor(rf_input$RAD_ID) %>% as.numeric(.)
rf_input$Diagnosis <- as.factor(rf_input$Diagnosis) %>% as.numeric(.)

write.csv(rf_input, 'RF_input.csv')

# load results from RF++ 
rf_out_sample <- read.csv('RF_output_sample_level.csv')
rf_out_subject <- read.csv('RF_output_subject_level.csv')

rf_out_sample <- rf_out_sample %>% 
  mutate(correct = ifelse(TRUE. == Predicted, 'yes', 'no'))

rf_out_subject <- rf_out_subject %>% 
  mutate(correct = ifelse(TRUE. == Predicted, 'yes', 'no'))

# classification error
table(rf_out_sample$correct)[[1]]/nrow(rf_out_sample)
table(rf_out_subject$correct)[[1]]/nrow(rf_out_subject)

# balanced error 
rf_conf_mat <- 
  rf_out_sample %>% # change to sample
  group_by(TRUE., correct) %>% summarise(n = n()) %>%
  as.data.frame()

0.5 * ((rf_conf_mat[1,3]/(rf_conf_mat[1,3] + rf_conf_mat[2,3])) +
  (rf_conf_mat[3,3]/(rf_conf_mat[3,3] + rf_conf_mat[4,3])))

#
#
#

# Longitudinal features without temporal awareness
data_lf <- data %>% pivot_wider(names_from = CoreVisit, values_from = !c(CoreVisit, Diagnosis, RAD_ID))
assay_lf <- data_lf[,-c(1:2)]

# PLS-DA (with NIPALS imputation)
lf_plsda <- plsda(X = assay_lf, Y = data_lf$Diagnosis, ncomp = 1, scale = TRUE)

perf_lf_plsda <- perf(lf_plsda, 
                        validation = 'Mfold',
                        folds = 6,
                        nrepet = 100)

perf_lf_plsda[["error.rate"]]
perf_lf_plsda[["error.rate.class"]]

dev.new()
plotLoadings(lf_plsda, ndisplay = 15)

# Multi-task learning
devtools::install_github("nguforche/MEml")
install_github("uyedaj/bayou")
library(MEml)
library(Formula)
library(gbm)
library(plyr)
library(caret)
library(lme4)
library(RRF)
library(inTrees)
library(ggplot2)
library(parallel)
library(class)
library(Hmisc)
library(kernlab)
require(flexmix)
require(gplots)
require(bayou)
require(PresenceAbsence)
require(DMwR)

data(heart.valve)
dat <- heart.valve

data$RAD_ID <- as.factor(data$RAD_ID)
data$CoreVisit <- as.numeric(as.factor(data$CoreVisit))
data$Diagnosis <- as.factor(data$Diagnosis)
data <- data %>% mutate(Diagnosis = ifelse(Diagnosis == 1, 0, 1))
resp.vars <- 'Diagnosis'
id <- 'RAD_ID'

seed = 123
para <- list(
  method = "cv", # internal cross-validation method for parameter tuning. See caret package
  tuneLength = 3, # grid size for parameter search 
  number = 3,  # number of internal cross-validation
  #n.trees = 100,   # number of trees in gbm 
  ntree = 400,   # number of trees in random forest
  interaction.depth = 5,
  shrinkage = 0.05,
  n.minobsinnode = 10,
  opt.para = TRUE, # perform parameter tuning through internal cross-validation 
  include.RE = TRUE,  ## to include estimated random effect as a predictor in the machine learning model
  max.iter = 10, ## maximum number of iterations for the "expectation maximization" like step  
  alpha = 0.05, 
  minsize = 20,
  maxdepth = 30,
  family = "binomial", 
  glmer.Control = glmerControl(optimizer = "bobyqa"), #glmerControl 
  likelihoodCheck = TRUE, 
  nAGQ = 0, 
  decay = 0.05, 
  K = 3, 
  tol = 1e-5,
  seed = seed
)

rhs.vars <- c(colnames(assay)) 
rand.vars  = 'CoreVisit'

data$dummy <- rep(1, 87)
data$dummy <- as.character(data$dummy)

# cannot set a random intercept for RAD_ID?
res <- MEml2(method = "MErf", 
             data = data, 
             id = 'RAD_ID',  
             resp.vars = resp.vars, 
             rhs.vars = rhs.vars,
             rand.vars = rand.vars, 
             para = para)



# does not accept factor as a response
res1 <- MErf(Diagnosis ~ Methylthioacetate + Methyl_thiocyanate, 
             dat = data,
             groups = id, 
             family = 'binomial',
             para = para)


