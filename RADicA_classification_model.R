## Predicitng asthma based on breath VOC samples

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(gridExtra)
library(vsn)
library(mixOmics)

# load summarised BG corrected dataset 
b1_imp_sum_corr_w <- read.csv('RADicA_BG_adjusted.csv')[,-1]

# w/o multivariate outliers
b1_corr_w1 <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b_corr_w1 <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]

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

# MINT
X <- b_corr_w1[,-c(1:4)]
Y <- b_corr_w1$Diagnosis
study <- b_corr_w1$CoreVisit

study[study == 'CV1'] <- 1
study[study == 'CV2'] <- 2

study <- as.factor(study)

base <- mint.plsda(X, Y, study = study, ncomp = 2, scale = TRUE)
base$prop_expl_var$X

# plot results
plotIndiv(base, legend = TRUE,
          study = 'global',
          title = 'MINT')

plotVar(base, cutoff = 0.6)

plotLoadings(base, ndisplay = 15, comp = 1, 
             study = 'global', 
             method = 'median',
             contrib = 'max')

# sparse MINT
sparse <- mint.splsda(X, Y, study = study,
                      keepX = c(5, 5))
sparse$prop_expl_var$X

plotIndiv(sparse, legend = TRUE,
          study = 'all.partial')
plotVar(sparse, cutoff = 0.6)

plotLoadings(sparse, ndisplay = 15, comp = 1, 
             study = 'global', 
             method = 'median',
             contrib = 'max')

# tuning
base_opt <- perf(base, dist = 'all')
plot(base_opt)

base_opt$global.error$BER
base_opt$global.error$error.rate.class

# compare with plsda
plsda <- plsda(X = X, Y = Y, ncomp = 2)
plotIndiv(plsda, pch = study, title = 'PLSDA')

plsda_opt <- perf(plsda, validation = 'Mfold', folds = 10, nrepeat = 100)
plot(plsda_opt)

plsda_opt$error.rate$BER
plsda_opt$error.rate.class

# tuning number of variables
# not working
base_tune <- tune(X = X, Y = Y, study = as.factor(study),
                  ncomp = 2,
                  test.keepX = seq(1,50,1),
                  method = 'mint.splsda',
                  measure = 'BER',
                  dist = 'centroids.dist')

plot(base_tune)
base_tune$choice.keepX

#
loads_glob <- base$loadings$X %>% as.data.frame() %>%
  mutate(comp = rownames(.))

loads_part <- base$loadings.partial$X %>% as.data.frame() %>%
  mutate(comp = rownames(.)) %>%
  rename(CV1_comp1 = X1.comp1,
         CV1_comp2 = X1.comp2,
         CV2_comp1 = X2.comp1,
         CV2_comp2 = X2.comp2)

#
loads_glob %>% 
  left_join(reg_valid_join %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>%
  arrange(desc(abs(comp1))) %>% 
  slice(1:20) %>% 
  ggplot(aes(x = comp1, y = fct_inorder(as.factor(comp)), fill = FC_filter)) + geom_col()

#
dev.new()
loads_part %>%
  left_join(reg_valid_join %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>%
  arrange(desc(abs(CV2_comp1))) %>% 
  slice(1:20) %>% 
  ggplot(aes(x = CV2_comp1, y = fct_inorder(as.factor(comp)), fill = FC_filter)) + geom_col()







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

data <- b_corr_w1

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

rda_samples_plot <- rda_scores %>% 
  ggplot(aes(x = RDA1, y = RDA2, colour = Diagnosis, shape = CoreVisit)) + 
  geom_point(size = 2) + theme_bw() +
  xlab(paste('RDA1:', rda_var_expl[1], '%')) +
  ylab(paste('RDA2:', rda_var_expl[2], '%')) +       
  scale_colour_manual(values = c('Asthma' = 'dodgerblue3',
                                 'Not Asthma' = 'darkorange1')) +
  scale_shape_manual(values = c('CV1' = 1,
                                'CV2' = 2)) +
  ggtitle('Redundancy Analysis sample score plot') +
  theme(plot.title = element_text(hjust = 0.5))

rda_samples_plot

ggsave('RDA_sample_scores_plot.tiff', units = 'mm', dpi = 300, width = 120, height = 80)

#
rda_loadings_voc <- scores(rda, choices = 1:3, display = 'species') %>%
  as.data.frame() %>%
  mutate(comp = rownames(.))

rda1 <- rda_loadings_voc %>%
  left_join(reg_valid_join %>% dplyr::select(comp, FC_filter)) %>%
  distinct() %>%
  arrange(desc(abs(RDA2))) %>%
  slice(1:15) %>%
  ggplot(aes(x = RDA2, y = fct_inorder(as.factor(comp)),
             fill = FC_filter)) + 
  geom_col() +
  theme_bw() +
  #ggtitle('Diagnosis') +
  ggtitle('Core Visit') +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank()) 

rda1

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

input <- b1_imp_sum_corr_w1 %>% dplyr::select(!Sample)

input1 <- input %>% mutate(Diagnosis = ifelse(Diagnosis  == 'Asthma', 1, 0))

input1 <- input1 %>% dplyr::select(!RAD_ID)

dup_id <- as.data.frame(table(input$RAD_ID)) %>% filter(Freq > 1)

input1 <- input %>% filter(RAD_ID %in% dup_id$Var1)

forest <- bimm_fit(data = input1,
                   formula = Diagnosis ~ . + (1 | CoreVisit:RAD_ID),
                   fun_model_ml = function (formula, data_train) 
                   {
                     ranger::ranger(formula = formula, 
                                    data = data_train, 
                                    classification = TRUE)
                   },
                   
                   n_iteration = 1,
                   verbose = TRUE)
summary(forest)

View(as.data.frame(forest[["model_ml"]][["variable.importance"]]))




# RF with ranger without random effect
library(ranger)
base_forest <- ranger(Diagnosis ~., 
                      data = input1 %>% dplyr::select(!c(CoreVisit, RAD_ID)),
                      probability = TRUE)

#
#
#

# MOFA



