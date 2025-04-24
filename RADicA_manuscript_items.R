## RADicA manuscript

library(dplyr)
library(tidyr)
library(foreign)
library(stringr)

'%ni%' <- Negate('%in%')
  
# Table 1
# load VOC and clinical data
b1 <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]
b2 <- read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN_summarised_BGfiltered.csv')[,-1]

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')
clin <- read.spss('RADicA Active 29052024.sav', to.data.frame = TRUE)

# format datasets
# keep clinical outcomes of interest
clin_core <- clin %>% dplyr::select(BBID, Ethnicity, MaleFemale, AgeInYears, SmokerPQ, 
                                    SmokerTypePQ, BMICV1, FeNOCV1, FeNOCV2, 
                                    FEV1PreCV1, FEV1PostCV1, FEV1PPrePredCV1, 
                                    FEV1PPostPredCV1, PD20CV2) 

# format clinical data 
clin_core$BBID <- str_sub(clin_core$BBID, end = -2)
colnames(clin_core)[1] <- 'RAD_ID'

clin_core <- clin_core %>% mutate(BDR = ((FEV1PostCV1-FEV1PreCV1)/FEV1PreCV1)*100,
                                BDRPPre = FEV1PPostPredCV1-FEV1PPrePredCV1)

# keep only patients with breath samples
clin_core <- clin_core %>% filter(RAD_ID %in% c(b1$RAD_ID, b2$RAD_ID))

# check for patients overlapping between datasets and remove from B2
b2 <- b2 %>% filter(RAD_ID %ni% c(intersect(b1$RAD_ID, b2$RAD_ID)))

#

# join VOC and clinical data
b1 <- b1 %>% dplyr::select(RAD_ID, CoreVisit) %>%
  left_join(meta %>% dplyr::select(RAD_ID, Diagnosis, age)) %>% distinct() 

b2 <- b2 %>% dplyr::select(RAD_ID, CoreVisit) %>%
  left_join(meta %>% dplyr::select(RAD_ID, Diagnosis, age))  %>% distinct()

#
b1_sum <- b1 %>% group_by(RAD_ID, Diagnosis) %>% summarise(n = n())

table(b1_sum$Diagnosis)
table(b1$CoreVisit)


b2_sum <- b2 %>% group_by(RAD_ID, Diagnosis) %>% summarise(n = n())

table(b2_sum$Diagnosis)
table(b2$CoreVisit)

b1 %>% group_by(Diagnosis, CoreVisit) %>% summarise(n = n())
b2 %>% group_by(Diagnosis, CoreVisit) %>% summarise(n = n())      

#
b1_sum <- b1_sum %>% left_join(clin_core)
b2_sum <- b2_sum %>% left_join(clin_core)


# Age
hist(b1_sum$AgeInYears, breaks = 5)
quantile(b1_sum$AgeInYears)

hist(b2_sum$AgeInYears, breaks = 5)
quantile(b2_sum$AgeInYears)

b1a <- b1_sum %>% filter(Diagnosis == 'Asthma')
b1na <- b1_sum %>% filter(Diagnosis == 'Not Asthma')

b2a <- b2_sum %>% filter(Diagnosis == 'Asthma')
b2na <- b2_sum %>% filter(Diagnosis == 'Not Asthma')

quantile(b1a$AgeInYears)
quantile(b1na$AgeInYears)

quantile(b2a$AgeInYears)
quantile(b2na$AgeInYears)

# smoking
table(b1_sum$SmokerPQ)
table(b1na$SmokerPQ)
table(b1a$SmokerPQ)

table(b2_sum$SmokerPQ)
table(b2na$SmokerPQ)
table(b2a$SmokerPQ)

#

table(b1_sum$SmokerTypePQ)
table(b1na$SmokerTypePQ)
table(b1a$SmokerTypePQ)

table(b2_sum$SmokerTypePQ)
table(b2na$SmokerTypePQ)
table(b2a$SmokerTypePQ)

# gender
table(b1_sum$MaleFemale)
table(b1a$MaleFemale)
table(b1na$MaleFemale)

table(b2_sum$MaleFemale)
table(b2a$MaleFemale)
table(b2na$MaleFemale)

# BMI
hist(b1_sum$BMICV1, breaks = 5)
quantile(b1_sum$BMICV1)
quantile(b1a$BMICV1)
quantile(b1na$BMICV1)

hist(b2_sum$BMICV1, breaks = 5)
quantile(b2_sum$BMICV1)
quantile(b2a$BMICV1)
quantile(b2na$BMICV1)

# asthma diagnostic tests
b1_sum$PD20CV2[b1_sum$PD20CV2 == '<0.015   '] <- 0.015
b2_sum$PD20CV2[b2_sum$PD20CV2 == '<0.015   '] <- 0.015

b1_sum$PD20CV2 <- as.numeric(b1_sum$PD20CV2)
b2_sum$PD20CV2 <- as.numeric(b2_sum$PD20CV2)

hist(b2_sum$FeNOCV1)
hist(b2_sum$BDR)
hist(b2_sum$PD20CV2)

# FeNO
quantile(b1_sum$FeNOCV1)
quantile(b1a$FeNOCV1)
quantile(b1na$FeNOCV1)

quantile(b2_sum$FeNOCV1, na.rm = TRUE)
quantile(b2a$FeNOCV1, na.rm = TRUE)
quantile(b2na$FeNOCV1)

# BDR
quantile(b1_sum$BDR)
quantile(b1a$BDR)
quantile(b1na$BDR)

quantile(b2_sum$BDR)
quantile(b2a$BDR)
quantile(b2na$BDR)

# BCT
quantile(b1_sum$PD20CV2, na.rm = TRUE)
quantile(b1a$PD20CV2, na.rm = TRUE)
quantile(b1na$PD20CV2, na.rm = TRUE)

quantile(b2_sum$PD20CV2, na.rm = TRUE)
quantile(b2a$PD20CV2, na.rm = TRUE)
quantile(b2na$PD20CV2, na.rm = TRUE)

# avaiability of test results
table(is.na(b2_sum$FeNOCV1) == TRUE)
table(is.na(b2_sum$PD20CV2) == TRUE)

table(is.na(b2na$PD20CV2) == TRUE)

# smoking

