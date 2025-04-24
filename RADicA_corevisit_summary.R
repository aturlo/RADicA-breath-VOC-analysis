## RADicA CoreVisits per patients summary

b1 <- read.csv('RADicA_B1_NAfiltered_imputed_CC2_PQN.csv')[,-1]
b2 <-  read.csv('RADicA_B2_NAfiltered_imputed_CC2_PQN.csv')[,-1]

meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

#

b1 <- b1 %>% filter(class == 'S1')
b2 <- b2 %>% filter(class == 'S1')

b1 <- b1 %>% left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, RAD_ID),
                       by = c('Sample' = 'Sample_ID'))

b2 <- b2 %>% left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, RAD_ID),
                       by = c('Sample' = 'Sample_ID')) 

#

sum_cv_b1 <- as.data.frame(table(b1$RAD_ID))

b1_cvs <- b1 %>% dplyr::select(RAD_ID, CoreVisit) %>%
  mutate(n = rep(1, nrow(.))) %>%
  pivot_wider(names_from = CoreVisit, values_from = n) 

b1_cvs[is.na(b1_cvs) == TRUE] <- 0

b1_cvs <- b1_cvs %>% mutate(pre_post = ifelse((CV1+CV2) >= 1 & (CV3+CV4) >= 1, 'Yes', 'No'))
  
#

sum_cv_b2 <- as.data.frame(table(b2$RAD_ID))

b2_cvs <- b2 %>% dplyr::select(RAD_ID, CoreVisit) %>%
  mutate(n = rep(1, nrow(.))) %>%
  pivot_wider(names_from = CoreVisit, values_from = n) 

b2_cvs[is.na(b2_cvs) == TRUE] <- 0

b2_cvs <- b2_cvs %>% mutate(pre_post = ifelse((CV1+CV2) >= 1 & (CV3+CV4) >= 1, 'Yes', 'No'))

#

b_cvs <- b1_cvs %>% filter(pre_post == 'Yes') %>%
  rbind(b2_cvs %>% filter(pre_post == 'Yes'))

b_cvs_ids <- b_cvs %>% dplyr::select(RAD_ID)

write.csv(b_cvs_ids, 'RADicA_RAD_IDs_pre_post_treatment.csv')
