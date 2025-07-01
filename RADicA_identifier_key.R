## Study sample name conversion

# author: Aggie Turlo
# project: RADicA
# date: 30/06/2025

#####################

library(tidyverse)

## load data
# load GC-MS peak area file
data <- read.csv('240508_Radica_VOC_peak_areas_v5.1.csv', check.names = FALSE)

# load patient metadata
meta <- read.csv('Metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

#

## replace spaces in compound names (columns) with underscores
colnames(data)[191] <- 'Beta phellandrene'
colnames(data) <- gsub(' ', '_', colnames(data))
colnames(data) <- gsub('-', '_', colnames(data))

# encode RAD IDs as consecutive integers
rad_ids <- data.frame(RAD_ID = unique(meta$RAD_ID)) %>%
  mutate(ID = sapply(1:length(unique(meta$RAD_ID)), function(id){
    paste('ID', id, sep = '')
  }))

# add new identifiers to metadata file
meta <- meta %>% left_join(rad_ids) 

# change file names swapping RAD ID for new identifier
split_ids <- str_split(meta$Sample_ID, pattern = '(?=RAD)', simplify = TRUE) %>%
  as.data.frame()
split_ids <- split_ids %>% cbind(meta %>% dplyr::select(ID, Sample, CoreVisit))
split_ids <- split_ids %>% 
  mutate(Sample_ID1 = ifelse(V3 == '', 
                         paste(V1, ID, str_sub(V2, start = 7), sep =''),
                         paste(V1, V2, ID, str_sub(V3, start = 7), sep = '')))
meta <- meta %>% left_join(split_ids %>% dplyr::select(ID, Sample_ID1, Sample, CoreVisit))

# remove one of each duplicated sample pairs (repeated due to analysis failure)
duplicates_b1 <- c('200320_RaDICA_RAD094_B1_434919_1', '200320_RaDICA_RAD074_B1_609366_1')

data <- data %>% filter(Sample %ni% duplicates_b1)

duplicates_b2 <- c('220318_24_RADicA_S1_RAD226_367150', '200417_RaDICA_RAD114_S2_367150_1')

data <- data %>% filter(Sample %ni% duplicates_b2)

# change file names in VOC peak area dataset
data <- data %>% dplyr::rename(Sample_ID = Sample) %>%
  left_join(meta %>% dplyr::select(Sample_ID, Sample_ID1))%>%
  relocate(Sample_ID1)

data <- data %>% mutate(Sample_ID1 = ifelse(is.na(Sample_ID1) == TRUE, Sample_ID, Sample_ID1)) 

# remove study sample names and RAD ID from analysis input files
data <- data %>% dplyr::select(!Sample_ID) %>%
  dplyr::rename(Sample = Sample_ID1)

meta <- meta %>% dplyr::select(!c(RAD_ID, Sample_ID)) 

meta <- meta %>% dplyr::rename(Sample_ID = Sample_ID1)

# remove RAD IDs from blank sample file names
split_snames <- str_split(data$Sample, pattern = '(?=RAD)', simplify = TRUE) %>%
  as.data.frame() %>% cbind(data %>% dplyr::select(Sample))
split_snames <- split_snames %>% mutate(Sample1  = ifelse(V2 == '', V1,
                                                          paste(V1, V2, str_sub(V3, start = 7), sep ='')))

data <- data %>% left_join(split_snames %>% dplyr::select(Sample, Sample1)) %>%
  relocate(Sample1)

data <- data %>% dplyr::select(!Sample) %>%
  rename(Sample = Sample1)

# save analysis input files
write.csv(data, 'RADicA_VOC_raw_peak_data.csv')
write.csv(meta, 'RADicA_VOC_metadata.csv')

#
#
#

# replace RAD identifiers with new IDs in backgrund-corrected dataset w/o outliers
b1_corr_out <- read.csv('RADicA_BG_adjusted_B1_outl_removed.csv')[,-1]
b2_corr_out <- read.csv('RADicA_BG_adjusted_B2_outl_removed.csv')[,-1]

key <- meta %>% dplyr::select(RAD_ID, ID, CoreVisit)

b1_corr_out <- b1_corr_out %>% left_join(key) %>%
  relocate(ID)
b2_corr_out <- b2_corr_out %>% left_join(key) %>%
  relocate(ID)

b1_corr_out <- b1_corr_out %>% mutate(Sample = paste(ID, CoreVisit, sep ='_')) %>%
  dplyr::select(!RAD_ID) 
b2_corr_out <- b2_corr_out %>% mutate(Sample = paste(ID, CoreVisit, sep ='_')) %>%
  dplyr::select(!RAD_ID)

write.csv(b1_corr_out, 'RADicA_BG_adjusted_B1_outl_removed.csv')
write.csv(b2_corr_out, 'RADicA_BG_adjusted_B2_outl_removed.csv')

#
#
#
