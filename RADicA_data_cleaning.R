## Data cleaning and subsetting

# author: Aggie Turlo
# project: RADicA
# date: 16/01/2025

#####################

library(dplyr)
library(tidyr)

## load data
# load GC-MS peak area file
data <- read.csv('240508_Radica_VOC_peak_areas_v5.1.csv', check.names = FALSE)

# load patient metadata
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

## replace spaces in compound names (columns) with underscores
colnames(data)[191] <- 'Beta phellandrene'
colnames(data) <- gsub(' ', '_', colnames(data))
colnames(data) <- gsub('-', '_', colnames(data))

## assign technical sample classes
# Blank - empty sorbent tube
# ES - external standards (empty sorbent tube with spiked-in standards)
# BG - room air sample
# S1 - breath sample replicate
# S2 - breath sample replicate
b_all <- data %>% 
  mutate(class = ifelse(grepl('tracker|Tr', Sample),'Blank',
                        ifelse(grepl('QC', Sample), 'ES',
                               ifelse(grepl('B1', Sample), 'BG',
                                      ifelse(grepl('S1', Sample), 'S1', 'S2'))))) %>%
  relocate(class)

b_all[b_all == 0] <- NA

## subset data into data collection campaigns (pre- and post-covid)
b_all <- b_all %>% separate(Acq_Date_Time, c('Date', 'Time'), sep = ' ')
b_all$Date <- as.Date(b_all$Date, format = "%d/%m/%Y")

b1_all <- b_all %>% filter(Batch %in% c(1,2,3,4))
b2_all <- b_all %>% filter(Batch %in% c(5, 6))

# remove one of each duplicated sample pairs (repeated due to analysis failure)
duplicates_b1 <- c('200320_RaDICA_RAD094_B1_434919_1', '200320_RaDICA_RAD074_B1_609366_1')

b1_all <- b1_all %>% filter(Sample %ni% duplicates_b1)

duplicates_b2 <- c('220318_24_RADicA_S1_RAD226_367150', '200417_RaDICA_RAD114_S2_367150_1')

b2_all <- b2_all %>% filter(Sample %ni% duplicates_b2)

## save outputs
write.csv(b1_all, 'RADicA_B1.csv')
write.csv(b2_all, 'RADicA_B2.csv')

#
#
#
