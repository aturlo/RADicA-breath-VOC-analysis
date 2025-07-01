## Data cleaning and subsetting

# author: Aggie Turlo
# project: RADicA
# date: 30/06/2025

#####################

library(tidyverse)

## load data
# load GC-MS peak area file
data <- read.csv('RADicA_VOC_raw_peak_data.csv')

# load patient metadata
meta <- read.csv('RADicA_VOC_metadata.csv')

# custom functions
'%ni%' <- Negate('%in%')

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


## save outputs
write.csv(b1_all, 'RADicA_B1.csv')
write.csv(b2_all, 'RADicA_B2.csv')

#
#
#
