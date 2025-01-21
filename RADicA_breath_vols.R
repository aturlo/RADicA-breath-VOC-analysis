# RADicA check breath volumes

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("rhdf5")

install.packages("hdf5r")

devtools::install_github("Azure/Microsoft365R")

library(rhdf5)
library(hdf5r)
library(Microsoft365R)
library(dplyr)
library(ggplot2)
library(stringr)
library(hms)
library(lubridate)

# load data from sharepoint
site <- get_sharepoint_site(site_url = 'https://livemanchesterac-my.sharepoint.com/personal/waqar_ahmed_manchester_ac_uk/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fwaqar%5Fahmed%5Fmanchester%5Fac%5Fuk%2FDocuments%2FFiles%2F2Projects%2FRADICA%2FRadica%20asthma%20vs%20not%20asthma%2FRadica%20VOCs%20shared%20folder&sortField=LinkFilename&isAscending=true&fromShare=true&ga=1')

drive <- get_business_onedrive()

# load collection dates from RADicA dataset
meta <- read.csv('Radica sample filenames aligned with clinical metadata.csv')[,-1]

b_meta <- read.csv('RaDICA Breath sample metadata.csv')

b2_all <- read.csv('RADicA_B2.csv')[,-1] %>% filter(class == 'S1') %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, Diagnosis),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

b1_all <- read.csv('RADicA_B1.csv')[,-1] %>% filter(class == 'S1') %>%
  left_join(meta %>% dplyr::select(Sample_ID, CoreVisit, Diagnosis),
            by = c('Sample' = 'Sample_ID')) %>%
  filter(CoreVisit %in% c('CV1', 'CV2'))

b_all <- rbind(b1_all, b2_all) %>%
  relocate(CoreVisit, Diagnosis)

# link with breath metadata
b_meta_date <- unlist(str_split(b_meta$Breath.sampling_date.time, ' ')) 
b_meta_date_df <- as.data.frame(matrix(b_meta_date , ncol = 2, byrow = T)) 
colnames(b_meta_date_df) <- c('Date', 'Time')
b_meta <- b_meta %>% mutate(CoreVisit = paste('CV', Visit.number, sep = ''),
                            RAD_ID = Patient.ID,
                            Collection_date = b_meta_date_df$Date,
                            Collection_time = b_meta_date_df$Time)

ids <- unlist(str_split(b_all$Sample, '_'))
ids_df <- as.data.frame(matrix(ids ,ncol = 6, byrow = T)) 
ids_df1 <- ids_df[1:104,]
ids_df2 <- ids_df[105:nrow(ids_df),]

rad_ids <- c(ids_df1[,3], ids_df2[,5])
rad_ids

b_all <- b_all %>% mutate(RAD_ID = rad_ids)  %>% relocate(RAD_ID)

b_all1 <- b_all %>% left_join(b_meta %>% 
                               dplyr::select(RAD_ID, CoreVisit, Collection_date, Collection_time)) %>%
  relocate(Collection_date, Collection_time)

b_all1$Collection_date <- as.Date(b_all1$Collection_date, format = '%d/%m/%Y')
b_all1$Collection_time <- hm(b_all1$Collection_time)

b_all1 <- b_all1 %>% dplyr::select(!c(Date, Time))

#
#
#

# save all the HDF5 file names in directory to a vector
listh5 <- dir(pattern = "*.h5")

dates <- sapply(listh5, function(fileN){
  date <- str_sub(fileN, start = 9, end = -8) 
})

dates <- as.Date(dates, format = '%d%m%Y')

setdiff(b_all1$Collection_date, dates)

datesInc <- dates[dates %in% b_all1$Collection_date] 

listh5 <- listh5[listh5 %in% names(datesInc)]

head(listh5)

h5ls("BCSdata_020720211235.h5")
file <- 'BCSdata_140220201434.h5'

# save breath volume and collection time for each HDF5 file
output <- bind_rows(lapply(listh5, function(file){
  
  file_name <- file
  groups <- h5ls(file_name)

  data <- h5read(file_name, 'Data') %>% as.data.frame.table() #%>%
    dplyr::select('Freq.Collection.time', 'Freq.Breathing.rate')
  data <- data %>% filter(Freq.Breathing.rate > 0)
  
  if ('Collection_info' %in% groups$name){
  info <- h5readAttributes(file_name, 'Collection_info')
  coll_info <-  data.frame(Vol_tube_L = info[['Collection per tube L']],
                           Vol_tube_R = info[['Collection per tube R']],
                           Coll_time = info[['Total collection time']],
                           dateTime = str_sub(file_name, start = 9, end = -4))
  
  coll_info <- coll_info %>% mutate(Sample = ifelse(Vol_tube_L > 0 &
                                                      Vol_tube_R > 0, 
                                                    'Breath', 'Blank'))
  
  coll_info <- coll_info %>% mutate(
    Time_breath = ifelse(
    Sample == 'Blank', Coll_time,
    Coll_time - data$Freq.Collection.time[1]),
    Time_no_breath = ifelse(Sample == 'Blank', 0,
                         Coll_time - Time_breath),
    Breaths = ifelse(Sample == 'Blank', 0,
                     mean(data$Freq.Breathing.rate[data$Freq.Breathing.rate > 0])*(Time_breath/60)))
  } else{
    coll_info <- data.frame(Vol_tube_L = NA, Vol_tube_R = NA, Coll_time = NA, date = NA,
                            Sample = NA, Time_no_breath = NA, Time_breath = NA, Breaths = NA)
  }
  
  }))

#

write.csv(output, 'ReCIVA_log_file_collection_info')

# distribution of breath sample volumes
output_b <- output %>% filter(Sample == 'Breath') %>%
  mutate(Collection_date = str_sub(dateTime, end = -5),
         Collection_time = paste(str_sub(dateTime, start = 9, end = -3), 
         str_sub(dateTime, start = 11), sep = ':'))

output_b$Collection_date <- as.Date(output_b$Collection_date, format = '%d%m%Y') 
output_b$Collection_time <- hm(output_b$Collection_time)

View(table(output_b$Collection_date))

# visualise distribution of breath sample volumes

tiff('Breath_sample_volumes.tiff', res = 300, unit = 'mm', width = 80, height = 130)

par(mfrow = c(2, 1)) # 21-by-2 grid of plots
par(oma = c(3, 3, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(1, 1, 4, 1))

barplot(table(round(output_b$Vol_tube_L,0)), main = 'Tube L', cex.axis = 0.6, cex.names = 0.6)
barplot(table(round(output_b$Vol_tube_R,0)), main = 'Tube R', cex.axis = 0.6, cex.names = 0.6)

mtext('Breath sample volume (mL)', side = 1, outer = TRUE, line = 2)
mtext('Count', side = 2, outer = TRUE, line = 2)

dev.off()

#
#
#

# filter collection information from HDF5 files by similarity of date and time with b_metadata
# match by date and hour?
b_all1$Collection_hour <- b_all1$Collection_time %>% hour()
b_all1 <- b_all1 %>% relocate(Collection_hour)

output_b$Collection_hour <- output_b$Collection_time %>% hour()

match_date_times <- output_b %>% 
  left_join(b_all1 %>% dplyr::select(Collection_date, Collection_hour,
                                     RAD_ID, CoreVisit))

table(is.na(match_date_times$RAD_ID))

# 196 samples matched between VOC and ReCiVA data
# check if samples < 500mL matched
# check if samples < 500mL had repeated measurement
rep_samp <- match_date_times %>% filter(Vol_tube_L < 500 | Vol_tube_R < 500)

rep_samp_date <- match_date_times %>% filter(Collection_date %in% rep_samp$Collection_date &
                                       Collection_hour %in% rep_samp$Collection_hour)

match_date_times1 <- match_date_times %>% dplyr::select(!date) %>% drop_na() 
match_date_times1 <- match_date_times1[-c(13, 96),]

match_date_times1$Sample = paste(match_date_times1$RAD_ID, match_date_times1$CoreVisit,
                                sep = '_')

setdiff(b_corr_w1$Sample, match_date_times1$Sample)
setdiff(b1_corr_w1$Sample, match_date_times1$Sample)


# visualise distribution of breath sample volumes
tiff('Breath_sample_volumes_VOC.tiff', res = 300, unit = 'mm', width = 80, height = 130)

par(mfrow = c(2, 1)) # 21-by-2 grid of plots
par(oma = c(3, 3, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar = c(1, 1, 4, 1))

barplot(table(round(match_date_times1$Vol_tube_L,0)), 
        main = 'Tube L', cex.axis = 0.6, cex.names = 0.6,
        space = 1)
barplot(table(round(match_date_times1$Vol_tube_R,0)), 
        main = 'Tube R', cex.axis = 0.6, cex.names = 0.6,
        width = c(0.1, 0.1, 0.1),
        space = 3)

mtext('Breath sample volume (mL)', side = 1, outer = TRUE, line = 2)
mtext('Count', side = 2, outer = TRUE, line = 2)

dev.off()

#

tiff('Breath_collection_times.tiff', res = 300 , unit = 'mm', width = 100, height = 80)
hist(match_date_times1$Time_breath/60, main = 'Breath collection time', 
     xlab = 'Time [min]', cex.axis = 0.8, breaks = 12)
dev.off()

tiff('Breath_collection_breaths.tiff', res = 300 , unit = 'mm', width = 100, height = 80)
hist(match_date_times1$Breaths, main = 'Total breaths collected', 
     xlab = 'Number of breaths', cex.axis = 0.8, breaks = 12)
dev.off()

quantile(output_b$Time_breath/60, na.rm = TRUE)
quantile(output_b$Breaths, na.rm = TRUE)

tiff('Breath_collection_timeVSbreaths.tiff', res = 300, unit = 'mm', width = 100, height = 80)
par(mar = c(4, 4, 1, 1))
plot(match_date_times1$Time_breath/60, match_date_times1$Breaths, xlab = 'Time [min]', ylab = 'Number of Breaths')
dev.off()

# Breath volume estimation?
match_date_times1$Mean_breath_vol <- (match_date_times1$Vol_tube_L + match_date_times1$Vol_tube_R)/match_date_times1$Breaths

tiff('Breath_collection_mean_breath_vol.tiff', res = 300, unit = 'mm', width = 100, height = 70)
par(mar = c(4, 4, 1, 1))
hist(match_date_times1$Mean_breath_vol,
     xlab = 'Mean breath volume (mL)', main = '')
dev.off()


#
match_date_times1 <- match_date_times1 %>% left_join(clin_dyn %>% dplyr::select(RAD_ID, CoreVisit, FVCPre,
                                                                                FEV1Pre))

tiff('Breath_collection_FVCvsBreathVol.tiff', res = 300, unit = 'mm', width = 100, height = 80)
par(mar = c(4, 4, 1, 1))
plot(match_date_times1$Mean_breath_vol, match_date_times1$FEV1Pre,
     xlab = 'Mean breath volume (mL)', ylab = 'FVCPre (L)')
dev.off()

cor(match_date_times1$Mean_breath_vol, match_date_times1$FEV1Pre, use = 'complete.obs')
cor(outl_rem$Breaths, outl_rem$FVCPre, use = 'complete.obs')

data %>% ggplot(aes(x = Freq.Collection.time, y = Freq.CO2stream)) + geom_line() + theme_bw() +
  xlab('Collection time (s)') + ylab('CO2 stream') +
  ggtitle('Recording of CO2 stream in one sample (010220211109)')
ggsave('CO2_stream_example.tiff', unit = 'mm', dpi = 300, width = 280, height = 100)

#
#
#

# 3-methylpentane
b_corr_w1 %>% ggplot(aes(x = CoreVisit, y = X3_methylpentane, fill = Diagnosis)) +
  geom_violin() + geom_boxplot()

b1_corr_w1 %>% ggplot(aes(x = CoreVisit, y = X3_methylpentane, fill = Diagnosis)) +
  geom_violin() + geom_boxplot()


