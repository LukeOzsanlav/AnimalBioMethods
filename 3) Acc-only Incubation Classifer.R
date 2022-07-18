##Luke Ozsanlav-Harris
## PhD incubation methods chapter

## Code initially adapted from scripts that identified incubation events from Acc and GPS data
## Turned those scripts into a loop that allowed sampling rates of Acc and GPS data to change and recalculate incubation lengths
## Finally removed the capacity for GPS data so only uses Acc data. This is for a comparison that comprises part of a methods chapter



## Created: 20/05/2020

## Updated: 05/05/2020
## Changes made:
## Updates was to include nsd in the classification process


## Updated: 09/03/2021
## Changes made:
## Need to automatically pick the lowest 24 days of ODBA then label the incubating, buffer and non-incubating points
## Added capacity to re-sample ODBA and loop through values
## Added capacity to re-sample GPS data and loop through values
## Added capacity to record the 95% quantities values so can look how they change with re-sampling
## Added capacity to change the length of the acc burst


## Updated: 07/06/2021
## Changes made:
## Adapted the script so that it just used the Accelerometer data (allow comparison with Acc + GPS approach)

## Updates: 13/06/2022
## Changed where the main outputs were saved to impovre the wrokflow





##loading Packages
library(amt)
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(geosphere)
library(DescTools)





##################
## ODBA prep 1. ##
################## Importing and preparing Ornitela ODBA data

## set working directory
setwd("~/Ornitela data/ALL_ORNI_DATA/Cleaned CSVs")

## read in data
orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))

## sub-setting out all of the acc data
orn_acc <- filter(orn_data, datatype == "SENSORS")
rm(orn_data)

## reduce data set to the columns required from here on
orn_acc <- subset(orn_acc, select = c("device_id", "UTC_datetime", "acc_x", "acc_y", "acc_z"))

## parse timestamp
#orn_acc$UTC_datetime <- as.character(ymd_hms(orn_acc$UTC_datetime))

## needed to prep data timestamp
orn_acc$UTC_datetime <- as.POSIXct(orn_acc$UTC_datetime, format= "%Y-%m-%d %H:%M:%S", tz= "GMT")
#orn_acc <- orn_acc %>% mutate(month=month(UTC_datetime), year = year(UTC_datetime))

## give it a tag year column 
orn_acc$tag_year <- paste0(orn_acc$device_id, "_", year(orn_acc$UTC_datetime))






##################
## Phen data 1. ##
################## Read in and clean the phenology data

#### trim GPS data to just Greenland period ####
## read in phenology data
setwd("~/Migration/Phenology")
Phenology <- fread("Full_phenology_fixes.csv")
Phenology$V1 <- NULL

##set a tag_year column
Phenology$tag_year <- paste(Phenology$ID, Phenology$year, sep = "_")

##set phenology columns to  give 1 day buffer before departure and after arrival date
Phenology$Green_dep <- (ymd_hms(paste0(Phenology$Green_dep, " 00:00:00"))) - 86400
Phenology$Green_arrive <- (ymd_hms(paste0(Phenology$Green_arrive, " 00:00:00"))) + 86400

##now join the two data sets
Phen <- subset(Phenology, select = c("tag_year", "Green_dep", "Green_arrive"))




##################
## ODBA Trim 1. ##
################## 

## Now trim ODBA data to Greenland period
orn_acc2 <- inner_join(orn_acc, Phen, by = "tag_year")
print("Should be equal"); length(unique(orn_acc2$tag_year)); length(unique(orn_acc$tag_year)) ## check haven't lost lots of tag_years
rm(orn_acc)

## filter out the period in Greenland
orn_Gr <- orn_acc2 %>%  filter(UTC_datetime > Green_arrive & UTC_datetime < Green_dep)
rm(orn_acc2)

## remove some of the extra rows that aren't needed anymore
orn_Gr$Green_dep <- NULL
orn_Gr$Green_arrive <- NULL




##################
## ODBA prep 2. ##
################## Labeling Acc data and checks for burst length and frequency

#### Labelling burst with unique code and by line #### 

## Set the required segment length for the acc bursts
freq<-10 # The Frequency of accelerometry data
secs<- 3 # the length of the burst in seconds
numrows<-freq*secs # the number of rows required to calculate metrics over the chosen period. 

## creating a time difference column
orn_Gr$timedif <- orn_Gr$UTC_datetime - lag(orn_Gr$UTC_datetime)

## changing all gaps between bursts to 1 and within burst gaps to 0 
#(problem: sometimes the time diffs go back in time where tags join together)
orn_Gr$timedif <- ifelse(orn_Gr$timedif < 0, (orn_Gr$timedif*-1) , orn_Gr$timedif) ##solve problem
orn_Gr$timedif <- ifelse(orn_Gr$timedif > 1, 1, 0) ##change to all 0 and 1s
orn_Gr$timedif[1] <- 1 ##first line is an NA so changing to 1 as it starts a burst

## now use cumsum to number the bursts, as there is a 1 at each burst boundary burst will be numbered increasingly in steps of 1
orn_Gr$burst_no <- cumsum(orn_Gr$timedif)



#### Check all burst are right length and frequency ####
## creating column with burst no and ID
orn_Gr$burst_no_ID <- paste0(orn_Gr$burst_no, "_", orn_Gr$device_id)
orn_Gr$dummy <- 1 ##dummy column

## creating data set with stats on each burst
Acc_data_checks <- orn_Gr %>% group_by(burst_no_ID) %>% summarise(
  Readings_in_burst = sum(dummy), end_time = max(UTC_datetime), 
  start_time =min(UTC_datetime), deviceID = min(device_id)) 
Acc_data_checks$burst_length_in_secs <- as.numeric(Acc_data_checks$end_time - Acc_data_checks$start_time) ##burst length(time)
Acc_data_checks$burst_freq <- Acc_data_checks$Readings_in_burst/Acc_data_checks$burst_length_in_secs ##burst_freq

## summary of checks
table(Acc_data_checks$Readings_in_burst) ##should all be 30
table(Acc_data_checks$burst_freq) ## should all be 10 and 15
sum(subset(Acc_data_checks, select = c("deviceID", "start_time")) %>% duplicated) ## should be 0
rm(Acc_data_checks) ## remove this from environment

## remove unwanted columns 
orn_Gr$dummy <- NULL




#############
## LOOP 1. ##
############# Set up matrix to loop through different sampling rates

## set the sampling frequencies fo the raw data
Acc_rate <- 6


## Create matrix or required sampling rates using expand.grid
sampling_matrix <- expand.grid(ODBA_sampling_rate = c(Acc_rate, Acc_rate*2, Acc_rate*4, Acc_rate*8, Acc_rate*12, Acc_rate*16, Acc_rate*20, Acc_rate*24), 
                               burst_length = c(30, 20, 10))

## Add empty column to add in training values
sampling_matrix$Inc_ODBA_q2.5 <- NA; sampling_matrix$Inc_ODBA_q97.5 <- NA; sampling_matrix$Inc_ODBA_mean <- NA 

sampling_matrix$Not_ODBA_q2.5 <- NA; sampling_matrix$Not_ODBA_q97.5 <- NA; sampling_matrix$Not_ODBA_mean <- NA




#####################
### Start of Loop ###
#####################

## loop through each round of sampling
for(B in 1:nrow(sampling_matrix)){
  
  message("Sampling round ", B, " out of ", nrow(sampling_matrix))
  
  
  
  #############
  ## ODBA 1. ##
  ############# Change the length of Acc bursts
  
  ## return the number of burst and label each line within a burst
  no_burst <- max(orn_Gr$burst_no)
  orn_Gr$burst_row <- rep(1:30, times = no_burst)
  
  ## Change the length of the burst
  if(sampling_matrix$burst_length[B] == 30){orn_acc_resamp <- orn_Gr} else{orn_acc_resamp <- orn_Gr %>% filter(burst_row <= sampling_matrix$burst_length[B])}
  
  
  ## need to convert the acc to g by dividing by 1000 (Only for Orni Data)
  orn_acc_resamp$X <- ((orn_acc_resamp$acc_x)/1000)
  orn_acc_resamp$Y <- ((orn_acc_resamp$acc_y)/1000)
  orn_acc_resamp$Z <- ((orn_acc_resamp$acc_z)/1000)
  
  ## streamline the orn_acc_resamp
  orn_acc_resamp <- subset(orn_acc_resamp, select = c("burst_no_ID", "X", "Y", "Z"))
  
  ## function to calculate ODBA per burst using data.table
  ODBA_calc_dt <- function(data){
    
    # change data to a data.table
    data <- as.data.table(data)
    
    # Minus each raw value from the satic acceleration value
    data <- data[,.(dynamic_x= abs(X - mean(X)),
                    dynamic_y= abs(Y - mean(Y)),
                    dynamic_z= abs(Z -mean(Z))),
                 by = .(burst_no_ID)]
    
    # calculate the dynamic acceleration value for each axis and sum all for ODBA
    output <- data[,.(dba_x = mean(dynamic_x),
                      dba_y = mean(dynamic_y),
                      dba_z = mean(dynamic_z),
                      odba = mean(dynamic_x)+ mean(dynamic_y)+ mean(dynamic_z)),
                   by = .(burst_no_ID)]
  }
  
  ## calculate ODBA for each burst using my function
  system.time(output <- ODBA_calc_dt(orn_acc_resamp))
  
  
  ## Combine ODBA values with burst info
  
  ## Extract only one line per burst
  dups <- subset(orn_Gr, select = c("burst_no_ID")) %>% duplicated
  sum(dups) ##how many rows are duplicates 
  ##now remove these duplicate rows
  orn_acc_short <- orn_Gr %>% filter(!dups)
  
  ## now join the data together
  orn_acc_short <- cbind(orn_acc_short, output$odba)
  setnames(orn_acc_short, old = c("V2"), new = c("ODBA"))
  
  ## remove some of the extra rows that aren't needed anymore
  ODBA_Gr <- subset(orn_acc_short, select = c("device_id", "UTC_datetime", "ODBA", "tag_year"))
  
  
  #############
  ## ODBA 2. ##
  ############# Sub sampling ODBA data
  
  message("Resampling ODBA")
  
  ## Re-sampling function from Liam
  ## adapted slightly to also loop through years as this makes the function quicker and I have 7+ years
  ## need to specify a time difference and a unit that you want to re-sample too
  ## TD is the tracking data set, dt is the new sampling rate and unit is the unit of dt
  trackSubSampyear <- function(TD = TD, dt=dt, unit=unit){
    
    TD <- TD[order(TD$UTC_datetime),] 
    years = unique(TD$year)
    
    for(j in 1:length(years)){
      
      message(j, " out of ", length(years))
      
      TD_sub = filter(TD, year == years[j]) # extract one year at a time
      
      # breakdown to datasets per bird
      unid = unique(TD_sub$tag_year) 
      nrid = length(unid)
      TDall = list(nrid)  
      TDsubred = list(nrid)
      timestep = paste(dt,unit)
      # create time sequence from min to max time with step sizes that are defined at the start of the function
      dati_start = min(TD_sub$UTC_datetime, na.rm= T)
      dati_end = max(TD_sub$UTC_datetime, na.rm= T)
      datiseq = seq(from=dati_start, to=dati_end, by=timestep)
      
      for (i in 1:nrid) {
        Dtemp = TD_sub[TD_sub$tag_year == unid[i],]
        idx = sapply(datiseq, function(x) which.min( abs( difftime( Dtemp$UTC_datetime, x, units='mins')))) # finds closest time in data to your created time series
        TDall[[i]] = Dtemp
        TDsubred[[i]] = unique(Dtemp[idx,]) # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
      }
      
      TDsubred2 <- do.call("rbind", TDsubred)
      if(j == 1){Resamp_Acc <- TDsubred2}
      else{Resamp_Acc <- rbind(Resamp_Acc, TDsubred2)}
      
    }
    
    return(Resamp_Acc)
    
  }
  
  
  
  ## prepare data for function, timestamp needs to be POSIX object and add year column
  ODBA_Gr$UTC_datetime <- as.POSIXct(ODBA_Gr$UTC_datetime, format= "%Y-%m-%d %H:%M:%OS", tz= "GMT")
  ODBA_Gr$year <- year(ODBA_Gr$UTC_datetime)
  
  ## re-sample the data if the required rate is different to that of the raw data
  if(sampling_matrix$ODBA_sampling_rate[B] == Acc_rate){print("Sampling rate already same as that requested")
    ODBA_Gr_resamp <-ODBA_Gr} else{## apply function to full data set
      ODBA_Gr_resamp <- trackSubSampyear(TD = ODBA_Gr, dt = sampling_matrix$ODBA_sampling_rate[B], unit = 'mins') }
  
  
  
  
  
  #############
  ## ODBA 3. ##
  ############# summarizing average daily ODBA
  
  ## create data set of just unique tag_day in the data set
  ODBA_Gr_resamp$tag_date <- paste0(ODBA_Gr_resamp$device_id, "_", as.Date(ODBA_Gr_resamp$UTC_datetime))
  ODBA_data <- distinct(ODBA_Gr_resamp, tag_date, .keep_all= TRUE)
  ODBA_data$ODBA <- NULL
  
  ##average daily ODBA for each tag then put back in original data set so retain all additional columns
  ODBA_Avg_pday <- ODBA_Gr_resamp %>% 
    group_by(tag_date) %>% 
    summarise(avg_odba = mean(ODBA)) %>% 
    full_join(ODBA_data, by = "tag_date") 
  
  
  
  #### filtering out tags that don't have data for the full breeding season or were dead ####
  ## create data set with tag years with little data removed
  no_of_days <- as.data.frame(table(ODBA_Avg_pday$tag_year)) %>% filter(Freq > 30)
  colnames(no_of_days)[1] <- "tag_year"
  no_of_days$tag_year <- as.character(no_of_days$tag_year)
  
  ##joining the above tag list with ODBA data so only have tag years that have data for most of a breeding season
  ODBA_Avg_pday <- inner_join(ODBA_Avg_pday, no_of_days, by = "tag_year")
  ODBA_Avg_pday$Freq <- NULL
  length(unique(ODBA_Avg_pday$tag_year)) ##should be the length as no_of_days data set
  
  
  
  
  

  #############
  ## Comb 1. ##
  ############# Refine ODBA data sets
  
  colnames(ODBA_Avg_pday)
  ##now bind to the ODBA data
  ##select the columns needed from ODBA data and rename so they match the movement data
  ODBA_Avg_pday$date <- as.Date(ODBA_Avg_pday$UTC_datetime)
  ODBA_Avg_pday <- subset(ODBA_Avg_pday, select = c("tag_date","avg_odba", "device_id", "year", "tag_year", "date"))
  ODBA_Avg_pday %>% setnames(old=c("device_id", "tag_year"), 
                             new=c("id", "Tag_year"))
  ODBA_Avg_pday$date <- as.character(ODBA_Avg_pday$date)
  ODBA_Avg_pday$id <- as.character(ODBA_Avg_pday$id)
  
  ##now full join the two together
  Inc <- ODBA_Avg_pday
  
  
  
  
  
  
  ###################
  ## Classifier 1. ##
  ################### Adding additional tag info
   
  ##importing dataset with breeding histories and adding info Inca data
  setwd("~/Additional data files")
  tag_info <- read.csv("Tagged bird summary data new.csv")
  
  ##creating ID column to bind everthing by
  tag_info$Bird.ID <- as.character(tag_info$Bird.ID)
  tag_info$S.N <- as.character(tag_info$S.N)
  tag_info$S.N <- ifelse(is.na(tag_info$S.N) == TRUE, tag_info$Bird.ID, tag_info$S.N)
  
  ##selecting important columns
  tag_info2 <- subset(tag_info, select = c("Age", "Mass", "Sex", "Ringing.location", 
                                           "S.N", "X2017.Gos","X2018.Gos", "X2019.Gos", "X2020.Gos"))
  colnames(tag_info2)[5] <- "id"
  
  ###joining together odba data and tag info
  Incex <- inner_join(Inc, tag_info2, by = "id")
  nrow(Inc); nrow(Incex)
  length(unique(Inc$Tag_year));length(unique(Incex$Tag_year)) ## two lengths should be the same or have lost tag_years
  
  
  
  
  
  ###################
  ## Classifier 2. ##
  ################### Extract training set
  
  ##filtering out the required data for each plot
  Breeders_2017 <- Incex %>% filter(year == 2017 & X2017.Gos %in% c(1:8) & Sex == "F")
  
  Breeders_2018 <- Incex %>% filter(year == 2018 & X2018.Gos %in% c(1:8) & Sex == "F")
  
  Breeders_2019 <- Incex %>% filter(year == 2019 & X2019.Gos %in% c(1:8) & Sex == "F")
  
  Breeders_2020 <- Incex %>% filter(year == 2020 & X2020.Gos %in% c(1:8) & Sex == "F")
  
  Breeders <- rbind(Breeders_2020, Breeders_2019, Breeders_2018, Breeders_2017)
  
  ## extract the unknown breeders, first just filter out adult females
  Females <- Incex %>% filter(Sex == "F" & Age == "A")
  ## now remove the Tag_years that are in the breeders training set
  training_birds <- unique(Breeders$Tag_year)
  Females <- filter(Females, !Tag_year %in% training_birds)
  Females <-  Females %>%  filter(!Tag_year == "17791_2020") # this tag has lots of missing data so going to take it now
  
  ##plotting ODBA through breeding season
  ##functions for plotting ODBA
  ODBA_plot <- function(h){
    ggplot(h, aes(x = date, y = avg_odba, col = Tag_year)) + geom_point() +
      facet_wrap(~Tag_year, scales = "free_x")
  }
  

  ODBA_plot2 <- function(h){
    ggplot(h, aes(x = date, y = avg_odba, col = Incubating.acc)) + geom_point() +
      facet_wrap(~Tag_year, scales = "free_x") + theme(legend.position = "none")
  }
  
  
  ###Plotting the various groups
  #ODBA_plot(Breeders)
  #ODBA_plot(Irish)
  #ODBA_plot(Scot)
  #ODBA_plot(Males)
  
  
  
  
  
  ###################
  ## Classifier 3. ##
  ################### Label training set
  
  Breeders <- rbind(Breeders_2020, Breeders_2019, Breeders_2018, Breeders_2017)
  ## First remove the 2020 observed successful breeders so they can be used for testing
  Breed_train <- Breeders 
  
  
  ## create 24 day rolling average of daily ODBA that is left centered
  Roll24_function <- function(x){ rollmean(x$avg_odba, k = 24, align = "left", fill = NA) } 
  
  Breed_train <- Breed_train %>% 
                 group_by(Tag_year) %>% 
                 nest() %>% 
                 mutate(Roll24_ODBA = map(data, Roll24_function)) %>% 
                 unnest(cols = c(data, Roll24_ODBA))
  
  ## Now Identify the lowest value 24 day period
  ## label that as incubating, the 2 day either side as buffers and everything else not incubating
  Tag_years <- unique(Breed_train$Tag_year)
  
  Inc_days <- 24 ## the length of minimum breeding period
  Buff <- 3 ## the buffer period
  
  ## loop through each tag_year
  for(j in 1:length(Tag_years)){
    
    message(j)
    
    Tag_sub <- filter(Breed_train, Tag_year == Tag_years[j])
    
    ## find the start of the 24 day period with the lowest ODBA
    low <- which.min(Tag_sub$Roll24_ODBA)
    
    Tag_sub$status <- NA
    
    Tag_sub$status[low:(low+Inc_days-1)] <- "Incubating"
    
    Tag_sub$status[(low-Buff):(low-1)] <- "Buffer"
    
    Tag_sub$status[(low+Inc_days):(low+Inc_days+Buff-1)] <- "Buffer"
    
    if(j ==1){Breed_train2 <- Tag_sub}
    else{Breed_train2 <- rbind(Breed_train2, Tag_sub)}
    
  }
  
  
  ## label all other days as not Incubating
  Breed_train2$status <- ifelse(is.na(Breed_train2$status) == T, "N_incubating", Breed_train2$status)
  
  
  
  
  
  ###################
  ## Classifier 4. ##
  ################### Calculate threshold from training set
  
  ##From the training data calculate quantities for ODBA, nsd and ddist of incubating and non incubating days
  # for incubating days calculate 95th quantile
  incubating_days <- Breed_train2 %>% filter(status == "Incubating")
  qinc_acc <- as.numeric(quantile(incubating_days$avg_odba, probs = 0.975)) 
  
  # for non-incubating days calculate 5th quantile
  non_incubating_days <- Breed_train2 %>% filter(status == "N_incubating")
  qnoinc_acc <- as.numeric(quantile(non_incubating_days$avg_odba, probs = 0.025, na.rm = TRUE))
  
  
  ## Save the 95% CIs from the training data
  error <- qt(0.975,df=length(incubating_days$avg_odba)-1)*sd(incubating_days$avg_odba)/sqrt(length(incubating_days$avg_odba))
  mean(incubating_days$avg_odba)
  
  
  
  
  ###################
  ## Classifier 5. ##
  ################### Use training data to label all other birds
  
  ## Use quantiles from acc to label days and accounts for if the two quantile values overlap
  ## No overlap: 2 = below qinc, 1 = in between qinc and qnoinc & 0 = above qnoinc
  ## Overlap: 2 = below qnoinc, 1 = in between qnoinc and qinc & 0 = above qinc
  if(qinc_acc < qnoinc_acc){
    Females$Incubating.acc <- ifelse(Females$avg_odba < qinc_acc, 2, 0)
    Females$Incubating.acc <- ifelse(Females$avg_odba > qinc_acc & Females$avg_odba < qnoinc_acc, 1, Females$Incubating.acc)
  }else{
    
    Females$Incubating.acc <- ifelse(Females$avg_odba < qnoinc_acc, 2, 0)
    Females$Incubating.acc <- ifelse(Females$avg_odba < qinc_acc & Females$avg_odba > qnoinc_acc, 1, Females$Incubating.acc)
  }
  
  ## plot the points labeled above
  ODBA_plot2(Females)

  
  
  ## loop to extend incubation periods based of days were I have labeled as 1 and 2 in the previous steps
  ## create a list of all the tag years
  tag_years <- (unique(Females$Tag_year))
  ## run the loop
  for (i in 1:length(tag_years)) {
    
    ## subset individual tag years for the loop
    single_tag <- filter(Females, Tag_year == tag_years[i])
    ## have as incubating if there is a 2 in the acc column
    single_tag$Incubating.cert <- ifelse(single_tag$Incubating.acc == 2, 1, 0)
    
    ## if there is a 1 in the acc column then only specifiy if incubating if it is near a 2
    for(j in 4:(nrow(single_tag)-3)) {
      
      if(single_tag$Incubating.cert[j] == 1){next} 
      
      for(z in 1:3) { 
        
        if(single_tag$Incubating.cert[j] == 1){next} 
        ## if we are below acc qinc threshold search up to 3 days ahead for certainly incubating day
        single_tag$Incubating.cert[j] <-  ifelse(single_tag$Incubating.acc[j] == 1 & single_tag$Incubating.acc[j+z] == 2, 
                                                 1, single_tag$Incubating.cert[j])
      }
      for(z in 1:3) { 
        
        if(single_tag$Incubating.cert[j] == 1){next} 
        ## if we are below acc qinc threshold search 1 day behind for certainly incubating day
        single_tag$Incubating.cert[j] <-  ifelse(single_tag$Incubating.acc[j] == 1 & single_tag$Incubating.acc[j-z] == 2,
                                                 1, single_tag$Incubating.cert[j])
      }
    
    }
    
    ## re join all the individual tag years into a new data frame
    if(i == 1){all_tags <- single_tag
    } else {all_tags <- rbind(all_tags, single_tag)} 
    
  }
  
  ## make sure there are no NAs in the Incubating.cert column
  all_tags$Incubating.cert <- ifelse(is.na(all_tags$Incubating.cert) == T, 0, all_tags$Incubating.cert)
  
  ## plot the outputs to check loop hasn't gone awol
  #ODBA_plot2(all_tags)
  
  
  ## Here is a verbal run through of the rules used to classify incubation
  ##  0. For each day in Greenland calculate average ODBA, net squared displacement (NSD) and distance between median locations of successive days (ddist)
  ##  1. Filter out adults females from the original data set (> 3rd breeding season as adult (17812 will be in 3rd BS in 2020))
  ##  2. Label known breeders based off low periods of ODBA. Label middle 24-26 days of low ODBA period as incubating and place buffers 1-2 days either side
  ##  3. Remove buffers from training data set sp just have the definitely incubating or not incubating days
  ##  4. For Incubating days calculate 95th quantile (qinc) for ODBA, nsd and ddist. For non Incubating days calculate 5th quantile (qnoinc) for ODBA, nsd and ddist.
  ##  5. Use quantiles from nsd and ddist to label days as; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
  ##  6. Next label days as certainly incubating if they are below the qinc threshold for ODBA, nsd and ddist
  ##  7. Label all days that have ODBA below max value from training data. 
  ##  8. Use a loop to extend the incubation periods based off days labeled as certainly incubating to those days that fall below the ODBA threshold
  ##  9. If acc runs out while incubating then use nsd and ddist to label the rest of incubation, must have a value of 2 in either the quantile labeling from nsd or ddist to be extended
  
  
  ## Changes slightly for the accelerometer only stuff
  ## Now have a incubating and a non-incubating threshold. Label as 2 if below the 97.5th quantile for incubating days
  ## and 1 if between the 97.5th quantile for incubating days and 2.5th quantile for non-incubating days
  ## everything else receives a 0. Then change any days labeled as 1 to a 2 if either 3 days ahead or 1 day behind there is a day
  ## labeled as a 2. All other values stay the same
  
  
  
  
  ###################
  ## Classifier 6. ##
  ################### Extract Attempt start and end dates
  
  ## Changed: 04/06/2020
  ## Changed method to include movement metrics so this bit of the code changed as well
  ## Changed: 29/09/2021
  ## At some of the higher sampling rates there are more than one possible incubation, need to come up with some rules to pick the most likely one
  
  
  ## label each incubation with it's own code
  ## Filter out the days that have been classified as incubation
  all_tags_inc <- filter(all_tags, Incubating.cert == 1)
  ## calculate number of days in between each row
  all_tags_inc$date <- as.Date(all_tags_inc$date)
  all_tags_inc$day_diff <- as.numeric(all_tags_inc$date -lag(all_tags_inc$date))
  all_tags_inc$change <- ifelse(all_tags_inc$day_diff == 1, 0, 1)
  all_tags_inc$change[1] <- 1
  
  ## use cumsum function to label each incubation
  all_tags_inc$Inc_ID <- cumsum(all_tags_inc$change)
  
  ## extract all the tag years used so i don't miss the birds that didn't have any incubation
  unique_tag_year <- subset(all_tags, select = c("Tag_year")) %>% unique()
  
  ##now use dplyr to create summary of each putative incubation
  All_attempts <- all_tags_inc %>% 
    group_by(Tag_year, Inc_ID) %>% 
    summarise(attempt_start = min(date), 
              attempt_end = max(date), 
              length = as.numeric(((attempt_end-attempt_start)+1)),
              score_tot = sum(Incubating.acc, na.rm = T),
              max_score = max(Incubating.acc, na.rm = T),
              ave_score = score_tot/length) %>% 
    full_join(unique_tag_year, by = "Tag_year")
  
  ##change NAs in the length column to 0's
  All_attempts$length <- ifelse(is.na(All_attempts$length) == TRUE, 0, All_attempts$length)
  
  ## add year day column for the start dates
  All_attempts$yday <- yday(All_attempts$attempt_start)
  
  ## Identify any tag-years with multiple putative incubation
  Dups22 <- duplicated(All_attempts$Tag_year)
  
  ## extract tag_years with multiple incubation attempts
  Dup_tags <- All_attempts[Dups22 == T, ]
  
  #### ADD THE IF() ELSE() STATEMENT IN HERE
  
  ## separate out the tags with multiple incubation and those without multiple incubation
  Dup_incs <- All_attempts %>% filter(Tag_year %in% c(Dup_tags$Tag_year))
  Non_dups <- All_attempts %>% filter(!Tag_year %in% c(Dup_tags$Tag_year))
  
  ## If there are no duplicates then skip the whole next section as it doesn't need doing
  if(nrow(Dup_incs) == 0){Attempt_length <- Non_dups}else{
    ## Rules for picking the best putative incubation
    ## 1) Remove incubation that don't have at least one day that was classified as a 6
    ## 2) Remove incubation during the molting period if all in the molting period then just pick the one with the best average daily score 
    ## 3) Pick the incubation with the best average score
    ## (Single days with a score of 6 might mess this up)
    
    ## work out the earliest start date and maximum incubation length for each tag_year in the duplicated data set
    Dup_sum <- Dup_incs %>% 
      group_by(Tag_year) %>% 
      summarise(earlist_start = min(yday, na.rm = T),
                max_length = max(length, na.rm = T)) %>% 
      full_join(Dup_incs, by = "Tag_year")
    
    ## first step remove any rows without a max score of 6
    ## this removes those 1/2 days of incubation that aren't connected to the main incubation
    ## Need to change this score to 4 or 2 in the GPS only and Acc only incubation
    Dup_sum2 <- filter(Dup_sum, max_score == 2)
    
    ## loop through each tag_year
    ## get unique tag_years in the data set
    unique_dups <- unique(Dup_sum2$Tag_year)
    
    ## run loop to select top incubation
    for(qq in 1:length(unique_dups)) {
      
      ## extract the tag year
      dup_TY <- filter(Dup_sum2, Tag_year == unique_dups[qq])
      
      ## if only one incubation left then just skip to the end
      if(nrow(dup_TY) == 1){dup_top <- dup_TY}
      else{
        ## take different paths if all putative incubation are in July
        ## If all in July then pick the one with the highest ave_score
        ## If ones before July then filter out for July
        if(dup_TY$earlist_start[1] <= 183){dup_TY <- filter(dup_TY, yday <= 183)}
        else{dup_TY <- filter(dup_TY, ave_score == max(dup_TY$ave_score))}
        
        ## remove 1 day incubation if there are longer ones
        if(max(dup_TY$length) == 1){dup_TY <- dup_TY}
        else{dup_TY <- filter(dup_TY, length > 1)}
        
        ## If any duplicates left then just pick the max score
        if(nrow(dup_TY) == 1){dup_top <- dup_TY}
        else{dup_top <- filter(dup_TY, ave_score == max(dup_TY$ave_score))
             dup_top <- filter(dup_top, length == max(dup_top$length))
             dup_top <- filter(dup_top, yday == min(dup_top$yday))}
      }
      
      ## join together all of the top incubation
      if(qq == 1){dup_all <- dup_top}
      else{dup_all <- rbind(dup_all, dup_top)}
      
    }
    
    ## Now join together the duplicate set and the non duplicate set
    dup_all <- dplyr::select(dup_all, c("Tag_year", "Inc_ID", "attempt_start", "attempt_end",   
                                 "length", "score_tot", "max_score", "ave_score", "yday"))
    Attempt_length <- rbind(as.data.frame(Non_dups), as.data.frame(dup_all))
    
    
  }
  
  
  ## add info on what the sampling rate were like for this round
  Attempt_length$ODBA_sampling <- sampling_matrix$ODBA_sampling_rate[B]
  Attempt_length$burst_length <- sampling_matrix$burst_length[B]
  Attempt_length$Inc_ID <- NULL
  
  ## add classifier values to the sampling matrix
  sampling_matrix$Inc_ODBA_q2.5[B]   <- as.numeric(quantile(incubating_days$avg_odba, probs = 0.025, na.rm = TRUE)) 
  sampling_matrix$Inc_ODBA_q97.5[B]  <- qinc_acc
  sampling_matrix$Inc_ODBA_mean[B]   <- mean(incubating_days$avg_odba, na.rm = TRUE)

  sampling_matrix$Not_ODBA_q2.5[B]   <- as.numeric(quantile(non_incubating_days$avg_odba, probs = 0.975, na.rm = TRUE))
  sampling_matrix$Not_ODBA_q97.5[B]  <- qnoinc_acc
  sampling_matrix$Not_ODBA_mean[B]   <- mean(non_incubating_days$avg_odba, na.rm = TRUE)

  
  if(B == 1){All_samples <- Attempt_length}else{All_samples <- rbind(All_samples, Attempt_length)}
  
  if(B == nrow(sampling_matrix)){classifier_values <- sampling_matrix
  message("Code finished")} else{print("Loop done")}
}



#################
## END OF LOOP ##
#################

write_csv(All_samples, file = "~/Incubation methods chapter/All classification outputs/Main MS outputs/Acc classifier outputM.csv")
write_csv(sampling_matrix, file = "~/Incubation methods chapter/All classification outputs/Main MS outputs/Acc classifier Sampling matrixM.csv")



#################
## Plotting 1. ##
################# Plot the training values from the loops



## Plot the ODBA training values
class_vals <- classifier_values[,1:5]
class_vals2 <- classifier_values[,c(1:2, 6:8)]
colnames(class_vals2)[3:ncol(class_vals2)] <- colnames(class_vals)[3:ncol(class_vals2)]
class_vals$Class <- "Incubating"
class_vals2$Class <- "Not Incubating"
ODBA_vals <- rbind(class_vals, class_vals2)

#colnames(class_vals)
ODBA_vals %>% 
  filter(burst_length == 30) %>% 
  ggplot(aes(x=ODBA_sampling_rate, y=Inc_ODBA_mean, group= Class, colour = Class)) + 
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin=Inc_ODBA_q2.5, ymax=Inc_ODBA_q97.5), width=1) +
    theme_light() +
    labs(y="Daily average ODBA", x = "ACC burst rate (mins)", title= "Training values as Acc burst rate increases")











