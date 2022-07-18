##Luke Ozsanlav-Harris
## PhD incubation methods chapter

## Code initially adapted from scripts that identified incubation events from Acc and GPS data
## Turned those scripts into a loop that allowed sampling rates of Acc and GPS data to change and recalculate incubation lengths
## Finally removed the capacity for Acc data so only uses GPS data. This is for a comparison that comprises part of a methods chapter



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
## Adapted the script so that it just used the GPS data (allow comparison with Acc + GPS approach)



### UPDATES as of 10/06/21:
## Need to change the precision/recall code to the new code in the Acc-only script
## Need to check that script is working properly as have some problems in the Acc script





##loading Packages
library(amt)
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(geosphere)
library(arm)
library(DescTools)


#################
## GPS prep 1. ##
################# Read in and clean GPS data

## set working directory
setwd("~/Ornitela data/ALL_ORNI_DATA/Cleaned CSVs")

## read in data
orn_data <-list.files(pattern = "*.csv") %>% map_df(~fread(.))

## Data Cleaning
## extracting just the GPS data
orn_gps <- subset(orn_data, orn_data$datatype == "GPS")
rm(orn_data)

## parse timestamp
orn_gps$UTC_datetime <- as.character(ymd_hms(orn_gps$UTC_datetime))

## Delete observations where missing lat or long or a timestamp. 
colnames(orn_gps)
ind <- complete.cases(orn_gps[, c("Latitude", "Longitude", "UTC_datetime")])

## The number of re locations with missing coordinates or timestamp (if any).
table(ind)

##removing rows with missing data (if any)
orn_gps <- orn_gps %>% filter(ind)

##Check for duplicated observations (ones with same lat, long, timestamp, and device ID)
ind2 <- orn_gps %>% 
  select(UTC_datetime, Longitude, Latitude, device_id, acc_x, acc_z, acc_y) %>%
  duplicated
sum(ind2) 
orn_gps <- orn_gps %>% filter(!ind2)




#################
## GPS prep 2. ##
################# Trim tracking data to Greenland period

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

##now join the two data sets
orn_gps$UTC_datetime<- ymd_hms(orn_gps$UTC_datetime)
orn_gps$tag_year <- paste(orn_gps$device_id, year(orn_gps$UTC_datetime), sep = "_") ## tag_year column
orn_gps2 <- inner_join(orn_gps, Phen, by = "tag_year")
length(unique(orn_gps$tag_year)); length(unique(orn_gps2$tag_year)) ## check haven't lost lots of tag_years
rm(orn_gps)


## filter out the period in Greenland
orn_Gr <- orn_gps2 %>%  filter(UTC_datetime >= Green_arrive & UTC_datetime <= Green_dep) 
print("number of tag years with GPS"); length(unique(orn_Gr$tag_year))
rm(orn_gps2)



#############
## LOOP 1. ##
############# Set up matrix to loop through different sampling rates

## set the sampling frequencies fo the raw data
GPS_rate <- 15

## Create matrix or required sampling rates using expand.grid
sampling_matrix <- expand.grid(GPS_fix_rate = c(GPS_rate, GPS_rate*2, GPS_rate*3, GPS_rate*4, GPS_rate*5, GPS_rate*6, GPS_rate*7, GPS_rate*8, GPS_rate*9, GPS_rate*10, GPS_rate*11, GPS_rate*12, GPS_rate*14, GPS_rate*16))

## Add empty column to add in training values
sampling_matrix$Inc_nsd_q2.5 <- NA; sampling_matrix$Inc_nsd_q97.5 <- NA; sampling_matrix$Inc_nsd_mean <- NA 
sampling_matrix$Inc_ddist_q2.5 <- NA; sampling_matrix$Inc_ddist_q97.5 <- NA; sampling_matrix$Inc_ddist_mean <- NA

sampling_matrix$Not_nsd_q2.5 <- NA; sampling_matrix$Not_nsd_q97.5 <- NA; sampling_matrix$Not_nsd_mean <- NA 
sampling_matrix$Not_ddist_q2.5 <- NA; sampling_matrix$Not_ddist_q97.5 <- NA; sampling_matrix$Not_ddist_mean <- NA




#####################
### Start of Loop ###
#####################

## loop through each round of sampling
for(B in 1:nrow(sampling_matrix)){
  
  message("Sampling round ", B, " out of ", nrow(sampling_matrix))
  
  
  
  ############
  ## GPS 1. ##
  ############ Re-sample the GPS data
  
  message("Resampling GPS")
  
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
      if(j == 1){Resamp_GPS <- TDsubred2}
      else{Resamp_GPS <- rbind(Resamp_GPS, TDsubred2)}
      
    }
    
    return(Resamp_GPS)
    
  }
  
  
  ## prepare data for function, timestamp needs to be POSIX object and add year column
  orn_Gr$UTC_datetime <- as.POSIXct(orn_Gr$UTC_datetime, format= "%Y-%m-%d %H:%M:%OS", tz= "GMT")
  orn_Gr$year <- year(orn_Gr$UTC_datetime)
  
  ## re-sample the data if the required rate is different to that of the raw data
  if(sampling_matrix$GPS_fix_rate[B] == GPS_rate){print("Sampling rate already same as that requested")
    orn_Gr_resamp <- orn_Gr}else{## apply function to full data set
      orn_Gr_resamp <- trackSubSampyear(TD = orn_Gr, dt = sampling_matrix$GPS_fix_rate[B], unit = 'mins') }
  
  
  
  
  ############
  ## GPS 2. ##
  ############ Use amt to calculate step length and net squared displacement
  
  ## setting timestamp as POSIX object
  orn_Gr_resamp$UTC_datetime <- as.POSIXct(orn_Gr_resamp$UTC_datetime, format= "%Y-%m-%d %H:%M:%S", tz= "GMT")
  
  ##need to convert the data to a track to work with AMT
  trk <- make_track(orn_Gr_resamp, .x = Longitude, .y = Latitude, 
                    .t = UTC_datetime, id = device_id, crs = CRS("+init=epsg:4326"))
  
  ##adding columns for month, year, day and hour
  trk <- trk %>% 
    mutate(day = mday(t_),
           month = month(t_), 
           year = year(t_),
           hour = hour(t_), 
           date = as.Date(t_),
           day_ID = paste(date, id, sep= "_"),
           julian_date = julian(date))
  
  ## create nsd calculator function
  nsd_calc <- function(x) { distGeo(cbind(x$x_[1], x$y_[1]), 
                                    cbind(x$x_, x$y_))/1000 } 
  
  ##calculating daily nsd, need a different grouping to above
  trk <-as.data.frame(trk)
  trk_nsd <- trk %>% 
    group_by(id, year, month, day) %>% 
    nest() %>% 
    mutate(nsd_ =map(data, nsd_calc)) %>% 
    unnest(cols = c(data, nsd_))
  
  
  
  ##Added: 02/06/2020
  ## Created my own metric of movement
  ## For each day it calculates a median lat and a long
  ## then within each tag day it calculates the distance between a fix and
  ## the median location. These distances are then averaged over each tag_day
  
  ## turn track into data frame
  track <- as.data.frame(trk)
  
  ##first work out median lat/long for each tag day
  track <- track %>% 
    group_by(day_ID) %>% 
    summarise(med_long = median(x_), med_lat = median(y_)) %>% 
    full_join(track, by = "day_ID") ##join back to original data set, each tag day has own med lat/long
  
  ## calculate distance between fix and med lat/long for that tag day
  track$distance <- distGeo(cbind(track$x_, track$y_), 
                            cbind(track$med_long, track$med_lat))/1000
  
  ## calculate the average distance to the daily median lat/long
  ave_dist_med_loc <- track %>%
    group_by(id, date) %>%
    summarise(ave_dist = mean(distance), n_fixes = length(x_)) %>% 
    mutate(year= year(ymd(date)), Tag_year = paste(id, year, sep = "_"))
  
  
  
  ## Added: 03/06/2020
  ## added another one of my own metrics
  ## calculate the distance between daily median locations
  Med_locs <-  track %>% 
    group_by(day_ID) %>% 
    summarise(med_long = median(x_), med_lat = median(y_),
              id = min(id), date = min(date), year = min(year))
  
  ##distance calculator function
  Diff_day_dist_calc <- function(x){
    
    distGeo(cbind(x$med_long, x$med_lat), 
            cbind(lag(x$med_long), lag(x$med_lat)))/1000
    
  }
  
  ##nest by tag and then calculate distances
  Med_locs <- Med_locs %>% group_by(id) %>% nest() %>% mutate(
    ddist = map(data, Diff_day_dist_calc)) %>% unnest(cols=c(data, ddist))
  
  ##formatting data set so it will join with other movement metrics
  Med_locs$med_long <- NULL
  Med_locs$med_lat <- NULL
  Med_locs$day_ID <- NULL
  Med_locs$Tag_year <- paste(Med_locs$id, Med_locs$year, sep = "_")
  
  
  
  
  
  ############
  ## GPS 3. ##
  ############ Create daily averages for movement metrics
  
  ##PLOTTING AVERAGE AND MEDIAN DAILY NSD FOR EACH TAG
  ##need to get rid of really high values as this is when the bird was migrating
  
  ##aggregate functions for mean and median nsd
  ave_daily_nsd <- aggregate(nsd_ ~ day_ID, FUN=mean, data=trk_nsd) 
  ave_daily_nsd <- aggregate(id ~ day_ID, FUN=min, data=trk_nsd) %>% inner_join(ave_daily_nsd)
  
  ##slightly different method so the plotting looks better
  median_daily_nsd <- aggregate(nsd_ ~ day_ID, FUN=median, data=trk_nsd) 
  median_daily_nsd <- separate(median_daily_nsd, day_ID,into = c("date", "id"),  sep= "_")
  median_daily_nsd$date <- as.Date.POSIXct(median_daily_nsd$date,  format= "%Y-%m-%d", tz= "GMT")
  median_daily_nsd <- median_daily_nsd %>% mutate(year=year(date))
  median_daily_nsd$Tag_year <- paste(median_daily_nsd$id, median_daily_nsd$year, sep = "_")
  
  ##avergae daily nsd plot
  #ggplot(ave_daily_nsd, aes(x = day_ID, y = log(nsd_))) + geom_point() + facet_wrap(~id, scales = "free_x")
  
  ###logging the nsd seems to make it easier to visalise, get a similar patter to the ODBA just less clear
  #ggplot(median_daily_nsd, aes(x = date, y = log(nsd_))) + geom_point() + scale_x_date(
  #  breaks = date_breaks("months"),labels = date_format("%b")) + facet_wrap(~Tag_year, nrow = 4, scales = "free_x")
  
  
  ## Now join together some of the info above
  ave_dist_med_loc$id <- as.character(ave_dist_med_loc$id)
  Med_locs$id <- as.character(Med_locs$id)
  Move_charactristics <- full_join(ave_dist_med_loc, median_daily_nsd) %>% full_join(Med_locs)
  length(unique(Move_charactristics$Tag_year))
  
  
  
  
  #############
  ## Comb 1. ##
  ############# Combine ODBA and GPS derived data sets
  
  ##create tag_date column in move data
  Move_charactristics$tag_date <- paste(Move_charactristics$id, Move_charactristics$date, sep = "_")
  Move_charactristics$Tag_year <- paste(Move_charactristics$id, Move_charactristics$year, sep = "_")
  Move_charactristics$date <- as.character(Move_charactristics$date)
  
  ##now full join the two together
  Inc <- Move_charactristics
  
  
  
  
  
  
  
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
  nsd_plot <- function(h){
    ggplot(h, aes(x = date, y = log(nsd_), col = Tag_year)) + geom_point() +
      facet_wrap(~Tag_year, scales = "free_x")
  }
  
  nsd_plot2 <- function(h){
    ggplot(h, aes(x = date, y = log(nsd_), col = Incubating.cert)) + geom_point() +
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
  
  ## now only using the GPS data need to come up with a new way of identifying the training data
  ## could run it with both and see what the difference is 
  ## alternatively could use a combination of the two? Could scale the two and then add the together and pick the lowest value
  
  ## extract the observed breeders
  Breeders <- rbind(Breeders_2020, Breeders_2019, Breeders_2018, Breeders_2017)
  
  ## Scale the nsd and ddist columns and then add together
  Breeders$nsd_ <- ifelse(Breeders$nsd_ == 0, 0.000001, Breeders$nsd_) ## bug-fix: in case there is a nsd value of 0. This arises due to a day only having 1 fix
  Breeders$scale_nsd <- rescale(log(Breeders$nsd_))
  Breeders$scale_ddist <- rescale(Breeders$ddist)
  
  ## Add the scaled variables together
  Breeders$scaled_move <- Breeders$scale_nsd + Breeders$scale_ddist
  
  ## First remove the 2020 observed successful breeders so they can be used for testing
  Breed_train <- Breeders  
  
  ## create 24 day rolling average of daily ODBA that is left centered
  Roll24_function <- function(x){ rollmean(x$scaled_move, k = 24, align = "left", fill = NA) } 
  
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
  qinc_nsd <- as.numeric(quantile(incubating_days$nsd_, probs = 0.975))
  qinc_ddist <- as.numeric(quantile(incubating_days$ddist, probs = 0.975))
  
  # for non-incubating days calculate 5th quantile
  non_incubating_days <- Breed_train2 %>% filter(status == "N_incubating")
  qnoinc_nsd <- as.numeric(quantile(non_incubating_days$nsd_, probs = 0.025, na.rm = TRUE))
  qnoinc_ddist <- as.numeric(quantile(non_incubating_days$ddist, probs = 0.025, na.rm = TRUE))
  
  
  
  
  ###################
  ## Classifier 5. ##
  ################### Use training data to label all other birds
  
  ## Again this needs changing to remove the acc data
  ## Have to use a combination of the two movement metrics to label then extend
  
  ##Use quantiles from nsd and ddist to label days; 2 = below qinc, 1 = in between qinc and qnoinc and 0 = above qnoinc
  ##label using ddist
  if(qinc_ddist < qnoinc_ddist){
    Females$Incubating1 <- ifelse(Females$ddist < qinc_ddist, 2, 0)
    Females$Incubating1 <- ifelse(Females$ddist > qinc_ddist & Females$ddist < qnoinc_ddist, 1, Females$Incubating1)
  }else{
    
    Females$Incubating1 <- ifelse(Females$ddist < qnoinc_ddist, 2, 0)
    Females$Incubating1 <- ifelse(Females$ddist < qinc_ddist & Females$ddist > qnoinc_ddist, 1, Females$Incubating1)
  }

  ## End up with NAs in this column which causes problems later
  Females$Incubating1 <- ifelse(is.na(Females$Incubating1) == T, 0, Females$Incubating1)
  
  
  ##label using nsd
  if(qinc_nsd < qnoinc_nsd){
    Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd, 2, 0)
    Females$Incubating2 <- ifelse(Females$nsd_ > qinc_nsd & Females$nsd_ < qnoinc_nsd, 1, Females$Incubating2)
  }else{
    
    Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd, 2, 0)
    Females$Incubating2 <- ifelse(Females$nsd_ < qinc_nsd & Females$nsd_ > qnoinc_nsd, 1, Females$Incubating2)
  }
  
  
  ##combine nsd and ddist labeling by adding
  Females$Incubating.comb <- Females$Incubating1 + Females$Incubating2
  
  ##label as definitely incubating if is labeled with 2 for nsd and ddist and below ODBA threshold
  Females$Incubating.cert <- ifelse(Females$Incubating.comb >= 4, 1, 0)
  
  
  #ODBA_plot2(Females)
  #nsd_plot2(Females)
  
  
  ## loop to extend incubation periods based of days were I have labeled as certainly incubating based off previous steps
  ##create a list of all the tag years
  tag_years <- (unique(Females$Tag_year))
  
  Females$Incubating.cert2 <- Females$Incubating.cert
  
  
  for (i in 1:length(tag_years)) {
    
    single_tag <- filter(Females, Tag_year == tag_years[i]) ## subset individual tag years for the loop
    
    for(j in 4:(nrow(single_tag)-3)) {
      
      if(single_tag$Incubating.cert2[j] == 1){next} 
      
      for(z in 1:3) { 
        
        if(single_tag$Incubating.cert2[j] == 1){next} 
        ## if we are below acc qinc threshold search up to 3 days ahead for certainly incubating day
        single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >=2 & single_tag$Incubating.cert[j+z] == 1, 
                                                 1, single_tag$Incubating.cert2[j])
      }
      for(z in 1:3) { 
        
        if(single_tag$Incubating.cert2[j] == 1){next} 
        ## if we are below acc qinc threshold search 1 day behind for certainly incubating day
        single_tag$Incubating.cert2[j] <-  ifelse(single_tag$Incubating.comb[j] >=2 & single_tag$Incubating.cert[j-z] == 1,
                                                  1, single_tag$Incubating.cert2[j])
      }

      
    }
    ## re join all the individual tag years into a new data frame
    if(i == 1){all_tags <- single_tag
    } else {all_tags <- rbind(all_tags, single_tag)} 
    
  }
  
  ## plot the outputs to check loop hasn't gone awol
  #nsd_plot2(all_tags)
  
  ## make sure there are no NAs in the Incubating.cert column
  all_tags$Incubating.cert2 <- ifelse(is.na(all_tags$Incubating.cert2) == T, 0, all_tags$Incubating.cert2)
  
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
  
  
  ## Changes slightly for the GPS only stuff
  ## Added in a section if the 2.5th quantile from the non-incubating days and 97.5th quantile of the incubating days overlap
  ## Used a 4 point tier system using the quantities of nsd and ddist, only days with 4 points were certainly incubating
  
  
  
  
  ###################
  ## Classifier 6. ##
  ################### Extract Attempt start and end dates
  
  ## Changed: 04/06/2020
  ## Changed method to include movement metrics so this bit of the code changed as well
  ## Changed: 29/09/2021
  ## At some of the higher sampling rates there are more than one possible incubation, need to come up with some rules to pick the most likely one
  
  
  ## label each incubation with it's own code
  ## Filter out the days that have been classified as incubation
  all_tags_inc <- filter(all_tags, Incubating.cert2 == 1)
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
              score_tot = sum(Incubating.comb, na.rm = T),
              max_score = max(Incubating.comb, na.rm = T),
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
    Dup_sum2 <- filter(Dup_sum, max_score == 4)
    
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
  
  ## Make sure there are no more duplicates
  stopifnot(duplicated(Attempt_length$Tag_year) == F)
    
  ## add info on what the sampling rate were like for this round
  Attempt_length$GPS_sampling <- sampling_matrix$GPS_fix_rate[B]
  Attempt_length$Inc_ID <- NULL
  
  ## add classifier values to the sampling matrix
  sampling_matrix$Inc_nsd_q2.5[B]    <- as.numeric(quantile(incubating_days$nsd_, probs = 0.025, na.rm = TRUE))
  sampling_matrix$Inc_nsd_q97.5[B]   <- qinc_nsd
  sampling_matrix$Inc_nsd_mean[B]    <- mean(incubating_days$nsd_, na.rm = TRUE)
  sampling_matrix$Inc_ddist_q2.5[B]  <- as.numeric(quantile(incubating_days$ddist, probs = 0.025, na.rm = TRUE))
  sampling_matrix$Inc_ddist_q97.5[B] <- qinc_ddist
  sampling_matrix$Inc_ddist_mean[B]  <- mean(incubating_days$ddist, na.rm = TRUE)

  sampling_matrix$Not_nsd_q2.5[B]    <- qnoinc_nsd
  sampling_matrix$Not_nsd_q97.5[B]   <- as.numeric(quantile(non_incubating_days$nsd_, probs = 0.975, na.rm = TRUE))
  sampling_matrix$Not_nsd_mean[B]    <- mean(non_incubating_days$nsd_, na.rm = TRUE)
  sampling_matrix$Not_ddist_q2.5[B]  <- qnoinc_ddist
  sampling_matrix$Not_ddist_q97.5[B] <- as.numeric(quantile(non_incubating_days$ddist, probs = 0.975, na.rm = TRUE))
  sampling_matrix$Not_ddist_mean[B]  <- mean(non_incubating_days$ddist, na.rm = TRUE) 
  

  
  if(B == 1){All_samples <- Attempt_length}else{All_samples <- rbind(All_samples, Attempt_length)}
  
  if(B == nrow(sampling_matrix)){classifier_values <- sampling_matrix
  message("Code finished")} else{print("Loop done")}
}



#################
## END OF LOOP ##
#################

write_csv(All_samples, file = "~/Incubation methods chapter/All classification outputs/Main MS outputs/GPS classifier outputM.csv")
write_csv(sampling_matrix, file = "~/Incubation methods chapter/All classification outputs/Main MS outputs/GPS classifier Sampling matrixM.csv")






#################
## Plotting 1. ##
################# Plot the training values from the loops


## Plot the movement metrics training data
index <- duplicated(classifier_values$GPS_fix_rate)
GPS_vals <- filter(classifier_values, !index)

move_vals <- GPS_vals[,1:7]
move_vals2 <- GPS_vals[,c(1, 8:13)]
colnames(move_vals2)[2:ncol(move_vals2)] <- colnames(move_vals)[2:ncol(move_vals)]
move_vals$Class <- "Incubating"
move_vals2$Class <- "Not Incubating"
Move_vals <- rbind(move_vals, move_vals2)
colnames(Move_vals)

ggplot(Move_vals, aes(x= GPS_fix_rate, y= Inc_nsd_mean, group= Class, colour = Class)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=Inc_nsd_q2.5, ymax=Inc_nsd_q97.5), width=6) +
  theme_light() + 
  coord_cartesian(ylim=c(0,0.5)) +
  labs(y="Median daily nsd (km)", x = "GPS fix rate (mins)", title= "Training values as GPS fix rate increases")


ggplot(Move_vals, aes(x= GPS_fix_rate, y= Inc_ddist_mean, group= Class, colour = Class)) + 
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=Inc_ddist_q2.5, ymax=Inc_ddist_q97.5), width=6) +
  theme_light() +
  labs(y="Distance from previous day (km)", x = "GPS fix rate (mins)", title= "Training values as GPS fix rate increases") +
  coord_cartesian(ylim=c(0,0.1)) 












