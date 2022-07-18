## Luke Ozsanlav-Harris

## Produce all the plots for my methods MS investigating how biologging can be used to ID incubation
## Three scripts proceed this one that run the three different classifiers, need to read in the outputs for each of these scripts 
## Plots need to be high quality for a MS sent to Animal Biotelemetry




## packages required
##loading Packages
library(amt)
library(tidyverse)
library(lubridate)
library(data.table)
library(zoo)
library(geosphere)
library(ggpubr)
library(DescTools)
library(RColorBrewer)






##
#### 1. Read in outputs from previous scripts ####
##

setwd("~/Incubation methods chapter/All classification outputs/Main MS outputs")
All_samples <- fread("GPS & Acc classifier outputM.csv")
GPS_samples <- fread("GPS classifier outputM.csv")
ACC_samples <- fread("Acc classifier outputM.csv")



# All_samples <- GPS_Acc_classifier_output4
# 
# Ref_samples2 <- Ref_samples

##
#### 2. Calculate precision and recall for all scripts ####
##

## extract the reference sampling regime
colnames(All_samples)
Ref_samples <- filter(All_samples, ODBA_sampling == 6, GPS_sampling == 15, burst_length == 30)
nrow(Ref_samples)





##
#### 2.1 JOINT-classifier precision and recall ####
##

## just check that there aren't any duplicates
duptest <- All_samples %>% dplyr::select(c("Tag_year", "GPS_sampling", "ODBA_sampling", "burst_length")) %>% duplicated()
table(duptest)

## extract the reference sampling regime
colnames(All_samples)
Ref_samples <- filter(All_samples, ODBA_sampling == 6, GPS_sampling == 15, burst_length == 30)


## Loop through all samples and determine the number of false positives and false negatives
for(j in 1:nrow(Ref_samples)) {
  print(j)
  
  ## extract one tag_year at a time
  Sub <- filter(All_samples, Tag_year == Ref_samples$Tag_year[j])
  
  ## initiate some columns
  Sub$false_pos <- 0
  Sub$false_neg <- 0
  
  
  ## If there was an incubation identified in the reference sample
  if (is.na(Ref_samples$attempt_start[j]) == F) {
    
    ## Calculate if there is overlap between the reference sequence and the query sequence. This also gives the number of true positives
    Ref_range <- c(yday(Ref_samples$attempt_start[j]), yday(Ref_samples$attempt_end[j]))
    Sample_range <- cbind(c(yday(Sub$attempt_start)), c(yday(Sub$attempt_end)))
    Sample_range[is.na(Sample_range) == T] <- 1000 ## need this in here in case at some sampling rates no incubation was found
    Sub$true_pos <- DescTools::Overlap(Ref_range, Sample_range)
    Sub$true_pos <- ifelse(Sub$true_pos > 0, Sub$true_pos + 1, Sub$true_pos)
    
    Sub$true_pos <- ifelse(Sub$length == 1 & Ref_samples$length[j] ==1 & Sub$attempt_start == Ref_samples$attempt_start[j], 1, Sub$true_pos)
    Sub$true_pos <- ifelse(Sub$length == 1 & Ref_samples$length[j] > 1 & Sub$attempt_start >= Ref_samples$attempt_start[j] & Sub$attempt_end <= Ref_samples$attempt_end[j], 1, Sub$true_pos)
    Sub$true_pos <- ifelse(Sub$length > 1 & Ref_samples$length[j] == 1 & Sub$attempt_start >= Ref_samples$attempt_start[j] & Sub$attempt_end <= Ref_samples$attempt_end[j], 1, Sub$true_pos)
    
    
    ## When the predicted incubation is shifted compared to the reference incubation but they still overlap 
    ## identify false positives and false negatives at the beginning end ends of incubation
    Sub$false_pos <- ifelse(Sub$true_pos > 0 & Ref_samples$attempt_start[j] > Sub$attempt_start, 
                            abs(as.numeric(as.Date(Sub$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub$false_pos)
    
    Sub$false_pos <- ifelse(Sub$true_pos > 0 & Ref_samples$attempt_end[j] < Sub$attempt_end, 
                            Sub$false_pos + abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub$false_pos)
    
    Sub$false_neg <- ifelse(Sub$true_pos > 0 & Ref_samples$attempt_start[j] < Sub$attempt_start, 
                            abs(as.numeric(as.Date(Sub$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub$false_neg)
    
    Sub$false_neg <- ifelse(Sub$true_pos > 0 & Ref_samples$attempt_end[j] > Sub$attempt_end, 
                            Sub$false_neg + abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub$false_neg)
    
    ## Need to account for if the predicted window doesn't overlap with the reference window
    ## In there is no overlap then the length of the reference incubation become the false negative 
    ## and the length fo the predicted incubation become the false positives
    Sub$false_neg <- ifelse((Sub$true_pos == 0), Ref_samples$length[j], Sub$false_neg)
    
    Sub$false_pos <- ifelse((Sub$true_pos == 0), Sub$length, Sub$false_pos)
    
    
    ## if the reference incubation lies completely within the predicted incubation
    Sub$false_pos <- ifelse((Ref_samples$attempt_start[j] > Sub$attempt_start & Ref_samples$attempt_end[j] < Sub$attempt_end),  
                            Sub$length - Ref_samples$length[j]  ,Sub$false_pos)
    
    Sub$false_neg <- ifelse((Ref_samples$attempt_start[j] > Sub$attempt_start & Ref_samples$attempt_end[j] < Sub$attempt_end),  
                            0 ,Sub$false_neg)
    
    
    ## If predicted incubation length is zero but reference incubation int then put the reference incubation length in the false neg column
    Sub$false_neg <- ifelse(is.na(Sub$attempt_end) == T, 
                            Ref_samples$length[j], Sub$false_neg)
    
    
  }else{
    
    ## if the reference sample finds no incubation and the current sampling regime does then add the length of the predicted incubation ot false pos column
    Sub$false_pos <- ifelse(is.na(Sub$attempt_end) == F, 
                            abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Sub$attempt_start)))+1,
                            Sub$false_pos)
    
    Sub$true_pos <- 0
  }
  
  ## bind together each sub section into one data frame
  if(j ==1){All_samples_class <- Sub}else{All_samples_class <- rbind(All_samples_class, Sub)}
  
}



## summarize the number of observations and false positives and false negatives per sampling regime
#Refs <- Ref_samples %>% dplyr::select(Tag_year, length)
#colnames(Refs)[2] <- "Ref_length"
#All_samples_class <- full_join(All_samples_class, Refs, by = "Tag_year")



##### CAN SKIP STRAIGHT TO THIS STAGE TO DO SOME OF THE PLOTTING ###
## write out the data set for the precision plot
#write_csv(All_samples_class, "~/Incubation methods chapter/Classifier_precision_across_sampling_regimes_and_tags.csv")
#setwd("~/Incubation methods chapter")
#All_samples_class <- fread("Classifier_precision_across_sampling_regimes_and_tags.csv")

## calculate the number of days mis-classified
All_samples_class$misclas <- All_samples_class$false_pos + All_samples_class$false_neg

## summarize by sampling regime
False_neg_pos <- All_samples_class %>% 
  group_by(ODBA_sampling, GPS_sampling, burst_length) %>% 
  summarise(no_obs = sum(length, na.rm = T), 
            no_Fpos = sum(false_pos, na.rm = T),
            no_Fneg = sum(false_neg, na.rm = T),
            no_Tpos = sum(true_pos, na.rm = T))

## Calculate precision and recall
False_neg_pos$precis <- False_neg_pos$no_Tpos/(False_neg_pos$no_Tpos + False_neg_pos$no_Fpos)
False_neg_pos$recall <- False_neg_pos$no_Tpos/(False_neg_pos$no_Tpos + False_neg_pos$no_Fneg)








##
#### 2.2 GPS-classifier precision and recall ####
##

## just check that there aren't any duplicates
duptest <- GPS_samples %>% dplyr::select(c("Tag_year", "GPS_sampling")) %>% duplicated()
table(duptest)


## Loop through all samples and determine the number of false positives and false negatives
for(j in 1:nrow(Ref_samples)) {
  print(j)
  
  ## extract one tag_year at a time
  Sub2 <- filter(GPS_samples, Tag_year == Ref_samples$Tag_year[j])
  
  ## initiate some columns
  Sub2$false_pos <- 0
  Sub2$false_neg <- 0
  
  
  ## If there was an incubation identified in the reference sample
  if (is.na(Ref_samples$attempt_start[j]) == F) {
    
    ## Calculate if there is overlap between the reference sequence and the query sequence. This also gives the number of true positives
    Ref_range <- c(yday(Ref_samples$attempt_start[j]), yday(Ref_samples$attempt_end[j]))
    Sample_range <- cbind(c(yday(Sub2$attempt_start)), c(yday(Sub2$attempt_end)))
    Sample_range[is.na(Sample_range) == T] <- 1000 ## need this in here in case at some sampling rates no incubation was found
    Sub2$true_pos <- DescTools::Overlap(Ref_range, Sample_range)
    Sub2$true_pos <- ifelse(Sub2$true_pos > 0, Sub2$true_pos + 1, Sub2$true_pos)
    
    ## When the predicted incubation is shifted compared to the reference incubation but they still overlap 
    ## identify false positives and false negatives at the beginning end ends of incubation
    Sub2$false_pos <- ifelse(Sub2$true_pos > 0 & Ref_samples$attempt_start[j] > Sub2$attempt_start, 
                            abs(as.numeric(as.Date(Sub2$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub2$false_pos)
    
    Sub2$false_pos <- ifelse(Sub2$true_pos > 0 & Ref_samples$attempt_end[j] < Sub2$attempt_end, 
                            Sub2$false_pos + abs(as.numeric(as.Date(Sub2$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub2$false_pos)
    
    Sub2$false_neg <- ifelse(Sub2$true_pos > 0 & Ref_samples$attempt_start[j] < Sub2$attempt_start, 
                            abs(as.numeric(as.Date(Sub2$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub2$false_neg)
    
    Sub2$false_neg <- ifelse(Sub2$true_pos > 0 & Ref_samples$attempt_end[j] > Sub2$attempt_end, 
                            Sub2$false_neg + abs(as.numeric(as.Date(Sub2$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub2$false_neg)
    
    ## Need to account for if the predicted window doesn't overlap with the reference window
    ## In there is no overlap then the length of the reference incubation become the false negative 
    ## and the length fo the predicted incubation become the false positives
    Sub2$false_neg <- ifelse((Sub2$true_pos == 0), Ref_samples$length[j], Sub2$false_neg)
    
    Sub2$false_pos <- ifelse((Sub2$true_pos == 0), Sub2$length, Sub2$false_pos)
    
    
    ## if the reference incubation lies completely within the predicted incubation
    Sub2$false_pos <- ifelse((Ref_samples$attempt_start[j] > Sub2$attempt_start & Ref_samples$attempt_end[j] < Sub2$attempt_end),  
                            Sub2$length - Ref_samples$length[j]  ,Sub2$false_pos)
    
    Sub2$false_neg <- ifelse((Ref_samples$attempt_start[j] > Sub2$attempt_start & Ref_samples$attempt_end[j] < Sub2$attempt_end),  
                            0 ,Sub2$false_neg)
    
    
    ## If predicted incubation length is zero but reference incubation int then put the reference incubation length in the false neg column
    Sub2$false_neg <- ifelse(is.na(Sub2$attempt_end) == T, 
                            Ref_samples$length[j], Sub2$false_neg)
    
    
  }else{
    
    ## if the reference sample finds no incubation and the current sampling regime does then add the length of the predicted incubation ot false pos column
    Sub2$false_pos <- ifelse(is.na(Sub2$attempt_end) == F, 
                            abs(as.numeric(as.Date(Sub2$attempt_end) - as.Date(Sub2$attempt_start)))+1,
                            Sub2$false_pos)
    
    Sub2$true_pos <- 0
  }
  
  ## bind together each sub section into one data frame
  if(j ==1){GPS_samples_class <- Sub2}else{GPS_samples_class <- rbind(GPS_samples_class, Sub2)}
  
}





## summarize accuracy by sampling regime
False_neg_posGPS <- GPS_samples_class %>% 
                  group_by(GPS_sampling) %>% 
                  summarise(no_obs = sum(length, na.rm = T), 
                            no_Fpos = sum(false_pos, na.rm = T),
                            no_Fneg = sum(false_neg, na.rm = T),
                            no_Tpos = sum(true_pos, na.rm = T))

## Calculate precision and recall
False_neg_posGPS$precis <- False_neg_posGPS$no_Tpos/(False_neg_posGPS$no_Tpos + False_neg_posGPS$no_Fpos)
False_neg_posGPS$recall <- False_neg_posGPS$no_Tpos/(False_neg_posGPS$no_Tpos + False_neg_posGPS$no_Fneg)







##
#### 2.3 ACC-classifier precision and recall ####
##

## just check that there aren't any duplicates
duptest <- ACC_samples %>% dplyr::select(c("Tag_year", "ODBA_sampling", "burst_length")) %>% duplicated()
table(duptest)


## Loop through all samples and determine the number of false positives and false negatives
for(j in 1:nrow(Ref_samples)) {
  print(j)
  
  ## extract one tag_year at a time
  Sub3 <- filter(ACC_samples, Tag_year == Ref_samples$Tag_year[j])
  
  ## initiate some columns
  Sub3$false_pos <- 0
  Sub3$false_neg <- 0
  
  
  ## If there was an incubation identified in the reference sample
  if (is.na(Ref_samples$attempt_start[j]) == F) {
    
    ## Calculate if there is overlap between the reference sequence and the query sequence. This also gives the number of true positives
    Ref_range <- c(yday(Ref_samples$attempt_start[j]), yday(Ref_samples$attempt_end[j]))
    Sample_range <- cbind(c(yday(Sub3$attempt_start)), c(yday(Sub3$attempt_end)))
    Sample_range[is.na(Sample_range) == T] <- 1000 ## need this in here in case at some sampling rates no incubation was found
    Sub3$true_pos <- DescTools::Overlap(Ref_range, Sample_range)
    Sub3$true_pos <- ifelse(Sub3$true_pos > 0, Sub3$true_pos + 1, Sub3$true_pos)
    
    ## When the predicted incubation is shifted compared to the reference incubation but they still overlap 
    ## identify false positives and false negatives at the beginning end ends of incubation
    Sub3$false_pos <- ifelse(Sub3$true_pos > 0 & Ref_samples$attempt_start[j] > Sub3$attempt_start, 
                            abs(as.numeric(as.Date(Sub3$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub3$false_pos)
    
    Sub3$false_pos <- ifelse(Sub3$true_pos > 0 & Ref_samples$attempt_end[j] < Sub3$attempt_end, 
                            Sub3$false_pos + abs(as.numeric(as.Date(Sub3$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub3$false_pos)
    
    Sub3$false_neg <- ifelse(Sub3$true_pos > 0 & Ref_samples$attempt_start[j] < Sub3$attempt_start, 
                            abs(as.numeric(as.Date(Sub3$attempt_start) - as.Date(Ref_samples$attempt_start[j]))),
                            Sub3$false_neg)
    
    Sub3$false_neg <- ifelse(Sub3$true_pos > 0 & Ref_samples$attempt_end[j] > Sub3$attempt_end, 
                            Sub3$false_neg + abs(as.numeric(as.Date(Sub3$attempt_end) - as.Date(Ref_samples$attempt_end[j]))),
                            Sub3$false_neg)
    
    ## Need to account for if the predicted window doesn't overlap with the reference window
    ## In there is no overlap then the length of the reference incubation become the false negative 
    ## and the length fo the predicted incubation become the false positives
    Sub3$false_neg <- ifelse((Sub3$true_pos == 0), Ref_samples$length[j], Sub3$false_neg)
    
    Sub3$false_pos <- ifelse((Sub3$true_pos == 0), Sub3$length, Sub3$false_pos)
    
    
    ## if the reference incubation lies completely within the predicted incubation
    Sub3$false_pos <- ifelse((Ref_samples$attempt_start[j] > Sub3$attempt_start & Ref_samples$attempt_end[j] < Sub3$attempt_end),  
                            Sub3$length - Ref_samples$length[j]  ,Sub3$false_pos)
    
    Sub3$false_neg <- ifelse((Ref_samples$attempt_start[j] > Sub3$attempt_start & Ref_samples$attempt_end[j] < Sub3$attempt_end),  
                            0 ,Sub3$false_neg)
    
    
    ## If predicted incubation length is zero but reference incubation int then put the reference incubation length in the false neg column
    Sub3$false_neg <- ifelse(is.na(Sub3$attempt_end) == T, 
                            Ref_samples$length[j], Sub3$false_neg)
    
    
  }else{
    
    ## if the reference sample finds no incubation and the current sampling regime does then add the length of the predicted incubation ot false pos column
    Sub3$false_pos <- ifelse(is.na(Sub3$attempt_end) == F, 
                            abs(as.numeric(as.Date(Sub3$attempt_end) - as.Date(Sub3$attempt_start)))+1,
                            Sub3$false_pos)
    
    Sub3$true_pos <- 0
  }
  
  ## bind together each sub section into one data frame
  if(j ==1){ACC_samples_class <- Sub3}else{ACC_samples_class <- rbind(ACC_samples_class, Sub3)}
  
}



## summarize the number of observations and false positives and false negatives per sampling regime
## calculate the number of true positives for each row
Refs <- Ref_samples %>% select(Tag_year, length)
colnames(Refs)[2] <- "Ref_length"
ACC_samples_class <- full_join(ACC_samples_class, Refs, by = "Tag_year")



## summarize by sampling regime
False_neg_posACC <- ACC_samples_class %>% 
  group_by(ODBA_sampling, burst_length) %>% 
  summarise(no_obs = sum(length, na.rm = T), 
            no_Fpos = sum(false_pos, na.rm = T),
            no_Fneg = sum(false_neg, na.rm = T),
            no_Tpos = sum(true_pos, na.rm = T))

## Calculate precision and recall
False_neg_posACC$precis <- False_neg_posACC$no_Tpos/(False_neg_posACC$no_Tpos + False_neg_posACC$no_Fpos)
False_neg_posACC$recall <- False_neg_posACC$no_Tpos/(False_neg_posACC$no_Tpos + False_neg_posACC$no_Fneg)







##
#### 3. Joint classifier only plots ####
##


##
#### 3.1 JOINT- Precision and accuracy of classifier at different sampling rates ####
##


## plot how the precision and recall vary with sampling regime
False_neg_pos$GPS_sampling2 <- as.factor(False_neg_pos$GPS_sampling) # need this as a factor to get the color right
False_neg_pos$burst_l <- ifelse(False_neg_pos$burst_length == 10, "1s burst",
                                ifelse(False_neg_pos$burst_length == 20, "2s burst", "3s burst"))


## Plot of classifier precision (GPS and Acc data)
prec <- ggplot(False_neg_pos, aes(x = ODBA_sampling, y= precis, group = GPS_sampling2, colour =GPS_sampling2)) + 
  geom_vline(xintercept = 24, colour = "black", linetype = 2) +
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.925, 1) +
  ylab("Precision") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS interval (mins)") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))


## Plot of classifier recall (GPS and Acc data)
rec <- ggplot(False_neg_pos, aes(x = ODBA_sampling, y= recall, group = GPS_sampling2, colour =GPS_sampling2)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.925, 1) +
  ylab("Recall") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS interval (mins)") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))


## save the plot output
ggarrange(prec, rec, labels = c("A", "B"), nrow = 2)
setwd("~/Incubation methods chapter/High quality MS plots/Main MS")
ggsave("JOINT-precision & recall plot.png", 
       width = 22, height = 24, units = "cm")






##
#### 3.2 JOINT- Number of days different to reference sampling rate for different incubation lengths ####
##

## For each reference sample extract the length of the identified incubation
Ref_lengths <- select(Ref_samples, c("Tag_year", "length"))
setnames(Ref_lengths, old = "length", new = "Orig_length")

## add on the original lengths to the main data set
All_samples_class2 <-  inner_join(All_samples_class, Ref_lengths, by = "Tag_year")

## create categories of different incubation lengths
All_samples_class2$categ <- ifelse(All_samples_class2$Orig_length == 0, "Deferal", 
                                   ifelse(All_samples_class2$Orig_length > 0 & All_samples_class2$Orig_length <=5, "1-5 days",
                                          ifelse(All_samples_class2$Orig_length > 5 & All_samples_class2$Orig_length <=10, "6-10 days",
                                                 ifelse(All_samples_class2$Orig_length > 10 & All_samples_class2$Orig_length <=24, "11-24 days",
                                                        ifelse(All_samples_class2$Orig_length > 24, "Hatched", "error")))))


## Group the average number of days incorrectly classified by different incubation length in the reference sample
All_samples_class2$dummy <- 1
Sample_errors_length <- All_samples_class2 %>% 
  filter(burst_length == 30) %>% 
  group_by(categ, ODBA_sampling, GPS_sampling) %>% 
  summarise(Tot_misclass = sum(misclas, na.rm = T),
            Tot_cases = sum(dummy),
            Ave_misclass = Tot_misclass/Tot_cases) 


## create plot of days mis-classified per group
## specify the levels of the factor for plotting
Sample_errors_length$categ <- as.factor(Sample_errors_length$categ)
levels(Sample_errors_length$categ)
Sample_errors_length$categ <- ordered(Sample_errors_length$categ, levels = c("Deferal", "1-5 days", "6-10 days", "11-24 days", "Hatched"))
Sample_errors_length$GPS_sampling <- as.factor(Sample_errors_length$GPS_sampling)
ggplot(Sample_errors_length, aes(x = ODBA_sampling, y= Ave_misclass, colour =GPS_sampling)) + 
  geom_point() +
  facet_wrap(~categ) +
  geom_line() +
  ylim(0, 1) +
  ylab("Days misclassifed per individual") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS interval (mins)") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))


## save the plot
setwd("~/Incubation methods chapter/High quality MS plots/Main MS")
ggsave("JOINT-Days missclassified per individ.png", 
       width = 26, height = 24, units = "cm")






##
#### 4. JOINT vs ACC- classifiers precison and recall ####
##

False_neg_posACC$GPS_sampling <- "Acc-only"
False_neg_posACC$burst_l <- ifelse(False_neg_posACC$burst_length == 10, "1s burst",
                                ifelse(False_neg_posACC$burst_length == 20, "2s burst", "3s burst"))

## bind together the two data sets
Presic_both <- plyr::rbind.fill(False_neg_pos, False_neg_posACC)

## plot the two data sets together
Presic_both$GPS_sampling2 <- as.factor(Presic_both$GPS_sampling) # need this as a factor to get the color right



## Plot of classifier precision (GPS + Acc comparison with Acc-only)
a <-ggplot(Presic_both, aes(x = ODBA_sampling, y = precis, group = GPS_sampling2, colour = GPS_sampling2)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.825, 1) +
  ylab("Precision") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS fix rate (mins)  ") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines")) +
  theme(legend.position = "none")


## Plot of classifier recall (GPS + Acc comparison with Acc-only)
b <- ggplot(Presic_both, aes(x = ODBA_sampling, y = recall, group = GPS_sampling2, colour = GPS_sampling2)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.825, 1) +
  ylab("Recall") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS fix rate (mins)  ") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))

## arrange the two plots, need the widths so the plots are all the same width
p1 <- ggarrange(a,b,labels = c("A", "B"), ncol = 2, widths = c(1, 1.32))






##
#### 5. JOINT vs GPS- classifiers precison and recall ####
##

False_neg_posGPS$ODBA_sampling <- "GPS-only"
False_neg_posGPS$burst_l <- "3s burst"

## bind together the two data sets
Presic_both2 <- plyr::rbind.fill(False_neg_pos, False_neg_posGPS)

## plot the two data sets together
## need to set ODBA sampling to a factor and then make sure it is ordered correctly
Presic_both2$ODBA_sampling <- as.factor(Presic_both2$ODBA_sampling) # need this as a factor to get the color right, 
Presic_both2$ODBA_sampling <- ordered(Presic_both2$ODBA_sampling, levels = c("6", "12", "24", "48", "72", "96", "120", "144", "GPS-only"))
Presic_both2$GPS_sampling <- as.numeric(Presic_both2$GPS_sampling) # needs to be numeric
Presic_both2 <- filter(Presic_both2, GPS_sampling <= 90)

## create colour palette so GPS only line is obvious
levs <- length(levels(Presic_both2$ODBA_sampling))-1
pallette <- brewer.pal(levs, "Paired")

## Plot of classifier precision (GPS + Acc comparison with GPS-only)
c <- Presic_both2 %>% 
  ggplot(aes(x = GPS_sampling, y= precis, group = ODBA_sampling, colour =ODBA_sampling)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.825, 1) +
  ylab("Precision") + xlab("GPS fix rate (mins)") + 
  labs(colour = "Burst interval (mins)") + theme_bw() +
  scale_color_manual(values=c(pallette, "#000000")) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=10), legend.title=element_text(size=12, face ="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines")) +
  theme(legend.position = "none")


## Plot of classifier precision (GPS + Acc comparison with GPS-only)
d <- Presic_both2 %>% 
  ggplot(aes(x = GPS_sampling, y= recall, group = ODBA_sampling, colour =ODBA_sampling)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0.825, 1) +
  ylab("Recall") + xlab("GPS fix rate (mins)") + 
  labs(colour = "Burst interval (mins)") + theme_bw() +
  scale_color_manual(values=c(pallette, "#000000")) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=10), legend.title=element_text(size=12, face ="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))


## comnine all the above plots
p2 <- ggarrange(c,d,labels = c("C", "D"), ncol = 2, widths = c(1, 1.32))
ggarrange(p1, p2, nrow=2)

## save the two compariosn plots
setwd("~/Incubation methods chapter/High quality MS plots/Main MS")
ggsave("JOINT vs ACC&GPS- precision & recall plot.png", 
       width = 34, height = 24, units = "cm")






##
#### 6. Calculate number of additional incubation for each classifier ####
##

##
#### 6.1 JOINT- additional incubation calculation ####
##

## Find where any incubation missed or any mistakenly classified
## Calculate the difference in days from reference start and end date
for(j in 1:nrow(Ref_samples)) {
  
  Sub <- filter(All_samples, Tag_year == Ref_samples$Tag_year[j])
  
  Sub$miss_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == F & is.na(Sub$attempt_start) == T, 1 , 0)
  
  Sub$extra_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == T & is.na(Sub$attempt_start) == F, 1 , 0)
  
  if (is.na(Ref_samples$attempt_start[j]) == T) {
    
    Sub$Start_off <- NA
    Sub$End_off <- NA
    
  }else{
    
    Sub$Start_off <- abs(as.numeric(as.Date(Sub$attempt_start) - as.Date(Ref_samples$attempt_start[j])))
    Sub$End_off <- abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Ref_samples$attempt_end[j])))
    
  }
  
  
  
  if(j ==1){All_samples2 <- Sub}else{All_samples2 <- rbind(All_samples2, Sub)}
  
}


## Summarize the accuracy measure by sampling regime
####**** CHECK IF 53 is the right number here ****####
colnames(All_samples2)
Sample_errors <- All_samples2 %>% 
                  group_by(ODBA_sampling, GPS_sampling, burst_length) %>% 
                  summarise(Miss_inc = sum(miss_inc, na.rm = T), 
                            Extra_inc = sum(extra_inc, na.rm = T),
                            Start_miss = sum(Start_off, na.rm = T), 
                            End_miss = sum(End_off, na.rm = T),
                            Start_ave = Start_miss/53,
                            End_ave = End_miss/53,
                            All_ave = (Start_miss + End_miss)/53)

Ref_samples

##
#### 6.2 GPS- additional incubation calculation ####
##

## Find where any incubation missed or any mistakenly classified
## Calculate the difference in days from reference start and end date

for(j in 1:nrow(Ref_samples)) {
  
  print(j)
  
  Sub <- filter(GPS_samples, Tag_year == Ref_samples$Tag_year[j])
  
  Sub$miss_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == F & is.na(Sub$attempt_start) == T, 1 , 0)
  
  Sub$extra_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == T & is.na(Sub$attempt_start) == F, 1 , 0)
  
  if (is.na(Ref_samples$attempt_start[j]) == T) {
    
    Sub$Start_off <- NA
    Sub$End_off <- NA
    
  }else{
    
    Sub$Start_off <- abs(as.numeric(as.Date(Sub$attempt_start) - as.Date(Ref_samples$attempt_start[j])))
    Sub$End_off <- abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Ref_samples$attempt_end[j])))
    
  }
  
  if(j ==1){GPS_samples2 <- Sub}else{GPS_samples2 <- rbind(GPS_samples2, Sub)}
  
}


## Summarize the accuracy measure by sampling regime
## NOTE there is 68 repeats for each sampling regime
colnames(GPS_samples2)
Sample_errorsGPS <- GPS_samples2 %>% 
                  group_by(GPS_sampling) %>% 
                  summarise(Miss_inc = sum(miss_inc, na.rm = T), 
                            Extra_inc = sum(extra_inc, na.rm = T),
                            Start_miss = sum(Start_off, na.rm = T), 
                            End_miss = sum(End_off, na.rm = T),
                            Start_ave = Start_miss/53,
                            End_ave = End_miss/53,
                            All_ave = (Start_miss + End_miss)/53)




##
#### 6.3 ACC- additional incubation calculation ####
##

## Find where any incubation missed or any mistakenly classified
## Calculate the difference in days from reference start and end date

for(j in 1:nrow(Ref_samples)) {
  
  Sub <- filter(ACC_samples, Tag_year == Ref_samples$Tag_year[j])
  
  Sub$miss_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == F & is.na(Sub$attempt_start) == T, 1 , 0)
  
  Sub$extra_inc <- ifelse(is.na(Ref_samples$attempt_start[j]) == T & is.na(Sub$attempt_start) == F, 1 , 0)
  
  if (is.na(Ref_samples$attempt_start[j]) == T) {
    
    Sub$Start_off <- NA
    Sub$End_off <- NA
    
  }else{
    
    Sub$Start_off <- abs(as.numeric(as.Date(Sub$attempt_start) - as.Date(Ref_samples$attempt_start[j])))
    Sub$End_off <- abs(as.numeric(as.Date(Sub$attempt_end) - as.Date(Ref_samples$attempt_end[j])))
    
  }
  
  
  
  if(j ==1){ACC_samples2 <- Sub}else{ACC_samples2 <- rbind(ACC_samples2, Sub)}
  
}


## Summarize the accuracy measure by sampling regime
## NOTE there is 68 repeats for each sampling regime
colnames(ACC_samples2)
Sample_errorsACC <- ACC_samples2 %>% 
                  group_by(ODBA_sampling, burst_length) %>% 
                  summarise(Miss_inc = sum(miss_inc, na.rm = T), 
                            Extra_inc = sum(extra_inc, na.rm = T),
                            Start_miss = sum(Start_off, na.rm = T), 
                            End_miss = sum(End_off, na.rm = T),
                            Start_ave = Start_miss/53,
                            End_ave = End_miss/53,
                            All_ave = (Start_miss + End_miss)/53)


##
#### 7. Plot number of additional incubation ####
##

GPS_error <- Sample_errorsGPS
Acc_error <- Sample_errorsACC


## add extra column to main data set from plotting
Sample_errors$burst_l <- ifelse(Sample_errors$burst_length == 10, "1s burst",
                                ifelse(Sample_errors$burst_length == 20, "2s burst", "3s burst"))


## Add a ODBA_sampling and burst_l column so that it will bind with the other data set and the plot will work
GPS_error$ODBA_sampling <- "GPS-only"
GPS_error$burst_l <- "3s burst"


## bind together the two data sets
Error_both <- plyr::rbind.fill(Sample_errors, GPS_error)


## plot the two data sets together
## need to set ODBA sampling to a factor and then make sure it is ordered correctly
Error_both$ODBA_sampling <- as.factor(Error_both$ODBA_sampling) # need this as a factor to get the color right, 
Error_both$ODBA_sampling <- ordered(Error_both$ODBA_sampling, levels = c("6", "12", "24", "48", "72", "96", "120", "144", "GPS-only"))
Error_both$GPS_sampling <- as.numeric(Error_both$GPS_sampling) # needs to be numeric
Error_both <- filter(Error_both, GPS_sampling <= 90)

## create colour palette so GPS only line is obvious
levs <- length(levels(Error_both$ODBA_sampling))-1
pallette <- brewer.pal(levs, "Paired")

## Plot of classifier precision (GPS + Acc comparison with GPS-only)
E1 <- Error_both %>% 
  ggplot(aes(x = GPS_sampling, y= Extra_inc, group = ODBA_sampling, colour =ODBA_sampling)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0, 15) +
  ylab("Additional Incubations") + xlab("GPS fix rate (mins)") + 
  labs(colour = "Burst interval (mins)") + theme_bw() +
  scale_color_manual(values=c(pallette, "#000000")) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=10), legend.title=element_text(size=12, face ="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))







## Add a GPS_sampling column so that it will bind with the other data set
Acc_error$GPS_sampling <- "Acc-only"

## bind together the two data sets
Error_both2 <- plyr::rbind.fill(Sample_errors, Acc_error)

## plot the two data sets together
Error_both2$GPS_sampling2 <- as.factor(Error_both2$GPS_sampling) # need this as a factor to get the color right

Error_both2$burst_l <- ifelse(Error_both2$burst_length == 10, "1s burst",
                              ifelse(Error_both2$burst_length == 20, "2s burst", "3s burst"))


## Plot of classifier precision (GPS + Acc comparison with Acc-only)
E2 <-ggplot(Error_both2, aes(x = ODBA_sampling, y = Extra_inc, group = GPS_sampling2, colour = GPS_sampling2)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~burst_l) + 
  ylim(0, 15) +
  ylab("Additional Incubations") + xlab("Burst interval (mins)") + 
  labs(colour = "GPS interval (mins)  ") + theme_bw() +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        axis.text=element_text(size=13, face = "bold"), axis.title=element_text(size=15, face = "bold"), 
        plot.title = element_text(size=14, face="bold"), legend.text=element_text(size=12), legend.title=element_text(size=12, face="bold"),
        strip.text.x = element_text(size = 13, face = "bold"), panel.spacing = unit(2, "lines"))


ggarrange(E1, E2, nrow=2)
setwd("~/Incubation methods chapter/High quality MS plots/Main MS")
ggsave("JOINT vs ACC&GPS- additional incubations plot.png", 
       width = 26, height = 24, units = "cm")

