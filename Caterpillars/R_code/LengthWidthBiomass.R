rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe #####


cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating.csv", 
                  col_types = cols(year = col_factor(levels = c("2014", 
                                                                "2015", "2016", "2017", "2018"))))

measured1718 <- read_csv("Documents/PhD/Caterpillar Measurements/CaterMeasurements1718.csv", 
                         col_types = cols(Year = col_factor(levels = c("2017", 
                                                                       "2018"))))

cater$STDY <- paste(cater$site, cater$tree, cater$date, cater$year)
measured1718$STDY <- paste(measured1718$Site, measured1718$Tree, measured1718$Date, measured1718$Year)
pmatch(measured1718$STDY, cater$STDY, duplicates.ok = TRUE)
cater <- rename(cater, Biomass="caterpillar mass")
measured1718$Biomass <- cater$Biomass
MeasuredBiomass1718 <- merge(measured1718, cater, by="STDY", duplicates.ok=TRUE)
#### Why have 9 rows been removed?!

MeasuredBiomass1718$year <- NULL
MeasuredBiomass1718$notes <- NULL
MeasuredBiomass1718$aphids <- NULL
MeasuredBiomass1718$other <- NULL
MeasuredBiomass1718$bibio <- NULL
MeasuredBiomass1718$beetles <- NULL
MeasuredBiomass1718$site <- NULL
MeasuredBiomass1718$tree <- NULL
MeasuredBiomass1718$`tree species` <- NULL
MeasuredBiomass1718$spiders <- NULL
MeasuredBiomass1718$recorder <- NULL
MeasuredBiomass1718$weather <- NULL
MeasuredBiomass1718$date <- NULL

###### So many strange things