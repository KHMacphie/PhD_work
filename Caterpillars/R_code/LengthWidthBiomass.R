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


cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")

measured1718 <- read.csv("Documents/PhD/Caterpillar Measurements/CaterMeasurements1718.csv")

cater$STDY <- paste(cater$site, cater$tree, cater$date, cater$year)
measured1718$STDY <- paste(measured1718$Site, measured1718$Tree, measured1718$Date, measured1718$Year)
pmatch(measured1718$STDY, cater$STDY, duplicates.ok = TRUE)
cater <- rename(cater, Biomass="caterpillar.mass")
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

## some sites have 2 entries for the same STDY
table(table(cater$STDY)) ## in whole dataframe
#  2x same STDY= 384, 3x same STDY= 8 
table(table(MeasuredBiomass1718$STDY)) ## in measured samples dataframe
#  2x same STDY= 7
which(table(MeasuredBiomass1718$STDY)==2) # which STDY's have >1 row
# AVI 13 141 2018 / AVI 14 141 2018 / GLF 6 147 2018 / PIT 10 147 2018 / RTH 4 168 2018 / STY 1 147 2018 / STY 8 147 2018
table(MeasuredBiomass1718$caterpillars)
# caterpillars: 0=10, 2=1, 6=1

### removing ones that dont have 1 caterpillar reported
onecaterpillar <- subset(MeasuredBiomass1718, caterpillars==1)
  # still 3 duplicate STDYs but leaving in for now   AVI 13 141 2018 / STY 1 147 2018 / STY 8 147 2018

hist(onecaterpillar$Width)
hist(onecaterpillar$Length)
#hist(onecaterpillar$Biomass)  doesnt work because of <0.01s

onecaterpillar.minustiny <- subset(MeasuredBiomass1718, Biomass!="<0.01")
onecaterpillar.minustiny <- subset(onecaterpillar.minustiny, Biomass!="0")
onecaterpillar.minustiny$Biomass <- as.numeric(as.character(onecaterpillar.minustiny$Biomass))
onecaterpillar.minustiny <- subset(onecaterpillar.minustiny, Biomass!="NA")
hist(onecaterpillar.minustiny$Biomass) ## definitely not Gaussian

# volume of a cylinder: V=π(r^2)h
onecaterpillar.minustiny$radiussqu <- (onecaterpillar.minustiny$Width/2)^2
onecaterpillar.minustiny$volume <- pi*onecaterpillar.minustiny$radiussqu*onecaterpillar.minustiny$Length

plot(onecaterpillar.minustiny$volume, onecaterpillar.minustiny$Biomass)
### pretty good line, a few rogue ones..

massvolumelm <- lm(Biomass~volume, data=onecaterpillar.minustiny)
summary(massvolumelm)
autoplot(massvolumelm, smooth.colour = NA) #not good..?

ggplot(onecaterpillar.minustiny, aes(volume, Biomass))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

ggplot(onecaterpillar.minustiny, aes(volume, Biomass, colour=Year))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()
## different betwen years..?
massvolumeyearlm <- lm(Biomass~volume*Year, data=onecaterpillar.minustiny)
summary(massvolumeyearlm)
# not significant

