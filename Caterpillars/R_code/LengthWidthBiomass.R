rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(plyr)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe #####


cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")

measured1718 <- read.csv("/Users/s1205615/Documents/PhD (stopped using 12:2:19)/Caterpillar Measurements/CaterMeasuredPlus.1718.csv")

cater$STDY <- paste(cater$site, cater$tree, cater$date, cater$year)
measured1718$STDY <- paste(measured1718$Site, measured1718$Tree, measured1718$Date, measured1718$Year)
pmatch(measured1718$STDY, cater$STDY, duplicates.ok = TRUE)
cater <- rename(cater, Biomass="caterpillar.mass")
MeasuredBiomass1718 <- merge(measured1718, cater, by="STDY", duplicates.ok=TRUE)
#### Why have 5 rows been removed?!

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
MeasuredBiomass1718$mass.uncertain <- NULL
MeasuredBiomass1718$cater.sample <- NULL

### removing ones that dont have 1 caterpillar reported or have no mass
onecaterpillar <- subset(MeasuredBiomass1718, caterpillars==1)
onecaterpillar <- subset(onecaterpillar, Biomass!="")
onecaterpillar$Year <- as.factor(onecaterpillar$Year)
#interval cencoring <0.01-0.02
onecaterpillar$Biomass1 <- revalue(onecaterpillar$Biomass, c("<0.01"="0.00", "0.01"="0.00", "0.02"="0.00"))
onecaterpillar$Biomass2 <- revalue(onecaterpillar$Biomass, c("<0.01"="0.02", "0.01"="0.02", "0.02"="0.02"))
onecaterpillar$Biomass1 <- as.numeric(as.character(onecaterpillar$Biomass1))
onecaterpillar$Biomass2 <- as.numeric(as.character(onecaterpillar$Biomass2))

# volume of a cylinder: V=π(r^2)h
onecaterpillar$radiussqu <- (onecaterpillar$Width/2)^2
onecaterpillar$Volume <- pi*onecaterpillar$radiussqu*onecaterpillar$Length

plot(onecaterpillar$Volume, onecaterpillar$Biomass2)
### pretty good line, a few rogue ones..

hist(onecaterpillar$Biomass1) ## definitely not Gaussian
hist(onecaterpillar$volume) ## definitely not Gaussian

#log transforming volume and biomass- problem with zero?!
onecaterpillar$logVolume <- log(onecaterpillar$Volume)
hist(onecaterpillar$logVolume) ## far more normally distributed

onecaterpillar$logBiomass1 <- log(onecaterpillar$Biomass1) ### zero becomes infinity, do I add 1?
onecaterpillar$logBiomass2 <- log(onecaterpillar$Biomass2)
hist(onecaterpillar$logBiomass1) # still not normal
hist(onecaterpillar$logBiomass2) # still not normal

# plots before logged
ggplot(onecaterpillar, aes(Volume, Biomass1))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()


ggplot(onecaterpillar, aes(Volume, Biomass2, colour=Year))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

#plot logged
ggplot(onecaterpillar, aes(logVolume, logBiomass2, colour=Year))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

## different between years..? cant put in logged as -inf in biomass 1
massvolumeyearlm <- lm(cbind(Biomass1,Biomass2)~Volume*Year, data=onecaterpillar)
summary(massvolumeyearlm)


