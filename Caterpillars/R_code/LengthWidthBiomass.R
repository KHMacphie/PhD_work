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
cater$Biomass <- cater$caterpillar.mass
MeasuredBiomass1718 <- merge(measured1718, cater, by="STDY", duplicates.ok=TRUE)
#### Why have 5 rows been removed?!- NAs when merged?

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
MeasuredBiomass1718$caterpillar.mass <- NULL

### removing ones that dont have 1 caterpillar reported or have no mass
onecaterpillar <- subset(MeasuredBiomass1718, caterpillars==1)
onecaterpillar <- subset(onecaterpillar, Biomass!="")
onecaterpillar$Year <- as.factor(onecaterpillar$Year)
#interval cencoring <0.01-0.02
onecaterpillar$Biomass1 <- revalue(onecaterpillar$Biomass, c("<0.01"="0.001", "0.01"="0.001"))
onecaterpillar$Biomass2 <- revalue(onecaterpillar$Biomass, c("<0.01"="0.01", "0.01"="0.01"))
onecaterpillar$Biomass1 <- as.numeric(as.character(onecaterpillar$Biomass1))
onecaterpillar$Biomass2 <- as.numeric(as.character(onecaterpillar$Biomass2))


# volume of a cylinder: V=π(r^2)h
onecaterpillar$radiussqu <- (onecaterpillar$Width/2)^2
onecaterpillar$Volume <- pi*onecaterpillar$radiussqu*onecaterpillar$Length

par(mfrow = c(1, 1))
plot(onecaterpillar$Volume, onecaterpillar$Biomass2)
### pretty good line, a few rogue ones..

hist(onecaterpillar$Biomass1) ## definitely not Gaussian
hist(onecaterpillar$volume) ## definitely not Gaussian

#log transforming volume and biomass- problem with zero?!
onecaterpillar$logVolume <- log(onecaterpillar$Volume)
hist(onecaterpillar$logVolume) ## far more normally distributed

onecaterpillar$logBiomass1 <- log(onecaterpillar$Biomass1)
onecaterpillar$logBiomass2 <- log(onecaterpillar$Biomass2)
hist(onecaterpillar$logBiomass1) # still not normal
hist(onecaterpillar$logBiomass2) # still not normal

# plots before logged
ggplot(onecaterpillar, aes(Volume, Biomass1, colour=y))+
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

#need priors...
prior<-list(R=list(V=diag(1), nu=0.002))
#BiomassVolume <- MCMCglmm(cbind(logBiomass1, logBiomass2) ~ logVolume, family= "cengaussian", data=onecaterpillar, prior=prior, nitt=400000, burnin=40000) 
#save(BiomassVolume, file = "~/Documents/Models/BiomassVolume.RData")
load("~/Documents/Models/BiomassVolume.RData")  

par(mfcol=c(1,1))
hist(BiomassVolume$VCV)
#### get ally to check histogram of residuals

# Put together equation for conversion including variance
# logbiomass = coeff*logvolume + intercept
  # coeff= 1.167569
  # intercept= -7.961385

coeff <- mean(BiomassVolume$Sol[,"logVolume"])
intercept <- mean(BiomassVolume$Sol[,"(Intercept)"])

predlogvol <- seq(-1.4, 7, 0.0001)
predlogmass <- coeff*predlogvol + intercept
predlogmassvar <- var(BiomassVolume$Sol[,"logVolume"])*predlogvol + var(BiomassVolume$Sol[,"(Intercept)"]) #not sure?!

plot(onecaterpillar$logVolume, onecaterpillar$logBiomass1, col=2)
points(predlogvol, predlogmass, type="l")
points(onecaterpillar$logVolume, onecaterpillar$logBiomass2, col=1)

plot(onecaterpillar$Volume, onecaterpillar$Biomass2, col=1)
points(exp(predlogvol), exp(predlogmass), type="l")
points(onecaterpillar$Volume, onecaterpillar$Biomass1, col=2)

#check prediction in onecaterpillar
onecaterpillar$predmass <- exp(coeff*onecaterpillar$logVolume + intercept)
# loop for var doesnt work but also dont know how i want to calc var or really what it means in terms of how we will use it
#for (i in length(onecaterpillar$logVolume)){
#  onecaterpillar$var[i] <- exp(var(BiomassVolume$Sol[,"logVolume"]*i + BiomassVolume$Sol[,"(Intercept)"]))
#}
onecaterpillar$difference <- onecaterpillar$predmass-onecaterpillar$Biomass1
hist(onecaterpillar$difference)
plot(onecaterpillar$Biomass2, onecaterpillar$predmass)
points(onecaterpillar$Biomass2, onecaterpillar$Biomass2, type="l")

# Calcutale predicted biomass for jacks caterpillars
jacks <- read.csv("/Users/s1205615/Dropbox/master_data/inverts/Caterpillars.csv")
jacks$radiussqu <- (jacks$width/2)^2
jacks$Volume <- pi*jacks$radiussqu*jacks$length
jacks$logVolume <- log(jacks$Volume)
jacks$STDY <- paste(jacks$site, jacks$tree, jacks$date, jacks$year)
jacks$logbiomass <- coeff*jacks$logVolume + intercept
jackscater <- data.frame(year=jacks$year, site=jacks$site, tree=jacks$tree, date=jacks$date, STDY=jacks$STDY, jbiomass=exp(jacks$logbiomass)) 
jackscaterSTDY <- aggregate(jbiomass~ STDY, jackscater, sum, na.action = na.pass)
jackscaterSTDY$jbiomass <- round(jackscaterSTDY$jbiomass, digits=2)
jackscaterSTDY$jbiomass <- as.factor(jackscaterSTDY$jbiomass)
jackscaterSTDY$jbiomass <- revalue(jackscaterSTDY$jbiomass, c("0"="0.01"))
jackscaterSTDY$jbiomass<- as.numeric(as.character(jackscaterSTDY$jbiomass))


# Bring all together and interval censor ones with caterpillar biomass making sure NAs stay as NA
pmatch(cater$STDY, jackscaterSTDY$STDY)
allbiomass<- merge(cater, jackscaterSTDY, by="STDY", duplicates.ok=TRUE, all.x=TRUE)
allbiomass$Biomass <- revalue(allbiomass$Biomass, c("<0.01"="0.01", "0.01"="0.01"))
allbiomass$Biomass <- as.numeric(as.character(allbiomass$Biomass))
allbiomass <- data.frame(STDY=allbiomass$STDY, caterpillars=allbiomass$caterpillars, Biomass=allbiomass$Biomass, jbiomass=allbiomass$jbiomass)
allbiomass <- subset(allbiomass, Biomass!="NA" | jbiomass!="NA") ## subset samples with caterpillar biomass so not applying zero mass to unknown ones
allbiomass$biomass2 <- rowSums(allbiomass[,c("Biomass", "jbiomass")], na.rm=TRUE)
allbiomass$biomass1 <- as.factor(allbiomass$biomass2)
allbiomass$biomass1 <- revalue(allbiomass$biomass1, c("0.01"="0.001"))
allbiomass$biomass1 <- as.numeric(as.character(allbiomass$biomass1))
allbiomass <- data.frame(STDY=allbiomass$STDY, biomass1=allbiomass$biomass1, biomass2=allbiomass$biomass2)

complete <- merge(cater, allbiomass, by="STDY", duplicates.ok=TRUE, all.x=TRUE)

#### !!!!!! IF EVER HAVE NA IN CATERPILLARS IT WILL BE ASSIGNED BIOMASS OF ZERO



#seeing if a quadratic fits better? aware that it shouldnt
#prior<-list(R=list(V=diag(1), nu=0.002))
#BiomassVolumeSq <- MCMCglmm(cbind(logBiomass1, logBiomass2) ~ logVolume+I(logVolume^2), family= "cengaussian", data=onecaterpillar, prior=prior, nitt=200000, burnin=20000) 

#predlogvol2 <- seq(-1.4, 7, 0.0001)
#predlogmass2 <- mean(BiomassVolumeSq$Sol[,"logVolume"])*predlogvol2 + mean(BiomassVolumeSq$Sol[,"(Intercept)"]) + mean(BiomassVolumeSq$Sol[,"I(logVolume^2)"])*predlogvol2^2

#plot(predlogvol2, predlogmass2, type="l")
#points(onecaterpillar$logVolume, onecaterpillar$logBiomass2, col=2)
#points(onecaterpillar$logVolume, onecaterpillar$logBiomass1, col=1)

#plot(onecaterpillar$Volume, onecaterpillar$Biomass2, col=2)
#points(exp(predlogvol2), exp(predlogmass2), type="l")
#points(onecaterpillar$Volume, onecaterpillar$Biomass1, col=1)


#### Try and see if larger interval censor improved fit?
