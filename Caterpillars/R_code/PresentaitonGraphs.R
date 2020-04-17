###################################
#### Extra Presentation Graphs ####
###################################

#### Abundance/peak height variation
pred_day <- seq(100,200,1) # different to all other graphs!
H.TempPoisson.2018 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2018.2 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.3) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2018.3 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.6) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.4 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.9) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]-0.3) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.6 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]-0.6) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

par(mfcol=c(1,1),mar=c(2,1,1,1), oma=c(1,3,0,0),cex=1.8)
plot(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=8, xlab="Date", ylab="Abundance", ylim=c(0,0.2)) #black
points(pred_day, exp(H.TempPoisson.2018.2), type="l",lwd=3, col=4) 
points(pred_day, exp(H.TempPoisson.2018.3), type="l",lwd=3,col=5) 
points(pred_day, exp(H.TempPoisson.2018.4), type="l", lwd=3, col=3)
points(pred_day, exp(H.TempPoisson.2018.5), type="l", lwd=3, col=6)
points(pred_day, exp(H.TempPoisson.2018.6), type="l", lwd=3, col=2)
title(ylab="Abundance", outer=TRUE, line = 1)
title( xlab="Date", outer=TRUE, line = 0)

# to make graph with different peak date
plot(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=2, xlab="Date", ylab="Abundance", ylim=c(0,0.12), xlim=c(100,200)) 


# elev/lat/temp graph for Apr temp
rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)

site <- read.csv("Dropbox/master_data/site/site_details.csv")

site<- rename(site, latitude="Mean.Lat")
site<- rename(site, longitude="Mean.Long")
site<- rename(site, elevation="Mean.Elev")
glimpse(site)

temp <- read.csv("Dropbox/master_data/site/temperatures.csv")
temp$year <- as.factor(temp$year)
temp$yearsite<- paste(temp$site, temp$year) #nests site and year?

##making data frame for the required data- site, year, latitude, elevation- need to add mean temp
mean_temps <- data.frame(site=temp$site, year=temp$year)
View(mean_temps)
mean_temps <- mean_temps[!duplicated(mean_temps), ] #remove duplicated rows
mean_temps$yearsite <- paste(mean_temps$site, mean_temps$year)
pmatch(mean_temps$yearsite, temp$yearsite)
mean_temps <- mean_temps %>% arrange(site) #arrange by site to match means order

# for april
Means <- data.frame(site=mean_temps$site, year=mean_temps$year, yearsite=mean_temps$yearsite, MeanTemp=tapply(apply(temp[, 1059:1778],1, mean), temp$yearsite, mean))
View(Means)
pmatch(Means$site, site$site)
temp_lat_elev <- merge(Means, site, by="site", duplicates.ok=TRUE)
View(temp_lat_elev)
temp_lat_elev <- select(temp_lat_elev, site, year, MeanTemp, latitude, elevation)

temp_lat_elev_wide <- spread(temp_lat_elev, year, MeanTemp)
View(temp_lat_elev_wide)

temp_lat_elev_wide <- arrange(temp_lat_elev_wide, latitude)

par(oma=c(1,1,1,4))
plot(temp_lat_elev_wide$latitude, temp_lat_elev_wide$elevation, type="l", col=1, ylim=c(0,850), xlab=expression("Latitude ("~degree~"N)"), ylab="Elevation (m)", cex.axis=1.5, cex.lab=1.5)
par(new = T)
plot(temp_lat_elev_wide$latitude, temp_lat_elev_wide$`2018`, col=2, type="l", axes=F, xlab=NA, ylab=NA, ylim=c(0,10), cex.axis=1.5)
axis(side = 4, cex.axis=1.5)
mtext(side = 4, line = 3, expression("Mean Temperature (°C)"), cex=1.5)
legend("topright", legend=c("Elevation","Temperature"), lty=c(1,1), pch=c(NA, NA), col=c("black", "red3"), cex=1.2)

#### Bar chart of site lat, elev and 2019 apr temp ####
temp_lat_elev_wide$site <- factor(temp_lat_elev_wide$site, levels = temp_lat_elev_wide$site[order(temp_lat_elev_wide$latitude)])
temp_lat_elev_wide$temp <- temp_lat_elev_wide$'2019'

ggplot(temp_lat_elev_wide, aes(site, elevation))+
  geom_bar(aes(fill=temp),stat = "identity")+
  theme_bw()+
  ylab("Elevation (m)")+
  xlab("Site by Latitude")+
  theme(text = element_text(size=20),axis.text.x= element_text(angle=90))+
  scale_fill_gradient(name = "Temperature (°C)", 
                        low = "gold", high = "red3") # saves as 7"x13"


transectAD <- subset(temp_lat_elev_wide, latitude < 56.9)
par(mfrow=c(2,2))
hist(temp_lat_elev_wide$elevation, breaks=100, xlim=c(10, 433), xlab=("Elevation (m)"), main="")
text(370, 2.8, "Full Transect", cex=1, col=1)
hist(transectAD$elevation, breaks=100, xlim=c(10, 433), ylim=c(0,3), xlab=("Elevation (m)"), main="")
text(370, 2.8, "Half Transect", cex=1, col=1)
hist(temp_lat_elev_wide$temp, breaks=200, xlim=c(6.6, 9.5), xlab=("Temperature (°C)"), main="")
text(9.1, 2.8, "Full Transect", cex=1, col=1)
hist(transectAD$temp, breaks=100, xlim=c(6.6, 9.5), ylim=c(0,3), xlab=("Temperature (°C)"), main="")
text(9.1, 2.8, "Half Transect", cex=1, col=1) # saved as 7"x8"
