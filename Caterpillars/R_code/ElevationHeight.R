rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(plyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)
library(forcats)
library(gridExtra)

################################
#### Elevation peak height  ####
################################

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)
cater[,6:9] <- NULL
cater[,7:15] <- NULL
cater$siteyear <- paste(cater$site, cater$year)
cater$siteday <- paste(cater$site, cater$date, cater$year) 
cater.date <- aggregate(cater$caterpillars~cater$siteday, FUN= mean)
names(cater.date)[c(1:2)] <-  c("siteday", "mean.cater")

store <- pmatch(cater.date$siteday, cater$siteday)
cater.date <- cbind(cater.date, cater$year[store], cater$site[store], cater$date[store])
names(cater.date)[c(3:5)] <-  c("year", "site", "date")

cater.date <- cater.date[order(cater.date$date),]
cater.date$siteday <- NULL

cater.spread <- spread(cater.date, date, mean.cater, fill=NA)

### need to get the mean from each 4 column combo
caterdays <- data.frame(year =cater.spread$year, site=cater.spread$site)
caterdays$d1178 <- rowSums(cater.spread[,3:4], na.rm=TRUE)
caterdays$d11920<- rowSums(cater.spread[,5:6], na.rm=TRUE)
caterdays$d1212 <- rowSums(cater.spread[,7:8], na.rm=TRUE)
caterdays$d1234 <- rowSums(cater.spread[,9:10], na.rm=TRUE)
caterdays$d1256 <- rowSums(cater.spread[,11:12], na.rm=TRUE)
caterdays$d1278 <- rowSums(cater.spread[,13:14], na.rm=TRUE)
caterdays$d12930 <- rowSums(cater.spread[,15:16], na.rm=TRUE)
caterdays$d1312 <- rowSums(cater.spread[,17:18], na.rm=TRUE)
caterdays$d1334 <- rowSums(cater.spread[,19:20], na.rm=TRUE)
caterdays$d1356 <- rowSums(cater.spread[,21:22], na.rm=TRUE)
caterdays$d1378 <- rowSums(cater.spread[,23:24], na.rm=TRUE)
caterdays$d13940 <- rowSums(cater.spread[,25:26], na.rm=TRUE)
caterdays$d1412 <- rowSums(cater.spread[,27:28], na.rm=TRUE)
caterdays$d1434 <- rowSums(cater.spread[,29:30], na.rm=TRUE)
caterdays$d1456 <- rowSums(cater.spread[,31:32], na.rm=TRUE)
caterdays$d1478 <- rowSums(cater.spread[,33:34], na.rm=TRUE)
caterdays$d14950 <- rowSums(cater.spread[,35:36], na.rm=TRUE)
caterdays$d1512 <- rowSums(cater.spread[,37:38], na.rm=TRUE)
caterdays$d1534 <- rowSums(cater.spread[,39:40], na.rm=TRUE)
caterdays$d1556 <- rowSums(cater.spread[,41:42], na.rm=TRUE)
caterdays$d1578 <- rowSums(cater.spread[,43:44], na.rm=TRUE)
caterdays$d15960 <- rowSums(cater.spread[,45:46], na.rm=TRUE)
caterdays$d1612 <- rowSums(cater.spread[,47:48], na.rm=TRUE)
caterdays$d1634 <- rowSums(cater.spread[,49:50], na.rm=TRUE)
caterdays$d1656 <- rowSums(cater.spread[,51:52], na.rm=TRUE)
caterdays$d1678 <- rowSums(cater.spread[,53:54], na.rm=TRUE)
caterdays$d16970 <- rowSums(cater.spread[,55:56], na.rm=TRUE)
caterdays$d1712 <- rowSums(cater.spread[,57:58], na.rm=TRUE)
caterdays$d1734 <- rowSums(cater.spread[,59:60], na.rm=TRUE)
caterdays$d175 <- cater.spread[,61]

cater4days <- data.frame(year =caterdays$year, site=caterdays$site)
cater4days$d1178920 <- rowMeans(caterdays[,3:4], na.rm=TRUE)
cater4days$d1192012 <- rowMeans(caterdays[,4:5], na.rm=TRUE)
cater4days$d121234 <- rowMeans(caterdays[,5:6], na.rm=TRUE)
cater4days$d123456 <- rowMeans(caterdays[,6:7], na.rm=TRUE)
cater4days$d125678 <- rowMeans(caterdays[,7:8], na.rm=TRUE)
cater4days$d1278930 <- rowMeans(caterdays[,8:9], na.rm=TRUE)
cater4days$d1293012 <- rowMeans(caterdays[,9:10], na.rm=TRUE)
cater4days$d131234 <- rowMeans(caterdays[,10:11], na.rm=TRUE)
cater4days$d133456 <- rowMeans(caterdays[,11:12], na.rm=TRUE)
cater4days$d135678 <- rowMeans(caterdays[,12:13], na.rm=TRUE)
cater4days$d1378940 <- rowMeans(caterdays[,13:14], na.rm=TRUE)
cater4days$d1394012 <- rowMeans(caterdays[,14:15], na.rm=TRUE)
cater4days$d141234 <- rowMeans(caterdays[,15:16], na.rm=TRUE)
cater4days$d143456 <- rowMeans(caterdays[,16:17], na.rm=TRUE)
cater4days$d145678 <- rowMeans(caterdays[,17:18], na.rm=TRUE)
cater4days$d1478950 <- rowMeans(caterdays[,18:19], na.rm=TRUE)
cater4days$d1495012 <- rowMeans(caterdays[,19:20], na.rm=TRUE)
cater4days$d151234 <- rowMeans(caterdays[,20:21], na.rm=TRUE)
cater4days$d153456 <- rowMeans(caterdays[,21:22], na.rm=TRUE)
cater4days$d155678 <- rowMeans(caterdays[,22:23], na.rm=TRUE)
cater4days$d1578960 <- rowMeans(caterdays[,23:24], na.rm=TRUE)
cater4days$d1596012 <- rowMeans(caterdays[,24:25], na.rm=TRUE)
cater4days$d161234 <- rowMeans(caterdays[,25:26], na.rm=TRUE)
cater4days$d163456 <- rowMeans(caterdays[,26:27], na.rm=TRUE)
cater4days$d165678 <- rowMeans(caterdays[,27:28], na.rm=TRUE)
cater4days$d1678970 <- rowMeans(caterdays[,28:29], na.rm=TRUE)
cater4days$d1697012 <- rowMeans(caterdays[,29:30], na.rm=TRUE)
cater4days$d171234 <- rowMeans(caterdays[,30:31], na.rm=TRUE)
cater4days$d17345 <- rowMeans(caterdays[,31:32], na.rm=TRUE)


cater4days$maxperiod <- colnames(cater4days[,3:31])[apply(cater4days[,3:31],1,which.max)]
cater4days$max <- apply(cater4days[,3:31],1,max)

site <- read.csv("Dropbox/master_data/site/site_details.csv")
store <- pmatch(cater4days$site, site$site, duplicates.ok=TRUE)
elevmax <- cbind(cater4days, site$Mean.Elev[store], site$Mean.Lat[store])
names(elevmax)[c(34:35)] <-  c("elevation", "latitude")
elevmax <- elevmax[-225,]  # removed pit19
ggplot(elevmax, aes(elevation, max, col=year))+
  geom_point()+
  theme_bw()

hist(elevmax$max, breaks=200)

elevmodel <- MCMCglmm(log(max)~elevation + year, data=elevmax, family="gaussian")
summary(elevmodel)

ggplot(elevmax, aes(latitude, max, col=year))+
  geom_point()+
  theme_bw()
