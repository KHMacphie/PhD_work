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
library(forcats)

#### plotting data with models ####

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)
cater$YSD <- paste(cater$year, cater$site, cater$date)
cater.date <- cater[,1:10]
cater.date$YSD <- paste(cater.date$year, cater.date$site, cater.date$date)
cater.date[,3:4] <- NULL
cater.date[,4:7] <- NULL
cater.date$plus1 <- cater.date$caterpillars+1
cater.date$logcater <- log(cater.date$plus1)
cater.date <- aggregate(cater.date$logcater~cater.date$YSD, FUN=mean)
names(cater.date)[c(1:2)] <-  c("YSD", "mean.logcaterpillars")

store <- pmatch(cater.date$YSD, cater$YSD)
cater.date <- cbind(cater.date, cater$year[store], cater$site[store], cater$date[store])
names(cater.date)[c(3:5)] <-  c("year", "site", "date")
cater.date$mean.caters <- exp(cater.date$mean.logcaterpillars)-1
plot(cater.date$date, cater.date$mean.caters)

site <- read.csv("Dropbox/master_data/site/site_details.csv")
store2 <- pmatch(cater.date$site, site$site, duplicates.ok = TRUE)
cater.date <- cbind(cater.date, site$Mean.Lat[store2], site$Mean.Elev[store2])
names(cater.date)[c(7:8)] <-  c("latitude", "elevation")

cater2018 <- cater.date[which(cater.date$year=='2018'), ]
ggplot(cater2018, aes(date, mean.caters, col=site))+
  geom_point()+
  theme_bw()

cater18 <- cater2018[order(cater2018$elevation),]
ggplot(cater18, aes(date, mean.caters, col=fct_inorder(site)))+
  geom_point()+
  theme_bw()

date18 <- aggregate(cater18$mean.caters~cater18$date, FUN=mean)
names(date18)[c(1:2)] <-  c("date", "mean.caterpillars")
plot(date18$date, date18$mean.caterpillars)

cater$YD <- paste(cater$year, cater$date)
cater$logcater <- log(cater$caterpillars+1)
allyears <- aggregate(cater$logcater~cater$YD, FUN=mean)
names(allyears)[c(1:2)] <-  c("YD", "meanlogcater")
store3 <- pmatch(allyears$YD, cater$YD)
allyears <- cbind(allyears, cater$year[store3], cater$date[store3])
names(allyears)[c(3:4)] <-  c("year", "date")
allyears$meancater <- exp(allyears$meanlogcater)-1

ggplot(allyears, aes(date, meancater, col=year))+
  geom_point()+
  theme_bw()

ggplot(allyears, aes(date, meancater, col=year))+
  geom_smooth()+
  theme_bw()

cater$sitetree <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$datecentred <- cater$date - mean(cater$date)

## core model to plot against yearday means
k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#basic<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2), random=~site+sitetree+siteday, family="poisson", data=cater, prior=prior2, nitt=200000, burnin=20000)
#save(basic, file = "~/Documents/Models/basic.RData")
load("~/Documents/Models/basic.RData")

predday <- seq(-26.4,28.6,0.5)
c14 <- mean(basic$Sol[,"(Intercept)"]) + mean(basic$Sol[,"datecentred"])*predday +
       mean(basic$Sol[,"I(datecentred^2)"])*predday^2
c15 <- mean(basic$Sol[,"(Intercept)"]+basic$Sol[,"year2015"]) + mean(basic$Sol[,"datecentred"]+basic$Sol[,"datecentred:year2015"])*predday +
  mean(basic$Sol[,"I(datecentred^2)"])*predday^2
c16 <- mean(basic$Sol[,"(Intercept)"]+basic$Sol[,"year2016"]) + mean(basic$Sol[,"datecentred"]+basic$Sol[,"datecentred:year2016"])*predday +
  mean(basic$Sol[,"I(datecentred^2)"])*predday^2
c17 <- mean(basic$Sol[,"(Intercept)"]+basic$Sol[,"year2017"]) + mean(basic$Sol[,"datecentred"]+basic$Sol[,"datecentred:year2017"])*predday +
  mean(basic$Sol[,"I(datecentred^2)"])*predday^2
c18 <- mean(basic$Sol[,"(Intercept)"]+basic$Sol[,"year2018"]) + mean(basic$Sol[,"datecentred"]+basic$Sol[,"datecentred:year2018"])*predday +
  mean(basic$Sol[,"I(datecentred^2)"])*predday^2
c19 <- mean(basic$Sol[,"(Intercept)"]+basic$Sol[,"year2019"]) + mean(basic$Sol[,"datecentred"]+basic$Sol[,"datecentred:year2019"])*predday +
  mean(basic$Sol[,"I(datecentred^2)"])*predday^2
predday2 <- seq(120,175,0.5)


plot(predday2, exp(c19), type="l", col=6, lwd=2, ylim=c(0,0.5))
points(predday2, exp(c18), type="l", col=5, lwd=2)
points(predday2, exp(c17), type="l", col=4, lwd=2)
points(predday2, exp(c16), type="l", col=3, lwd=2)
points(predday2, exp(c15), type="l", col=2, lwd=2)
points(predday2, exp(c14), type="l", col=1, lwd=2)
points(allyears$date, allyears$meancater, col=allyears$year)
legend(120,0.15,unique(allyears$year),col=1:length(allyears$year),pch=1)

# year * date^2
k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#basic2<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2)*year, random=~site+sitetree+siteday, family="poisson", data=cater, prior=prior2, nitt=200000, burnin=20000)
#save(basic2, file = "~/Documents/Models/basic2.RData")
load("~/Documents/Models/basic2.RData")

predday <- seq(-26.4,28.6,0.5)
c142 <- mean(basic2$Sol[,"(Intercept)"]) + mean(basic2$Sol[,"datecentred"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"])*predday^2
c152 <- mean(basic2$Sol[,"(Intercept)"]+basic2$Sol[,"year2015"]) + mean(basic2$Sol[,"datecentred"]+basic2$Sol[,"datecentred:year2015"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"]+basic2$Sol[,"year2015:I(datecentred^2)"])*predday^2
c162 <- mean(basic2$Sol[,"(Intercept)"]+basic2$Sol[,"year2016"]) + mean(basic2$Sol[,"datecentred"]+basic2$Sol[,"datecentred:year2016"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"]+basic2$Sol[,"year2016:I(datecentred^2)"])*predday^2
c172 <- mean(basic2$Sol[,"(Intercept)"]+basic2$Sol[,"year2017"]) + mean(basic2$Sol[,"datecentred"]+basic2$Sol[,"datecentred:year2017"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"]+basic2$Sol[,"year2017:I(datecentred^2)"])*predday^2
c182 <- mean(basic2$Sol[,"(Intercept)"]+basic2$Sol[,"year2018"]) + mean(basic2$Sol[,"datecentred"]+basic2$Sol[,"datecentred:year2018"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"]+basic2$Sol[,"year2018:I(datecentred^2)"])*predday^2
c192 <- mean(basic2$Sol[,"(Intercept)"]+basic2$Sol[,"year2019"]) + mean(basic2$Sol[,"datecentred"]+basic2$Sol[,"datecentred:year2019"])*predday +
  mean(basic2$Sol[,"I(datecentred^2)"]+basic2$Sol[,"year2019:I(datecentred^2)"])*predday^2
predday2 <- seq(120,175,0.5)


plot(predday2, exp(c192), type="l", col=6, lwd=2, ylim=c(0,0.5))
points(predday2, exp(c182), type="l", col=5, lwd=2)
points(predday2, exp(c172), type="l", col=4, lwd=2)
points(predday2, exp(c162), type="l", col=3, lwd=2)
points(predday2, exp(c152), type="l", col=2, lwd=2)
points(predday2, exp(c142), type="l", col=1, lwd=2)
points(allyears$date, allyears$meancater, col=allyears$year)
legend(120,0.15,unique(allyears$year),col=1:length(allyears$year),pch=1)

#### elevation peak height ####

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)
cater$YSD <- paste(cater$year, cater$site, cater$date)
cater.date <- cater[,1:10]
cater.date$YSD <- paste(cater.date$year, cater.date$site, cater.date$date)
cater.date[,3:4] <- NULL
cater.date[,4:7] <- NULL
cater.date$logcater <- log(cater.date$caterpillars+1)
cater.date <- aggregate(cater.date$logcater~cater.date$YSD, FUN=mean)
names(cater.date)[c(1:2)] <-  c("YSD", "mean.logcaterpillars")

store <- pmatch(cater.date$YSD, cater$YSD)
cater.date <- cbind(cater.date, cater$year[store], cater$site[store], cater$date[store])
names(cater.date)[c(3:5)] <-  c("year", "site", "date")
cater.date$mean.caters <- exp(cater.date$mean.logcaterpillars)-1
plot(cater.date$date, cater.date$mean.caters)

site <- read.csv("Dropbox/master_data/site/site_details.csv")
store2 <- pmatch(cater.date$site, site$site, duplicates.ok = TRUE)
cater.date <- cbind(cater.date, site$Mean.Lat[store2], site$Mean.Elev[store2])
names(cater.date)[c(7:8)] <-  c("latitude", "elevation")
cater.date$SY <- paste(cater.date$site, cater.date$year)

df.agg <- aggregate(mean.caters ~ SY, cater.date, max)
# then simply merge with the original
df.max <- merge(df.agg, cater.date)
plot(df.max$elevation, df.max$mean.caters, ylim=c(0,10), col=df.max$year)
legend(0,8,unique(df.max$year),col=1:length(df.max$year),pch=1)
plot(df.max$elevation, df.max$date, col=df.max$year)

library(lme4)
lm(df.max$mean.caters~df.max$elevation)
lm(df.max$date~df.max$elevation)
