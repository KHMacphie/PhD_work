rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)

cater$siteday <- paste(cater$site, cater$day, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$datecent <- cater$date-mean(cater$date)

a<-1000
prior<-list(R=list(V=diag(1), nu=0.002), 
               G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

DateCubed<- MCMCglmm(caterpillars~datecent*year+I(datecent^2)+I(datecent^3), random=~recorder+siteday+treeID+site, family="poisson", data=cater, prior=prior, nitt=250000, burnin=25000)
save(DateCubed, file = "~/Documents/Models/DateCubed.RData")
load("~/Documents/Models/DateCubed.RData") # DIC 20664.23 

DateSquared<- MCMCglmm(caterpillars~datecent*year+I(datecent^2), 
 random=~recorder+siteday+treeID+site, family="poisson", data=cater, prior=prior, nitt=250000, burnin=25000)
save(DateSquared, file = "~/Documents/Models/DateSquared.RData")
load("~/Documents/Models/DateSquared.RData") # DIC 20689.56

predday <- seq(-30,30, 0.1)
cubic2014 <-mean(DateCubed$Sol[,1])+
			mean(DateCubed$Sol[,2])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3
			
cubic2015 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,3])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,10])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2016 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,4])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,11])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2017 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,5])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,12])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2018 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,6])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,13])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2019 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,7])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,14])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

			
sqrd2014 <- mean(DateSquared$Sol[,1])+
			mean(DateSquared$Sol[,2])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2015 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,3])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,9])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2016 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,4])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,10])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2017 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,5])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,11])*predday+
			mean(DateSquared$Sol[,8])*predday^2

sqrd2018 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,6])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,12])*predday+
			mean(DateSquared$Sol[,8])*predday^2

sqrd2019 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,7])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,13])*predday+
			mean(DateSquared$Sol[,8])*predday^2
	
	
par(mfcol=c(1,2))			
plot(predday,exp(sqrd2014), type="l", col=2, ylim=c(0,0.12))
points(predday,exp(sqrd2015), type="l", col=2)
points(predday,exp(sqrd2016), type="l", col=2)
points(predday,exp(sqrd2017), type="l", col=2)
points(predday,exp(sqrd2018), type="l", col=2)
points(predday,exp(sqrd2019), type="l", col=2)

plot(predday, exp(cubic2014), type="l", ylim=c(0,0.12))
points(predday,exp(cubic2015), type="l", col=1)
points(predday,exp(cubic2016), type="l", col=1)
points(predday,exp(cubic2017), type="l", col=1)
points(predday,exp(cubic2018), type="l", col=1)
points(predday,exp(cubic2019), type="l", col=1)

#### Peak dates for cubic
a <- 3*mean(DateCubed$Sol[,9])
b <- 2*mean(DateCubed$Sol[,8])
c <- mean(DateCubed$Sol[,2])

abline(v=((-b - sqrt(b^2-4*a*c))/(2*a)), col=2)
peak14c <- 
 

