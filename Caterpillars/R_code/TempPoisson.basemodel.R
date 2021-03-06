rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe #####

site <- read_csv("Dropbox/master_data/site/site_details.csv", 
                 col_types = cols(`Mean Elev` = col_double()))

#site <- read_csv("~/Dropbox/master_data/site/site_details.csv", 
#col_types = cols(`Mean Elev` = col_double()))
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
temp <- read.csv("Dropbox/master_data/site/temperatures.csv")

cater$year <- as.factor(cater$year)
temp$year <- as.factor(temp$year)

#cater<- mutate(cater, catbinom=caterpillars)
#cater$catbinom[which(cater$caterpillars>1)]<-1
#cater<-cater[-which(is.na(cater$catbinom)==TRUE),]
pmatch(cater$site,site$site,duplicates.ok=TRUE)
all_data<- merge(cater, site, by="site", duplicates.ok=TRUE)
all_data<- rename(all_data, latitude="Mean Lat")
all_data<- rename(all_data, longitude="Mean Long")
all_data<- rename(all_data, elevation="Mean Elev")
all_data$sitetree <- paste(all_data$tree, all_data$site)

temp$yearsite<- paste(temp$site, temp$year) #nests site and year?

mean_temps <- data.frame(site=temp$site, year=temp$year)
#View(mean_temps)
mean_temps <- mean_temps[!duplicated(mean_temps), ] #remove duplicated rows
mean_temps$yearsite <- paste(mean_temps$site, mean_temps$year)
pmatch(mean_temps$yearsite, temp$yearsite)
mean_temps <- mean_temps %>% arrange(site) #arrange by site to match means order

Means <- data.frame(site=mean_temps$site, year=mean_temps$year, yearsite=mean_temps$yearsite, Mar1=tapply(apply(temp[, 339:698],1, mean), temp$yearsite, mean))  #gives mean temperature between the two loggers per sie between the specified window
Means$Mar2 <-  tapply(apply(temp[, 699:1058],1, mean), temp$yearsite, mean)
Means$Apr1 <-  tapply(apply(temp[, 1059:1418],1, mean), temp$yearsite, mean)
Means$Apr2 <-  tapply(apply(temp[, 1419:1778],1, mean), temp$yearsite, mean)
Means$Mar  <-  tapply(apply(temp[, 339:1058],1, mean), temp$yearsite, mean)
Means$M2A1 <-  tapply(apply(temp[, 699:1418],1, mean), temp$yearsite, mean)
Means$Apr <-   tapply(apply(temp[, 1059:1778],1, mean), temp$yearsite, mean)
Means$MarA1 <- tapply(apply(temp[, 339:1418],1, mean), temp$yearsite, mean)
Means$M2Apr <- tapply(apply(temp[, 699:1778],1, mean), temp$yearsite, mean)
Means$MA <-    tapply(apply(temp[, 339:1778],1, mean), temp$yearsite, mean)

all_data$yearsite<- paste(all_data$site, all_data$year)
pmatch(all_data$yearsite, Means$yearsite)
Means <- select(Means, -site, -year)
cater.temp<- merge(all_data, Means, by="yearsite", duplicates.ok=TRUE)
cater.temp$sitetree <- paste(cater.temp$tree, cater.temp$site)
cater.temp$siteday <- paste(cater.temp$site, cater.temp$date, cater.temp$year)
cater.temp$obs<-as.factor(seq(1,length(cater.temp[,1])))


#### TempMCMCglmm with Apr*date ####

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#TempPoisson<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2), random=~site+sitetree+siteday, family="poisson", data=cater.temp, prior=prior2, nitt=200000, burnin=20000)
#save(TempPoisson, file = "~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempPoisson.RData")
load("~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempPoisson.RData")
load("/Users/s1205615/Documents/PhD (stopped using 12:2:19)/GitHub/R/Caterpillar analysis/results/TempPoisson.RData")
summary(TempPoisson)

plot(TempPoisson$VCV) #random effects
plot(TempPoisson$Sol) #fixedeffects
autocorr(TempPoisson$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
TempPoisson.Sim<-simulate(TempPoisson,nsim=100)
sum(cater.temp$caterpillars)
par(mfcol=c(1,1))
hist(apply(TempPoisson.Sim,2,sum))
abline(v=sum(cater.temp$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(TempPoisson.Sim,2,propzero))
abline(v=propzero(cater.temp$caterpillars), col="red")

#### B ####
mean(-TempPoisson$Sol[,"date:Apr"]/(2*TempPoisson$Sol[,"I(date^2)"])) # -4.688788
HPDinterval(-TempPoisson$Sol[,"date:Apr"]/(2*TempPoisson$Sol[,"I(date^2)"])) # -6.620825 to -2.688257

#### Figure at 2 temperatures TempPoisson ####

# Q1 for mean Apr temps input for lower temperature curves by year                                      
pred_day<-seq(120,175,1)
L.TempPoisson.2014 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.817*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.817*TempPoisson$Sol[,"Apr"]) #Apr


L.TempPoisson.2015 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.817*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.817*TempPoisson$Sol[,"Apr"])  #Apr

L.TempPoisson.2016 <-
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.817*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.817*TempPoisson$Sol[,"Apr"])  #Apr

L.TempPoisson.2017 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.817*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.817*TempPoisson$Sol[,"Apr"])  #Apr

L.TempPoisson.2018 <- 
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.817*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.817*TempPoisson$Sol[,"Apr"])  #Apr


M.TempPoisson.2014 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*TempPoisson$Sol[,"Apr"]) #Apr

M.TempPoisson.2015 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*TempPoisson$Sol[,"Apr"])  #Apr

M.TempPoisson.2016 <-
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*TempPoisson$Sol[,"Apr"])  #Apr

M.TempPoisson.2017 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*TempPoisson$Sol[,"Apr"])  #Apr

M.TempPoisson.2018 <- 
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*TempPoisson$Sol[,"Apr"])  #Apr


#Q3 for mean Apr temps input for higher temperature curves by year

H.TempPoisson.2014 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2015 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2016 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2017 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2018 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

par(mfcol=c(1,1),mar=c(2,1,1,1), oma=c(1,3,0,0),cex=1.8)
plot(pred_day, exp(H.TempPoisson.2014), type="l", lwd=3, xlab="Date", ylab="Abundance", ylim=c(0,0.17)) #black
points(pred_day, exp(H.TempPoisson.2015), type="l",lwd=3, col=2) #red
points(pred_day, exp(H.TempPoisson.2016), type="l",lwd=3,col=3) #green
points(pred_day, exp(H.TempPoisson.2017), type="l", lwd=3, col=4)
points(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.95, seg.len=0.8)
title(ylab="Abundance", outer=TRUE, line = 1)
title( xlab="Date", outer=TRUE, line = 0)
legend("topright", legend="8.4°C", bty="n")


par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(pred_day, exp(L.TempPoisson.2014), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.3), axes=F) #black
points(pred_day, exp(L.TempPoisson.2015), type="l",lwd=3, col=2) #red
points(pred_day, exp(L.TempPoisson.2016), type="l",lwd=3,col=3) #green
points(pred_day, exp(L.TempPoisson.2017), type="l", lwd=3, col=4)
points(pred_day, exp(L.TempPoisson.2018), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="6.8°C", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(pred_day, exp(H.TempPoisson.2014), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.3), axes=F) #black
points(pred_day, exp(H.TempPoisson.2015), type="l",lwd=3, col=2) #red
points(pred_day, exp(H.TempPoisson.2016), type="l",lwd=3,col=3) #green
points(pred_day, exp(H.TempPoisson.2017), type="l", lwd=3, col=4)
points(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=5)
legend("topright", legend="8.4°C", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
title(ylab="Caterpillar Abundance", outer=TRUE, line = 2)
title( xlab="Date", outer=TRUE, line = 0)

plot(pred_day, exp(H.TempPoisson.2014), type="l", lwd=3, xlab="Date", ylab="Abundance", ylim=c(0,0.27))
points(pred_day, exp(M.TempPoisson.2014), type="l", lwd=3, col=2)
points(pred_day, exp(L.TempPoisson.2014), type="l", lwd=3, col=3)
#### Calculating peak dates #### 

#Peak for each year at Apr Q1, mean and Q3 with CI

low.peak.14 <- mean(-((TempPoisson$Sol[,"date"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
low.peak.14.CI <- HPDinterval(-((TempPoisson$Sol[,"date"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

low.peak.15 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
low.peak.15.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

low.peak.16 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
low.peak.16.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

low.peak.17 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
low.peak.17.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

low.peak.18 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
low.peak.18.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(6.935*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))


hi.peak.14 <- mean(-((TempPoisson$Sol[,"date"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                     (2*TempPoisson$Sol[,"I(date^2)"]))
hi.peak.14.CI <- HPDinterval(-((TempPoisson$Sol[,"date"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                               (2*TempPoisson$Sol[,"I(date^2)"]))

hi.peak.15 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                     (2*TempPoisson$Sol[,"I(date^2)"]))
hi.peak.15.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                               (2*TempPoisson$Sol[,"I(date^2)"]))

hi.peak.16 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                     (2*TempPoisson$Sol[,"I(date^2)"]))
hi.peak.16.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                               (2*TempPoisson$Sol[,"I(date^2)"]))

hi.peak.17 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                     (2*TempPoisson$Sol[,"I(date^2)"]))
hi.peak.17.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                               (2*TempPoisson$Sol[,"I(date^2)"]))

hi.peak.18 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                     (2*TempPoisson$Sol[,"I(date^2)"]))
hi.peak.18.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(8.268*TempPoisson$Sol[,"date:Apr"]))/
                               (2*TempPoisson$Sol[,"I(date^2)"]))


mid.peak.14 <- mean(-((TempPoisson$Sol[,"date"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
mid.peak.14.CI <- HPDinterval(-((TempPoisson$Sol[,"date"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

mid.peak.15 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
mid.peak.15.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2015"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

mid.peak.16 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
mid.peak.16.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2016"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

mid.peak.17 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
mid.peak.17.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2017"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

mid.peak.18 <- mean(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                      (2*TempPoisson$Sol[,"I(date^2)"]))
mid.peak.18.CI <- HPDinterval(-((TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])+(7.545*TempPoisson$Sol[,"date:Apr"]))/
                                (2*TempPoisson$Sol[,"I(date^2)"]))

TempPoisson.caterpeaks <- data_frame(year=as.factor (c(2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018)), 
                                      peakdate= c(low.peak.14, low.peak.15, low.peak.16, low.peak.17, low.peak.18, 
                                                  hi.peak.14, hi.peak.15, hi.peak.16, hi.peak.17, hi.peak.18,
                                                  mid.peak.14, mid.peak.15, mid.peak.16, mid.peak.17, mid.peak.18),
                                      temp= as.double(c(6.935,6.935,6.935,6.935,6.935,
                                                        8.268,8.268,8.268,8.268,8.268,
                                                        7.545,7.545,7.545,7.545,7.545)),
                                      lowerCI= as.double(c(low.peak.14.CI["var1","lower"], low.peak.15.CI["var1","lower"], low.peak.16.CI["var1","lower"], low.peak.17.CI["var1","lower"], low.peak.18.CI["var1","lower"], 
                                                           hi.peak.14.CI["var1","lower"], hi.peak.15.CI["var1","lower"], hi.peak.16.CI["var1","lower"], hi.peak.17.CI["var1","lower"], hi.peak.18.CI["var1","lower"],
                                                           mid.peak.14.CI["var1","lower"], mid.peak.15.CI["var1","lower"], mid.peak.16.CI["var1","lower"], mid.peak.17.CI["var1","lower"], mid.peak.18.CI["var1","lower"])),
                                      upperCI= as.double(c(low.peak.14.CI["var1","upper"], low.peak.15.CI["var1","upper"], low.peak.16.CI["var1","upper"], low.peak.17.CI["var1","upper"], low.peak.18.CI["var1","upper"], 
                                                           hi.peak.14.CI["var1","upper"], hi.peak.15.CI["var1","upper"], hi.peak.16.CI["var1","upper"], hi.peak.17.CI["var1","upper"], hi.peak.18.CI["var1","upper"],
                                                           mid.peak.14.CI["var1","upper"], mid.peak.15.CI["var1","upper"], mid.peak.16.CI["var1","upper"], mid.peak.17.CI["var1","upper"], mid.peak.18.CI["var1","upper"])))



ggplot(TempPoisson.caterpeaks, aes(year, peakdate, colour= as.factor(temp)))+
  geom_point(size=3, alpha=0.5)+
  #geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, width=0.3))+
  theme_bw()

TempPoisson14 <- subset(TempPoisson.caterpeaks, year=="2014")

TempTiming <- ggplot(TempPoisson14, aes(x=temp, y=peakdate))+
  geom_point()+
  geom_smooth(stat="identity")+
  geom_errorbar(aes(ymax=upperCI, ymin=lowerCI, width=0.05))+
  theme_bw()+
  xlab("Temperature (°C)")+
  ylab("Peak Date")+
  theme(text = element_text(size=30))
TempTiming


# B
-mean(TempPoisson$Sol[,"date:Apr"]/(2*TempPoisson$Sol[,"I(date^2)"])) #-4.688788
-HPDinterval(TempPoisson$Sol[,"date:Apr"]/(2*TempPoisson$Sol[,"I(date^2)"])) # -2.688257 -6.620825
# graph gradient
(TempPoisson14[1,2]-TempPoisson14[2,2])/(TempPoisson14[1,3]-TempPoisson14[2,3]) #-4.688788


#### loop for posterior distribution of peak date 2014 lower temp ####
pred_daydetailed <- seq(120,175, 0.001)
#peakdateposteriorloop <- for(x in 1:18000) {
#  mean(TempPoisson$Sol[x,"(Intercept)"]) +  #Intercept+yearchange
#    mean(TempPoisson$Sol[x,"date"])*pred_daydetailed +  #date+yearchange
#    mean(TempPoisson$Sol[x,"I(date^2)"])*pred_daydetailed^2 +  #I(date^2)
#    mean(6.935*TempPoisson$Sol[x,"date:Apr"])*pred_daydetailed +  #date*Apr
#    mean(6.935*TempPoisson$Sol[x,"Apr"])
#}

peakheightposteriorloop<-c()
for(x in 1:18000) {
  peakheightposteriorloop[x] <- TempPoisson$Sol[x,"(Intercept)"] +  
  TempPoisson$Sol[x,"date"]*162.9719 +  
  TempPoisson$Sol[x,"I(date^2)"]*162.9719^2 +  
  6.935*TempPoisson$Sol[x,"date:Apr"]*162.9719 +  
  6.935*TempPoisson$Sol[x,"Apr"]
  }


#### Posterior distribution for peak date and height 2014 low temp####

postdistributions <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions$date <- TempPoisson$Sol[,"date"]
postdistributions$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions$peakdate <- -((postdistributions$date)+(6.935*postdistributions$dateApr))/(2*postdistributions$datesq)

plot(postdistributions$peakdate)


lowtemppeakheight2014 <-  
  exp(postdistributions$intercept +  #Intercept+yearchange
  postdistributions$date*postdistributions$peakdate +  #date+yearchange
  postdistributions$datesq*postdistributions$peakdate^2 +  #I(date^2)
  6.935*postdistributions$dateApr*postdistributions$peakdate +  #date*Apr
  6.935*postdistributions$Apr) #Apr

postdistributions$peakheight <- exp(postdistributions$intercept +  #Intercept+yearchange
                                      postdistributions$date*postdistributions$peakdate +  #date+yearchange
                                      postdistributions$datesq*postdistributions$peakdate^2 +  #I(date^2)
                                      6.935*postdistributions$dateApr*postdistributions$peakdate +  #date*Apr
                                      6.935*postdistributions$Apr) #Apr
plot(postdistributions$peakheight)

mean(postdistributions$peakdate) # 162.9719
median(postdistributions$peakdate) # 161.376
HPDinterval(postdistributions$peakdate) # lower 155.3744  upper 170.901
mean(postdistributions$peakheight) # ridiculous answer..???
median(postdistributions$peakheight) # 0.2515833
HPDinterval(postdistributions$peakheight) # lower 0.09813638 upper 0.4848575

##lots of ridiculously extreme values in the post distributions...



##### trying a 3d plot   temp 6.5-8.5 ####
df3d <- data.frame(date=seq(120,175,1))
df3d$'6.5' <-  
  exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
  mean(6.5*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
  mean(6.5*TempPoisson$Sol[,"Apr"]))  #Apr
df3d$'8.5' <-  
  exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
        mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
        mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
        mean(8.5*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
        mean(8.5*TempPoisson$Sol[,"Apr"]))  #Apr

df3d$'6.6' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(6.6*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(6.6*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'6.7' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(6.7*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(6.7*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'6.8' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(6.8*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(6.8*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'6.9' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(6.9*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(6.9*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.0' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.0*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.0*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.1' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.1*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.1*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.2' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.2*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.2*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.3' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.3*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.3*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.4' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.4*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.4*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.5' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.5*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.5*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.6' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.6*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.6*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.7' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.7*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.7*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.8' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.8*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.8*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'7.9' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(7.9*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(7.9*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.0' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.0*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.0*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.1' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.1*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.1*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.2' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.2*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.2*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.3' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.3*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.3*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.4' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.4*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.4*TempPoisson$Sol[,"Apr"]))  #Apr
 df3d$'8.5' <-  
  +   exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
            +         mean(TempPoisson$Sol[,"date"])*df3d$date +  #date+yearchange
            +         mean(TempPoisson$Sol[,"I(date^2)"])*df3d$date^2 +  #I(date^2)
            +         mean(8.5*TempPoisson$Sol[,"date:Apr"])*df3d$date +  #date*Apr
            +         mean(8.5*TempPoisson$Sol[,"Apr"]))  #Apr

library(tidyr)
df3dlong <- gather(df3d, temp, cater,2:22)
df3dlong$temp <- as.numeric(as.character(df3dlong$temp))
library(lattice)
cloud(df3dlong$cater ~ df3dlong$date*df3dlong$temp)
library(rgl)
library(car)
scatter3d(df3dlong$cater,df3dlong$date,df3dlong$temp) # does work but not in r studio


#### Model with temp^2 ####

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#TempSquPoisson<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2)+I(Apr^2), random=~site+sitetree+siteday, family="poisson", data=cater.temp, prior=prior2, nitt=250000, burnin=25000)
save(TempSquPoisson, file = "~/Documents/Models/TempSquPoisson.RData")
load("~/Documents/Models/TempSquPoisson.RData")
summary(TempSquPoisson)

pred_day<-seq(120,175,1)
plot(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                       +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                       +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                       +                        mean(6.817*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                       +                        mean(6.817*TempSquPoisson$Sol[,"Apr"]) +
                                                mean(6.817^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(7.317*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(7.317*TempSquPoisson$Sol[,"Apr"]) +
                       mean(7.317^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(7.817*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(7.817*TempSquPoisson$Sol[,"Apr"]) +
                       mean(7.817^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(8.217*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(8.217*TempSquPoisson$Sol[,"Apr"]) +
                       mean(8.217^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(8.517*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(8.517*TempSquPoisson$Sol[,"Apr"]) +
                       mean(8.517^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(8.817*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(8.817*TempSquPoisson$Sol[,"Apr"]) +
                       mean(8.817^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(9.217*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(9.217*TempSquPoisson$Sol[,"Apr"]) +
                       mean(9.217^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(10.17*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(10.17*TempSquPoisson$Sol[,"Apr"]) +
                       mean(10.17^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")
points(pred_day, exp(mean(TempSquPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                         +                        mean(TempSquPoisson$Sol[,"date"])*pred_day +  #date+yearchange
                         +                        mean(TempSquPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
                         +                        mean(10.87*TempSquPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
                         +                        mean(10.87*TempSquPoisson$Sol[,"Apr"]) +
                       mean(10.87^2*TempSquPoisson$Sol[,"I(Apr^2)"])), type="l", ylab="")

#### still quadratic but broader

#### plot peak height against temp
