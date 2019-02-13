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
cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating.csv", 
                  col_types = cols(year = col_factor(levels = c("2014", 
                                                                "2015", "2016", "2017", "2018"))))
#cater <- read_csv("~/Dropbox/master_data/inverts/Branch_Beating.csv", 
#col_types = cols(year = col_factor(levels = c("2014", 
#                                             "2015", "2016", "2017", "2018"))))

temp <- read_csv("Dropbox/master_data/site/temperatures.csv", 
                 col_types = cols(year = col_factor(levels = c("2014", 
                                                               "2015", "2016", "2017", "2018"))))

cater<- mutate(cater, catbinom=caterpillars)
cater$catbinom[which(cater$caterpillars>1)]<-1
cater<-cater[-which(is.na(cater$catbinom)==TRUE),]
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


#### TempMCMCglmm with Apr*year ####

#k<-1000
#prior1<-list(R=list(V=1,fix=1),
#             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#TempYear<- MCMCglmm(catbinom~date*year+Apr*date+Apr*year+I(date^2), random=~site+sitetree+siteday, family="categorical", data=cater.temp, prior=prior1, nitt=200000)
#save(TempYear, file = "~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempYear.RData")
load("~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempYear.RData")


#### looking at model output ####

#summary(TempYear)
#plot(TempYear$VCV) #random effects
#plot(TempYear$Sol) #fixedeffects
#autocorr(TempYear$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
#TempYear.Sim<-simulate(TempYear,nsim=100)
#sum(cater.temp$catbinom)
#par(mfcol=c(1,1))
#hist(apply(TempYear.Sim,2,sum))
#abline(v=sum(cater.temp$catbinom),col=2)

#### Figure at 2 temperatures ####

##TempYear                 intercept               date*year                        I(date^2)                 date*Apr                    Apr*year
# Q1                                       
pred_day<-seq(120,175,1)
L.TempYear.2014 <-       -1.815e+02            +   1.957e+00*pred_day            + -4.942e-03*pred_day^2 + (6.935*-5.048e-02)*pred_day + 6.935*7.287e+00
L.TempYear.2015 <-      (-1.815e+02-1.549e+01) +  (1.957e+00+3.519e-02)*pred_day + -4.942e-03*pred_day^2 + (6.935*-5.048e-02)*pred_day + 6.935*(7.287e+00+9.600e-01)
L.TempYear.2016 <-      (-1.815e+02+6.748e+00) +  (1.957e+00-7.187e-02)*pred_day + -4.942e-03*pred_day^2 + (6.935*-5.048e-02)*pred_day + 6.935*(7.287e+00+2.168e-01)
L.TempYear.2017 <-      (-1.815e+02+7.720e+00) +  (1.957e+00-6.867e-02)*pred_day + -4.942e-03*pred_day^2 + (6.935*-5.048e-02)*pred_day + 6.935*(7.287e+00+1.118e-01)
L.TempYear.2018 <-      (-1.815e+02+3.225e+00) +  (1.957e+00-6.862e-02)*pred_day + -4.942e-03*pred_day^2 + (6.935*-5.048e-02)*pred_day + 6.935*(7.287e+00+7.997e-01)

# Q3
H.TempYear.2014 <-       -1.815e+02            +   1.957e+00*pred_day            + -4.942e-03*pred_day^2 + (8.268*-5.048e-02)*pred_day + 8.268*7.287e+00
H.TempYear.2015 <-      (-1.815e+02-1.549e+01) +  (1.957e+00+3.519e-02)*pred_day + -4.942e-03*pred_day^2 + (8.268*-5.048e-02)*pred_day + 8.268*(7.287e+00+9.600e-01)
H.TempYear.2016 <-      (-1.815e+02+6.748e+00) +  (1.957e+00-7.187e-02)*pred_day + -4.942e-03*pred_day^2 + (8.268*-5.048e-02)*pred_day + 8.268*(7.287e+00+2.168e-01)
H.TempYear.2017 <-      (-1.815e+02+7.720e+00) +  (1.957e+00-6.867e-02)*pred_day + -4.942e-03*pred_day^2 + (8.268*-5.048e-02)*pred_day + 8.268*(7.287e+00+1.118e-01)
H.TempYear.2018 <-      (-1.815e+02+3.225e+00) +  (1.957e+00-6.862e-02)*pred_day + -4.942e-03*pred_day^2 + (8.268*-5.048e-02)*pred_day + 8.268*(7.287e+00+7.997e-01)

par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(pred_day, plogis(L.TempYear.2014), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.45), axes=F) #black
points(pred_day, plogis(L.TempYear.2015), type="l",lwd=3, col=2) #red
points(pred_day, plogis(L.TempYear.2016), type="l",lwd=3,col=3) #green
points(pred_day, plogis(L.TempYear.2017), type="l", lwd=3, col=4)
points(pred_day, plogis(L.TempYear.2018), type="l", lwd=3, col=5)
#legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5))
legend("topright", legend="Q1", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(pred_day, plogis(H.TempYear.2014), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.45), axes=F) #black
points(pred_day, plogis(H.TempYear.2015), type="l",lwd=3, col=2) #red
points(pred_day, plogis(H.TempYear.2016), type="l",lwd=3,col=3) #green
points(pred_day, plogis(H.TempYear.2017), type="l", lwd=3, col=4)
points(pred_day, plogis(H.TempYear.2018), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="Q3", bty="n")
box()


#### Calculating peak dates ####  NOT RIGHT

#Peak for each year at Apr Q1, mean and Q3 with CI

low.peak.14 <- mean(-((TempYear$Sol[,"date"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
low.peak.14.CI <- HPDinterval(-((TempYear$Sol[,"date"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

low.peak.15 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
low.peak.15.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

low.peak.16 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
low.peak.16.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

low.peak.17 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
low.peak.17.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

low.peak.18 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
low.peak.18.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(6.935*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))


hi.peak.14 <- mean(-((TempYear$Sol[,"date"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                     (2*TempYear$Sol[,"I(date^2)"]))
hi.peak.14.CI <- HPDinterval(-((TempYear$Sol[,"date"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                               (2*TempYear$Sol[,"I(date^2)"]))

hi.peak.15 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                     (2*TempYear$Sol[,"I(date^2)"]))
hi.peak.15.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                               (2*TempYear$Sol[,"I(date^2)"]))

hi.peak.16 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                     (2*TempYear$Sol[,"I(date^2)"]))
hi.peak.16.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                               (2*TempYear$Sol[,"I(date^2)"]))

hi.peak.17 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                     (2*TempYear$Sol[,"I(date^2)"]))
hi.peak.17.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                               (2*TempYear$Sol[,"I(date^2)"]))

hi.peak.18 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                     (2*TempYear$Sol[,"I(date^2)"]))
hi.peak.18.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(8.268*TempYear$Sol[,"date:Apr"]))/
                               (2*TempYear$Sol[,"I(date^2)"]))


mid.peak.14 <- mean(-((TempYear$Sol[,"date"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
mid.peak.14.CI <- HPDinterval(-((TempYear$Sol[,"date"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

mid.peak.15 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
mid.peak.15.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2015"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

mid.peak.16 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
mid.peak.16.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2016"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

mid.peak.17 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
mid.peak.17.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2017"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

mid.peak.18 <- mean(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                      (2*TempYear$Sol[,"I(date^2)"]))
mid.peak.18.CI <- HPDinterval(-((TempYear$Sol[,"date"]+TempYear$Sol[,"date:year2018"])+(7.545*TempYear$Sol[,"date:Apr"]))/
                                (2*TempYear$Sol[,"I(date^2)"]))

TempYear.caterpeaks <- data_frame(year=as.factor (c(2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018)), 
                                 peakdate= c(low.peak.14, low.peak.15, low.peak.16, low.peak.17, low.peak.18, 
                                             hi.peak.14, hi.peak.15, hi.peak.16, hi.peak.17, hi.peak.18,
                                             mid.peak.14, mid.peak.15, mid.peak.16, mid.peak.17, mid.peak.18),
                                 temp= as.double(c(6.935,6.935,6.935,6.935,6.935,
                                                   8.268,8.268,8.268,8.268,8.268,
                                                   7.545,7.545,7.545,7.545,7.545)))

TempYear.caterpeaksCI <- data.frame(lowerCI= as.double(c(low.peak.14.CI["var1","lower"], low.peak.15.CI["var1","lower"], low.peak.16.CI["var1","lower"], low.peak.17.CI["var1","lower"], low.peak.18.CI["var1","lower"], 
                                                        hi.peak.14.CI["var1","lower"], hi.peak.15.CI["var1","lower"], hi.peak.16.CI["var1","lower"], hi.peak.17.CI["var1","lower"], hi.peak.18.CI["var1","lower"],
                                                        mid.peak.14.CI["var1","lower"], mid.peak.15.CI["var1","lower"], mid.peak.16.CI["var1","lower"], mid.peak.17.CI["var1","lower"], mid.peak.18.CI["var1","lower"])),
                                   upperCI= as.double(c(low.peak.14.CI["var1","upper"], low.peak.15.CI["var1","upper"], low.peak.16.CI["var1","upper"], low.peak.17.CI["var1","upper"], low.peak.18.CI["var1","upper"], 
                                                        hi.peak.14.CI["var1","upper"], hi.peak.15.CI["var1","upper"], hi.peak.16.CI["var1","upper"], hi.peak.17.CI["var1","upper"], hi.peak.18.CI["var1","upper"],
                                                        mid.peak.14.CI["var1","upper"], mid.peak.15.CI["var1","upper"], mid.peak.16.CI["var1","upper"], mid.peak.17.CI["var1","upper"], mid.peak.18.CI["var1","upper"])))

TempYear.caterpeaks$lowerCI <- TempYear.caterpeaksCI$lowerCI  
TempYear.caterpeaks$upperCI <- TempYear.caterpeaksCI$upperCI

ggplot(TempYear.caterpeaks, aes(temp, peakdate, colour= year))+
  geom_point(size=3, alpha=0.5)+
  geom_smooth(aes(ymin=lowerCI, ymax=upperCI, fill=year), stat='identity')+
  theme_bw()

#### Calculate B by year #### rthghkejfjwebkjrileutmvwelj;e 

# Calc B and CI  -date:Apr/2*quad
mean(-TempYear$Sol[,"date:Apr"]/(2*TempYear$Sol[,"I(date^2)"]))
HPDinterval(-MCMC.AprTempNQ$Sol[,"date:Apr"]/(2*MCMC.AprTempNQ$Sol[,"I(date^2)"]))


#using peak dates to compare..  NOT RIGHT
(hi.peak.14-low.peak.14)/(8.268-6.935) # = -5.123913 
(hi.peak.15-low.peak.15)/(8.268-6.935)
(hi.peak.16-low.peak.16)/(8.268-6.935)
(hi.peak.17-low.peak.17)/(8.268-6.935)
(hi.peak.18-low.peak.18)/(8.268-6.935)



##### Adding date*year*temp #####

#TempYearDate<- MCMCglmm(catbinom~date*year*Apr +I(date^2), random=~site+sitetree+siteday, family="categorical", data=cater.temp, prior=prior1, nitt=200000)
#save(TempYearDate, file = "~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempYearDate.RData")
load("~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempYearDate.RData")

#### looking at model output ####

summary(TempYearDate)
plot(TempYearDate$VCV) #random effects
plot(TempYearDate$Sol) #fixedeffects
autocorr(TempYearDate$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
TempYearDate.Sim<-simulate(TempYearDate,nsim=100)
sum(cater.temp$catbinom)
par(mfcol=c(1,1))
hist(apply(TempYearDate.Sim,2,sum))
abline(v=sum(cater.temp$catbinom),col=2)

#### Calculate B by year ####

# Calc B and CI  -date:Apr/2*quad
B2014 <- mean(-TempYearDate$Sol[,"date:Apr"]/(2*TempYearDate$Sol[,"I(date^2)"]))
B2014CI <- HPDinterval(-TempYearDate$Sol[,"date:Apr"]/(2*TempYearDate$Sol[,"I(date^2)"]))

B2015 <- mean(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))
B2015CI <- HPDinterval(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))

B2016 <- mean(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))
B2016CI <- HPDinterval(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))

B2017 <- mean(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))
B2017CI <- HPDinterval(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))

B2018 <- mean(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))
B2018CI <- HPDinterval(-(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])/(2*TempYearDate$Sol[,"I(date^2)"]))

BbyYear <- data.frame(year=as.factor(c(2014, 2015, 2016, 2017, 2018)),
                      B=c(B2014, B2015, B2016, B2017, B2018),
                      lowerCI=as.double(c(B2014CI["var1","lower"], B2015CI["var1","lower"], B2016CI["var1","lower"], B2017CI["var1","lower"], B2018CI["var1","lower"])),
                      upperCI=as.double(c(B2014CI["var1","upper"], B2015CI["var1","upper"], B2016CI["var1","upper"], B2017CI["var1","upper"], B2018CI["var1","upper"])))

View(BbyYear)

#### Figure at 2 temperatures TempYearDate ####

##TempYear                      intercept               date*year                        I(date^2)                 date*Apr                    Apr*year                   
# Q1                                       
pred_day<-seq(120,175,1)
L.TempYearDate.2014 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.935*TempYearDate$Sol[,"date:Apr"])*pred_day +  #date*Apr+yearchange
  mean(6.935*TempYearDate$Sol[,"Apr"])  #Apr+yearchange

L.TempYearDate.2015 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(6.935*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2015:Apr"]))  #Apr+yearchange

L.TempYearDate.2016 <-
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(6.935*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2016:Apr"]))  #Apr+yearchange

L.TempYearDate.2017 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(6.935*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2017:Apr"]))  #Apr+yearchange

L.TempYearDate.2018 <- 
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(6.935*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2018:Apr"]))  #Apr+yearchange

H.TempYearDate.2014 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.268*TempYearDate$Sol[,"date:Apr"])*pred_day +  #date*Apr+yearchange
  mean(8.268*TempYearDate$Sol[,"Apr"])  #Apr+yearchange

H.TempYearDate.2015 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(8.268*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2015:Apr"]))  #Apr+yearchange

H.TempYearDate.2016 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(8.268*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2016:Apr"]))  #Apr+yearchange

H.TempYearDate.2017 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(8.268*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2017:Apr"]))  #Apr+yearchange

H.TempYearDate.2018 <-  
  mean(TempYearDate$Sol[,"(Intercept)"]+TempYearDate$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempYearDate$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"]))*pred_day +  #date*Apr+yearchange
  mean(8.268*(TempYearDate$Sol[,"Apr"]+TempYearDate$Sol[,"year2018:Apr"]))  #Apr+yearchange


par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(pred_day, plogis(L.TempYearDate.2014), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.45), axes=F) #black
points(pred_day, plogis(L.TempYearDate.2015), type="l",lwd=3, col=2) #red
points(pred_day, plogis(L.TempYearDate.2016), type="l",lwd=3,col=3) #green
points(pred_day, plogis(L.TempYearDate.2017), type="l", lwd=3, col=4)
points(pred_day, plogis(L.TempYearDate.2018), type="l", lwd=3, col=5)
#legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5))
legend("topright", legend="Q1", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(pred_day, plogis(H.TempYearDate.2014), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.45), axes=F) #black
points(pred_day, plogis(H.TempYearDate.2015), type="l",lwd=3, col=2) #red
points(pred_day, plogis(H.TempYearDate.2016), type="l",lwd=3,col=3) #green
points(pred_day, plogis(H.TempYearDate.2017), type="l", lwd=3, col=4)
points(pred_day, plogis(H.TempYearDate.2018), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="Q3", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
title(ylab="Probability of caterpillar occurrence", outer=TRUE, line = 2)
title( xlab="Day of the year", outer=TRUE, line = 0)


#### Calculating peak dates TempYearDate####  

#Peak for each year at Apr Q1, mean and Q3 with CI

low.peak.14 <- mean(-((TempYearDate$Sol[,"date"])+(6.935*TempYearDate$Sol[,"date:Apr"]))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
low.peak.14.CI <- HPDinterval(-((TempYearDate$Sol[,"date"])+(6.935*TempYearDate$Sol[,"date:Apr"]))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

low.peak.15 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
low.peak.15.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

low.peak.16 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
low.peak.16.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

low.peak.17 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
low.peak.17.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

low.peak.18 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
low.peak.18.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(6.935*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))



hi.peak.14 <- mean(-((TempYearDate$Sol[,"date"])+(8.268*TempYearDate$Sol[,"date:Apr"]))/
                     (2*TempYearDate$Sol[,"I(date^2)"]))
hi.peak.14.CI <- HPDinterval(-((TempYearDate$Sol[,"date"])+(8.268*TempYearDate$Sol[,"date:Apr"]))/
                               (2*TempYearDate$Sol[,"I(date^2)"]))

hi.peak.15 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                     (2*TempYearDate$Sol[,"I(date^2)"]))
hi.peak.15.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                               (2*TempYearDate$Sol[,"I(date^2)"]))

hi.peak.16 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                     (2*TempYearDate$Sol[,"I(date^2)"]))
hi.peak.16.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                               (2*TempYearDate$Sol[,"I(date^2)"]))

hi.peak.17 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                     (2*TempYearDate$Sol[,"I(date^2)"]))
hi.peak.17.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                               (2*TempYearDate$Sol[,"I(date^2)"]))

hi.peak.18 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                     (2*TempYearDate$Sol[,"I(date^2)"]))
hi.peak.18.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(8.268*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                               (2*TempYearDate$Sol[,"I(date^2)"]))



mid.peak.14 <- mean(-((TempYearDate$Sol[,"date"])+(7.545*TempYearDate$Sol[,"date:Apr"]))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
mid.peak.14.CI <- HPDinterval(-((TempYearDate$Sol[,"date"])+(7.545*TempYearDate$Sol[,"date:Apr"]))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

mid.peak.15 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
mid.peak.15.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2015"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2015:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

mid.peak.16 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
mid.peak.16.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2016"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2016:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

mid.peak.17 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
mid.peak.17.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2017"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2017:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))

mid.peak.18 <- mean(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                      (2*TempYearDate$Sol[,"I(date^2)"]))
mid.peak.18.CI <- HPDinterval(-((TempYearDate$Sol[,"date"]+TempYearDate$Sol[,"date:year2018"])+(7.545*(TempYearDate$Sol[,"date:Apr"]+TempYearDate$Sol[,"date:year2018:Apr"])))/
                                (2*TempYearDate$Sol[,"I(date^2)"]))


TempYearDate.caterpeaks <- data_frame(year=as.factor (c(2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018, 2014, 2015, 2016, 2017, 2018)), 
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



ggplot(TempYearDate.caterpeaks, aes(temp, peakdate, colour= year))+
  geom_point(size=3, alpha=0.5)+
  geom_smooth(aes(ymin=lowerCI, ymax=upperCI, fill=year), stat='identity')+
  theme_bw()

#using peak dates to compare..
(hi.peak.14-low.peak.14)/(8.268-6.935)  #-4.483048
(hi.peak.15-low.peak.15)/(8.268-6.935)  #-3.624830
(hi.peak.16-low.peak.16)/(8.268-6.935)  #-4.512576
(hi.peak.17-low.peak.17)/(8.268-6.935)  #-4.512576
(hi.peak.18-low.peak.18)/(8.268-6.935)  #-5.007452


#### Without 2014 B's ####

TempYearDate.caterpeaks.No2014 <- filter(TempYearDate.caterpeaks, year == "2015"|year == "2016"|year == "2017"|year == "2018")
ggplot(TempYearDate.caterpeaks.No2014, aes(temp, peakdate, colour= year))+
  geom_point(size=3, alpha=0.5)+
  geom_smooth(aes(ymin=lowerCI, ymax=upperCI, fill=year), stat='identity')+
  theme_bw()

#scatter with CI bars for B
ggplot(BbyYear, aes(year, B))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, width=0.5))+
  theme_bw()
  
#without 2014 for closer look
BbyYear.No2014 <- filter(BbyYear, year == "2015"|year == "2016"|year == "2017"|year == "2018")
ggplot(BbyYear.No2014, aes(year, B))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=lowerCI, ymax=upperCI, width=0.5))+
  theme_bw()
  