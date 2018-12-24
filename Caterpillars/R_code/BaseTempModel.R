#################################
#### Visser paper dates temp ####
#################################


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

#cater<- mutate(cater, catbinom=caterpillars)
#cater$catbinom[which(cater$caterpillars>1)]<-1
#cater<-cater[-which(is.na(cater$catbinom)==TRUE),]
pmatch(cater$site,site$site,duplicates.ok=TRUE)
all_data<- merge(cater, site, by="site", duplicates.ok=TRUE)
all_data<- rename(all_data, latitude="Mean Lat")
all_data<- rename(all_data, longitude="Mean Long")
all_data<- rename(all_data, elevation="Mean Elev")
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data<-all_data[-which(is.na(all_data$caterpillars)==TRUE),] ## not used in VisserTempBaseModel

temp$yearsite<- paste(temp$site, temp$year) #nests site and year

#### Temperature data frame ####
mean_temps <- data.frame(site=temp$site, year=temp$year)
#View(mean_temps)
mean_temps <- mean_temps[!duplicated(mean_temps), ] #remove duplicated rows
mean_temps$yearsite <- paste(mean_temps$site, mean_temps$year)
pmatch(mean_temps$yearsite, temp$yearsite)
mean_temps <- mean_temps %>% arrange(site) #arrange by site to match means order

## Visser et al dates for temp cue of cater peak= 8th March (67) - 17th May (137)= columns 507 to 2209

mean_temps$VisserTemp <- tapply(apply(temp[, 507:2209],1, mean), temp$yearsite, mean) #calculating mean temp within time window for each logger then the mean of the two values per site

## Putting into full dataframe with all beating and site data
all_data$yearsite<- paste(all_data$site, all_data$year)
pmatch(all_data$yearsite, mean_temps$yearsite)
mean_temps <- select(mean_temps, -site, -year)
all_data<- merge(all_data, mean_temps, by="yearsite", duplicates.ok=TRUE)
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data$siteday <- paste(all_data$site, all_data$date, all_data$year)
all_data$obs<-as.factor(seq(1,length(all_data[,1])))

###############################################################################
##### MCMCglmm for base temp model with Visser paper date range mean temp #####

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#VisserTempBaseModel<- MCMCglmm(caterpillars~date*year+VisserTemp*date+I(date^2), random=~site+sitetree+siteday, family="poisson", data=all_data, prior=prior2, nitt=200000, burnin=20000)
#save(VisserTempBaseModel, file = "~/Dropbox/KirstyMacPhie/caterpillar analysis/results/VisserTempBaseModel.RData")
load("~/Dropbox/KirstyMacPhie/caterpillar analysis/results/VisserTempBaseModel.RData")
summary(VisserTempBaseModel)
## Higher DIC than Apr temp in model

plot(VisserTempBaseModel$VCV) #random effects
plot(VisserTempBaseModel$Sol) #fixedeffects
autocorr(VisserTempBaseModel$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
VisserTempBaseModel.Sim<-simulate(VisserTempBaseModel,nsim=100)
sum(all_data$caterpillars)
par(mfcol=c(1,1))
hist(apply(VisserTempBaseModel.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(VisserTempBaseModel.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")

summary(all_data$VisserTemp) ## looking at temp range in VisserTemp
q1 <- 6.786
q2 <- 7.344 #median
q3 <- 8.094

#### B for VisserTempBaseModel ####
mean(-VisserTempBaseModel$Sol[,"date:VisserTemp"]/(2*VisserTempBaseModel$Sol[,"I(date^2)"])) # -4.390252
HPDinterval(-VisserTempBaseModel$Sol[,"date:VisserTemp"]/(2*VisserTempBaseModel$Sol[,"I(date^2)"])) # lower:-7.186703, upper:-2.709484

######################################################################
#### Calc peak date and height with CI for each year at Q1/Med/Q3 ####

#### Posterior distribution for peak date and height 2014 ####

PD14 <- data_frame(intercept= VisserTempBaseModel$Sol[,"(Intercept)"])
PD14$date <- VisserTempBaseModel$Sol[,"date"]
PD14$datesq <- VisserTempBaseModel$Sol[,"I(date^2)"]
PD14$dateVisserTemp <- VisserTempBaseModel$Sol[,"date:VisserTemp"]
PD14$VisserTemp <-  VisserTempBaseModel$Sol[,"VisserTemp"]
PD14$Q1peakdate <- -((PD14$date)+(q1*PD14$dateVisserTemp))/(2*PD14$datesq)
PD14$Q1peakheight <- exp(PD14$intercept +  #Intercept+yearchange
                                          PD14$date*PD14$Q1peakdate +  #date+yearchange
                                          PD14$datesq*PD14$Q1peakdate^2 +  #I(date^2)
                                          q1*PD14$dateVisserTemp*PD14$Q1peakdate +  #date*VisserTemp
                                          q1*PD14$VisserTemp) #VisserTemp
PD14$Q2peakdate <- -((PD14$date)+(q2*PD14$dateVisserTemp))/(2*PD14$datesq)
PD14$Q2peakheight <- exp(PD14$intercept +  #Intercept+yearchange
                                          PD14$date*PD14$Q2peakdate +  #date+yearchange
                                          PD14$datesq*PD14$Q2peakdate^2 +  #I(date^2)
                                          q2*PD14$dateVisserTemp*PD14$Q2peakdate +  #date*VisserTemp
                                          q2*PD14$VisserTemp) #VisserTemp
PD14$Q3peakdate <- -((PD14$date)+(q3*PD14$dateVisserTemp))/(2*PD14$datesq)
PD14$Q3peakheight <- exp(PD14$intercept +  #Intercept+yearchange
                                          PD14$date*PD14$Q3peakdate +  #date+yearchange
                                          PD14$datesq*PD14$Q3peakdate^2 +  #I(date^2)
                                          q3*PD14$dateVisserTemp*PD14$Q3peakdate +  #date*VisserTemp
                                          q3*PD14$VisserTemp) #VisserTemp

Q1peakdate14 <- median(PD14$Q1peakdate) 
Q1peakdate14CI <- HPDinterval(PD14$Q1peakdate) 
Q1peakheight14 <- median(PD14$Q1peakheight) 
Q1peakheight14CI <- HPDinterval(PD14$Q1peakheight) 

Q2peakdate14 <- median(PD14$Q2peakdate) 
Q2peakdate14CI <- HPDinterval(PD14$Q2peakdate) 
Q2peakheight14 <- median(PD14$Q2peakheight) 
Q2peakheight14CI <- HPDinterval(PD14$Q2peakheight) 

Q3peakdate14 <- median(PD14$Q3peakdate) 
Q3peakdate14CI <- HPDinterval(PD14$Q3peakdate) 
Q3peakheight14 <- median(PD14$Q3peakheight) 
Q3peakheight14CI <- HPDinterval(PD14$Q3peakheight) 


#### Posterior distribution for peak date and height 2015 ####

PD15 <- data_frame(intercept= VisserTempBaseModel$Sol[,"(Intercept)"])
PD15$date <- VisserTempBaseModel$Sol[,"date"]
PD15$datesq <- VisserTempBaseModel$Sol[,"I(date^2)"]
PD15$dateVisserTemp <- VisserTempBaseModel$Sol[,"date:VisserTemp"]
PD15$VisserTemp <-  VisserTempBaseModel$Sol[,"VisserTemp"]
PD15$intercept15 <- VisserTempBaseModel$Sol[,"year2015"]
PD15$date15 <- VisserTempBaseModel$Sol[,"date:year2015"]
PD15$Q1peakdate <- -((PD15$date)+(q1*PD15$dateVisserTemp)+(PD15$date15))/(2*PD15$datesq)
PD15$Q1peakheight <- exp((PD15$intercept+PD15$intercept15) +  #Intercept+yearchange
                                          (PD15$date+PD15$date15)*PD15$Q1peakdate +  #date+yearchange
                                          PD15$datesq*PD15$Q1peakdate^2 +  #I(date^2)
                                          q1*PD15$dateVisserTemp*PD15$Q1peakdate +  #date*VisserTemp
                                          q1*PD15$VisserTemp) #VisserTemp
PD15$Q2peakdate <- -((PD15$date)+(q2*PD15$dateVisserTemp)+(PD15$date15))/(2*PD15$datesq)
PD15$Q2peakheight <- exp((PD15$intercept+PD15$intercept15) +  #Intercept+yearchange
                                          (PD15$date+PD15$date15)*PD15$Q2peakdate +  #date+yearchange
                                          PD15$datesq*PD15$Q2peakdate^2 +  #I(date^2)
                                          q2*PD15$dateVisserTemp*PD15$Q2peakdate +  #date*VisserTemp
                                          q2*PD15$VisserTemp) #VisserTemp
PD15$Q3peakdate <- -((PD15$date)+(q3*PD15$dateVisserTemp)+(PD15$date15))/(2*PD15$datesq)
PD15$Q3peakheight <- exp((PD15$intercept+PD15$intercept15) +  #Intercept+yearchange
                                          (PD15$date+PD15$date15)*PD15$Q3peakdate +  #date+yearchange
                                          PD15$datesq*PD15$Q3peakdate^2 +  #I(date^2)
                                          q3*PD15$dateVisserTemp*PD15$Q3peakdate +  #date*VisserTemp
                                          q3*PD15$VisserTemp) #VisserTemp

Q1peakdate15 <- median(PD15$Q1peakdate) 
Q1peakdate15CI <- HPDinterval(PD15$Q1peakdate) 
Q1peakheight15 <- median(PD15$Q1peakheight) 
Q1peakheight15CI <- HPDinterval(PD15$Q1peakheight) 

Q2peakdate15 <- median(PD15$Q2peakdate) 
Q2peakdate15CI <- HPDinterval(PD15$Q2peakdate) 
Q2peakheight15 <- median(PD15$Q2peakheight) 
Q2peakheight15CI <- HPDinterval(PD15$Q2peakheight) 

Q3peakdate15 <- median(PD15$Q3peakdate) 
Q3peakdate15CI <- HPDinterval(PD15$Q3peakdate) 
Q3peakheight15 <- median(PD15$Q3peakheight) 
Q3peakheight15CI <- HPDinterval(PD15$Q3peakheight) 

#### Posterior distribution for peak date and height 2016 ####

PD16 <- data_frame(intercept= VisserTempBaseModel$Sol[,"(Intercept)"])
PD16$date <- VisserTempBaseModel$Sol[,"date"]
PD16$datesq <- VisserTempBaseModel$Sol[,"I(date^2)"]
PD16$dateVisserTemp <- VisserTempBaseModel$Sol[,"date:VisserTemp"]
PD16$VisserTemp <-  VisserTempBaseModel$Sol[,"VisserTemp"]
PD16$intercept16 <- VisserTempBaseModel$Sol[,"year2016"]
PD16$date16 <- VisserTempBaseModel$Sol[,"date:year2016"]
PD16$Q1peakdate <- -((PD16$date)+(q1*PD16$dateVisserTemp)+(PD16$date16))/(2*PD16$datesq)
PD16$Q1peakheight <- exp((PD16$intercept+PD16$intercept16) +  #Intercept+yearchange
                                          (PD16$date+PD16$date16)*PD16$Q1peakdate +  #date+yearchange
                                          PD16$datesq*PD16$Q1peakdate^2 +  #I(date^2)
                                          q1*PD16$dateVisserTemp*PD16$Q1peakdate +  #date*VisserTemp
                                          q1*PD16$VisserTemp) #VisserTemp
PD16$Q2peakdate <- -((PD16$date)+(q2*PD16$dateVisserTemp)+(PD16$date16))/(2*PD16$datesq)
PD16$Q2peakheight <- exp((PD16$intercept+PD16$intercept16) +  #Intercept+yearchange
                                          (PD16$date+PD16$date16)*PD16$Q2peakdate +  #date+yearchange
                                          PD16$datesq*PD16$Q2peakdate^2 +  #I(date^2)
                                          q2*PD16$dateVisserTemp*PD16$Q2peakdate +  #date*VisserTemp
                                          q2*PD16$VisserTemp) #VisserTemp
PD16$Q3peakdate <- -((PD16$date)+(q3*PD16$dateVisserTemp)+(PD16$date16))/(2*PD16$datesq)
PD16$Q3peakheight <- exp((PD16$intercept+PD16$intercept16) +  #Intercept+yearchange
                                          (PD16$date+PD16$date16)*PD16$Q3peakdate +  #date+yearchange
                                          PD16$datesq*PD16$Q3peakdate^2 +  #I(date^2)
                                          q3*PD16$dateVisserTemp*PD16$Q3peakdate +  #date*VisserTemp
                                          q3*PD16$VisserTemp) #VisserTemp

Q1peakdate16 <- median(PD16$Q1peakdate) 
Q1peakdate16CI <- HPDinterval(PD16$Q1peakdate) 
Q1peakheight16 <- median(PD16$Q1peakheight) 
Q1peakheight16CI <- HPDinterval(PD16$Q1peakheight) 

Q2peakdate16 <- median(PD16$Q2peakdate) 
Q2peakdate16CI <- HPDinterval(PD16$Q2peakdate) 
Q2peakheight16 <- median(PD16$Q2peakheight) 
Q2peakheight16CI <- HPDinterval(PD16$Q2peakheight) 

Q3peakdate16 <- median(PD16$Q3peakdate) 
Q3peakdate16CI <- HPDinterval(PD16$Q3peakdate) 
Q3peakheight16 <- median(PD16$Q3peakheight) 
Q3peakheight16CI <- HPDinterval(PD16$Q3peakheight) 

#### Posterior distribution for peak date and height 2017 ####

PD17 <- data_frame(intercept= VisserTempBaseModel$Sol[,"(Intercept)"])
PD17$date <- VisserTempBaseModel$Sol[,"date"]
PD17$datesq <- VisserTempBaseModel$Sol[,"I(date^2)"]
PD17$dateVisserTemp <- VisserTempBaseModel$Sol[,"date:VisserTemp"]
PD17$VisserTemp <-  VisserTempBaseModel$Sol[,"VisserTemp"]
PD17$intercept17 <- VisserTempBaseModel$Sol[,"year2017"]
PD17$date17 <- VisserTempBaseModel$Sol[,"date:year2017"]
PD17$Q1peakdate <- -((PD17$date)+(q1*PD17$dateVisserTemp)+(PD17$date17))/(2*PD17$datesq)
PD17$Q1peakheight <- exp((PD17$intercept+PD17$intercept17) +  #Intercept+yearchange
                                          (PD17$date+PD17$date17)*PD17$Q1peakdate +  #date+yearchange
                                          PD17$datesq*PD17$Q1peakdate^2 +  #I(date^2)
                                          q1*PD17$dateVisserTemp*PD17$Q1peakdate +  #date*VisserTemp
                                          q1*PD17$VisserTemp) #VisserTemp
PD17$Q2peakdate <- -((PD17$date)+(q2*PD17$dateVisserTemp)+(PD17$date17))/(2*PD17$datesq)
PD17$Q2peakheight <- exp((PD17$intercept+PD17$intercept17) +  #Intercept+yearchange
                                          (PD17$date+PD17$date17)*PD17$Q2peakdate +  #date+yearchange
                                          PD17$datesq*PD17$Q2peakdate^2 +  #I(date^2)
                                          q2*PD17$dateVisserTemp*PD17$Q2peakdate +  #date*VisserTemp
                                          q2*PD17$VisserTemp) #VisserTemp
PD17$Q3peakdate <- -((PD17$date)+(q3*PD17$dateVisserTemp)+(PD17$date17))/(2*PD17$datesq)
PD17$Q3peakheight <- exp((PD17$intercept+PD17$intercept17) +  #Intercept+yearchange
                                          (PD17$date+PD17$date17)*PD17$Q3peakdate +  #date+yearchange
                                          PD17$datesq*PD17$Q3peakdate^2 +  #I(date^2)
                                          q3*PD17$dateVisserTemp*PD17$Q3peakdate +  #date*VisserTemp
                                          q3*PD17$VisserTemp) #VisserTemp

Q1peakdate17 <- median(PD17$Q1peakdate) 
Q1peakdate17CI <- HPDinterval(PD17$Q1peakdate) 
Q1peakheight17 <- median(PD17$Q1peakheight) 
Q1peakheight17CI <- HPDinterval(PD17$Q1peakheight) 

Q2peakdate17 <- median(PD17$Q2peakdate) 
Q2peakdate17CI <- HPDinterval(PD17$Q2peakdate) 
Q2peakheight17 <- median(PD17$Q2peakheight) 
Q2peakheight17CI <- HPDinterval(PD17$Q2peakheight) 

Q3peakdate17 <- median(PD17$Q3peakdate) 
Q3peakdate17CI <- HPDinterval(PD17$Q3peakdate) 
Q3peakheight17 <- median(PD17$Q3peakheight) 
Q3peakheight17CI <- HPDinterval(PD17$Q3peakheight) 

#### Posterior distribution for peak date and height 2018 ####

PD18 <- data_frame(intercept= VisserTempBaseModel$Sol[,"(Intercept)"])
PD18$date <- VisserTempBaseModel$Sol[,"date"]
PD18$datesq <- VisserTempBaseModel$Sol[,"I(date^2)"]
PD18$dateVisserTemp <- VisserTempBaseModel$Sol[,"date:VisserTemp"]
PD18$VisserTemp <-  VisserTempBaseModel$Sol[,"VisserTemp"]
PD18$intercept18 <- VisserTempBaseModel$Sol[,"year2018"]
PD18$date18 <- VisserTempBaseModel$Sol[,"date:year2018"]
PD18$Q1peakdate <- -((PD18$date)+(q1*PD18$dateVisserTemp)+(PD18$date18))/(2*PD18$datesq)
PD18$Q1peakheight <- exp((PD18$intercept+PD18$intercept18) +  #Intercept+yearchange
                                          (PD18$date+PD18$date18)*PD18$Q1peakdate +  #date+yearchange
                                          PD18$datesq*PD18$Q1peakdate^2 +  #I(date^2)
                                          q1*PD18$dateVisserTemp*PD18$Q1peakdate +  #date*VisserTemp
                                          q1*PD18$VisserTemp) #VisserTemp
PD18$Q2peakdate <- -((PD18$date)+(q2*PD18$dateVisserTemp)+(PD18$date18))/(2*PD18$datesq)
PD18$Q2peakheight <- exp((PD18$intercept+PD18$intercept18) +  #Intercept+yearchange
                                          (PD18$date+PD18$date18)*PD18$Q2peakdate +  #date+yearchange
                                          PD18$datesq*PD18$Q2peakdate^2 +  #I(date^2)
                                          q2*PD18$dateVisserTemp*PD18$Q2peakdate +  #date*VisserTemp
                                          q2*PD18$VisserTemp) #VisserTemp
PD18$Q3peakdate <- -((PD18$date)+(q3*PD18$dateVisserTemp)+(PD18$date18))/(2*PD18$datesq)
PD18$Q3peakheight <- exp((PD18$intercept+PD18$intercept18) +  #Intercept+yearchange
                                          (PD18$date+PD18$date18)*PD18$Q3peakdate +  #date+yearchange
                                          PD18$datesq*PD18$Q3peakdate^2 +  #I(date^2)
                                          q3*PD18$dateVisserTemp*PD18$Q3peakdate +  #date*VisserTemp
                                          q3*PD18$VisserTemp) #VisserTemp

Q1peakdate18 <- median(PD18$Q1peakdate) 
Q1peakdate18CI <- HPDinterval(PD18$Q1peakdate) 
Q1peakheight18 <- median(PD18$Q1peakheight) 
Q1peakheight18CI <- HPDinterval(PD18$Q1peakheight) 

Q2peakdate18 <- median(PD18$Q2peakdate) 
Q2peakdate18CI <- HPDinterval(PD18$Q2peakdate) 
Q2peakheight18 <- median(PD18$Q2peakheight) 
Q2peakheight18CI <- HPDinterval(PD18$Q2peakheight) 

Q3peakdate18 <- median(PD18$Q3peakdate) 
Q3peakdate18CI <- HPDinterval(PD18$Q3peakdate) 
Q3peakheight18 <- median(PD18$Q3peakheight) 
Q3peakheight18CI <- HPDinterval(PD18$Q3peakheight) 

#### Dataframe for peak date/height and CIs for each year at temps:Q1, Q2 & Q3 ####

VisserTempBaseModel.PeakInfo <- data_frame(year=as.factor(c(2014, 2014, 2014, 2015, 2015, 2015, 2016, 2016, 2016, 2017, 2017, 2017, 2018, 2018, 2018)),
                                   temp=as.double(c(6.786, 7.344, 8.094,
                                                    6.786, 7.344, 8.094,
                                                    6.786, 7.344, 8.094,
                                                    6.786, 7.344, 8.094,
                                                    6.786, 7.344, 8.094)),
                                   peakdate= c(Q1peakdate14, Q2peakdate14, Q3peakdate14,
                                               Q1peakdate15, Q2peakdate15, Q3peakdate15,
                                               Q1peakdate16, Q2peakdate16, Q3peakdate16,
                                               Q1peakdate17, Q2peakdate17, Q3peakdate17,
                                               Q1peakdate18, Q2peakdate18, Q3peakdate18),
                                   peakdatelowerCI= as.double(c(Q1peakdate14CI["var1","lower"], Q2peakdate14CI["var1","lower"], Q3peakdate14CI["var1","lower"],
                                                                Q1peakdate15CI["var1","lower"], Q2peakdate15CI["var1","lower"], Q3peakdate15CI["var1","lower"],
                                                                Q1peakdate16CI["var1","lower"], Q2peakdate16CI["var1","lower"], Q3peakdate16CI["var1","lower"],
                                                                Q1peakdate17CI["var1","lower"], Q2peakdate17CI["var1","lower"], Q3peakdate17CI["var1","lower"],
                                                                Q1peakdate18CI["var1","lower"], Q2peakdate18CI["var1","lower"], Q3peakdate18CI["var1","lower"])),
                                   peakdateupperCI= as.double(c(Q1peakdate14CI["var1","upper"], Q2peakdate14CI["var1","upper"], Q3peakdate14CI["var1","upper"],
                                                                Q1peakdate15CI["var1","upper"], Q2peakdate15CI["var1","upper"], Q3peakdate15CI["var1","upper"],
                                                                Q1peakdate16CI["var1","upper"], Q2peakdate16CI["var1","upper"], Q3peakdate16CI["var1","upper"],
                                                                Q1peakdate17CI["var1","upper"], Q2peakdate17CI["var1","upper"], Q3peakdate17CI["var1","upper"],
                                                                Q1peakdate18CI["var1","upper"], Q2peakdate18CI["var1","upper"], Q3peakdate18CI["var1","upper"])),
                                   peakheight= c(Q1peakheight14, Q2peakheight14, Q3peakheight14,
                                                 Q1peakheight15, Q2peakheight15, Q3peakheight15,
                                                 Q1peakheight16, Q2peakheight16, Q3peakheight16,
                                                 Q1peakheight17, Q2peakheight17, Q3peakheight17,
                                                 Q1peakheight18, Q2peakheight18, Q3peakheight18),
                                   peakheightlowerCI= as.double(c(Q1peakheight14CI["var1","lower"], Q2peakheight14CI["var1","lower"], Q3peakheight14CI["var1","lower"],
                                                                  Q1peakheight15CI["var1","lower"], Q2peakheight15CI["var1","lower"], Q3peakheight15CI["var1","lower"],
                                                                  Q1peakheight16CI["var1","lower"], Q2peakheight16CI["var1","lower"], Q3peakheight16CI["var1","lower"],
                                                                  Q1peakheight17CI["var1","lower"], Q2peakheight17CI["var1","lower"], Q3peakheight17CI["var1","lower"],
                                                                  Q1peakheight18CI["var1","lower"], Q2peakheight18CI["var1","lower"], Q3peakheight18CI["var1","lower"])),
                                   peakheightupperCI= as.double(c(Q1peakheight14CI["var1","upper"], Q2peakheight14CI["var1","upper"], Q3peakheight14CI["var1","upper"],
                                                                  Q1peakheight15CI["var1","upper"], Q2peakheight15CI["var1","upper"], Q3peakheight15CI["var1","upper"],
                                                                  Q1peakheight16CI["var1","upper"], Q2peakheight16CI["var1","upper"], Q3peakheight16CI["var1","upper"],
                                                                  Q1peakheight17CI["var1","upper"], Q2peakheight17CI["var1","upper"], Q3peakheight17CI["var1","upper"],
                                                                  Q1peakheight18CI["var1","upper"], Q2peakheight18CI["var1","upper"], Q3peakheight18CI["var1","upper"])))
View(VisserTempBaseModel.PeakInfo)



#############################################################################
#### Plot peak date/height with CI error bars for each year at Q1/Med/Q3 ####  

ggplot(VisserTempBaseModel.PeakInfo, aes(peakdate, peakheight, colour=year))+
  geom_point(size=1, alpha=0.5)+
  geom_errorbar(aes(ymin=peakheightlowerCI, ymax=peakheightupperCI))+
  geom_errorbarh(aes(xmin=peakdatelowerCI, xmax=peakdateupperCI))+
  facet_grid(. ~ temp)+
  xlab("Date")+
  ylab("Peak Abundance")+
  theme_bw()  
#ggsave("~/Dropbox/KirstyMacPhie/PhD_work_KM/Caterpillars/Results/VisserTemp_PeakDHCIsQ123.png")

##plot of peak height ad date by year and temp
ggplot(VisserTempBaseModel.PeakInfo, aes(peakdate, peakheight))+
  geom_point(aes(colour=as.factor(temp), shape=year))+
  theme_bw()

##plot of peak date at different temperatures by year
ggplot(VisserTempBaseModel.PeakInfo, aes(temp, peakdate))+
  geom_point(aes(colour=year))+
  theme_bw()

##plot of peak height at different temperatures by year
ggplot(VisserTempBaseModel.PeakInfo, aes(temp, peakheight, colour=year))+
  geom_point(aes(colour=year))+
  #geom_smooth(aes(ymin=peakheightlowerCI, ymax=peakheightupperCI, fill=year), stat="identity", width=0.1)+
  theme_bw()
##?? so temperature is having different effects on height of the peak in different years.. could that be because there are different temperatures in the different years

#peak height by year colour=temp
ggplot(VisserTempBaseModel.PeakInfo, aes(year, peakheight, colour=temp))+
  geom_point()+
  #geom_smooth(aes(ymin=peakheightlowerCI, ymax=peakheightupperCI, fill=year), stat="identity", width=0.1)+
  theme_bw()

#### figure of peak date by temp in 2018 ####
VisserTempBaseModel.PeakInfo18 <- VisserTempBaseModel.PeakInfo
VisserTempBaseModel.PeakInfo18 <- filter(VisserTempBaseModel.PeakInfo18, year=="2018")
ggplot(VisserTempBaseModel.PeakInfo18, aes(temp, peakdate))+
  geom_point(size=3, alpha=0.5)+
  geom_smooth(stat='identity')+
  geom_errorbar(aes(ymin=peakdatelowerCI, ymax=peakdateupperCI), width=0.1)+
  ylab("Peak Date")+
  xlab("Temperature ("~degree~"C)")+
  theme_bw()





#######################################################
#### Plot abundance over time by year at Q1 and Q3 ####

# Q1 from mean VisserTemp for lower temperature curves by year                                      
pred_day<-seq(120,175,0.5)
L.VisserTempBaseModel.2014 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q1*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q1*VisserTempBaseModel$Sol[,"VisserTemp"]) #VisserTemp


L.VisserTempBaseModel.2015 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q1*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q1*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

L.VisserTempBaseModel.2016 <-
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q1*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q1*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

L.VisserTempBaseModel.2017 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q1*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q1*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

L.VisserTempBaseModel.2018 <- 
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q1*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q1*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

#Q3 from mean VisserTemp for higher temperature curves by year
H.VisserTempBaseModel.2014 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q3*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q3*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

H.VisserTempBaseModel.2015 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2015"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2015"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q3*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q3*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

H.VisserTempBaseModel.2016 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2016"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2016"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q3*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q3*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

H.VisserTempBaseModel.2017 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2017"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2017"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q3*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q3*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp

H.VisserTempBaseModel.2018 <-  
  mean(VisserTempBaseModel$Sol[,"(Intercept)"]+VisserTempBaseModel$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(VisserTempBaseModel$Sol[,"date"]+VisserTempBaseModel$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(VisserTempBaseModel$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(q3*VisserTempBaseModel$Sol[,"date:VisserTemp"])*pred_day +  #date*VisserTemp
  mean(q3*VisserTempBaseModel$Sol[,"VisserTemp"])  #VisserTemp


par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(pred_day, exp(L.VisserTempBaseModel.2014), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.27), axes=F) #black
points(pred_day, exp(L.VisserTempBaseModel.2015), type="l",lwd=3, col=2) #red
points(pred_day, exp(L.VisserTempBaseModel.2016), type="l",lwd=3,col=3) #green
points(pred_day, exp(L.VisserTempBaseModel.2017), type="l", lwd=3, col=4)
points(pred_day, exp(L.VisserTempBaseModel.2018), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="V Q1", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(pred_day, exp(H.VisserTempBaseModel.2014), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.27), axes=F) #black
points(pred_day, exp(H.VisserTempBaseModel.2015), type="l",lwd=3, col=2) #red
points(pred_day, exp(H.VisserTempBaseModel.2016), type="l",lwd=3,col=3) #green
points(pred_day, exp(H.VisserTempBaseModel.2017), type="l", lwd=3, col=4)
points(pred_day, exp(H.VisserTempBaseModel.2018), type="l", lwd=3, col=5)
legend("topright", legend="V Q3", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
title(ylab="Caterpillar Abundance", outer=TRUE, line = 2)
title( xlab="Date", outer=TRUE, line = 0)

