rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe for beating, site info and temps #####

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

#### Base model for temp using Poisson ####

#TempPoisson<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2), random=~site+sitetree+siteday, family="poisson", data=cater.temp, prior=prior2, nitt=200000, burnin=20000)
#save(TempPoisson, file = "~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempPoisson.RData")
load("~/Dropbox/KirstyMacPhie/caterpillar analysis/results/TempPoisson.RData")
summary(TempPoisson)

mean(-postdistributions14$dateApr/(2*postdistributions14$datesq))

#### Posterior distribution for peak date and height 2014 ####

postdistributions14 <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions14$date <- TempPoisson$Sol[,"date"]
postdistributions14$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions14$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions14$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions14$Q1peakdate <- -((postdistributions14$date)+(6.817*postdistributions14$dateApr))/(2*postdistributions14$datesq)
postdistributions14$Q1peakheight <- exp(postdistributions14$intercept +  #Intercept+yearchange
                                      postdistributions14$date*postdistributions14$Q1peakdate +  #date+yearchange
                                      postdistributions14$datesq*postdistributions14$Q1peakdate^2 +  #I(date^2)
                                      6.817*postdistributions14$dateApr*postdistributions14$Q1peakdate +  #date*Apr
                                      6.817*postdistributions14$Apr) #Apr
postdistributions14$Q2peakdate <- -((postdistributions14$date)+(7.562*postdistributions14$dateApr))/(2*postdistributions14$datesq)
postdistributions14$Q2peakheight <- exp(postdistributions14$intercept +  #Intercept+yearchange
                                          postdistributions14$date*postdistributions14$Q2peakdate +  #date+yearchange
                                          postdistributions14$datesq*postdistributions14$Q2peakdate^2 +  #I(date^2)
                                          7.562*postdistributions14$dateApr*postdistributions14$Q2peakdate +  #date*Apr
                                          7.562*postdistributions14$Apr) #Apr
postdistributions14$Q3peakdate <- -((postdistributions14$date)+(8.442*postdistributions14$dateApr))/(2*postdistributions14$datesq)
postdistributions14$Q3peakheight <- exp(postdistributions14$intercept +  #Intercept+yearchange
                                          postdistributions14$date*postdistributions14$Q3peakdate +  #date+yearchange
                                          postdistributions14$datesq*postdistributions14$Q3peakdate^2 +  #I(date^2)
                                          8.442*postdistributions14$dateApr*postdistributions14$Q3peakdate +  #date*Apr
                                          8.442*postdistributions14$Apr) #Apr

Q1peakdate14 <- median(postdistributions14$Q1peakdate) 
Q1peakdate14CI <- HPDinterval(postdistributions14$Q1peakdate) 
Q1peakheight14 <- median(postdistributions14$Q1peakheight) 
Q1peakheight14CI <- HPDinterval(postdistributions14$Q1peakheight) 

Q2peakdate14 <- median(postdistributions14$Q2peakdate) 
Q2peakdate14CI <- HPDinterval(postdistributions14$Q2peakdate) 
Q2peakheight14 <- median(postdistributions14$Q2peakheight) 
Q2peakheight14CI <- HPDinterval(postdistributions14$Q2peakheight) 

Q3peakdate14 <- median(postdistributions14$Q3peakdate) 
Q3peakdate14CI <- HPDinterval(postdistributions14$Q3peakdate) 
Q3peakheight14 <- median(postdistributions14$Q3peakheight) 
Q3peakheight14CI <- HPDinterval(postdistributions14$Q3peakheight) 


#### Posterior distribution for peak date and height 2015 ####

postdistributions15 <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions15$date <- TempPoisson$Sol[,"date"]
postdistributions15$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions15$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions15$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions15$intercept15 <- TempPoisson$Sol[,"year2015"]
postdistributions15$date15 <- TempPoisson$Sol[,"date:year2015"]
postdistributions15$Q1peakdate <- -((postdistributions15$date)+(6.817*postdistributions15$dateApr)+(postdistributions15$date15))/(2*postdistributions15$datesq)
postdistributions15$Q1peakheight <- exp((postdistributions15$intercept+postdistributions15$intercept15) +  #Intercept+yearchange
                                          (postdistributions15$date+postdistributions15$date15)*postdistributions15$Q1peakdate +  #date+yearchange
                                          postdistributions15$datesq*postdistributions15$Q1peakdate^2 +  #I(date^2)
                                          6.817*postdistributions15$dateApr*postdistributions15$Q1peakdate +  #date*Apr
                                          6.817*postdistributions15$Apr) #Apr
postdistributions15$Q2peakdate <- -((postdistributions15$date)+(7.562*postdistributions15$dateApr)+(postdistributions15$date15))/(2*postdistributions15$datesq)
postdistributions15$Q2peakheight <- exp((postdistributions15$intercept+postdistributions15$intercept15) +  #Intercept+yearchange
                                          (postdistributions15$date+postdistributions15$date15)*postdistributions15$Q2peakdate +  #date+yearchange
                                          postdistributions15$datesq*postdistributions15$Q2peakdate^2 +  #I(date^2)
                                          7.562*postdistributions15$dateApr*postdistributions15$Q2peakdate +  #date*Apr
                                          7.562*postdistributions15$Apr) #Apr
postdistributions15$Q3peakdate <- -((postdistributions15$date)+(8.442*postdistributions15$dateApr)+(postdistributions15$date15))/(2*postdistributions15$datesq)
postdistributions15$Q3peakheight <- exp((postdistributions15$intercept+postdistributions15$intercept15) +  #Intercept+yearchange
                                          (postdistributions15$date+postdistributions15$date15)*postdistributions15$Q3peakdate +  #date+yearchange
                                          postdistributions15$datesq*postdistributions15$Q3peakdate^2 +  #I(date^2)
                                          8.442*postdistributions15$dateApr*postdistributions15$Q3peakdate +  #date*Apr
                                          8.442*postdistributions15$Apr) #Apr

Q1peakdate15 <- median(postdistributions15$Q1peakdate) 
Q1peakdate15CI <- HPDinterval(postdistributions15$Q1peakdate) 
Q1peakheight15 <- median(postdistributions15$Q1peakheight) 
Q1peakheight15CI <- HPDinterval(postdistributions15$Q1peakheight) 

Q2peakdate15 <- median(postdistributions15$Q2peakdate) 
Q2peakdate15CI <- HPDinterval(postdistributions15$Q2peakdate) 
Q2peakheight15 <- median(postdistributions15$Q2peakheight) 
Q2peakheight15CI <- HPDinterval(postdistributions15$Q2peakheight) 

Q3peakdate15 <- median(postdistributions15$Q3peakdate) 
Q3peakdate15CI <- HPDinterval(postdistributions15$Q3peakdate) 
Q3peakheight15 <- median(postdistributions15$Q3peakheight) 
Q3peakheight15CI <- HPDinterval(postdistributions15$Q3peakheight) 

#### Posterior distribution for peak date and height 2016 ####

postdistributions16 <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions16$date <- TempPoisson$Sol[,"date"]
postdistributions16$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions16$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions16$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions16$intercept16 <- TempPoisson$Sol[,"year2016"]
postdistributions16$date16 <- TempPoisson$Sol[,"date:year2016"]
postdistributions16$Q1peakdate <- -((postdistributions16$date)+(6.817*postdistributions16$dateApr)+(postdistributions16$date16))/(2*postdistributions16$datesq)
postdistributions16$Q1peakheight <- exp((postdistributions16$intercept+postdistributions16$intercept16) +  #Intercept+yearchange
                                          (postdistributions16$date+postdistributions16$date16)*postdistributions16$Q1peakdate +  #date+yearchange
                                          postdistributions16$datesq*postdistributions16$Q1peakdate^2 +  #I(date^2)
                                          6.817*postdistributions16$dateApr*postdistributions16$Q1peakdate +  #date*Apr
                                          6.817*postdistributions16$Apr) #Apr
postdistributions16$Q2peakdate <- -((postdistributions16$date)+(7.562*postdistributions16$dateApr)+(postdistributions16$date16))/(2*postdistributions16$datesq)
postdistributions16$Q2peakheight <- exp((postdistributions16$intercept+postdistributions16$intercept16) +  #Intercept+yearchange
                                          (postdistributions16$date+postdistributions16$date16)*postdistributions16$Q2peakdate +  #date+yearchange
                                          postdistributions16$datesq*postdistributions16$Q2peakdate^2 +  #I(date^2)
                                          7.562*postdistributions16$dateApr*postdistributions16$Q2peakdate +  #date*Apr
                                          7.562*postdistributions16$Apr) #Apr
postdistributions16$Q3peakdate <- -((postdistributions16$date)+(8.442*postdistributions16$dateApr)+(postdistributions16$date16))/(2*postdistributions16$datesq)
postdistributions16$Q3peakheight <- exp((postdistributions16$intercept+postdistributions16$intercept16) +  #Intercept+yearchange
                                          (postdistributions16$date+postdistributions16$date16)*postdistributions16$Q3peakdate +  #date+yearchange
                                          postdistributions16$datesq*postdistributions16$Q3peakdate^2 +  #I(date^2)
                                          8.442*postdistributions16$dateApr*postdistributions16$Q3peakdate +  #date*Apr
                                          8.442*postdistributions16$Apr) #Apr

Q1peakdate16 <- median(postdistributions16$Q1peakdate) 
Q1peakdate16CI <- HPDinterval(postdistributions16$Q1peakdate) 
Q1peakheight16 <- median(postdistributions16$Q1peakheight) 
Q1peakheight16CI <- HPDinterval(postdistributions16$Q1peakheight) 

Q2peakdate16 <- median(postdistributions16$Q2peakdate) 
Q2peakdate16CI <- HPDinterval(postdistributions16$Q2peakdate) 
Q2peakheight16 <- median(postdistributions16$Q2peakheight) 
Q2peakheight16CI <- HPDinterval(postdistributions16$Q2peakheight) 

Q3peakdate16 <- median(postdistributions16$Q3peakdate) 
Q3peakdate16CI <- HPDinterval(postdistributions16$Q3peakdate) 
Q3peakheight16 <- median(postdistributions16$Q3peakheight) 
Q3peakheight16CI <- HPDinterval(postdistributions16$Q3peakheight) 

#### Posterior distribution for peak date and height 2017 ####

postdistributions17 <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions17$date <- TempPoisson$Sol[,"date"]
postdistributions17$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions17$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions17$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions17$intercept17 <- TempPoisson$Sol[,"year2017"]
postdistributions17$date17 <- TempPoisson$Sol[,"date:year2017"]
postdistributions17$Q1peakdate <- -((postdistributions17$date)+(6.817*postdistributions17$dateApr)+(postdistributions17$date17))/(2*postdistributions17$datesq)
postdistributions17$Q1peakheight <- exp((postdistributions17$intercept+postdistributions17$intercept17) +  #Intercept+yearchange
                                          (postdistributions17$date+postdistributions17$date17)*postdistributions17$Q1peakdate +  #date+yearchange
                                          postdistributions17$datesq*postdistributions17$Q1peakdate^2 +  #I(date^2)
                                          6.817*postdistributions17$dateApr*postdistributions17$Q1peakdate +  #date*Apr
                                          6.817*postdistributions17$Apr) #Apr
postdistributions17$Q2peakdate <- -((postdistributions17$date)+(7.562*postdistributions17$dateApr)+(postdistributions17$date17))/(2*postdistributions17$datesq)
postdistributions17$Q2peakheight <- exp((postdistributions17$intercept+postdistributions17$intercept17) +  #Intercept+yearchange
                                          (postdistributions17$date+postdistributions17$date17)*postdistributions17$Q2peakdate +  #date+yearchange
                                          postdistributions17$datesq*postdistributions17$Q2peakdate^2 +  #I(date^2)
                                          7.562*postdistributions17$dateApr*postdistributions17$Q2peakdate +  #date*Apr
                                          7.562*postdistributions17$Apr) #Apr
postdistributions17$Q3peakdate <- -((postdistributions17$date)+(8.442*postdistributions17$dateApr)+(postdistributions17$date17))/(2*postdistributions17$datesq)
postdistributions17$Q3peakheight <- exp((postdistributions17$intercept+postdistributions17$intercept17) +  #Intercept+yearchange
                                          (postdistributions17$date+postdistributions17$date17)*postdistributions17$Q3peakdate +  #date+yearchange
                                          postdistributions17$datesq*postdistributions17$Q3peakdate^2 +  #I(date^2)
                                          8.442*postdistributions17$dateApr*postdistributions17$Q3peakdate +  #date*Apr
                                          8.442*postdistributions17$Apr) #Apr

Q1peakdate17 <- median(postdistributions17$Q1peakdate) 
Q1peakdate17CI <- HPDinterval(postdistributions17$Q1peakdate) 
Q1peakheight17 <- median(postdistributions17$Q1peakheight) 
Q1peakheight17CI <- HPDinterval(postdistributions17$Q1peakheight) 

Q2peakdate17 <- median(postdistributions17$Q2peakdate) 
Q2peakdate17CI <- HPDinterval(postdistributions17$Q2peakdate) 
Q2peakheight17 <- median(postdistributions17$Q2peakheight) 
Q2peakheight17CI <- HPDinterval(postdistributions17$Q2peakheight) 

Q3peakdate17 <- median(postdistributions17$Q3peakdate) 
Q3peakdate17CI <- HPDinterval(postdistributions17$Q3peakdate) 
Q3peakheight17 <- median(postdistributions17$Q3peakheight) 
Q3peakheight17CI <- HPDinterval(postdistributions17$Q3peakheight) 

#### Posterior distribution for peak date and height 2018 ####

postdistributions18 <- data_frame(intercept= TempPoisson$Sol[,"(Intercept)"])
postdistributions18$date <- TempPoisson$Sol[,"date"]
postdistributions18$datesq <- TempPoisson$Sol[,"I(date^2)"]
postdistributions18$dateApr <- TempPoisson$Sol[,"date:Apr"]
postdistributions18$Apr <-  TempPoisson$Sol[,"Apr"]
postdistributions18$intercept18 <- TempPoisson$Sol[,"year2018"]
postdistributions18$date18 <- TempPoisson$Sol[,"date:year2018"]
postdistributions18$Q1peakdate <- -((postdistributions18$date)+(6.817*postdistributions18$dateApr)+(postdistributions18$date18))/(2*postdistributions18$datesq)
postdistributions18$Q1peakheight <- exp((postdistributions18$intercept+postdistributions18$intercept18) +  #Intercept+yearchange
                                          (postdistributions18$date+postdistributions18$date18)*postdistributions18$Q1peakdate +  #date+yearchange
                                          postdistributions18$datesq*postdistributions18$Q1peakdate^2 +  #I(date^2)
                                          6.817*postdistributions18$dateApr*postdistributions18$Q1peakdate +  #date*Apr
                                          6.817*postdistributions18$Apr) #Apr
postdistributions18$Q2peakdate <- -((postdistributions18$date)+(7.562*postdistributions18$dateApr)+(postdistributions18$date18))/(2*postdistributions18$datesq)
postdistributions18$Q2peakheight <- exp((postdistributions18$intercept+postdistributions18$intercept18) +  #Intercept+yearchange
                                          (postdistributions18$date+postdistributions18$date18)*postdistributions18$Q2peakdate +  #date+yearchange
                                          postdistributions18$datesq*postdistributions18$Q2peakdate^2 +  #I(date^2)
                                          7.562*postdistributions18$dateApr*postdistributions18$Q2peakdate +  #date*Apr
                                          7.562*postdistributions18$Apr) #Apr
postdistributions18$Q3peakdate <- -((postdistributions18$date)+(8.442*postdistributions18$dateApr)+(postdistributions18$date18))/(2*postdistributions18$datesq)
postdistributions18$Q3peakheight <- exp((postdistributions18$intercept+postdistributions18$intercept18) +  #Intercept+yearchange
                                          (postdistributions18$date+postdistributions18$date18)*postdistributions18$Q3peakdate +  #date+yearchange
                                          postdistributions18$datesq*postdistributions18$Q3peakdate^2 +  #I(date^2)
                                          8.442*postdistributions18$dateApr*postdistributions18$Q3peakdate +  #date*Apr
                                          8.442*postdistributions18$Apr) #Apr

Q1peakdate18 <- median(postdistributions18$Q1peakdate) 
Q1peakdate18CI <- HPDinterval(postdistributions18$Q1peakdate) 
Q1peakheight18 <- median(postdistributions18$Q1peakheight) 
Q1peakheight18CI <- HPDinterval(postdistributions18$Q1peakheight) 

Q2peakdate18 <- median(postdistributions18$Q2peakdate) 
Q2peakdate18CI <- HPDinterval(postdistributions18$Q2peakdate) 
Q2peakheight18 <- median(postdistributions18$Q2peakheight) 
Q2peakheight18CI <- HPDinterval(postdistributions18$Q2peakheight) 

Q3peakdate18 <- median(postdistributions18$Q3peakdate) 
Q3peakdate18CI <- HPDinterval(postdistributions18$Q3peakdate) 
Q3peakheight18 <- median(postdistributions18$Q3peakheight) 
Q3peakheight18CI <- HPDinterval(postdistributions18$Q3peakheight) 

#### Dataframe for peak date/height and CIs for each year at temps:Q1, Q2 & Q3 ####

TempPoisson.PeakInfo <- data_frame(year=as.factor(c(2014, 2014, 2014, 2015, 2015, 2015, 2016, 2016, 2016, 2017, 2017, 2017, 2018, 2018, 2018)),
                                   temp=as.double(c(6.817, 7.562, 8.442,
                                                    6.817, 7.562, 8.442,
                                                    6.817, 7.562, 8.442,
                                                    6.817, 7.562, 8.442,
                                                    6.817, 7.562, 8.442)),
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
View(TempPoisson.PeakInfo)


##plot of peak height ad date by year and temp
ggplot(TempPoisson.PeakInfo, aes(peakdate, peakheight))+
  geom_point(aes(colour=as.factor(temp), shape=year))+
  theme_bw()

##plot of peak date at different temperatures by year
ggplot(TempPoisson.PeakInfo, aes(temp, peakdate))+
  geom_point(aes(colour=year))+
  theme_bw()

#### figure of peak height by date between years with CIs at median temp ####
TempPoisson.PeakInfoQ2 <- TempPoisson.PeakInfo
TempPoisson.PeakInfoQ2 <- as.factor(TempPoisson.PeakInfoQ2$temp)
TempPoisson.PeakInfoQ2 <- filter(TempPoisson.PeakInfoQ2, temp=="7.562")
TempPoisson.PeakInfoQ2
ggplot(TempPoisson.PeakInfoQ2, aes(peakdate, peakheight, colour=year))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=peakheightlowerCI, ymax=peakheightupperCI))+
  geom_errorbarh(aes(xmin=peakdatelowerCI, xmax=peakdateupperCI))+
  xlab("Date")+
  ylab("Abundance")+
  theme_bw()
  
TempPoisson.PeakInfoQ1Q2Q3 <- TempPoisson.PeakInfo
TempPoisson.PeakInfoQ1Q2Q3$Q <- c("Q1", "Q2", "Q3","Q1", "Q2", "Q3","Q1", "Q2", "Q3","Q1", "Q2", "Q3","Q1", "Q2", "Q3")
ggplot(TempPoisson.PeakInfoQ1Q2Q3, aes(peakdate, peakheight, colour=year))+
  geom_point(size=1, alpha=0.5)+
  geom_errorbar(aes(ymin=peakheightlowerCI, ymax=peakheightupperCI))+
  geom_errorbarh(aes(xmin=peakdatelowerCI, xmax=peakdateupperCI))+
  facet_grid(. ~ Q)+
  xlab("Date")+
  ylab("Abundance")+
  theme_bw()  
  
#### figure of peak date by temp in 2018 ####
TempPoisson.PeakInfo18 <- TempPoisson.PeakInfo
TempPoisson.PeakInfo18 <- filter(TempPoisson.PeakInfo18, year=="2018")
ggplot(TempPoisson.PeakInfo18, aes(temp, peakdate))+
  geom_point(size=3, alpha=0.5)+
  geom_smooth(stat='identity')+
  geom_errorbar(aes(ymin=peakdatelowerCI, ymax=peakdateupperCI), width=0.1)+
  ylab("Peak Date")+
  xlab("Temperature ("~degree~"C)")+
  theme_bw()

