#################################################
## Random effect output in base temp Apr model ##
#################################################


rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
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
all_data<-all_data[-which(is.na(all_data$caterpillars)==TRUE),] 

temp$yearsite<- paste(temp$site, temp$year) #nests site and year

#### Temperature data frame ####
mean_temps <- data.frame(site=temp$site, year=temp$year)
#View(mean_temps)
mean_temps <- mean_temps[!duplicated(mean_temps), ] #remove duplicated rows
mean_temps$yearsite <- paste(mean_temps$site, mean_temps$year)
pmatch(mean_temps$yearsite, temp$yearsite)
mean_temps <- mean_temps %>% arrange(site) #arrange by site to match means order

## Using mean temp through Apr

mean_temps$Apr <-   tapply(apply(temp[, 1059:1778],1, mean), temp$yearsite, mean) #calculating mean temp within time window for each logger then the mean of the two values per site

## Putting into full dataframe with all beating and site data
all_data$yearsite<- paste(all_data$site, all_data$year)
pmatch(all_data$yearsite, mean_temps$yearsite)
mean_temps <- select(mean_temps, -site, -year)
all_data<- merge(all_data, mean_temps, by="yearsite", duplicates.ok=TRUE)
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data$siteday <- paste(all_data$site, all_data$date, all_data$year)
all_data$obs<-as.factor(seq(1,length(all_data[,1])))

########################################################################################
##### MCMCglmm for base temp model with Apr temp and saving random effect variance #####

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#######################
## Added yearsite as a random effect with pr=TRUE also up to 300k iterations and 30k burnin
#AprBaseModel_RE has all previous REs, AprBaseModel_RE_2 has site day removed to reduce size

#AprBaseModel_RE<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2), random=~site+sitetree+siteday+yearsite, family="poisson", data=all_data, prior=prior3, nitt=300000, burnin=30000, pr=TRUE)
#AprBaseModel_RE_2<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2), random=~site+sitetree+yearsite, family="poisson", data=all_data, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)
#AprBaseModel_RE_3<- MCMCglmm(caterpillars~date*year+Apr*date+I(date^2), random=~site+sitetree+siteday, family="poisson", data=all_data, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)

#save(AprBaseModel_RE, file = "~/Documents/PhD/R/Caterpillar analysis/Models/AprBaseModel_RE.RData")
#save(AprBaseModel_RE_2, file = "~/Documents/PhD/R/Caterpillar analysis/Models/AprBaseModel_RE_2.RData")
#save(AprBaseModel_RE_3, file = "~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/AprBaseModel_RE_3.RData")
load("~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/AprBaseModel_RE.RData")
load("~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/AprBaseModel_RE_2.RData")
load("~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/AprBaseModel_RE_3.RData")

summary(AprBaseModel_RE) # DIC: 11394.1 
summary(AprBaseModel_RE_2) # DIC: 11541.33
summary(AprBaseModel_RE_3) # DIC: 11407.09

plot(AprBaseModel_RE$VCV) #random effects
#plot(AprBaseModel_RE$Sol) #fixedeffects   too many.. R shut down..
autocorr(AprBaseModel_RE$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
# AprBaseModel_RE
AprBaseModel_RE.Sim<-simulate(AprBaseModel_RE,nsim=100)
sum(all_data$caterpillars)
par(mfcol=c(1,1))
hist(apply(AprBaseModel_RE.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AprBaseModel_RE.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")

# AprBaseModel_RE_2
AprBaseModel_RE_2.Sim<-simulate(AprBaseModel_RE_2,nsim=100)
par(mfcol=c(1,1))
hist(apply(AprBaseModel_RE_2.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AprBaseModel_RE_2.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")

# AprBaseModel_RE_3
AprBaseModel_RE_3.Sim<-simulate(AprBaseModel_RE_3,nsim=100)
par(mfcol=c(1,1))
hist(apply(AprBaseModel_RE_3.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AprBaseModel_RE_3.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")


###################################################
#### Graphs showing different curves for sites ####

## looking at no.3 because original random effects
# plot 2014 as example with different sites

## base r graph all lines black
pred_day <- seq(120,175,0.5)
average2014 <- mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day +  #date+yearchange
  mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]) #Apr 

ALN <- mean(AprBaseModel_RE_3$Sol[,"site.ALN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
ART <- mean(AprBaseModel_RE_3$Sol[,"site.ART"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
AVI <- mean(AprBaseModel_RE_3$Sol[,"site.AVI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
AVN <- mean(AprBaseModel_RE_3$Sol[,"site.AVN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
BAD <- mean(AprBaseModel_RE_3$Sol[,"site.BAD"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
BIR <- mean(AprBaseModel_RE_3$Sol[,"site.BIR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
BLA <- mean(AprBaseModel_RE_3$Sol[,"site.BLA"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
BLG <- mean(AprBaseModel_RE_3$Sol[,"site.BLG"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
CAL <- mean(AprBaseModel_RE_3$Sol[,"site.CAL"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
CAR <- mean(AprBaseModel_RE_3$Sol[,"site.CAR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
CRU <- mean(AprBaseModel_RE_3$Sol[,"site.CRU"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DAV <- mean(AprBaseModel_RE_3$Sol[,"site.DAV"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DEL <- mean(AprBaseModel_RE_3$Sol[,"site.DEL"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DLW <- mean(AprBaseModel_RE_3$Sol[,"site.DLW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DNC <- mean(AprBaseModel_RE_3$Sol[,"site.DNC"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DNM <- mean(AprBaseModel_RE_3$Sol[,"site.DNM"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DNS <- mean(AprBaseModel_RE_3$Sol[,"site.DNS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DOR <- mean(AprBaseModel_RE_3$Sol[,"site.DOR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DOW <- mean(AprBaseModel_RE_3$Sol[,"site.DOW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
DUN <- mean(AprBaseModel_RE_3$Sol[,"site.DUN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
EDI <- mean(AprBaseModel_RE_3$Sol[,"site.EDI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
FOF <- mean(AprBaseModel_RE_3$Sol[,"site.FOF"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
FOU <- mean(AprBaseModel_RE_3$Sol[,"site.FOU"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
FSH <- mean(AprBaseModel_RE_3$Sol[,"site.FSH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
GLF <- mean(AprBaseModel_RE_3$Sol[,"site.GLF"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
HWP <- mean(AprBaseModel_RE_3$Sol[,"site.HWP"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
INS <- mean(AprBaseModel_RE_3$Sol[,"site.INS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
KCK <- mean(AprBaseModel_RE_3$Sol[,"site.KCK"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
KCZ <- mean(AprBaseModel_RE_3$Sol[,"site.KCZ"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
LVN <- mean(AprBaseModel_RE_3$Sol[,"site.LVN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
MCH <- mean(AprBaseModel_RE_3$Sol[,"site.MCH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
MUN <- mean(AprBaseModel_RE_3$Sol[,"site.MUN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
NEW <- mean(AprBaseModel_RE_3$Sol[,"site.NEW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
OSP <- mean(AprBaseModel_RE_3$Sol[,"site.OSP"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
PIT <- mean(AprBaseModel_RE_3$Sol[,"site.PIT"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
PTH <- mean(AprBaseModel_RE_3$Sol[,"site.PTH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
RSY <- mean(AprBaseModel_RE_3$Sol[,"site.RSY"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
RTH <- mean(AprBaseModel_RE_3$Sol[,"site.RTH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
SER <- mean(AprBaseModel_RE_3$Sol[,"site.SER"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
SLS <- mean(AprBaseModel_RE_3$Sol[,"site.SLS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
SPD <- mean(AprBaseModel_RE_3$Sol[,"site.SPD"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
STY <- mean(AprBaseModel_RE_3$Sol[,"site.STY"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
TAI <- mean(AprBaseModel_RE_3$Sol[,"site.TAI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])
TOM <- mean(AprBaseModel_RE_3$Sol[,"site.TOM"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*pred_day + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*pred_day^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*pred_day + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"])


plot(pred_day, exp(ALN), type="l", ylim=c(0,1.2), ylab="Caterpillars", xlab = "Date")
points(pred_day, exp(ART), type="l")
points(pred_day, exp(AVI), type="l")
points(pred_day, exp(AVN), type="l")
points(pred_day, exp(BAD), type="l")
points(pred_day, exp(BIR), type="l")
points(pred_day, exp(BLA), type="l")
points(pred_day, exp(BLG), type="l")
points(pred_day, exp(CAL), type="l")
points(pred_day, exp(CAR), type="l")
points(pred_day, exp(CRU), type="l")
points(pred_day, exp(DAV), type="l")
points(pred_day, exp(DEL), type="l")
points(pred_day, exp(DLW), type="l")
points(pred_day, exp(DNC), type="l")
points(pred_day, exp(DNM), type="l")
points(pred_day, exp(DNS), type="l")
points(pred_day, exp(DOR), type="l")
points(pred_day, exp(DOW), type="l")
points(pred_day, exp(DUN), type="l")
points(pred_day, exp(EDI), type="l")
points(pred_day, exp(FOF), type="l")
points(pred_day, exp(FOU), type="l")
points(pred_day, exp(FSH), type="l")
points(pred_day, exp(GLF), type="l")
points(pred_day, exp(HWP), type="l")
points(pred_day, exp(INS), type="l")
points(pred_day, exp(KCK), type="l")
points(pred_day, exp(KCZ), type="l")
points(pred_day, exp(LVN), type="l")
points(pred_day, exp(MCH), type="l")
points(pred_day, exp(MUN), type="l")
points(pred_day, exp(NEW), type="l")
points(pred_day, exp(OSP), type="l")
points(pred_day, exp(PIT), type="l")
points(pred_day, exp(PTH), type="l")
points(pred_day, exp(RSY), type="l")
points(pred_day, exp(RTH), type="l")
points(pred_day, exp(SER), type="l")
points(pred_day, exp(SLS), type="l")
points(pred_day, exp(SPD), type="l")
points(pred_day, exp(STY), type="l")
points(pred_day, exp(TAI), type="l")
points(pred_day, exp(TOM), type="l")

## ggplot
BySite <- data.frame(Date=seq(120,175,0.5))
BySite$ALN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.ALN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$ART <- exp(mean(AprBaseModel_RE_3$Sol[,"site.ART"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$AVI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.AVI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$AVN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.AVN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$BAD <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BAD"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$BIR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BIR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$BLA <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BLA"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$BLG <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BLG"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$CAL <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CAL"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$CAR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CAR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$CRU <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CRU"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DAV <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DAV"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DEL <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DEL"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DLW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DLW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DNC <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNC"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DNM <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNM"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DNS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DOR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DOR"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DOW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DOW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$DUN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DUN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$EDI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.EDI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$FOF <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FOF"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$FOU <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FOU"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$FSH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FSH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$GLF <- exp(mean(AprBaseModel_RE_3$Sol[,"site.GLF"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$HWP <- exp(mean(AprBaseModel_RE_3$Sol[,"site.HWP"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$INS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.INS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$KCK <- exp(mean(AprBaseModel_RE_3$Sol[,"site.KCK"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$KCZ <- exp(mean(AprBaseModel_RE_3$Sol[,"site.KCZ"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$LVN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.LVN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$MCH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.MCH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$MUN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.MUN"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$NEW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.NEW"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$OSP <- exp(mean(AprBaseModel_RE_3$Sol[,"site.OSP"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$PIT <- exp(mean(AprBaseModel_RE_3$Sol[,"site.PIT"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$PTH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.PTH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$RSY <- exp(mean(AprBaseModel_RE_3$Sol[,"site.RSY"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$RTH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.RTH"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$SER <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SER"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$SLS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SLS"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$SPD <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SPD"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$STY <- exp(mean(AprBaseModel_RE_3$Sol[,"site.STY"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$TAI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.TAI"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))
BySite$TOM <- exp(mean(AprBaseModel_RE_3$Sol[,"site.TOM"]) + mean(AprBaseModel_RE_3$Sol[,"(Intercept)"]) + mean(AprBaseModel_RE_3$Sol[,"date"])*BySite$Date + mean(AprBaseModel_RE_3$Sol[,"I(date^2)"])*BySite$Date^2 + mean(7.545*AprBaseModel_RE_3$Sol[,"date:Apr"])*BySite$Date + mean(7.545*AprBaseModel_RE_3$Sol[,"Apr"]))

BySiteLong <- gather(BySite, key= "Site", value="Caterpillars", select=2:45)

ggplot(BySiteLong, aes(Date, Caterpillars, colour=Site))+
  geom_smooth(aes(Date, Caterpillars, colour=Site), stat="identity", size=0.8)+
  theme_bw()

###########################################################
#### Calculating the icient for each site with CI #### 

#ALN ART AVI AVN BAD BIR BLA BLG CAL CAR CRU DAV DEL DLW DNC DNM DNS DOR DOW DUN EDI FOF FOU FSH GLF HWP INS KCK KCZ LVN MCH MUN NEW OSP PIT PTH RSY RTH SER SLS SPD STY TAI TOM )

### Further down with out exp() shows positive and negative this shows difference from mean intercept if multipled
ALN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.ALN"]))  
ART <- exp(mean(AprBaseModel_RE_3$Sol[,"site.ART"]))  
AVI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.AVI"]))  
AVN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.AVN"]))  
BAD <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BAD"]))  
BIR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BIR"]))  
BLA <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BLA"]))  
BLG <- exp(mean(AprBaseModel_RE_3$Sol[,"site.BLG"]))  
CAL <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CAL"]))  
CAR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CAR"]))  
CRU <- exp(mean(AprBaseModel_RE_3$Sol[,"site.CRU"]))  
DAV <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DAV"]))  
DEL <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DEL"]))  
DLW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DLW"]))  
DNC <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNC"]))  
DNM <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNM"]))  
DNS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DNS"]))  
DOR <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DOR"]))  
DOW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DOW"]))  
DUN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.DUN"]))  
EDI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.EDI"]))  
FOF <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FOF"]))  
FOU <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FOU"]))  
FSH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.FSH"]))  
GLF <- exp(mean(AprBaseModel_RE_3$Sol[,"site.GLF"]))  
HWP <- exp(mean(AprBaseModel_RE_3$Sol[,"site.HWP"]))  
INS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.INS"]))  
KCK <- exp(mean(AprBaseModel_RE_3$Sol[,"site.KCK"]))  
KCZ <- exp(mean(AprBaseModel_RE_3$Sol[,"site.KCZ"]))  
LVN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.LVN"]))  
MCH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.MCH"]))  
MUN <- exp(mean(AprBaseModel_RE_3$Sol[,"site.MUN"]))  
NEW <- exp(mean(AprBaseModel_RE_3$Sol[,"site.NEW"])) 
OSP <- exp(mean(AprBaseModel_RE_3$Sol[,"site.OSP"])) 
PIT <- exp(mean(AprBaseModel_RE_3$Sol[,"site.PIT"])) 
PTH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.PTH"])) 
RSY <- exp(mean(AprBaseModel_RE_3$Sol[,"site.RSY"])) 
RTH <- exp(mean(AprBaseModel_RE_3$Sol[,"site.RTH"]))
SER <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SER"])) 
SLS <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SLS"])) 
SPD <- exp(mean(AprBaseModel_RE_3$Sol[,"site.SPD"]))  
STY <- exp(mean(AprBaseModel_RE_3$Sol[,"site.STY"]))  
TAI <- exp(mean(AprBaseModel_RE_3$Sol[,"site.TAI"]))
TOM <- exp(mean(AprBaseModel_RE_3$Sol[,"site.TOM"])) 

ALNCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.ALN"]))  
ARTCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.ART"]))  
AVICIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.AVI"]))  
AVNCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.AVN"]))  
BADCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.BAD"]))  
BIRCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.BIR"]))  
BLACIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.BLA"]))  
BLGCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.BLG"]))  
CALCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.CAL"]))  
CARCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.CAR"]))  
CRUCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.CRU"]))  
DAVCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DAV"]))  
DELCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DEL"]))  
DLWCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DLW"]))  
DNCCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNC"]))  
DNMCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNM"]))  
DNSCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNS"]))  
DORCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DOR"]))  
DOWCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DOW"]))  
DUNCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.DUN"]))  
EDICIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.EDI"]))  
FOFCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.FOF"]))  
FOUCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.FOU"]))  
FSHCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.FSH"]))  
GLFCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.GLF"]))  
HWPCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.HWP"]))  
INSCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.INS"]))  
KCKCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.KCK"]))  
KCZCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.KCZ"]))  
LVNCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.LVN"]))  
MCHCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.MCH"]))  
MUNCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.MUN"]))  
NEWCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.NEW"])) 
OSPCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.OSP"])) 
PITCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.PIT"])) 
PTHCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.PTH"])) 
RSYCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.RSY"])) 
RTHCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.RTH"]))
SERCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.SER"])) 
SLSCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.SLS"])) 
SPDCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.SPD"]))  
STYCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.STY"]))  
TAICIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.TAI"]))
TOMCIs <- exp( HPDinterval(AprBaseModel_RE_3$Sol[,"site.TOM"])) 

Coeff <- data.frame(Site= c("ALN", "ART", "AVI", "AVN", "BAD", "BIR", "BLA", "BLG", "CAL", "CAR", "CRU", "DAV", "DEL", "DLW", "DNC", "DNM", "DNS", "DOR", "DOW", "DUN", "EDI", "FOF", "FOU", "FSH", "GLF", "HWP", "INS", "KCK", "KCZ", "LVN", "MCH", "MUN", "NEW", "OSP", "PIT", "PTH", "RSY", "RTH", "SER", "SLS", "SPD", "STY", "TAI", "TOM"),
                    MeanCoeff=c(ALN, ART, AVI, AVN, BAD, BIR, BLA, BLG, CAL, CAR, CRU, DAV, DEL, DLW, DNC, DNM, DNS, DOR, DOW, DUN, EDI, FOF, FOU, FSH, GLF, HWP, INS, KCK, KCZ, LVN, MCH, MUN, NEW, OSP, PIT, PTH, RSY, RTH, SER, SLS, SPD, STY, TAI, TOM),
                    LCI=c(ALNCIs["var1","lower"], ARTCIs["var1","lower"], AVICIs["var1","lower"], AVNCIs["var1","lower"], BADCIs["var1","lower"], BIRCIs["var1","lower"], BLACIs["var1","lower"], BLGCIs["var1","lower"], CALCIs["var1","lower"], CARCIs["var1","lower"], CRUCIs["var1","lower"], DAVCIs["var1","lower"], DELCIs["var1","lower"], DLWCIs["var1","lower"], DNCCIs["var1","lower"], DNMCIs["var1","lower"], DNSCIs["var1","lower"], DORCIs["var1","lower"], DOWCIs["var1","lower"], DUNCIs["var1","lower"], EDICIs["var1","lower"], FOFCIs["var1","lower"], FOUCIs["var1","lower"], FSHCIs["var1","lower"], GLFCIs["var1","lower"], HWPCIs["var1","lower"], INSCIs["var1","lower"], KCKCIs["var1","lower"], KCZCIs["var1","lower"], LVNCIs["var1","lower"], MCHCIs["var1","lower"], MUNCIs["var1","lower"], NEWCIs["var1","lower"], OSPCIs["var1","lower"], PITCIs["var1","lower"], PTHCIs["var1","lower"], RSYCIs["var1","lower"], RTHCIs["var1","lower"], SERCIs["var1","lower"], SLSCIs["var1","lower"], SPDCIs["var1","lower"], STYCIs["var1","lower"], TAICIs["var1","lower"], TOMCIs["var1","lower"]),
                    UCI=c(ALNCIs["var1","upper"], ARTCIs["var1","upper"], AVICIs["var1","upper"], AVNCIs["var1","upper"], BADCIs["var1","upper"], BIRCIs["var1","upper"], BLACIs["var1","upper"], BLGCIs["var1","upper"], CALCIs["var1","upper"], CARCIs["var1","upper"], CRUCIs["var1","upper"], DAVCIs["var1","upper"], DELCIs["var1","upper"], DLWCIs["var1","upper"], DNCCIs["var1","upper"], DNMCIs["var1","upper"], DNSCIs["var1","upper"], DORCIs["var1","upper"], DOWCIs["var1","upper"], DUNCIs["var1","upper"], EDICIs["var1","upper"], FOFCIs["var1","upper"], FOUCIs["var1","upper"], FSHCIs["var1","upper"], GLFCIs["var1","upper"], HWPCIs["var1","upper"], INSCIs["var1","upper"], KCKCIs["var1","upper"], KCZCIs["var1","upper"], LVNCIs["var1","upper"], MCHCIs["var1","upper"], MUNCIs["var1","upper"], NEWCIs["var1","upper"], OSPCIs["var1","upper"], PITCIs["var1","upper"], PTHCIs["var1","upper"], RSYCIs["var1","upper"], RTHCIs["var1","upper"], SERCIs["var1","upper"], SLSCIs["var1","upper"], SPDCIs["var1","upper"], STYCIs["var1","upper"], TAICIs["var1","upper"], TOMCIs["var1","upper"]))

ggplot(Coeff, aes(Site, MeanCoeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymin=LCI, ymax=UCI, width=0.5))+
  theme_bw()


## removing exp so positive and negative?!

ALN <- mean(AprBaseModel_RE_3$Sol[,"site.ALN"])  
ART <- mean(AprBaseModel_RE_3$Sol[,"site.ART"])  
AVI <- mean(AprBaseModel_RE_3$Sol[,"site.AVI"])  
AVN <- mean(AprBaseModel_RE_3$Sol[,"site.AVN"])  
BAD <- mean(AprBaseModel_RE_3$Sol[,"site.BAD"])  
BIR <- mean(AprBaseModel_RE_3$Sol[,"site.BIR"])  
BLA <- mean(AprBaseModel_RE_3$Sol[,"site.BLA"])  
BLG <- mean(AprBaseModel_RE_3$Sol[,"site.BLG"])  
CAL <- mean(AprBaseModel_RE_3$Sol[,"site.CAL"])  
CAR <- mean(AprBaseModel_RE_3$Sol[,"site.CAR"])  
CRU <- mean(AprBaseModel_RE_3$Sol[,"site.CRU"])  
DAV <- mean(AprBaseModel_RE_3$Sol[,"site.DAV"])  
DEL <- mean(AprBaseModel_RE_3$Sol[,"site.DEL"])  
DLW <- mean(AprBaseModel_RE_3$Sol[,"site.DLW"])  
DNC <- mean(AprBaseModel_RE_3$Sol[,"site.DNC"])  
DNM <- mean(AprBaseModel_RE_3$Sol[,"site.DNM"])  
DNS <- mean(AprBaseModel_RE_3$Sol[,"site.DNS"])  
DOR <- mean(AprBaseModel_RE_3$Sol[,"site.DOR"])  
DOW <- mean(AprBaseModel_RE_3$Sol[,"site.DOW"])  
DUN <- mean(AprBaseModel_RE_3$Sol[,"site.DUN"])  
EDI <- mean(AprBaseModel_RE_3$Sol[,"site.EDI"])  
FOF <- mean(AprBaseModel_RE_3$Sol[,"site.FOF"])  
FOU <- mean(AprBaseModel_RE_3$Sol[,"site.FOU"])  
FSH <- mean(AprBaseModel_RE_3$Sol[,"site.FSH"])  
GLF <- mean(AprBaseModel_RE_3$Sol[,"site.GLF"])  
HWP <- mean(AprBaseModel_RE_3$Sol[,"site.HWP"])  
INS <- mean(AprBaseModel_RE_3$Sol[,"site.INS"])  
KCK <- mean(AprBaseModel_RE_3$Sol[,"site.KCK"])  
KCZ <- mean(AprBaseModel_RE_3$Sol[,"site.KCZ"])  
LVN <- mean(AprBaseModel_RE_3$Sol[,"site.LVN"])  
MCH <- mean(AprBaseModel_RE_3$Sol[,"site.MCH"])  
MUN <- mean(AprBaseModel_RE_3$Sol[,"site.MUN"])  
NEW <- mean(AprBaseModel_RE_3$Sol[,"site.NEW"]) 
OSP <- mean(AprBaseModel_RE_3$Sol[,"site.OSP"]) 
PIT <- mean(AprBaseModel_RE_3$Sol[,"site.PIT"]) 
PTH <- mean(AprBaseModel_RE_3$Sol[,"site.PTH"]) 
RSY <- mean(AprBaseModel_RE_3$Sol[,"site.RSY"]) 
RTH <- mean(AprBaseModel_RE_3$Sol[,"site.RTH"])
SER <- mean(AprBaseModel_RE_3$Sol[,"site.SER"]) 
SLS <- mean(AprBaseModel_RE_3$Sol[,"site.SLS"]) 
SPD <- mean(AprBaseModel_RE_3$Sol[,"site.SPD"])  
STY <- mean(AprBaseModel_RE_3$Sol[,"site.STY"])  
TAI <- mean(AprBaseModel_RE_3$Sol[,"site.TAI"])
TOM <- mean(AprBaseModel_RE_3$Sol[,"site.TOM"]) 

ALNCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.ALN"])  
ARTCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.ART"])  
AVICIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.AVI"])  
AVNCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.AVN"])  
BADCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.BAD"])  
BIRCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.BIR"])  
BLACIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.BLA"])  
BLGCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.BLG"])  
CALCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.CAL"])  
CARCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.CAR"])  
CRUCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.CRU"])  
DAVCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DAV"])  
DELCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DEL"])  
DLWCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DLW"])  
DNCCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNC"])  
DNMCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNM"])  
DNSCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DNS"])  
DORCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DOR"])  
DOWCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DOW"])  
DUNCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.DUN"])  
EDICIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.EDI"])  
FOFCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.FOF"])  
FOUCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.FOU"])  
FSHCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.FSH"])  
GLFCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.GLF"])  
HWPCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.HWP"])  
INSCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.INS"])  
KCKCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.KCK"])  
KCZCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.KCZ"])  
LVNCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.LVN"])  
MCHCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.MCH"])  
MUNCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.MUN"])  
NEWCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.NEW"]) 
OSPCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.OSP"]) 
PITCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.PIT"]) 
PTHCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.PTH"]) 
RSYCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.RSY"]) 
RTHCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.RTH"])
SERCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.SER"]) 
SLSCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.SLS"]) 
SPDCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.SPD"])  
STYCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.STY"])
TAICIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.TAI"])
TOMCIs <-  HPDinterval(AprBaseModel_RE_3$Sol[,"site.TOM"]) 

Coeff <- data.frame(Site= c("ALN", "ART", "AVI", "AVN", "BAD", "BIR", "BLA", "BLG", "CAL", "CAR", "CRU", "DAV", "DEL", "DLW", "DNC", "DNM", "DNS", "DOR", "DOW", "DUN", "EDI", "FOF", "FOU", "FSH", "GLF", "HWP", "INS", "KCK", "KCZ", "LVN", "MCH", "MUN", "NEW", "OSP", "PIT", "PTH", "RSY", "RTH", "SER", "SLS", "SPD", "STY", "TAI", "TOM"),
                    MeanCoeff=c(ALN, ART, AVI, AVN, BAD, BIR, BLA, BLG, CAL, CAR, CRU, DAV, DEL, DLW, DNC, DNM, DNS, DOR, DOW, DUN, EDI, FOF, FOU, FSH, GLF, HWP, INS, KCK, KCZ, LVN, MCH, MUN, NEW, OSP, PIT, PTH, RSY, RTH, SER, SLS, SPD, STY, TAI, TOM),
                    LCI=c(ALNCIs["var1","lower"], ARTCIs["var1","lower"], AVICIs["var1","lower"], AVNCIs["var1","lower"], BADCIs["var1","lower"], BIRCIs["var1","lower"], BLACIs["var1","lower"], BLGCIs["var1","lower"], CALCIs["var1","lower"], CARCIs["var1","lower"], CRUCIs["var1","lower"], DAVCIs["var1","lower"], DELCIs["var1","lower"], DLWCIs["var1","lower"], DNCCIs["var1","lower"], DNMCIs["var1","lower"], DNSCIs["var1","lower"], DORCIs["var1","lower"], DOWCIs["var1","lower"], DUNCIs["var1","lower"], EDICIs["var1","lower"], FOFCIs["var1","lower"], FOUCIs["var1","lower"], FSHCIs["var1","lower"], GLFCIs["var1","lower"], HWPCIs["var1","lower"], INSCIs["var1","lower"], KCKCIs["var1","lower"], KCZCIs["var1","lower"], LVNCIs["var1","lower"], MCHCIs["var1","lower"], MUNCIs["var1","lower"], NEWCIs["var1","lower"], OSPCIs["var1","lower"], PITCIs["var1","lower"], PTHCIs["var1","lower"], RSYCIs["var1","lower"], RTHCIs["var1","lower"], SERCIs["var1","lower"], SLSCIs["var1","lower"], SPDCIs["var1","lower"], STYCIs["var1","lower"], TAICIs["var1","lower"], TOMCIs["var1","lower"]),
                    UCI=c(ALNCIs["var1","upper"], ARTCIs["var1","upper"], AVICIs["var1","upper"], AVNCIs["var1","upper"], BADCIs["var1","upper"], BIRCIs["var1","upper"], BLACIs["var1","upper"], BLGCIs["var1","upper"], CALCIs["var1","upper"], CARCIs["var1","upper"], CRUCIs["var1","upper"], DAVCIs["var1","upper"], DELCIs["var1","upper"], DLWCIs["var1","upper"], DNCCIs["var1","upper"], DNMCIs["var1","upper"], DNSCIs["var1","upper"], DORCIs["var1","upper"], DOWCIs["var1","upper"], DUNCIs["var1","upper"], EDICIs["var1","upper"], FOFCIs["var1","upper"], FOUCIs["var1","upper"], FSHCIs["var1","upper"], GLFCIs["var1","upper"], HWPCIs["var1","upper"], INSCIs["var1","upper"], KCKCIs["var1","upper"], KCZCIs["var1","upper"], LVNCIs["var1","upper"], MCHCIs["var1","upper"], MUNCIs["var1","upper"], NEWCIs["var1","upper"], OSPCIs["var1","upper"], PITCIs["var1","upper"], PTHCIs["var1","upper"], RSYCIs["var1","upper"], RTHCIs["var1","upper"], SERCIs["var1","upper"], SLSCIs["var1","upper"], SPDCIs["var1","upper"], STYCIs["var1","upper"], TAICIs["var1","upper"], TOMCIs["var1","upper"]))

ggplot(Coeff, aes(Site, MeanCoeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=UCI, ymin=LCI, width=0.5))+
  theme_bw()

