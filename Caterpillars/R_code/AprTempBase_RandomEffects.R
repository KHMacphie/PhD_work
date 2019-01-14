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
#### Calculating the coefficient for each site with CI #### 

# AprBaseModel_RE_3 Site coefficients 
which(colnames(AprBaseModel_RE_3$Sol)=="site.ALN") #= 14
which(colnames(AprBaseModel_RE_3$Sol)=="site.TOM") #= 57

siteREcropped <- AprBaseModel_RE_3$Sol[,14:57] # crop to just the columns wanted
site.df <- data.frame(site=c(colnames(siteREcropped))) #column for yearsite 
site.df$coeff <- apply(siteREcropped,2, mean) # mean 
for(i in 1:length(site.df$site)) {   # loop for CIs
  A <- HPDinterval(siteREcropped[,i])
  site.df$lowci[i] <- A["var1","lower"] 
  site.df$upci[i] <- A["var1","upper"] 
} 
site.df$site <- gsub("site.","", site.df$site)

## sites in alphabetical order
ggplot(site.df, aes(site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))

pmatch(site.df$site, site$site)
site.df<- merge(site.df, site, by="site", duplicates.ok=TRUE)
site.df<- rename(site.df, latitude="Mean Lat")
site.df<- rename(site.df, elevation="Mean Elev")

## sites ordered by latitude
site.df.lat <- site.df
site.df.lat$Site <- site.df.lat$site
site.df.lat$Site <- order(site.df$latitude) 
ggplot(site.df.lat, aes(Site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))


## sites ordered by elevation
site.df.el <- site.df
site.df.el$Site <- site.df.el$site
site.df.el$Site <- order(site.df$elevation) 
ggplot(site.df.el, aes(Site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))


#####################################
#### YearSite in AprBaseModel_RE ####


# Trying to look at yearsite random effect in AprBaseModel_RE
summary(AprBaseModel_RE)
head(AprBaseModel_RE$Sol)
which( colnames(AprBaseModel_RE$Sol)=="yearsite.ALN 2014")
(14*13)+7
4982+188
which(colnames(AprBaseModel_RE$Sol)=="yearsite.TOM 2018")
## yearsite in columns 4982:5170


## Data frame with mean and Cis for yearsite RE
cropped <- AprBaseModel_RE$Sol[,4982:5170] # crop to just the columns wanted
yearsite.df <- data.frame(ys=c(colnames(cropped))) #column for yearsite 
yearsite.df$coeff <- apply(cropped,2, mean) # mean 
for(i in 1:length(yearsite.df$ys)) {   # loop for CIs
A <- HPDinterval(cropped[,i])
yearsite.df$lowci[i] <- A["var1","lower"]
yearsite.df$upci[i] <- A["var1","upper"]
} 

## putting in year and site seperately
yearsite.df$yearsite <- (unique(all_data$yearsite))
yearsite.df$site <- yearsite.df$yearsite
yearsite.df$year <- yearsite.df$yearsite

numbers <- c(0,1,2,3,4,5,6,7,8,9)
for(i in c(numbers)){
  yearsite.df$site <- gsub(i, "",yearsite.df$site)
}

letters <- c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
for(i in c(letters)){
  yearsite.df$year <- gsub(i, "",yearsite.df$year)
}


#### Plotting yearsite REs by site    !!!!!!! should this also have year (and site) coef added?
ggplot(yearsite.df, aes(year, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  facet_wrap(.~site)+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))
  
#### Plotting yearsite REs by year    !!!!!! should this also have site (and year) added?
ggplot(yearsite.df, aes(site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  facet_grid(year~.)+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))
  