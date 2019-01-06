#################################################
## Random effect output in base temp Apr model ##
#################################################


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

