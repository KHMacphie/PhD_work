###########################################################
#### Abundance peak date/height and variances for each site:year ####
###########################################################

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

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")

#site <- read.csv("Dropbox/master_data/site/site_details.csv")
#temp <- read.csv("Dropbox/master_data/site/temperatures.csv")

cater$year <- as.factor(cater$year)
cater$siteyear <- paste(cater$site, cater$year)
cater$sitetree <- paste(cater$site, cater$tree)

#SiteYear <- glmer(caterpillars~date+I(date^2)+(date|year)+(date|siteyear)+(date|site), family=poisson , data=cater)
#save(SiteYear, file = "~/Documents/Models/SiteYear.RData")
load("~/Documents/Models/SiteYear.RData")

summary(SiteYear)
coef(SiteYear)$site
coef(SiteYear)$siteyear
coef(SiteYear)$year






## bayesian?
k<-1000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

### is idh() different variances and us() uniform variances?
SiteCurves<- MCMCglmm(caterpillars~date*year+I(date^2), 
                       random=~us(1+date):siteyear+us(1+date):site+sitetree+siteday, 
                       family="poisson", data=cater_habitat_condensed, prior=prior, nitt=250000, burnin=25000, pr=TRUE)
save(MultiMembRE, file = "~/Documents/Models/SiteCurves.RData")
load("~/Documents/Models/SiteCurves.RData")