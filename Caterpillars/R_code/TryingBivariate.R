##########################
#### Trying bivariate ####
##########################
rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(plyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)
library(forcats)
library(gridExtra)

##########################
#### Sorting caterpillar Mass ####
##########################

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)

# set all =<0.02 to 0.02 ready for interval censoring
cater$caterpillar.mass <- revalue(cater$caterpillar.mass, c("0"="")) #where there is no caterpillar mass= NA. MCMCglmmm treats this as missing at random. 
cater$caterpillar.mass2 <- revalue(cater$caterpillar.mass, c("<0.01"="0.02", "0.01"="0.02"))
cater$caterpillar.mass2 <- as.numeric(as.character(cater$caterpillar.mass2))

#mass per caterpillar and interval censoring
cater$mpc2 <- cater$caterpillar.mass2/cater$caterpillars
cater$mpc1 <- as.character(cater$caterpillar.mass2)
cater$mpc1 <- revalue(cater$mpc1, c("0.02"="0.001"))
cater$mpc1 <- as.numeric(cater$mpc1)
cater$mpc1 <- cater$mpc1/cater$caterpillars
cater$mpc1 <- ifelse(cater$mpc1 < 0.001,0.001,cater$mpc1)

# Sort rows with mass.uncertain
which(cater$mass.uncertain=="1") # 11694 21835 22385 25555 25743 25876 26665 27967 28137
cater$notes[11694]
cater$mpc1[11694] <- NA
cater$mpc2[11694] <- NA
cater[21835,]
cater$mpc2[21835] <- cater$caterpillar.mass2[21835]
cater[22385,]
cater$mpc1[22385] <- cater$caterpillar.mass2[22385]/4
cater$mpc2[22385] <- cater$caterpillar.mass2[22385]/4
cater[25555,]
cater$mpc2[25555] <- cater$caterpillar.mass2[25555]
cater[25743,]
cater$mpc2[25743] <- cater$caterpillar.mass2[25743]
cater[25876,]
cater$mpc1[25876] <- cater$caterpillar.mass2[25876]/3
cater$mpc2[25876] <- cater$caterpillar.mass2[25876]/3
cater[26665,]
cater$mpc1[26665] <- cater$caterpillar.mass2[26665]/7
cater$mpc2[26665] <- cater$caterpillar.mass2[26665]/7
cater[27967,]
cater$mpc1[27967] <- cater$caterpillar.mass2[27967]
cater$mpc2[27967] <- cater$caterpillar.mass2[27967]
cater[28137,]
cater$mpc1[28137] <- cater$caterpillar.mass2[28137]/3
cater$mpc2[28137] <- cater$caterpillar.mass2[28137]/3


#log
cater$logmpc1 <- log(cater$mpc1)
cater$logmpc2 <- log(cater$mpc2)

# extra variables
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$siteyear <- paste(cater$site, cater$year)
cater$datecent <- cater$date - mean(cater$date)
cater$datescaled <- cater$date/max(cater$date)

#### bivariate model ####
#trait-1 =- fits a different intercept for the two columns. Same for any contiobnuous variable
#trait latitude fits a seperate latitide for each
#us(trait):site estimates variance across sites in caterpillar abundance and mass. Also the covariance between the two aceross sites.
#for fixed effects you can use at.level(trait,1):covariate to fit a predictor (e.g., date^2) to just one of the reponse variables.
#e.g. at.level(trait,1):date:elevation
k<-1000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

bivariate<-MCMCglmm(cbind(caterpillars, logmpc1,logmpc2)~trait-1 + trait:datescaled + at.level(trait,1):datescaled^2,
                    random=~us(trait):tree.species + us(trait):site + us(trait):year + us(trait):treeID + us(trait):recorder + us(trait):siteday,
                    rcov=~us(trait):units, family=c("poisson","cengaussian"), data=cater, prior=prior, pr=TRUE, thin=100, nitt=300000, burnin=30000)
save(bivariate, file = "~/Documents/Models/bivariate.RData")
