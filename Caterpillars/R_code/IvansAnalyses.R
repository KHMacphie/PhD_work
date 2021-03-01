rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(plyr)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)
library(forcats)
library(gridExtra)

#### Working out what to expect of Ivvan's analyses ####

#### Dataframe ####
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater$treeID <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$siteyear <- paste(cater$site, cater$year)
cater$datescale <- scale(cater$date) #unscale(x, center= 146.77, scale=14.04305)   now from -2.11991 to 2.01024
cater$tree.species<- revalue(cater$tree.species, c("?"="OthDecid", "Cherry"="OthDecid", "Aspen"="OthDecid", "Chestnut"="OthDecid", "Lime"="OthDecid", "Field maple"="OthDecid", "Damson"="OthDecid", "Whitebeam"="OthDecid"))
cater<- subset(cater, tree.species!="OthDecid")
cater$resid <- seq(1,length(cater$year),1)
cater <- subset(cater, year!="2014")

mod <- glmer(beetles~datescale+I(datescale^2)+(1+datescale+I(datescale^2)|site)+(1+datescale+I(datescale^2)|tree.species)+(1+datescale+I(datescale^2)|year)+(1|siteday)+(1|siteyear)+(1|recorder)+(1|treeID)+(1|resid),data=cater,family=poisson,nAGQ=0) 

dummydayofyearscaled <- seq(-2.12, 2.01, 0.01)

mean <- summary(mod)$coef[1,1]+summary(mod)$coef[2,1]*dummydayofyearscaled+summary(mod)$coef[3,1]*dummydayofyearscaled^2

cioverallbeetles<-confint.merMod(mod,method="Wald")

lower<-cioverallbeetles[26,1]* dummydayofyearscaled ^2+cioverallbeetles[25,1]* dummydayofyearscaled +cioverallbeetles[24,1] 

upper<-cioverallbeetles[26,2]* dummydayofyearscaled ^2+cioverallbeetles[25,2]* dummydayofyearscaled +cioverallbeetles[24,2] 

plot(dummydayofyearscaled,exp(mean), ylim=c(0,0.3))
lines(dummydayofyearscaled,exp(lower),lty=4,col="blue")
lines(dummydayofyearscaled,exp(upper),lty=4,col="blue")

pred <- data.frame(datescale=seq(-2.12, 2.01, 0.01))
predict(mod, newdata=pred, interval="confidence")

sim <- simulate(mod, 1000)
hist(apply(sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
#normally not too many rogue values but reduce x axis a bit to see main distribution relative to observed value 
hist(apply(sim,2,sum), breaks=100000, xlim=c(0,100000))
abline(v=sum(cater$beetles),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(sim,2,propzero), breaks=1000) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater$beetles), col="red")

#### Ivan loop problem for all year peak figure data points
data3 <- read.csv("Downloads/data3.csv")
dateyear <- unique(paste(data3$date, data3$year)) #list of unique year date combinations
store <- pmatch(dateyear, data3$dateyear) #matching those to the year and date
dateyear <- cbind(dateyear, data3$year[store], data3$date[store]) 
table(dateyear[,3]) # shows how many years each date was sampled in 
#117, 118, 119, 120, 121, 122, 174, 175 have only 1 year which stops the model running 
#maybe remove from plot as it isnt going to be an average of multiple years? 
#just say in the legend that in G data points only for dates with >1 year's data
#or can keep them in if you prefer but make sure to mention it is just the mean of one year.
#can see in your previous plot on the googledoc that the first few days are the same as 2019.
#Here I'm removing them but you have those points if you want to add them back in.

data3 <-  data3[data3$date>122&data3$date<174, ] #keeping rows with dates between 123 and 173

#start of ally loop code for daily invert abundance across all years
store.peakallyearsmm2<-c()
day<- min(data3$date)
for (date in min(data3$date):max(data3$date)){
  dailydataallyearsmm2<-data3[which(data3$date==date),]
  if(sum(dailydataallyearsmm2$beetles)>0){  
    modallyearsmm2<-glmer(beetles~1+(1|site)+(1|resid)+(1|tree.species)+(1|year)+(1|siteyear),family=poisson,data=dailydataallyearsmm2,nAGQ=0)
    store.peakallyearsmm2[date-122]<-exp(fixef(modallyearsmm2))
  }
  if(sum(dailydataallyearsmm2$beetles)==0){  
    store.peakallyearsmm2[date-122]<-0
  }
}
#end of ally loop code for daily invert abundance 

#the list in store.peakallyearsmm2 is now for dates 123:173



#estimates, ending in l is link scale ending in d data
birchl <- -1.841
birchd <- 0.159
sel <- 0.169
lod <- 0.114
upd <- 0.221

bob <- rnorm(1000,-1.841,0.169) #random numbers drawn from distribution of link scale parameter estimate

#plotting link scale hist with mean and +/- sd
abline(v=-1.841, col=2)
abline(v=-1.841+0.169, col=3)
abline(v=-1.841-0.169, col=3)

#plotting data scale hist with exp(mean and +/- sd) 
abline(v=exp(-1.841), col=2)
abline(v=exp(-1.841+0.169), col=3)
abline(v=exp(-1.841-0.169), col=3)

# so exp range is are no longer even to either side of the mean

####Repeat in mcmc glmm
beetspid <- subset(cater_habitat, year!="2014")
beetspid$resid <- seq(1,nrow(beetspid),1) 
#modspiderpeaknewfull<-glmer(spiders~datescale+I(datescale^2)+(1+datescale+I(datescale^2)|site)+(1+datescale+I(datescale^2)|tree.species)+(1+datescale+I(datescale^2)|year)+(1|datesiteyear)+(1|siteyear)+(1|recorder)+(1|treeID)+(1|resid),data=data3,family=poisson,nAGQ=0)
modbeetlepeaknewfull<-glmer(beetles~datescaled+I(datescaled^2)+(1+datescaled+I(datescaled^2)|site)+(1+datescaled+I(datescaled^2)|tree.species)+(1+datescaled+I(datescaled^2)|year)+(1|siteday)+(1|siteyear)+(1|recorder)+(1|treeID)+(1|resid),data=beetspid,family=poisson,nAGQ=0) #14:42-15:00



k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

beetles<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
                   random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + siteyear + treeID + siteday + recorder, 
                   family="poisson", data=beetspid, prior=prior3, nitt=200000, burnin=100000, thin=100)


pred <- data.frame(dscal= seq(min(beetspid$datescaled), max(beetspid$datescaled), 0.001))
                   
for(i in 1:length(pred$dscal)){
       A <- beetles$Sol[,1]+beetles$Sol[,2]*pred$dscal[i]+beetles$Sol[,3]*pred$dscal[i]^2
       pred$mean[i] <- exp(mean(A))
       pred$lowci[i] <- exp(HPDinterval(A)[1])
       pred$upci[i] <- exp(HPDinterval(A)[2])
  }
                   
plot(pred$dscal, pred$mean, ylim=c(0,0.4), type="l", xlab="date scaled", ylab="abundance")
points(pred$dscal, pred$lowci, type="l", col=2)
points(pred$dscal, pred$upci, type="l", col=2)