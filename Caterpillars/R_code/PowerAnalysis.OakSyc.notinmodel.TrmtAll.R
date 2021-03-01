rm(list=ls())
setwd('/Users/s1205615/')
#library(ggplot2)
#library(plyr)
#library(dplyr)
#library(ggfortify)
#library(tidyr)
#library(readr)
#library(doBy)
#library(lme4)
library(MCMCglmm)
#library(forcats)
#library(gridExtra)

################################################
#### Testing Branch Beating- Power Analysis ####
################################################



#rpois(n,lambda)
#first beating will be removed because no treatment, so sampling begins before what is included here

#response <- #abundance of caterpillars
#treatment <- #5 level factor: duration between beatings: 2,4,6,8,10  use 4 days as reference level  
#labelled: d2, d4, d6, d8, d10
#treesp <- # oak and sycamore
#date <- #ordinal date as factor, help control for phenology, weather etc
#branch <- #treeID + branch no.
#tree <-  #treeID
#treedate <- paste(tree, date) #help control for tree/invert phenology
#branchdate <- paste(branch,date) #help control for branch size/general differences and also tree/invert phenology

# based around caterpillars, all need to be on log scale

#1: treatments spread between trees- would need too many trees to be possible
#2: all samples on each tree, keep treesp in model 
#3: all samples on each tree, remove treesp from model ********** This script
#4: all samples on each tree, slope not categorical

########################################
#### Parts that can be manipulated: ####
########################################
mintrees <- 15 #of each sp (so 10 = 20 trees total)
maxtrees <- 30  #of each sp
intervals <- 5 #increments of tree number to be tested
treecount <- rep(seq(mintrees,maxtrees,intervals),2) # number of trees in each simulation                        
difference <- c(rep(0.333,(length(treecount)/2)),rep(0.5,length(treecount)/2)) # proportional differences tested
repeats <- 50  # how many simulations to run                                     
########################################

for(x in 1:length(treecount)){

n.oak <- treecount[x]  # how many oak trees               
n.syc <- treecount[x]  # how many sycamore trees          

#5 treatments on 1 tree, staggering start date so every treatment should be on every day
startd2 <- c(rep(120,n.oak)[1:n.oak],rep(120,n.syc)[1:n.syc])
startd4 <- c(rep(seq(120,122,2),ceiling(n.oak/2))[1:n.oak],rep(seq(120,122,2),ceiling(n.syc/2))[1:n.syc])
startd6 <- c(rep(seq(120,124,2),ceiling(n.oak/3))[1:n.oak],rep(seq(120,124,2),ceiling(n.syc/3))[1:n.syc])
startd8 <- c(rep(seq(120,126,2),ceiling(n.oak/4))[1:n.oak],rep(seq(120,126,2),ceiling(n.syc/4))[1:n.syc])
startd10 <- c(rep(seq(120,128,2),ceiling(n.oak/5))[1:n.oak],rep(seq(120,128,2),ceiling(n.syc/5))[1:n.syc])
n.trees <- n.oak+n.syc
trsp <- c(rep("oak",n.oak),rep("syc",n.syc))

df <- data.frame(date=NA,treatment=NA, tree=NA)
df <- na.omit(df)

for(k in 1:n.trees){
d2 <- data.frame(date=as.factor(seq(startd2[k],170,2)),
                 treatment=rep("d2",length(seq(startd2[k],170,2))))
d4 <- data.frame(date=as.factor(seq(startd4[k],170,4)),
                 treatment=rep("d4",length(seq(startd4[k],170,4))))
d6 <- data.frame(date=as.factor(seq(startd6[k],170,6)),
                 treatment=rep("d6",length(seq(startd6[k],170,6))))
d8 <- data.frame(date=as.factor(seq(startd8[k],170,8)),
                 treatment=rep("d8",length(seq(startd8[k],170,8))))
d10 <- data.frame(date=as.factor(seq(startd10[k],170,10)),
                  treatment=rep("d10",length(seq(startd10[k],170,10)))) 

onetree <- rbind(d2,d4,d6,d8,d10)
onetree$treesp <- rep(trsp[k],length(onetree$date))
onetree$tree <- paste(onetree$treesp, rep(k,length(onetree$date)))

df <- rbind(df,onetree)
}

df$branch <- paste(df$tree, df$treatment)
df$treedate <- paste(df$tree, df$date)
df$branchdate <- paste(df$branch, df$date)

# number of levels to each random term
n.date <- length(unique(df$date))
n.tree <- length(unique(df$tree))
n.branch <- length(unique(df$branch))
n.treedate <- length(unique(df$treedate))
n.branchdate <- length(unique(df$branchdate))

# sd's from variances for each random term from TTHA model because has no date and date^2 in fixed
sd.tree <- sqrt(0.257/2) #treeID var split between tree and branch
sd.date <- sqrt(1.065) #from siteday (closest thing from that model.. site and site year are in model too)
sd.branch <- sqrt(0.257/2) #treeID var split between tree and branch
sd.treedate <- sqrt(0.666/2) #residual var split between treedate and branchdate
sd.branchdate <- sqrt(0.666/2) #residual var split between treedate and branchdate

#data frame for model coefficients and p-values

modelcoeffs <- data.frame(rep=seq(1,repeats,1)) #adjust here for number of repeats

#loop for simulating the response, running model and storing relevant outputs
for(i in 1:length(modelcoeffs$rep)){   
  
# assigning effect magnitude to each random effect
dates <- data.frame(date=unique(df$date))
dates$eff.date <- rnorm(n.date, mean=0, sd=sd.date)
trees <- data.frame(tree=unique(df$tree))
trees$eff.tree <- rnorm(n.tree, mean=0, sd=sd.tree)
branches <- data.frame(branch=unique(df$branch))
branches$eff.branch <- rnorm(n.branch, mean=0, sd=sd.branch)
treedates <- data.frame(treedate=unique(df$treedate))
treedates$eff.treedate <- rnorm(n.treedate, mean=0, sd=sd.treedate)
branchdates <- data.frame(branchdate=unique(df$branchdate))
branchdates$eff.branchdate <- rnorm(n.branchdate, mean=0, sd=sd.branchdate)

store1 <- pmatch(df$date, dates$date, duplicates.ok = TRUE)
store2 <- pmatch(df$tree, trees$tree, duplicates.ok = TRUE)
store3 <- pmatch(df$branch, branches$branch, duplicates.ok = TRUE)
store4 <- pmatch(df$treedate, treedates$treedate, duplicates.ok = TRUE)
store5 <- pmatch(df$branchdate, branchdates$branchdate,duplicates.ok = TRUE)

df$eff.date <- dates$eff.date[store1]
df$eff.tree <- trees$eff.tree[store2]
df$eff.branch <- branches$eff.branch[store3]
df$eff.treedate <- treedates$eff.treedate[store4]
df$eff.branchdate <- branchdates$eff.branchdate[store5]

df$treesptrmt <- paste(df$treesp, df$treatment)
df.d2.oak <- subset(df, treesptrmt=="oak d2")
df.d2.syc <- subset(df, treesptrmt=="syc d2")
df.d4.oak <- subset(df, treesptrmt=="oak d4")
df.d4.syc <- subset(df, treesptrmt=="syc d4")
df.d6.oak <- subset(df, treesptrmt=="oak d6")
df.d6.syc <- subset(df, treesptrmt=="syc d6")
df.d8.oak <- subset(df, treesptrmt=="oak d8")
df.d8.syc <- subset(df, treesptrmt=="syc d8")
df.d10.oak <- subset(df, treesptrmt=="oak d10")
df.d10.syc <- subset(df, treesptrmt=="syc d10")

# sample size for each treatment on each treesp
n.d2.oak <- length(df.d2.oak$branchdate)
n.d2.syc <- length(df.d2.syc$branchdate)
n.d4.oak <- length(df.d4.oak$branchdate)
n.d4.syc <- length(df.d4.syc$branchdate)
n.d6.oak <- length(df.d6.oak$branchdate)
n.d6.syc <- length(df.d6.syc$branchdate)
n.d8.oak <- length(df.d8.oak$branchdate)
n.d8.syc <- length(df.d8.syc$branchdate)
n.d10.oak <- length(df.d10.oak$branchdate)
n.d10.syc <- length(df.d10.syc$branchdate)

df <- rbind(df.d2.oak,df.d2.syc,df.d4.oak,df.d4.syc,df.d6.oak,df.d6.syc,df.d8.oak,df.d8.syc,df.d10.oak,df.d10.syc)
df$treatment<- relevel(df$treatment, ref="d4") #d4 as reference level, oak automatic because alphabetical

#Difference between each treatment: going for 2 days = 50% change in abundance
intercept <- -4.09 #intercept from TTHA
propdif <- difference[x]
d4.oak.mean <- intercept+0.45 #oak random effect in TTHA
d4.syc.mean <- intercept+0.12 #syc random effect in TTHA
d2.oak.mean <- d4.oak.mean+log(1-propdif)   
d2.syc.mean <- d4.syc.mean+log(1-propdif)
d6.oak.mean <- d4.oak.mean+log(1+propdif)
d6.syc.mean <- d4.syc.mean+log(1+propdif)
d8.oak.mean <- d4.oak.mean+log(1+2*propdif)
d8.syc.mean <- d4.syc.mean+log(1+2*propdif)
d10.oak.mean <- d4.oak.mean+log(1+3*propdif)
d10.syc.mean <- d4.syc.mean+log(1+3*propdif)


#simulating the response, running model and storing relevant outputs  
df$poisresponse <- c(
    rpois(n.d2.oak,  exp(d2.oak.mean + df.d2.oak$eff.date + df.d2.oak$eff.branch + df.d2.oak$eff.tree + df.d2.oak$eff.treedate + df.d2.oak$eff.branchdate)), 
    rpois(n.d2.syc,  exp(d2.syc.mean + df.d2.syc$eff.date + df.d2.syc$eff.branch + df.d2.syc$eff.tree + df.d2.syc$eff.treedate + df.d2.syc$eff.branchdate)), 
    rpois(n.d4.oak,  exp(d4.oak.mean + df.d4.oak$eff.date + df.d4.oak$eff.branch + df.d4.oak$eff.tree + df.d4.oak$eff.treedate + df.d4.oak$eff.branchdate)),
    rpois(n.d4.syc,  exp(d4.syc.mean + df.d4.syc$eff.date + df.d4.syc$eff.branch + df.d4.syc$eff.tree + df.d4.syc$eff.treedate + df.d4.syc$eff.branchdate)),
    rpois(n.d6.oak,  exp(d6.oak.mean + df.d6.oak$eff.date + df.d6.oak$eff.branch + df.d6.oak$eff.tree + df.d6.oak$eff.treedate + df.d6.oak$eff.branchdate)), 
    rpois(n.d6.syc,  exp(d6.syc.mean + df.d6.syc$eff.date + df.d6.syc$eff.branch + df.d6.syc$eff.tree + df.d6.syc$eff.treedate + df.d6.syc$eff.branchdate)),
    rpois(n.d8.oak,  exp(d8.oak.mean + df.d8.oak$eff.date + df.d8.oak$eff.branch + df.d8.oak$eff.tree + df.d8.oak$eff.treedate + df.d8.oak$eff.branchdate)),
    rpois(n.d8.syc,  exp(d8.syc.mean + df.d8.syc$eff.date + df.d8.syc$eff.branch + df.d8.syc$eff.tree + df.d8.syc$eff.treedate + df.d8.syc$eff.branchdate)),
    rpois(n.d10.oak, exp(d10.oak.mean + df.d10.oak$eff.date + df.d10.oak$eff.branch + df.d10.oak$eff.tree + df.d10.oak$eff.treedate + df.d10.oak$eff.branchdate)),
    rpois(n.d10.syc, exp(d10.syc.mean + df.d10.syc$eff.date + df.d10.syc$eff.branch + df.d10.syc$eff.tree + df.d10.syc$eff.treedate + df.d10.syc$eff.branchdate)))
  
  k<-10000
  prior<-list(R=list(V=1,nu=0.02),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
  
  poismodel <- MCMCglmm(poisresponse~treatment, random=~date+branch+tree+treedate, family = "poisson", data = df, nitt=500000, burnin=50000, thin=20, prior=prior)
  
  modelcoeffs$d4.mean[i] <- summary(poismodel)$solutions[1,1]
  modelcoeffs$d4.pval[i] <- summary(poismodel)$solutions[1,5]
  modelcoeffs$d2.mean[i] <- summary(poismodel)$solutions[2,1]
  modelcoeffs$d2.pval[i] <- summary(poismodel)$solutions[2,5]
  modelcoeffs$d6.mean[i] <- summary(poismodel)$solutions[3,1]
  modelcoeffs$d6.pval[i] <- summary(poismodel)$solutions[3,5]
  modelcoeffs$d8.mean[i] <- summary(poismodel)$solutions[4,1]
  modelcoeffs$d8.pval[i] <- summary(poismodel)$solutions[4,5]
  modelcoeffs$d10.mean[i] <- summary(poismodel)$solutions[5,1]
  modelcoeffs$d10.pval[i] <- summary(poismodel)$solutions[5,5]
  #modelcoeffs$syc.mean[i] <- summary(poismodel)$solutions[6,1]
  #modelcoeffs$syc.pval[i] <- summary(poismodel)$solutions[6,5]
}

#what the coefficients should be
true.d4.oak.mean <- intercept+0.45 #oak random effect in TTHA
#true.syc.dif <- d4.syc.mean-d4.oak.mean #syc random effect in TTHA
true.d2.dif <- log(1-propdif)   
true.d6.dif <- log(1+propdif)
true.d8.dif <- log(1+2*propdif)
true.d10.dif <- log(1+3*propdif)

# Checking mean prediction and p-values
# histogram of coefficients, red line= mean from models, green line= true value
# text= % of models with significant result
names <- c(paste(treecount,"Oak",treecount,"SycNIM_",difference,"dif"))
mypath <- file.path(paste("~/Documents/Models/PowerAnalysis/OakSyc_notinmodel_TrmtAll_",names[x],".pdf"))

pdf(file = mypath, width=8, height=6)

par(mfrow=c(3,2),ask=F, oma=c(0,0,2,0))
hist(modelcoeffs$d4.mean, breaks=100)
abline(v=mean(modelcoeffs$d4.mean), col=2)
abline(v=true.d4.oak.mean, col=3)
legend("topleft", legend=(length(which(modelcoeffs$d4.pval<0.05))/repeats)*100, bty = "n")
mtext(paste("Trees=",(n.oak+n.syc),"   Prop difference (2days)=",propdif,"   Nsim=",repeats),side = 3, line = -1, outer = TRUE)


hist(modelcoeffs$d2.mean, breaks=100)
abline(v=mean(modelcoeffs$d2.mean), col=2)
abline(v=true.d2.dif, col=3)
legend("topleft", legend=paste((length(which(modelcoeffs$d2.pval<0.05))/repeats)*100,"%"), bty = "n")

hist(modelcoeffs$d6.mean, breaks=100)
abline(v=mean(modelcoeffs$d6.mean), col=2)
abline(v=true.d6.dif, col=3)

legend("topleft", legend=paste((length(which(modelcoeffs$d6.pval<0.05))/repeats)*100,"%"), bty = "n")

hist(modelcoeffs$d8.mean, breaks=100)
abline(v=mean(modelcoeffs$d8.mean), col=2)
abline(v=true.d8.dif, col=3)

legend("topleft", legend=paste((length(which(modelcoeffs$d8.pval<0.05))/repeats)*100,"%"), bty = "n")

hist(modelcoeffs$d10.mean, breaks=100)
abline(v=mean(modelcoeffs$d10.mean), col=2)
abline(v=true.d10.dif, col=3)

legend("topleft", legend=paste((length(which(modelcoeffs$d10.pval<0.05))/repeats)*100,"%"), bty = "n")

#hist(modelcoeffs$syc.mean, breaks=100)
#abline(v=mean(modelcoeffs$syc.mean), col=2)
#abline(v=true.syc.dif, col=3)
#legend("topleft", legend=paste((length(which(modelcoeffs$syc.pval<0.05))/repeats)*100,"%"), bty = "n")

dev.off()

}
#save as 6"x8"