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

# LINES OF CODE THAT CAN BE CHANGED INDICATED BY **********

#rpois(n,lambda)
#first beating will be removed because no treatment (0 in spreadsheet and not included here)

#response <- #abundance of caterpillars
#treatment <- #5 level factor: duration between beatings: 2,4,6,8,10  use 4 days as reference level  
#labelled: d2, d4, d6, d8, d10
#treesp <- # oak and sycamore
#date <- #ordinal date as factor, help control for phenology, weather etc
#branch <- #tree code + branch no.
#tree <- # 3 or 4 branches per tree
#treedate <- paste(tree, date) #help control for tree/invert phenology
#branchdate <- paste(branch,date) #help control for branch size/general differences and also tree/invert phenology

# based around caterpillars, all need to be on log scale

#1: as it is now but just oak
#2: all samples on each tree, keep treesp in model
#3: all samples on each tree, remove treesp from model
#4: all samples on each tree, just oak

########################################
#### Parts that can be manipulated: ####
########################################
# multiples of 5 trees per species to spread treatments over dates and trees
reps.oak <- 8 # (so reps=2 = 10 trees)                                            **********
#reps.syc <- 3 #                                                                   **********
propdif <- 0.5 #proportional difference per 2 day difference in treatment         **********
repeats <- 50  # how many simulations to run                                     **********
########################################


#number of samples per treatment per 5 trees- 5 trees allows even spread of treatments among dates and trees
n.d2 <- 50
n.d4 <- 24
n.d6 <- 24
n.d8 <- 24
n.d10 <- 25

#setting up dataframe for each treatment across 5 trees the actual dates dont matter just example
df.d2 <- data.frame(date=as.factor(rep(seq(114,162,2),2)), 
                    tree=c(rep("A",(n.d2/2)),rep("C",(n.d2/2))), 
                    branch=c(rep("a",(n.d2/2)),rep("d",(n.d2/2))), 
                    treatment=rep("d2",n.d2)) 

df.d4 <- data.frame(date=as.factor(c(seq(116,160,4),seq(114,158,4))), 
                    tree=c(rep("B",(n.d4/2)),rep("D",(n.d4/2))), 
                    branch=rep("a",(n.d4)), 
                    treatment=rep("d4",n.d4))

df.d6 <- data.frame(date=as.factor(c(seq(116,158,6),seq(118,160,6),seq(120,162,6))), 
                    tree=c(rep("A",(n.d6/3)),rep("C",(n.d6/3)),rep("E",(n.d6/3))), 
                    branch=c(rep("b",(n.d6/3)),rep("a",(n.d6/3)),rep("a",(n.d6/3))), 
                    treatment=rep("d6",n.d6))

df.d8 <- data.frame(date=as.factor(c(seq(118,158,8),seq(124,164,8),seq(122,162,8),seq(120,160,8))), 
                    tree=c(rep("B",(n.d8/4)),rep("C",(n.d8/4)),rep("D",(n.d8/4)),rep("E",(n.d8/4))), 
                    branch=c(rep("b",(n.d8/4)),rep("b",(n.d8/4)),rep("b",(n.d8/4)),rep("b",(n.d8/4))), 
                    treatment=rep("d8",n.d8))

df.d10 <- data.frame(date=as.factor(c(seq(126,166,10),seq(122,162,10),seq(124,164,10),seq(118,158,10),seq(120,160,10))), 
                     tree=c(rep("A",(n.d10/5)),rep("B",(n.d10/5)),rep("C",(n.d10/5)),rep("D",(n.d10/5)),rep("E",(n.d10/5))), 
                     branch=c(rep("c",(n.d10/5)),rep("c",(n.d10/5)),rep("c",(n.d10/5)),rep("c",(n.d10/5)),rep("c",(n.d10/5))), 
                     treatment=rep("d10",n.d10))


# dataframe for samples collected for every 5 trees
df.5trees <- rbind(df.d2,df.d4,df.d6,df.d8,df.d10)

#Sort oak dataframe
df.oak <- df.5trees
df.oak$treesp <- rep("oak", length(df.oak$date)) #column for treesp
df.oak$tree <- paste(df.oak$treesp, df.oak$tree) #tree ID by species and tree letter

df.oak <- do.call("rbind", replicate(reps.oak, df.oak, simplify = FALSE)) #repeat the oak dataframe as many times as its being replicated..
df.oak$rep <- rep(1:reps.oak, each=(length(df.oak$date)/reps.oak)) #new column with repeat number for each time the df was repeated 

df.oak$tree <- paste(df.oak$tree, df.oak$rep) #Tree ID individual since df replication
df.oak$branch <- paste(df.oak$tree, df.oak$branch) #branch ID by species, tree and branch letters
df.oak$rep <- NULL

#Sort syc dataframe
#df.syc <- df.5trees
#df.syc$treesp <- rep("syc", length(df.syc$date)) #column for treesp
#df.syc$tree <- paste(df.syc$treesp, df.syc$tree) #tree ID by species and tree letter

#df.syc <- do.call("rbind", replicate(reps.syc, df.syc, simplify = FALSE)) #repeat the syc dataframe as many times as its being replicated..
#df.syc$rep <- rep(1:reps.syc, each=(length(df.syc$date)/reps.syc))#new column with repeat number for each time the df was repeated 

#df.syc$tree <- paste(df.syc$tree, df.syc$rep) #Tree ID individual since df replication
#df.syc$branch <- paste(df.syc$tree, df.syc$branch) #branch ID by species, tree and branch letters
#df.syc$rep <- NULL

# one dataframe for all samples
df <- df.oak #rbind(df.oak, df.syc)
df$treedate <- paste(df$tree, df$date) #treedate for random term
df$branchdate <- paste(df$branch, df$date) #branch date for random term (residual for simulation not in model)

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
#df.d2.syc <- subset(df, treesptrmt=="syc d2")
df.d4.oak <- subset(df, treesptrmt=="oak d4")
#df.d4.syc <- subset(df, treesptrmt=="syc d4")
df.d6.oak <- subset(df, treesptrmt=="oak d6")
#df.d6.syc <- subset(df, treesptrmt=="syc d6")
df.d8.oak <- subset(df, treesptrmt=="oak d8")
#df.d8.syc <- subset(df, treesptrmt=="syc d8")
df.d10.oak <- subset(df, treesptrmt=="oak d10")
#df.d10.syc <- subset(df, treesptrmt=="syc d10")

# sample size for each treatment on each treesp
n.d2.oak <- length(df.d2.oak$branchdate)
#n.d2.syc <- length(df.d2.syc$branchdate)
n.d4.oak <- length(df.d4.oak$branchdate)
#n.d4.syc <- length(df.d4.syc$branchdate)
n.d6.oak <- length(df.d6.oak$branchdate)
#n.d6.syc <- length(df.d6.syc$branchdate)
n.d8.oak <- length(df.d8.oak$branchdate)
#n.d8.syc <- length(df.d8.syc$branchdate)
n.d10.oak <- length(df.d10.oak$branchdate)
#n.d10.syc <- length(df.d10.syc$branchdate)

df <- rbind(df.d2.oak,df.d4.oak,df.d6.oak,df.d8.oak,df.d10.oak) #rbind(df.d2.oak,df.d2.syc,df.d4.oak,df.d4.syc,df.d6.oak,df.d6.syc,df.d8.oak,df.d8.syc,df.d10.oak,df.d10.syc)
df$treatment<- relevel(df$treatment, ref="d4") #d4 as reference level, oak automatic because alphabetical

#Difference between each treatment: going for 2 days = 50% change in abundance
intercept <- -4.09 #intercept from TTHA
d4.oak.mean <- intercept+0.45 #oak random effect in TTHA
#d4.syc.mean <- intercept+0.12 #syc random effect in TTHA
d2.oak.mean <- d4.oak.mean+log(1-propdif)   
#d2.syc.mean <- d4.syc.mean+log(1-propdif)
d6.oak.mean <- d4.oak.mean+log(1+propdif)
#d6.syc.mean <- d4.syc.mean+log(1+propdif)
d8.oak.mean <- d4.oak.mean+log(1+2*propdif)
#d8.syc.mean <- d4.syc.mean+log(1+2*propdif)
d10.oak.mean <- d4.oak.mean+log(1+3*propdif)
#d10.syc.mean <- d4.syc.mean+log(1+3*propdif)


#simulating the response, running model and storing relevant outputs  
df$poisresponse <- c(
    rpois(n.d2.oak,  exp(d2.oak.mean + df.d2.oak$eff.date + df.d2.oak$eff.branch + df.d2.oak$eff.tree + df.d2.oak$eff.treedate + df.d2.oak$eff.branchdate)), 
    #rpois(n.d2.syc,  exp(d2.syc.mean + df.d2.syc$eff.date + df.d2.syc$eff.branch + df.d2.syc$eff.tree + df.d2.syc$eff.treedate + df.d2.syc$eff.branchdate)), 
    rpois(n.d4.oak,  exp(d4.oak.mean + df.d4.oak$eff.date + df.d4.oak$eff.branch + df.d4.oak$eff.tree + df.d4.oak$eff.treedate + df.d4.oak$eff.branchdate)),
    #rpois(n.d4.syc,  exp(d4.syc.mean + df.d4.syc$eff.date + df.d4.syc$eff.branch + df.d4.syc$eff.tree + df.d4.syc$eff.treedate + df.d4.syc$eff.branchdate)),
    rpois(n.d6.oak,  exp(d6.oak.mean + df.d6.oak$eff.date + df.d6.oak$eff.branch + df.d6.oak$eff.tree + df.d6.oak$eff.treedate + df.d6.oak$eff.branchdate)), 
    #rpois(n.d6.syc,  exp(d6.syc.mean + df.d6.syc$eff.date + df.d6.syc$eff.branch + df.d6.syc$eff.tree + df.d6.syc$eff.treedate + df.d6.syc$eff.branchdate)),
    rpois(n.d8.oak,  exp(d8.oak.mean + df.d8.oak$eff.date + df.d8.oak$eff.branch + df.d8.oak$eff.tree + df.d8.oak$eff.treedate + df.d8.oak$eff.branchdate)),
    #rpois(n.d8.syc,  exp(d8.syc.mean + df.d8.syc$eff.date + df.d8.syc$eff.branch + df.d8.syc$eff.tree + df.d8.syc$eff.treedate + df.d8.syc$eff.branchdate)),
    rpois(n.d10.oak, exp(d10.oak.mean + df.d10.oak$eff.date + df.d10.oak$eff.branch + df.d10.oak$eff.tree + df.d10.oak$eff.treedate + df.d10.oak$eff.branchdate)))#,
    #rpois(n.d10.syc, exp(d10.syc.mean + df.d10.syc$eff.date + df.d10.syc$eff.branch + df.d10.syc$eff.tree + df.d10.syc$eff.treedate + df.d10.syc$eff.branchdate)))
  
  k<-10000
  prior<-list(R=list(V=1,nu=0.02),
              G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                     G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
  
  poismodel <- MCMCglmm(poisresponse~treatment, random=~date+branch+tree+treedate, family = "poisson", data = df, nitt=500000, burnin=50000, prior=prior)
  
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

par(mfrow=c(3,2),ask=F, oma=c(0,0,2,0))
hist(modelcoeffs$d4.mean, breaks=100)
abline(v=mean(modelcoeffs$d4.mean), col=2)
abline(v=true.d4.oak.mean, col=3)
legend("topleft", legend=(length(which(modelcoeffs$d4.pval<0.05))/repeats)*100, bty = "n")
#mtext(paste("Oak trees=",reps.oak*5,"   Sycamore trees=",reps.syc*5,"   Prop difference (2days)=",propdif,"   Nsim=",repeats),side = 3, line = -1, outer = TRUE)
mtext(paste("Oak trees=",reps.oak*5,"   Prop difference (2days)=",propdif,"   Nsim=",repeats),side = 3, line = -1, outer = TRUE)


hist(modelcoeffs$d2.mean, breaks=100)
abline(v=mean(modelcoeffs$d2.mean), col=2)
abline(v=true.d2.dif, col=3)
legend("topleft", legend=(length(which(modelcoeffs$d2.pval<0.05))/repeats)*100, bty = "n")

hist(modelcoeffs$d6.mean, breaks=100)
abline(v=mean(modelcoeffs$d6.mean), col=2)
abline(v=true.d6.dif, col=3)

legend("topleft", legend=(length(which(modelcoeffs$d6.pval<0.05))/repeats)*100, bty = "n")

hist(modelcoeffs$d8.mean, breaks=100)
abline(v=mean(modelcoeffs$d8.mean), col=2)
abline(v=true.d8.dif, col=3)

legend("topleft", legend=(length(which(modelcoeffs$d8.pval<0.05))/repeats)*100, bty = "n")

hist(modelcoeffs$d10.mean, breaks=100)
abline(v=mean(modelcoeffs$d10.mean), col=2)
abline(v=true.d10.dif, col=3)

legend("topleft", legend=(length(which(modelcoeffs$d10.pval<0.05))/repeats)*100, bty = "n")

#hist(modelcoeffs$syc.mean, breaks=100)
#abline(v=mean(modelcoeffs$syc.mean), col=2)
#abline(v=true.syc.dif, col=3)
#legend("topleft", legend=(length(which(modelcoeffs$syc.pval<0.05))/repeats)*100, bty = "n")

#save as 6"x8"