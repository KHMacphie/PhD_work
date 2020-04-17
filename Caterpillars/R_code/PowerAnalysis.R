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

################################################
#### Testing Branch Beating- Power Analysis ####
################################################

rpois(n,lambda)
#first beating will be removed because no treatment (0 in spreadsheet and not included here)

response <- #abundance of caterpillars
treatment <- #5 level factor: duration between beatings: 2,4,6,8,10  use 4 days as reference level
  
treesp <- # oak and sycamore
date <- #ordinal date as factor, help control for phenology, weather etc
branch <- #tree code + branch no.
tree <- #factor (use letters maybe for clarity)
treedate <- paste(tree, date) #help control for tree/invert phenology
branchdate <- paste(branch,date) #help control for branch size/general differences and also tree/invert phenology

# based around caterpillars, all need to be on log scale
# variances from TTHA model becasue no date and date^2 in fixed
n.tree <- 10
sd.tree <- sqrt(0.257/2) #treeID var split between tree and branch

n.date <- 27
sd.date <- sqrt(1.065) #from siteday (closest thing from that model.. site and site year are in model too)

n.branch <- 3*n.tree
sd.branch <- sqrt(0.257/2) #treeID var split between tree and branch

n.treedate <- 89*(n.tree/5)
sd.treedate <- sqrt(0.666/2) #residual var split between treedate and branchdate
  
n.branchdate <- 122*(n.tree/5) #sum of number of samples for each treatment*(ntree/5)
sd.branchdate <- sqrt(0.666/2) #residual var split between treedate and branchdate

intercept <- -4.09 #intercept from TTHA
n.oak <- ceiling(n.tree/2) #will round up if decimal
n.syc <- as.integer(n.tree/2) #will round down if decimal

#going for 2 days =50% change in abundance
propdif <- 0.5 #proportional difference per 2 days 
n.d4.oak <- 24*(n.oak/5)
d4.oak.mean <- intercept+0.45 #oak random effect in TTHA

n.d4.syc <- 24*(n.syc/5)
d4.syc.mean <- intercept+0.12 #syc random effect in TTHA

n.d2.oak <- 25*(n.oak/5)
d2.oak.mean <- d4.oak.mean+log(1-propdif)   

n.d2.syc <- 25*(n.syc/5)
d2.syc.mean <- d4.syc.mean+log(1-propdif)

n.d6.oak <- 24*(n.oak/5)
d6.oak.mean <- d4.oak.mean+log(1+propdif)

n.d6.syc <- 24*(n.syc/5)
d6.syc.mean <- d4.syc.mean+log(1+propdif)

n.d8.oak <- 24*(n.oak/5)
d8.oak.mean <- d4.oak.mean+log(1+2*propdif)

n.d8.syc <- 24*(n.syc/5)
d8.syc.mean <- d4.syc.mean+log(1+2*propdif)

n.d10.oak <- 25*(n.oak/5)
d10.oak.mean <- d4.oak.mean+log(1+3*propdif)

n.d10.syc <- 25*(n.syc/5)
d10.syc.mean <- d4.syc.mean+log(1+3*propdif)

df.d2.oak <- data.frame(branchdate=seq(1,n.d2.oak,1), 
                        date=as.factor(seq(114,162,2)), 
                        tree=rep("A",n.d2.oak), 
                        branch=rep("a",n.d2.oak), 
                        treatment=rep("d2",n.d2.oak), 
                        treesp=rep("oak",n.d2.oak))
df.d2.oak$branch <- paste(df.d2.oak$tree, df.d2.oak$branch)
df.d2.oak$branchdate <- paste(df.d2.oak$branch, df.d2.oak$date)
df.d2.oak$treedate <- paste(df.d2.oak$tree, df.d2.oak$date)



df.d4.oak <- data.frame(branchdate=seq(1,n.d4.oak,1), 
                        date=as.factor(c(seq(116,160,4),seq(114,158,4))), 
                        tree=c(rep("B",(n.d4.oak/2)),rep("D",(n.d4.oak/2))), 
                        branch=rep("a",(n.d4.oak)), 
                        treatment=rep("d4",n.d4.oak), 
                        treesp=rep("oak",n.d4.oak))
df.d4.oak$branch <- paste(df.d4.oak$tree, df.d4.oak$branch)
df.d4.oak$branchdate <- paste(df.d4.oak$branch, df.d4.oak$date)
df.d4.oak$treedate <- paste(df.d4.oak$tree, df.d4.oak$date)

df.d6.oak <- data.frame(branchdate=seq(1,n.d6.oak,1), 
                        date=as.factor(c(seq(116,158,6),seq(118,160,6),seq(120,162,6))), 
                        tree=c(rep("A",(n.d6.oak/3)),rep("C",(n.d6.oak/3)),rep("E",(n.d6.oak/3))), 
                        branch=c(rep("b",(n.d6.oak/3)),rep("a",(n.d6.oak/3)),rep("a",(n.d6.oak/3))), 
                        treatment=rep("d6",n.d6.oak), 
                        treesp=rep("oak",n.d6.oak))
df.d6.oak$branch <- paste(df.d6.oak$tree, df.d6.oak$branch)
df.d6.oak$branchdate <- paste(df.d6.oak$branch, df.d6.oak$date)
df.d6.oak$treedate <- paste(df.d6.oak$tree, df.d6.oak$date)

df.d8.oak <- data.frame(branchdate=seq(1,n.d8.oak,1), 
                        date=as.factor(c(seq(118,158,8),seq(124,164,8),seq(122,162,8),seq(120,160,8))), 
                        tree=c(rep("B",(n.d8.oak/4)),rep("C",(n.d8.oak/4)),rep("D",(n.d8.oak/4)),rep("E",(n.d8.oak/4))), 
                        branch=c(rep("b",(n.d8.oak/4)),rep("b",(n.d8.oak/4)),rep("b",(n.d8.oak/4)),rep("b",(n.d8.oak/4))), 
                        treatment=rep("d8",n.d8.oak), 
                        treesp=rep("oak",n.d8.oak))
df.d8.oak$branch <- paste(df.d8.oak$tree, df.d8.oak$branch)
df.d8.oak$branchdate <- paste(df.d8.oak$branch, df.d8.oak$date)
df.d8.oak$treedate <- paste(df.d8.oak$tree, df.d8.oak$date)

df.d10.oak <- data.frame(branchdate=seq(1,n.d10.oak,1), 
                        date=as.factor(c(seq(126,166,10),seq(122,162,10),seq(124,164,10),seq(118,158,10),seq(120,160,10))), 
                        tree=c(rep("A",(n.d10.oak/5)),rep("B",(n.d10.oak/5)),rep("C",(n.d10.oak/5)),rep("D",(n.d10.oak/5)),rep("E",(n.d10.oak/5))), 
                        branch=c(rep("c",(n.d10.oak/5)),rep("c",(n.d10.oak/5)),rep("c",(n.d10.oak/5)),rep("c",(n.d10.oak/5)),rep("c",(n.d10.oak/5))), 
                        treatment=rep("d10",n.d10.oak), 
                        treesp=rep("oak",n.d10.oak))
df.d10.oak$branch <- paste(df.d10.oak$tree, df.d10.oak$branch)
df.d10.oak$branchdate <- paste(df.d10.oak$branch, df.d10.oak$date)
df.d10.oak$treedate <- paste(df.d10.oak$tree, df.d10.oak$date)

df.d2.syc <- data.frame(branchdate=seq(1,n.d2.syc,1), 
                        date=as.factor(seq(114,162,2)), 
                        tree=rep("F",n.d2.syc), 
                        branch=rep("a",n.d2.syc), 
                        treatment=rep("d2",n.d2.syc), 
                        treesp=rep("syc",n.d2.syc))
df.d2.syc$branch <- paste(df.d2.syc$tree, df.d2.syc$branch)
df.d2.syc$branchdate <- paste(df.d2.syc$branch, df.d2.syc$date)
df.d2.syc$treedate <- paste(df.d2.syc$tree, df.d2.syc$date)

df.d4.syc <- data.frame(branchdate=seq(1,n.d4.syc,1), 
                        date=as.factor(c(seq(116,160,4),seq(114,158,4))), 
                        tree=c(rep("G",(n.d4.syc/2)),rep("I",(n.d4.syc/2))), 
                        branch=rep("a",(n.d4.syc)), 
                        treatment=rep("d4",n.d4.syc), 
                        treesp=rep("syc",n.d4.syc))
df.d4.syc$branch <- paste(df.d4.syc$tree, df.d4.syc$branch)
df.d4.syc$branchdate <- paste(df.d4.syc$branch, df.d4.syc$date)
df.d4.syc$treedate <- paste(df.d4.syc$tree, df.d4.syc$date)

df.d6.syc <- data.frame(branchdate=seq(1,n.d6.syc,1), 
                        date=as.factor(c(seq(116,158,6),seq(118,160,6),seq(120,162,6))), 
                        tree=c(rep("F",(n.d6.syc/3)),rep("H",(n.d6.syc/3)),rep("J",(n.d6.syc/3))), 
                        branch=c(rep("b",(n.d6.syc/3)),rep("a",(n.d6.syc/3)),rep("a",(n.d6.syc/3))), 
                        treatment=rep("d6",n.d6.syc), 
                        treesp=rep("syc",n.d6.syc))
df.d6.syc$branch <- paste(df.d6.syc$tree, df.d6.syc$branch)
df.d6.syc$branchdate <- paste(df.d6.syc$branch, df.d6.syc$date)
df.d6.syc$treedate <- paste(df.d6.syc$tree, df.d6.syc$date)

df.d8.syc <- data.frame(branchdate=seq(1,n.d8.syc,1), 
                        date=as.factor(c(seq(118,158,8),seq(124,164,8),seq(122,162,8),seq(120,160,8))), 
                        tree=c(rep("G",(n.d8.syc/4)),rep("H",(n.d8.syc/4)),rep("I",(n.d8.syc/4)),rep("J",(n.d8.syc/4))), 
                        branch=c(rep("b",(n.d8.syc/4)),rep("b",(n.d8.syc/4)),rep("b",(n.d8.syc/4)),rep("b",(n.d8.syc/4))), 
                        treatment=rep("d8",n.d8.syc), 
                        treesp=rep("syc",n.d8.syc))
df.d8.syc$branch <- paste(df.d8.syc$tree, df.d8.syc$branch)
df.d8.syc$branchdate <- paste(df.d8.syc$branch, df.d8.syc$date)
df.d8.syc$treedate <- paste(df.d8.syc$tree, df.d8.syc$date)

df.d10.syc <- data.frame(branchdate=seq(1,n.d10.syc,1), 
                         date=as.factor(c(seq(126,166,10),seq(122,162,10),seq(124,164,10),seq(118,158,10),seq(120,160,10))), 
                         tree=c(rep("F",(n.d10.syc/5)),rep("G",(n.d10.syc/5)),rep("H",(n.d10.syc/5)),rep("I",(n.d10.syc/5)),rep("J",(n.d10.syc/5))), 
                         branch=c(rep("c",(n.d10.syc/5)),rep("c",(n.d10.syc/5)),rep("c",(n.d10.syc/5)),rep("c",(n.d10.syc/5)),rep("c",(n.d10.syc/5))), 
                         treatment=rep("d10",n.d10.syc), 
                         treesp=rep("syc",n.d10.syc))
df.d10.syc$branch <- paste(df.d10.syc$tree, df.d10.syc$branch)
df.d10.syc$branchdate <- paste(df.d10.syc$branch, df.d10.syc$date)
df.d10.syc$treedate <- paste(df.d10.syc$tree, df.d10.syc$date)

df <- rbind(df.d2.oak,df.d2.syc,df.d4.oak,df.d4.syc,df.d6.oak,df.d6.syc,df.d8.oak,df.d8.syc,df.d10.oak,df.d10.syc)


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

store1 <- pmatch(df$date,dates$date, duplicates.ok = TRUE)
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

df <- rbind(df.d2.oak,df.d2.syc,df.d4.oak,df.d4.syc,df.d6.oak,df.d6.syc,df.d8.oak,df.d8.syc,df.d10.oak,df.d10.syc)
df$treesp<- relevel(df$treesp, ref="oak")
df$treatment<- relevel(df$treatment, ref="d4")

modelcoeffs <- data.frame(rep=seq(1,5,1))
#missing tree taxa

for(i in 1:length(modelcoeffs$rep)){   
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
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

poismodel <- MCMCglmm(poisresponse~treatment+treesp, random=~date+branch+tree+treedate, family = "poisson", data = df, nitt=100000, burnin=35000, prior=prior)

modelcoeffs$d4.mean[i] <- summary(poismodel)$solutions[1,1]
modelcoeffs$d2.mean[i] <- summary(poismodel)$solutions[2,1]
modelcoeffs$d2.pval[i] <- summary(poismodel)$solutions[2,5]
modelcoeffs$d6.mean[i] <- summary(poismodel)$solutions[3,1]
modelcoeffs$d6.pval[i] <- summary(poismodel)$solutions[3,5]
modelcoeffs$d8.mean[i] <- summary(poismodel)$solutions[4,1]
modelcoeffs$d8.pval[i] <- summary(poismodel)$solutions[4,5]
modelcoeffs$d10.mean[i] <- summary(poismodel)$solutions[5,1]
modelcoeffs$d10.pval[i] <- summary(poismodel)$solutions[5,5]
modelcoeffs$syc.mean[i] <- summary(poismodel)$solutions[6,1]
modelcoeffs$syc.pval[i] <- summary(poismodel)$solutions[6,5]
}


# Checking mean prediction and p-values
mean(modelcoeffs$d4.mean)
mean(modelcoeffs$d2.mean)
length(which(modelcoeffs$d2.pval<0.05))/10000
mean(modelcoeffs$d6.mean)
length(which(modelcoeffs$d6.pval<0.05))/10000
mean(modelcoeffs$d8.mean)
length(which(modelcoeffs$d8.pval<0.05))/10000
mean(modelcoeffs$d10.mean)
length(which(modelcoeffs$d10.pval<0.05))/10000
mean(modelcoeffs$syc.mean)
length(which(modelcoeffs$syc.pval<0.05))/10000






 