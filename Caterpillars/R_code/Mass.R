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
#### Caterpillar Mass ####
##########################

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)
summary(cater)
unique(cater$notes)
# set all =<0.02 to 0.02 ready for interval censoring
cater$caterpillar.mass <- revalue(cater$caterpillar.mass, c("0"=""))
cater$caterpillar.mass2 <- revalue(cater$caterpillar.mass, c("<0.01"="0.02", "0.01"="0.02"))
cater$caterpillar.mass2 <- as.numeric(as.character(cater$caterpillar.mass2))

#mass per caterpillar
cater$mpc2 <- cater$caterpillar.mass2/cater$caterpillars
cater$mpc1 <- as.character(cater$caterpillar.mass2)
cater$mpc1 <- revalue(cater$mpc1, c("0.02"="0.001"))
cater$mpc1 <- as.numeric(cater$mpc1)
cater$mpc1 <- cater$mpc1/cater$caterpillars
cater$mpc1 <- ifelse(cater$mpc1 < 0.001,0.001,cater$mpc1)

#cater$mpc1 <- round(cater$mpc1, digits=3)
#cater$mpc1 <- as.character(cater$mpc1)
#cater$mpc1 <- revalue(cater$mpc1, c("0"="0.002","0.001"="0.002"))
#cater$mpc1 <- as.numeric(cater$mpc1)

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


#Sort weather
cater$weather <- gsub(" ","", cater$weather)
cater$weather <- gsub("damp","wet", cater$weather)
cater$weather <- gsub("DRY","dry", cater$weather)
cater$weather <- gsub("dry/wind","dry", cater$weather)
cater$weather <- gsub("dry/v.windy","dry", cater$weather)
cater$weather <- gsub("WET","wet", cater$weather)

#interval censor and log
hist(cater$mpc1, breaks=100)
hist(cater$mpc2, breaks=100)
cater$logmpc1 <- log(cater$mpc1)
cater$logmpc2 <- log(cater$mpc2)
hist(cater$logmpc1, breaks=100)
hist(cater$logmpc2, breaks=100)

plot(cater$date, cater$logmpc1)
points(cater$date, cater$logmpc2, col=2)

# Checking normality using random value between 0.001 and 0.02
#cater$mass_rand<-cater$mpc1
#cater$mass_rand[which(cater$mass_rand==0.001)]<-runif(length(cater$mass_rand[which(cater$mass_rand==0.001)]),0.001,0.02)
#hist(log(cater$mass_rand))

# extra variables
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$siteyear <- paste(cater$site, cater$year)
cater$datecent <- cater$date - mean(cater$date)
# Model priors
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


Mass<- MCMCglmm(cbind(logmpc1, logmpc2)~ datecent, 
                       random=~us(1+datecent):year + us(1+datecent):site + treeID + recorder + siteday, 
                       family="cengaussian", data=cater, prior=prior, nitt=300000, burnin=30000, pr=TRUE)
save(Mass, file = "~/Documents/Models/Mass.RData")
load("~/Documents/Models/Mass.RData")

prior2<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


Mass2<- MCMCglmm(cbind(logmpc1, logmpc2)~ datecent, 
                random=~us(1+datecent):tree.species + us(1+datecent):year + us(1+datecent):site + treeID + recorder + siteday, 
                family="cengaussian", data=cater, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)
save(Mass2, file = "~/Documents/Models/Mass2.RData")
load("~/Documents/Models/Mass2.RData")

###########################################################
#### Calculating the coefficient for each tree sp. with CI #### 

#  tree species intercept coefficients : columns 3:20 
#  tree species slope coefficients : columns 21:38
TTintercept.cropped <- Mass2$Sol[,3:20] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- gsub(".tree.species.","", TTintercept$treetaxa)
TTintercept$treetaxa <- gsub("(Intercept)","", TTintercept$treetaxa)
TTintercept$treetaxa <- gsub("()","", TTintercept$treetaxa) # not working

par(mfcol=c(1,1), cex=1.5)
hist(Mass2$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,0.2))
abline(v=0, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
#legend("topright", legend="B", bty="n") 

par(mfcol=c(1,1), cex=0.5)
ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

TTslope.cropped <- Mass2$Sol[,21:38] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datecent.tree.species.","", TTslope$treetaxa)

par(mfcol=c(1,1), cex=1.5)
hist(Mass2$VCV[,"datecent:datecent.tree.species"], breaks=500)
abline(v=0, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
#legend("topright", legend="B", bty="n") 

par(mfcol=c(1,1), cex=0.5)
ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))


#### Try to expand grid so one row per caterpillar
cater.expanded <- cater[rep(row.names(cater), cater$caterpillars), 1:28]

#k<-10000
#prior<-list(R=list(V=1,nu=0.002),
#            G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
#                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#Mass3<- MCMCglmm(cbind(logmpc1, logmpc2)~ datecent, 
#                random=~us(1+datecent):year + us(1+datecent):site + treeID + recorder + siteday, 
#                family="cengaussian", data=cater.expanded, prior=prior, nitt=300000, burnin=30000, pr=TRUE)
#save(Mass3, file = "~/Documents/Models/Mass3.RData")
#load("~/Documents/Models/Mass3.RData")
k<-10000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#### Mass4 uses lower interval of 0.002 not 0.001
#Mass4<- MCMCglmm(cbind(logmpc1, logmpc2)~ datecent, 
#                 random=~us(1+datecent):tree.species + us(1+datecent):year + us(1+datecent):site + treeID + recorder + siteday, 
#                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)
#save(Mass4, file = "~/Documents/Models/Mass4.RData")
load("~/Documents/Models/Mass4.RData")
summary(Mass4)
par(mfcol=c(2,1), cex=1.5)
hist(Mass4$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,0.2), main="TreeTaxa intercept")
abline(v=0.001513915, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
hist(Mass4$VCV[,"datecent:datecent.tree.species"], breaks=500, xlim=c(0,0.002), main="TreeTaxa slope")
abline(v=5.986041e-05, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
#legend("topright", legend="B", bty="n") 

TTslope.cropped <- Mass4$Sol[,19:34] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datecent.tree.species.","", TTslope$treetaxa)

ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- Mass4$Sol[,3:18] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))



### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddaycent <- seq(-29.4,28.6,0.1)
predday <- preddaycent+mean(cater$date)
meanslope <- mean(Mass4$Sol[,1])+mean(Mass4$Sol[,2])*preddaycent
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(Mass4$Sol[,1])+TaxaCatGrowth[1,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[1,3])*preddaycent
ash <- mean(Mass4$Sol[,1])+TaxaCatGrowth[2,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[2,3])*preddaycent
aspen <- mean(Mass4$Sol[,1])+TaxaCatGrowth[3,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[3,3])*preddaycent
beech <- mean(Mass4$Sol[,1])+TaxaCatGrowth[4,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[4,3])*preddaycent
birch <- mean(Mass4$Sol[,1])+TaxaCatGrowth[5,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[5,3])*preddaycent
cherry <- mean(Mass4$Sol[,1])+TaxaCatGrowth[6,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[6,3])*preddaycent
chestnut <- mean(Mass4$Sol[,1])+TaxaCatGrowth[7,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[7,3])*preddaycent
damson <- mean(Mass4$Sol[,1])+TaxaCatGrowth[8,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[8,3])*preddaycent
elm <- mean(Mass4$Sol[,1])+TaxaCatGrowth[9,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[9,3])*preddaycent
fieldmaple <- mean(Mass4$Sol[,1])+TaxaCatGrowth[10,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[10,3])*preddaycent
hazel <- mean(Mass4$Sol[,1])+TaxaCatGrowth[11,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[11,3])*preddaycent
lime <- mean(Mass4$Sol[,1])+TaxaCatGrowth[12,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[12,3])*preddaycent
oak <- mean(Mass4$Sol[,1])+TaxaCatGrowth[13,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[13,3])*preddaycent
rowan <- mean(Mass4$Sol[,1])+TaxaCatGrowth[14,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[14,3])*preddaycent
sycamore <- mean(Mass4$Sol[,1])+TaxaCatGrowth[15,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[15,3])*preddaycent
willow <- mean(Mass4$Sol[,1])+TaxaCatGrowth[16,2]+mean(Mass4$Sol[,2]+TaxaCatGrowth[16,3])*preddaycent

plot(cater$date, cater$logmpc1)

points(cater$date, cater$logmpc2, col=2)
points(predday,alder, type="l", lwd=2, col=4)
points(predday,ash, type="l", lwd=2, col=4)
points(predday,aspen, type="l", lwd=2, col=4)
points(predday,beech, type="l", lwd=2, col=4)
points(predday,birch, type="l", lwd=2, col=4)
points(predday,cherry, type="l", lwd=2, col=4)
points(predday,chestnut, type="l", lwd=2, col=4)
points(predday,damson, type="l", lwd=2, col=4)
points(predday,elm, type="l", lwd=2, col=4)
points(predday,fieldmaple, type="l", lwd=2, col=4)
points(predday,hazel, type="l", lwd=2, col=4)
points(predday,lime, type="l", lwd=2, col=4)
points(predday,oak, type="l", lwd=2, col=4)
points(predday,rowan, type="l", lwd=2, col=4)
points(predday,sycamore, type="l", lwd=2, col=4)
points(predday,willow, type="l", lwd=2, col=4)
points(predday,meanslope, type="l", lwd=2, col=3)

plot(cater$date, exp(cater$logmpc1))
points(cater$date, exp(cater$logmpc2), col=2)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col=3)

#### same again with lower bound at 0.001 and date scaled not centred, also with weather
cater.expanded$datescaled <- cater.expanded$date/max(cater.expanded$date)

k<-10000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


Mass5<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + weather, 
                 random=~us(1+datescaled):tree.species + us(1+datescaled):year + us(1+datescaled):site + treeID + recorder + siteday, 
                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=250000, burnin=25000, pr=TRUE)
save(Mass5, file = "~/Documents/Models/Mass5.RData")
load("~/Documents/Models/Mass5.RData")

summary(Mass5)
par(mfcol=c(2,1), cex=1.5)
hist(Mass5$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,20), main="TreeTaxa intercept",xlab="Variance")
abline(v=0.5545415, col=2,type="l", lty=2)
hist(Mass5$VCV[,"datescaled:datescaled.tree.species"], breaks=500, xlim=c(0,30), main="TreeTaxa slope",xlab="Variance")
abline(v=1.018952, col=2,type="l", lty=2)


TTslope.cropped <- Mass5$Sol[,22:37] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- Mass5$Sol[,6:21] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))



### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass5$Sol[,1])+mean(Mass5$Sol[,2])*preddayscaled
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(Mass5$Sol[,1])+TaxaCatGrowth[1,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[1,3])*preddayscaled
ash <- mean(Mass5$Sol[,1])+TaxaCatGrowth[2,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[2,3])*preddayscaled
aspen <- mean(Mass5$Sol[,1])+TaxaCatGrowth[3,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[3,3])*preddayscaled
beech <- mean(Mass5$Sol[,1])+TaxaCatGrowth[4,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[4,3])*preddayscaled
birch <- mean(Mass5$Sol[,1])+TaxaCatGrowth[5,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[5,3])*preddayscaled
cherry <- mean(Mass5$Sol[,1])+TaxaCatGrowth[6,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[6,3])*preddayscaled
chestnut <- mean(Mass5$Sol[,1])+TaxaCatGrowth[7,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[7,3])*preddayscaled
damson <- mean(Mass5$Sol[,1])+TaxaCatGrowth[8,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[8,3])*preddayscaled
elm <- mean(Mass5$Sol[,1])+TaxaCatGrowth[9,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[9,3])*preddayscaled
fieldmaple <- mean(Mass5$Sol[,1])+TaxaCatGrowth[10,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[10,3])*preddayscaled
hazel <- mean(Mass5$Sol[,1])+TaxaCatGrowth[11,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[11,3])*preddayscaled
lime <- mean(Mass5$Sol[,1])+TaxaCatGrowth[12,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[12,3])*preddayscaled
oak <- mean(Mass5$Sol[,1])+TaxaCatGrowth[13,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[13,3])*preddayscaled
rowan <- mean(Mass5$Sol[,1])+TaxaCatGrowth[14,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[14,3])*preddayscaled
sycamore <- mean(Mass5$Sol[,1])+TaxaCatGrowth[15,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[15,3])*preddayscaled
willow <- mean(Mass5$Sol[,1])+TaxaCatGrowth[16,2]+mean(Mass5$Sol[,2]+TaxaCatGrowth[16,3])*preddayscaled

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, cater.expanded$logmpc1, xlab="Date", ylab="log(Mass)", pch=20, col="grey")
points(cater.expanded$date, cater.expanded$logmpc2, pch=20, col=1)
points(predday,alder, type="l", lwd=2, col=4)
points(predday,ash, type="l", lwd=2, col=4)
points(predday,aspen, type="l", lwd=2, col=4)
points(predday,beech, type="l", lwd=2, col=4)
points(predday,birch, type="l", lwd=2, col=4)
points(predday,cherry, type="l", lwd=2, col=4)
points(predday,chestnut, type="l", lwd=2, col=4)
points(predday,damson, type="l", lwd=2, col=4)
points(predday,elm, type="l", lwd=2, col=4)
points(predday,fieldmaple, type="l", lwd=2, col=4)
points(predday,hazel, type="l", lwd=2, col=4)
points(predday,lime, type="l", lwd=2, col=4)
points(predday,oak, type="l", lwd=2, col=4)
points(predday,rowan, type="l", lwd=2, col=4)
points(predday,sycamore, type="l", lwd=2, col=4)
points(predday,willow, type="l", lwd=2, col=4)
points(predday,meanslope, type="l", lwd=3, col="red")

plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

###########################################
#### Mass 6: Year in fixed, no weather ####

k<-10000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


Mass6<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled*year, 
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + treeID + recorder + siteday, 
                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=250000, burnin=25000, pr=TRUE)
save(Mass6, file = "~/Documents/Models/Mass6.RData")
load("~/Documents/Models/Mass6.RData")

par(mfcol=c(2,2), cex=1.5)
hist(Mass6$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,20), main="TreeTaxa intercept",xlab="Variance")
abline(v=0.5235048, col=2, lty=2)
hist(Mass6$VCV[,"datescaled:datescaled.tree.species"], breaks=500, xlim=c(0,30), main="TreeTaxa slope",xlab="Variance")
abline(v=1.031537, col=2, lty=2)
hist(Mass6$VCV[,"(Intercept):(Intercept).site"], breaks=5000, xlim=c(0,0.4), main="Site intercept",xlab="Variance")
abline(v=4.344364e-10, col=2, lty=2)
hist(Mass6$VCV[,"datescaled:datescaled.site"], breaks=5000, xlim=c(0,0.6), main="Site slope",xlab="Variance")
abline(v=2.834175e-10, col=2, lty=2)


TTslope.cropped <- Mass6$Sol[,23:38] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- Mass6$Sol[,7:22] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass6$Sol[,1])+mean(Mass6$Sol[,2])*preddayscaled
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(Mass6$Sol[,1])+TaxaCatGrowth[1,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[1,3])*preddayscaled
ash <- mean(Mass6$Sol[,1])+TaxaCatGrowth[2,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[2,3])*preddayscaled
aspen <- mean(Mass6$Sol[,1])+TaxaCatGrowth[3,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[3,3])*preddayscaled
beech <- mean(Mass6$Sol[,1])+TaxaCatGrowth[4,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[4,3])*preddayscaled
birch <- mean(Mass6$Sol[,1])+TaxaCatGrowth[5,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[5,3])*preddayscaled
cherry <- mean(Mass6$Sol[,1])+TaxaCatGrowth[6,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[6,3])*preddayscaled
chestnut <- mean(Mass6$Sol[,1])+TaxaCatGrowth[7,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[7,3])*preddayscaled
damson <- mean(Mass6$Sol[,1])+TaxaCatGrowth[8,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[8,3])*preddayscaled
elm <- mean(Mass6$Sol[,1])+TaxaCatGrowth[9,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[9,3])*preddayscaled
fieldmaple <- mean(Mass6$Sol[,1])+TaxaCatGrowth[10,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[10,3])*preddayscaled
hazel <- mean(Mass6$Sol[,1])+TaxaCatGrowth[11,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[11,3])*preddayscaled
lime <- mean(Mass6$Sol[,1])+TaxaCatGrowth[12,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[12,3])*preddayscaled
oak <- mean(Mass6$Sol[,1])+TaxaCatGrowth[13,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[13,3])*preddayscaled
rowan <- mean(Mass6$Sol[,1])+TaxaCatGrowth[14,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[14,3])*preddayscaled
sycamore <- mean(Mass6$Sol[,1])+TaxaCatGrowth[15,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[15,3])*preddayscaled
willow <- mean(Mass6$Sol[,1])+TaxaCatGrowth[16,2]+mean(Mass6$Sol[,2]+TaxaCatGrowth[16,3])*preddayscaled

### barchart of lower interval by date
intsamples <- subset(cater.expanded, mpc1 == 0.001, 
                     select=c(date, mpc1))
hist(intsamples$date, breaks=100)

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, cater.expanded$logmpc1, xlab="Date", ylab="log(Mass)", pch=20, col="grey")
points(cater.expanded$date, cater.expanded$logmpc2, pch=20, col=1)
points(predday,alder, type="l", lwd=2, col=4)
points(predday,ash, type="l", lwd=2, col=4)
points(predday,aspen, type="l", lwd=2, col=4)
points(predday,beech, type="l", lwd=2, col=4)
points(predday,birch, type="l", lwd=2, col=4)
points(predday,cherry, type="l", lwd=2, col=4)
points(predday,chestnut, type="l", lwd=2, col=4)
points(predday,damson, type="l", lwd=2, col=4)
points(predday,elm, type="l", lwd=2, col=4)
points(predday,fieldmaple, type="l", lwd=2, col=4)
points(predday,hazel, type="l", lwd=2, col=4)
points(predday,lime, type="l", lwd=2, col=4)
points(predday,oak, type="l", lwd=2, col=4)
points(predday,rowan, type="l", lwd=2, col=4)
points(predday,sycamore, type="l", lwd=2, col=4)
points(predday,willow, type="l", lwd=2, col=4)
points(predday,meanslope, type="l", lwd=3, col="red")

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, exp(cater.expanded$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

Siteslope.cropped <- Mass6$Sol[,83:126] # crop to just the columns wanted
Siteslope <- data.frame(site=c(colnames(Siteslope.cropped))) #dataframe with column for yearsite 
Siteslope$coeff <- apply(Siteslope.cropped,2, mean) # mean 
for(i in 1:length(Siteslope$site)) {   # loop for CIs
  A <- HPDinterval(Siteslope.cropped[,i])
  Siteslope$lowci[i] <- A["var1","lower"] 
  Siteslope$upci[i] <- A["var1","upper"] 
} 
Siteslope$site <- gsub("datescaled.site.","", Siteslope$site)

Siteslopecoeff <- ggplot(Siteslope, aes(site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

Siteintercept.cropped <- Mass6$Sol[,39:82] # crop to just the columns wanted
Siteintercept <- data.frame(site=c(colnames(Siteintercept.cropped))) #dataframe with column for yearsite 
Siteintercept$coeff <- apply(Siteintercept.cropped,2, mean) # mean 
for(i in 1:length(Siteintercept$site)) {   # loop for CIs
  A <- HPDinterval(Siteintercept.cropped[,i])
  Siteintercept$lowci[i] <- A["var1","lower"] 
  Siteintercept$upci[i] <- A["var1","upper"] 
} 
Siteintercept$site <- Siteslope$site

Siteinterceptcoeff <- ggplot(Siteintercept, aes(site, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

library(gridExtra)

row3 <- grid.arrange(Siteinterceptcoeff, ncol = 1, widths = 1)
row4 <- grid.arrange(Siteslopecoeff, ncol = 1, widths = 1)
SiteCoeff <- grid.arrange(row3, row4, nrow = 2, heights = c(1,1))
    
###########################################
#### Mass 6: Year in fixed, no weather, date^2, no 14-16 in data frame, extra random treeday ####

cater.expanded <- subset(cater.expanded, year == "2017"|year=="2018"|year=="2019")
cater.expanded$treeday <- paste(cater.expanded$treeID, cater.expanded$date, cater.expanded$year)
k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


Mass7<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + treeID + siteday, 
                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=250000, burnin=25000, pr=TRUE)
save(Mass7, file = "~/Documents/Models/Mass7.RData")
load("~/Documents/Models/Mass7.RData")

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass7$Sol[,1])+mean(Mass7$Sol[,2])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, exp(cater.expanded$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

par(mfcol=c(2,2), cex=1.5)
hist(Mass7$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,20), main="TreeTaxa intercept",xlab="Variance")
abline(v=0.5235048, col=2, lty=2)
hist(Mass7$VCV[,"datescaled:datescaled.tree.species"], breaks=500, xlim=c(0,30), main="TreeTaxa slope",xlab="Variance")
abline(v=1.031537, col=2, lty=2)

TTslope.cropped <- Mass7$Sol[,20:35] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- Mass7$Sol[,4:19] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass7$Sol[,1])+mean(Mass7$Sol[,2])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(Mass7$Sol[,1])+TaxaCatGrowth[1,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[1,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
ash <- mean(Mass7$Sol[,1])+TaxaCatGrowth[2,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[2,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
aspen <- mean(Mass7$Sol[,1])+TaxaCatGrowth[3,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[3,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
beech <- mean(Mass7$Sol[,1])+TaxaCatGrowth[4,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[4,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
birch <- mean(Mass7$Sol[,1])+TaxaCatGrowth[5,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[5,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
cherry <- mean(Mass7$Sol[,1])+TaxaCatGrowth[6,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[6,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
chestnut <- mean(Mass7$Sol[,1])+TaxaCatGrowth[7,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[7,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
damson <- mean(Mass7$Sol[,1])+TaxaCatGrowth[8,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[8,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
elm <- mean(Mass7$Sol[,1])+TaxaCatGrowth[9,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[9,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
fieldmaple <- mean(Mass7$Sol[,1])+TaxaCatGrowth[10,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[10,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
hazel <- mean(Mass7$Sol[,1])+TaxaCatGrowth[11,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[11,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
lime <- mean(Mass7$Sol[,1])+TaxaCatGrowth[12,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[12,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
oak <- mean(Mass7$Sol[,1])+TaxaCatGrowth[13,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[13,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
rowan <- mean(Mass7$Sol[,1])+TaxaCatGrowth[14,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[14,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
sycamore <- mean(Mass7$Sol[,1])+TaxaCatGrowth[15,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[15,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2
willow <- mean(Mass7$Sol[,1])+TaxaCatGrowth[16,2]+mean(Mass7$Sol[,2]+TaxaCatGrowth[16,3])*preddayscaled+mean(Mass7$Sol[,3])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater.expanded, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, exp(cater.expanded$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,172), main="")
plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

####################################
#### Mass 7, treespecies.datesquared ####
k<-10000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#Mass8<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled*year + I(datescaled^2), 
#                 random=~us(1+datescaled+I(datescaled^2)):tree.species  + us(1+datescaled):site + treeID + siteday, 
#                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=250000, burnin=25000, pr=TRUE)
#save(Mass8, file = "~/Documents/Models/Mass8.RData")
load("~/Documents/Models/Mass8.RData")

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass8$Sol[,1])+mean(Mass8$Sol[,2])*preddayscaled+mean(Mass8$Sol[,5])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, exp(cater.expanded$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- Mass8$Sol[,24:39] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- Mass8$Sol[,8:23] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

TTquad.cropped <- Mass8$Sol[,40:55] # crop to just the columns wanted
TTquad <- data.frame(treetaxa=c(colnames(TTquad.cropped))) #dataframe with column for yearsite 
TTquad$coeff <- apply(TTquad.cropped,2, mean) # mean 
for(i in 1:length(TTquad$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTquad.cropped[,i])
  TTquad$lowci[i] <- A["var1","lower"] 
  TTquad$upci[i] <- A["var1","upper"] 
} 
TTquad$treetaxa <- TTslope$treetaxa

TTquadcoeff <- ggplot(TTquad, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Quadratic")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
row3 <- grid.arrange(TTquadcoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, row3, nrow = 3, heights = c(1,1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2], Quad=TTquad[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater.expanded$date)
meanslope <- mean(Mass8$Sol[,1])+mean(Mass8$Sol[,2])*preddayscaled+mean(Mass8$Sol[,5])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(Mass8$Sol[,1])+TaxaCatGrowth[1,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[1,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[1,4])*preddayscaled^2
ash <- mean(Mass8$Sol[,1])+TaxaCatGrowth[2,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[2,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[2,4])*preddayscaled^2
aspen <- mean(Mass8$Sol[,1])+TaxaCatGrowth[3,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[3,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[3,4])*preddayscaled^2
beech <- mean(Mass8$Sol[,1])+TaxaCatGrowth[4,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[4,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[4,4])*preddayscaled^2
birch <- mean(Mass8$Sol[,1])+TaxaCatGrowth[5,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[5,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[5,4])*preddayscaled^2
cherry <- mean(Mass8$Sol[,1])+TaxaCatGrowth[6,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[6,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[6,4])*preddayscaled^2
chestnut <- mean(Mass8$Sol[,1])+TaxaCatGrowth[7,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[7,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[7,4])*preddayscaled^2
damson <- mean(Mass8$Sol[,1])+TaxaCatGrowth[8,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[8,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[8,4])*preddayscaled^2
elm <- mean(Mass8$Sol[,1])+TaxaCatGrowth[9,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[9,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[9,4])*preddayscaled^2
fieldmaple <- mean(Mass8$Sol[,1])+TaxaCatGrowth[10,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[10,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[10,4])*preddayscaled^2
hazel <- mean(Mass8$Sol[,1])+TaxaCatGrowth[11,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[11,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[11,4])*preddayscaled^2
lime <- mean(Mass8$Sol[,1])+TaxaCatGrowth[12,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[12,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[12,4])*preddayscaled^2
oak <- mean(Mass8$Sol[,1])+TaxaCatGrowth[13,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[13,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[13,4])*preddayscaled^2
rowan <- mean(Mass8$Sol[,1])+TaxaCatGrowth[14,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[14,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[14,4])*preddayscaled^2
sycamore <- mean(Mass8$Sol[,1])+TaxaCatGrowth[15,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[15,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[15,4])*preddayscaled^2
willow <- mean(Mass8$Sol[,1])+TaxaCatGrowth[16,2]+mean(Mass8$Sol[,2]+TaxaCatGrowth[16,3])*preddayscaled+mean(Mass8$Sol[,5]+TaxaCatGrowth[16,4])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater.expanded, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater.expanded$date, exp(cater.expanded$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,172), main="")
plot(cater.expanded$date, exp(cater.expanded$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater.expanded$date, exp(cater.expanded$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

# mass on 165
TTmass <- data.frame(mean=(exp(Mass8$Sol[,1]+Mass8$Sol[,2]*0.9+Mass8$Sol[,5]*0.9^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,8]+(Mass8$Sol[,2]+Mass8$Sol[,24])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,40])*0.9^2))-TTmass$mean
TTmass$ash <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,9]+(Mass8$Sol[,2]+Mass8$Sol[,25])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,41])*0.9^2))-TTmass$mean
TTmass$aspen <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,10]+(Mass8$Sol[,2]+Mass8$Sol[,26])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,42])*0.9^2))-TTmass$mean
TTmass$beech <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,11]+(Mass8$Sol[,2]+Mass8$Sol[,27])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,43])*0.9^2))-TTmass$mean
TTmass$birch <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,12]+(Mass8$Sol[,2]+Mass8$Sol[,28])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,44])*0.9^2))-TTmass$mean
TTmass$cherry <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,13]+(Mass8$Sol[,2]+Mass8$Sol[,29])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,45])*0.9^2))-TTmass$mean
TTmass$chestnut <- (exp(Mass8$Sol[,1]+Mass8$Sol[,14]+(Mass8$Sol[,2]+Mass8$Sol[,30])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,46])*0.9^2))-TTmass$mean
TTmass$damson <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,15]+(Mass8$Sol[,2]+Mass8$Sol[,31])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,47])*0.9^2))-TTmass$mean
TTmass$elm <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,16]+(Mass8$Sol[,2]+Mass8$Sol[,32])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,48])*0.9^2))-TTmass$mean
TTmass$fieldmaple<-(exp(Mass8$Sol[,1]+Mass8$Sol[,17]+(Mass8$Sol[,2]+Mass8$Sol[,33])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,49])*0.9^2))-TTmass$mean
TTmass$hazel <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,18]+(Mass8$Sol[,2]+Mass8$Sol[,34])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,50])*0.9^2))-TTmass$mean
TTmass$lime <-     (exp(Mass8$Sol[,1]+Mass8$Sol[,19]+(Mass8$Sol[,2]+Mass8$Sol[,35])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,51])*0.9^2))-TTmass$mean
TTmass$oak <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,20]+(Mass8$Sol[,2]+Mass8$Sol[,36])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,52])*0.9^2))-TTmass$mean
TTmass$rowan <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,21]+(Mass8$Sol[,2]+Mass8$Sol[,37])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,53])*0.9^2))-TTmass$mean
TTmass$sycamore <- (exp(Mass8$Sol[,1]+Mass8$Sol[,22]+(Mass8$Sol[,2]+Mass8$Sol[,38])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,54])*0.9^2))-TTmass$mean
TTmass$willow <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,23]+(Mass8$Sol[,2]+Mass8$Sol[,39])*0.9+(Mass8$Sol[,5]+Mass8$Sol[,55])*0.9^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:17] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  ylim(-0.05,0.15)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day165")+
  theme(text = element_text(size=15))

# mass on 170
TTmass <- data.frame(mean=(exp(Mass8$Sol[,1]+Mass8$Sol[,2]*1+Mass8$Sol[,5]*1^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,8]+(Mass8$Sol[,2]+Mass8$Sol[,24])*1+(Mass8$Sol[,5]+Mass8$Sol[,40])*1^2))-TTmass$mean
TTmass$ash <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,9]+(Mass8$Sol[,2]+Mass8$Sol[,25])*1+(Mass8$Sol[,5]+Mass8$Sol[,41])*1^2))-TTmass$mean
TTmass$aspen <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,10]+(Mass8$Sol[,2]+Mass8$Sol[,26])*1+(Mass8$Sol[,5]+Mass8$Sol[,42])*1^2))-TTmass$mean
TTmass$beech <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,11]+(Mass8$Sol[,2]+Mass8$Sol[,27])*1+(Mass8$Sol[,5]+Mass8$Sol[,43])*1^2))-TTmass$mean
TTmass$birch <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,12]+(Mass8$Sol[,2]+Mass8$Sol[,28])*1+(Mass8$Sol[,5]+Mass8$Sol[,44])*1^2))-TTmass$mean
TTmass$cherry <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,13]+(Mass8$Sol[,2]+Mass8$Sol[,29])*1+(Mass8$Sol[,5]+Mass8$Sol[,45])*1^2))-TTmass$mean
TTmass$chestnut <- (exp(Mass8$Sol[,1]+Mass8$Sol[,14]+(Mass8$Sol[,2]+Mass8$Sol[,30])*1+(Mass8$Sol[,5]+Mass8$Sol[,46])*1^2))-TTmass$mean
TTmass$damson <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,15]+(Mass8$Sol[,2]+Mass8$Sol[,31])*1+(Mass8$Sol[,5]+Mass8$Sol[,47])*1^2))-TTmass$mean
TTmass$elm <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,16]+(Mass8$Sol[,2]+Mass8$Sol[,32])*1+(Mass8$Sol[,5]+Mass8$Sol[,48])*1^2))-TTmass$mean
TTmass$fieldmaple<-(exp(Mass8$Sol[,1]+Mass8$Sol[,17]+(Mass8$Sol[,2]+Mass8$Sol[,33])*1+(Mass8$Sol[,5]+Mass8$Sol[,49])*1^2))-TTmass$mean
TTmass$hazel <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,18]+(Mass8$Sol[,2]+Mass8$Sol[,34])*1+(Mass8$Sol[,5]+Mass8$Sol[,50])*1^2))-TTmass$mean
TTmass$lime <-     (exp(Mass8$Sol[,1]+Mass8$Sol[,19]+(Mass8$Sol[,2]+Mass8$Sol[,35])*1+(Mass8$Sol[,5]+Mass8$Sol[,51])*1^2))-TTmass$mean
TTmass$oak <-      (exp(Mass8$Sol[,1]+Mass8$Sol[,20]+(Mass8$Sol[,2]+Mass8$Sol[,36])*1+(Mass8$Sol[,5]+Mass8$Sol[,52])*1^2))-TTmass$mean
TTmass$rowan <-    (exp(Mass8$Sol[,1]+Mass8$Sol[,21]+(Mass8$Sol[,2]+Mass8$Sol[,37])*1+(Mass8$Sol[,5]+Mass8$Sol[,53])*1^2))-TTmass$mean
TTmass$sycamore <- (exp(Mass8$Sol[,1]+Mass8$Sol[,22]+(Mass8$Sol[,2]+Mass8$Sol[,38])*1+(Mass8$Sol[,5]+Mass8$Sol[,54])*1^2))-TTmass$mean
TTmass$willow <-   (exp(Mass8$Sol[,1]+Mass8$Sol[,23]+(Mass8$Sol[,2]+Mass8$Sol[,39])*1+(Mass8$Sol[,5]+Mass8$Sol[,55])*1^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:17] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  ylim(-0.05,0.15)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day170")+
  theme(text = element_text(size=15))


###################
#### Weighting ####
###################

# use HabitatTreeTaxaCategories script
k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#MassWeighted<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
#                 random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + year + treeID + siteday + recorder + us(1/sqrt(caterpillars)):units, 
#                 family="cengaussian", data=cater_habitat, prior=prior3, nitt=250000, burnin=25000, pr=TRUE, thin=100)
#save(MassWeighted, file = "~/Documents/Models/MassWeighted.RData")
load("~/Documents/Models/MassWeighted.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MassWeighted2$Sol[,1])+mean(MassWeighted2$Sol[,2])*preddayscaled+mean(MassWeighted2$Sol[,3])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- MassWeighted2$Sol[,14:23] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- MassWeighted2$Sol[,4:13] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

TTquad.cropped <- MassWeighted2$Sol[,24:33] # crop to just the columns wanted
TTquad <- data.frame(treetaxa=c(colnames(TTquad.cropped))) #dataframe with column for yearsite 
TTquad$coeff <- apply(TTquad.cropped,2, mean) # mean 
for(i in 1:length(TTquad$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTquad.cropped[,i])
  TTquad$lowci[i] <- A["var1","lower"] 
  TTquad$upci[i] <- A["var1","upper"] 
} 
TTquad$treetaxa <- TTslope$treetaxa

TTquadcoeff <- ggplot(TTquad, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Quadratic")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
row3 <- grid.arrange(TTquadcoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, row3, nrow = 3, heights = c(1,1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2], Quad=TTquad[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MassWeighted2$Sol[,1])+mean(MassWeighted2$Sol[,2])*preddayscaled+mean(MassWeighted2$Sol[,3])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[1,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[1,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[1,4])*preddayscaled^2
ash <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[2,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[2,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[2,4])*preddayscaled^2
beech <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[3,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[3,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[3,4])*preddayscaled^2
birch <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[4,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[4,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[4,4])*preddayscaled^2
elm <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[5,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[5,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[5,4])*preddayscaled^2
hazel <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[6,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[6,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[6,4])*preddayscaled^2
oak <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[7,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[7,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[7,4])*preddayscaled^2
rowan <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[8,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[8,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[8,4])*preddayscaled^2
sycamore <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[9,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[9,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[9,4])*preddayscaled^2
willow <- mean(MassWeighted2$Sol[,1])+TaxaCatGrowth[10,2]+mean(MassWeighted2$Sol[,2]+TaxaCatGrowth[10,3])*preddayscaled+mean(MassWeighted2$Sol[,3]+TaxaCatGrowth[10,4])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

# mass on 165
TTmass <- data.frame(mean=(exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,2]*0.9+MassWeighted2$Sol[,3]*0.9^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,4]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,14])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,24])*0.9^2))-TTmass$mean
TTmass$ash <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,5]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,15])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,25])*0.9^2))-TTmass$mean
TTmass$beech <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,6]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,16])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,25])*0.9^2))-TTmass$mean
TTmass$birch <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,7]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,17])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,26])*0.9^2))-TTmass$mean
TTmass$elm <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,8]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,18])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,27])*0.9^2))-TTmass$mean
TTmass$hazel <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,9]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,19])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,28])*0.9^2))-TTmass$mean
TTmass$oak <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,10]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,20])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,29])*0.9^2))-TTmass$mean
TTmass$rowan <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,11]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,21])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,30])*0.9^2))-TTmass$mean
TTmass$sycamore <- (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,12]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,22])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,31])*0.9^2))-TTmass$mean
TTmass$willow <-   (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,13]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,23])*0.9+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,32])*0.9^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day165")+
  theme(text = element_text(size=15))

# mass on 170
TTmass <- data.frame(mean=(exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,2]*1+MassWeighted2$Sol[,3]*1^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,4]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,14])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,24])*1^2))-TTmass$mean
TTmass$ash <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,5]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,15])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,25])*1^2))-TTmass$mean
TTmass$beech <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,6]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,16])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,26])*1^2))-TTmass$mean
TTmass$birch <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,7]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,17])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,27])*1^2))-TTmass$mean
TTmass$elm <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,8]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,18])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,28])*1^2))-TTmass$mean
TTmass$hazel <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,9]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,19])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,29])*1^2))-TTmass$mean
TTmass$oak <-      (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,10]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,20])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,30])*1^2))-TTmass$mean
TTmass$rowan <-    (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,11]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,21])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,31])*1^2))-TTmass$mean
TTmass$sycamore <- (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,12]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,22])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,32])*1^2))-TTmass$mean
TTmass$willow <-   (exp(MassWeighted2$Sol[,1]+MassWeighted2$Sol[,13]+(MassWeighted2$Sol[,2]+MassWeighted2$Sol[,23])*1+(MassWeighted2$Sol[,3]+MassWeighted2$Sol[,33])*1^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day170")+
  theme(text = element_text(size=15))


#######################################################
#### Corrected mass weighting  and only 17-19 data ####
#######################################################

# Use HabitatTreeTaxaCategories
load("~/Documents/Models/MassWeighted1719.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MassWeighted1719$Sol[,1])+mean(MassWeighted1719$Sol[,4])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- MassWeighted1719$Sol[,16:25] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- MassWeighted1719$Sol[,6:15] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2,  nrow = 2, heights = c(1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MassWeighted1719$Sol[,1])+mean(MassWeighted1719$Sol[,4])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[1,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[1,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
ash <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[2,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[2,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
beech <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[3,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[3,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
birch <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[4,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[4,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
elm <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[5,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[5,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
hazel <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[6,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[6,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
oak <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[7,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[7,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
rowan <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[8,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[8,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
sycamore <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[9,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[9,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2
willow <- mean(MassWeighted1719$Sol[,1])+TaxaCatGrowth[10,2]+mean(MassWeighted1719$Sol[,4]+TaxaCatGrowth[10,3])*preddayscaled+mean(MassWeighted1719$Sol[,5])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

# mass on 168- last day on which all tree taxa have had a caterpillar
TTmass <- data.frame(mean=(exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,4]*0.96+MassWeighted1719$Sol[,5]*0.96^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,6]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,16])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$ash <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,7]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,17])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$beech <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,8]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,18])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$birch <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,9]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,19])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$elm <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,10]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,20])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$hazel <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,11]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,21])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$oak <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,12]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,22])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$rowan <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,13]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,23])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$sycamore <- (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,14]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,24])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean
TTmass$willow <-   (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,15]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,25])*0.96+(MassWeighted1719$Sol[,5])*0.96^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day168")+
  theme(text = element_text(size=15))

# day 160
TTmass <- data.frame(mean=(exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,4]*0.91+MassWeighted1719$Sol[,5]*0.91^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,6]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,16])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$ash <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,7]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,17])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$beech <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,8]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,18])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$birch <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,9]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,19])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$elm <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,10]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,20])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$hazel <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,11]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,21])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$oak <-      (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,12]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,22])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$rowan <-    (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,13]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,23])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$sycamore <- (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,14]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,24])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean
TTmass$willow <-   (exp(MassWeighted1719$Sol[,1]+MassWeighted1719$Sol[,15]+(MassWeighted1719$Sol[,4]+MassWeighted1719$Sol[,25])*0.91+(MassWeighted1719$Sol[,5])*0.91^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day160")+
  theme(text = element_text(size=15))

################################################
#### With date^2:tree and site interactions ####

load("~/Documents/Models/MWd2int.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWd2int$Sol[,1])+mean(MWd2int$Sol[,4])*preddayscaled+mean(MWd2int$Sol[,5])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- MWd2int$Sol[,16:25] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- MWd2int$Sol[,6:15] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

TTquad.cropped <- MWd2int$Sol[,26:35] # crop to just the columns wanted
TTquad <- data.frame(treetaxa=c(colnames(TTquad.cropped))) #dataframe with column for yearsite 
TTquad$coeff <- apply(TTquad.cropped,2, mean) # mean 
for(i in 1:length(TTquad$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTquad.cropped[,i])
  TTquad$lowci[i] <- A["var1","lower"] 
  TTquad$upci[i] <- A["var1","upper"] 
} 
TTquad$treetaxa <- TTslope$treetaxa

TTquadcoeff <- ggplot(TTquad, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Quadratic")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
row3 <- grid.arrange(TTquadcoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, row3,  nrow = 3, heights = c(1,1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2], Quad=TTquad[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWd2int$Sol[,1])+mean(MWd2int$Sol[,4])*preddayscaled+mean(MWd2int$Sol[,5])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[1,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[1,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[1,4])*preddayscaled^2
ash <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[2,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[2,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[2,4])*preddayscaled^2
beech <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[3,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[3,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[3,4])*preddayscaled^2
birch <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[4,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[4,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[4,4])*preddayscaled^2
elm <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[5,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[5,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[5,4])*preddayscaled^2
hazel <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[6,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[6,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[6,4])*preddayscaled^2
oak <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[7,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[7,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[7,4])*preddayscaled^2
rowan <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[8,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[8,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[8,4])*preddayscaled^2
sycamore <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[9,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[9,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[9,4])*preddayscaled^2
willow <- mean(MWd2int$Sol[,1])+TaxaCatGrowth[10,2]+mean(MWd2int$Sol[,4]+TaxaCatGrowth[10,3])*preddayscaled+mean(MWd2int$Sol[,5]+TaxaCatGrowth[10,4])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

# mass on 168- last day on which all tree taxa have had a caterpillar
TTmass <- data.frame(mean=(exp(MWd2int$Sol[,1]+MWd2int$Sol[,4]*0.96+MWd2int$Sol[,5]*0.96^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,6]+(MWd2int$Sol[,4]+MWd2int$Sol[,16])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,26])*0.96^2))-TTmass$mean
TTmass$ash <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,7]+(MWd2int$Sol[,4]+MWd2int$Sol[,17])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,27])*0.96^2))-TTmass$mean
TTmass$beech <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,8]+(MWd2int$Sol[,4]+MWd2int$Sol[,18])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,28])*0.96^2))-TTmass$mean
TTmass$birch <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,9]+(MWd2int$Sol[,4]+MWd2int$Sol[,19])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,29])*0.96^2))-TTmass$mean
TTmass$elm <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,10]+(MWd2int$Sol[,4]+MWd2int$Sol[,20])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,30])*0.96^2))-TTmass$mean
TTmass$hazel <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,11]+(MWd2int$Sol[,4]+MWd2int$Sol[,21])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,31])*0.96^2))-TTmass$mean
TTmass$oak <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,12]+(MWd2int$Sol[,4]+MWd2int$Sol[,22])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,32])*0.96^2))-TTmass$mean
TTmass$rowan <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,13]+(MWd2int$Sol[,4]+MWd2int$Sol[,23])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,33])*0.96^2))-TTmass$mean
TTmass$sycamore <- (exp(MWd2int$Sol[,1]+MWd2int$Sol[,14]+(MWd2int$Sol[,4]+MWd2int$Sol[,24])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,34])*0.96^2))-TTmass$mean
TTmass$willow <-   (exp(MWd2int$Sol[,1]+MWd2int$Sol[,15]+(MWd2int$Sol[,4]+MWd2int$Sol[,25])*0.96+(MWd2int$Sol[,5]+MWd2int$Sol[,35])*0.96^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day168")+
  theme(text = element_text(size=15))

# day 160
TTmass <- data.frame(mean=(exp(MWd2int$Sol[,1]+MWd2int$Sol[,4]*0.91+MWd2int$Sol[,5]*0.91^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,6]+(MWd2int$Sol[,4]+MWd2int$Sol[,16])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,26])*0.91^2))-TTmass$mean
TTmass$ash <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,7]+(MWd2int$Sol[,4]+MWd2int$Sol[,17])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,27])*0.91^2))-TTmass$mean
TTmass$beech <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,8]+(MWd2int$Sol[,4]+MWd2int$Sol[,18])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,28])*0.91^2))-TTmass$mean
TTmass$birch <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,9]+(MWd2int$Sol[,4]+MWd2int$Sol[,19])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,29])*0.91^2))-TTmass$mean
TTmass$elm <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,10]+(MWd2int$Sol[,4]+MWd2int$Sol[,20])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,30])*0.91^2))-TTmass$mean
TTmass$hazel <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,11]+(MWd2int$Sol[,4]+MWd2int$Sol[,21])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,31])*0.91^2))-TTmass$mean
TTmass$oak <-      (exp(MWd2int$Sol[,1]+MWd2int$Sol[,12]+(MWd2int$Sol[,4]+MWd2int$Sol[,22])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,32])*0.91^2))-TTmass$mean
TTmass$rowan <-    (exp(MWd2int$Sol[,1]+MWd2int$Sol[,13]+(MWd2int$Sol[,4]+MWd2int$Sol[,23])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,33])*0.91^2))-TTmass$mean
TTmass$sycamore <- (exp(MWd2int$Sol[,1]+MWd2int$Sol[,14]+(MWd2int$Sol[,4]+MWd2int$Sol[,24])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,34])*0.91^2))-TTmass$mean
TTmass$willow <-   (exp(MWd2int$Sol[,1]+MWd2int$Sol[,15]+(MWd2int$Sol[,4]+MWd2int$Sol[,25])*0.91+(MWd2int$Sol[,5]+MWd2int$Sol[,35])*0.91^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day160")+
  theme(text = element_text(size=15))

#### MWd2intAll : MWd2int with all tree taxa categories ####   Make sure correct dataset!! ALL BEATEN TREE TAXA
load("~/Documents/Models/MWd2intAll.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWd2intAll$Sol[,1])+mean(MWd2intAll$Sol[,4])*preddayscaled+mean(MWd2intAll$Sol[,5])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- MWd2intAll$Sol[,23:39] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- MWd2intAll$Sol[,6:22] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))

TTquad.cropped <- MWd2intAll$Sol[,40:56] # crop to just the columns wanted
TTquad <- data.frame(treetaxa=c(colnames(TTquad.cropped))) #dataframe with column for yearsite 
TTquad$coeff <- apply(TTquad.cropped,2, mean) # mean 
for(i in 1:length(TTquad$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTquad.cropped[,i])
  TTquad$lowci[i] <- A["var1","lower"] 
  TTquad$upci[i] <- A["var1","upper"] 
} 
TTquad$treetaxa <- TTslope$treetaxa

TTquadcoeff <- ggplot(TTquad, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Quadratic")+
  theme(text = element_text(size=10))

library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
row3 <- grid.arrange(TTquadcoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, row3,  nrow = 3, heights = c(1,1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2], Quad=TTquad[,2])   

preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWd2intAll$Sol[,1])+mean(MWd2intAll$Sol[,4])*preddayscaled+mean(MWd2intAll$Sol[,5])*preddayscaled^2
alder <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[2,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[2,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[2,4])*preddayscaled^2
ash <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[3,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[3,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[3,4])*preddayscaled^2
aspen <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[4,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[4,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[4,4])*preddayscaled^2
beech <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[5,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[5,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[5,4])*preddayscaled^2
birch <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[6,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[6,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[6,4])*preddayscaled^2
cherry <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[7,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[7,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[7,4])*preddayscaled^2
chestnut <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[8,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[8,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[8,4])*preddayscaled^2
damson <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[9,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[9,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[9,4])*preddayscaled^2
elm <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[10,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[10,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[10,4])*preddayscaled^2
fieldmaple <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[11,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[11,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[11,4])*preddayscaled^2
hazel <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[12,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[12,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[12,4])*preddayscaled^2
lime <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[13,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[13,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[13,4])*preddayscaled^2
oak <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[14,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[14,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[14,4])*preddayscaled^2
rowan <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[15,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[15,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[15,4])*preddayscaled^2
sycamore <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[16,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[16,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[16,4])*preddayscaled^2
willow <- mean(MWd2intAll$Sol[,1])+TaxaCatGrowth[17,2]+mean(MWd2intAll$Sol[,4]+TaxaCatGrowth[17,3])*preddayscaled+mean(MWd2intAll$Sol[,5]+TaxaCatGrowth[17,4])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(aspen), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(cherry), type="l", lwd=2, col=4)
points(predday,exp(chestnut), type="l", lwd=2, col=4)
points(predday,exp(damson), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(fieldmaple), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(lime), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

# mass on 168- last day on which all tree taxa have had a caterpillar
TTmass <- data.frame(mean=(exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,4]*0.96+MWd2intAll$Sol[,5]*0.96^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,7]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,24])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,41])*0.96^2))-TTmass$mean
TTmass$ash <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,8]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,25])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,42])*0.96^2))-TTmass$mean
TTmass$beech <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,10]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,27])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,44])*0.96^2))-TTmass$mean
TTmass$birch <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,11]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,28])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,45])*0.96^2))-TTmass$mean
TTmass$elm <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,15]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,32])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,49])*0.96^2))-TTmass$mean
TTmass$hazel <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,17]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,34])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,51])*0.96^2))-TTmass$mean
TTmass$oak <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,19]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,36])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,53])*0.96^2))-TTmass$mean
TTmass$rowan <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,20]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,37])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,54])*0.96^2))-TTmass$mean
TTmass$sycamore <- (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,21]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,38])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,55])*0.96^2))-TTmass$mean
TTmass$willow <-   (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,22]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,39])*0.96+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,56])*0.96^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day168")+
  theme(text = element_text(size=15))

# day 160
TTmass <- data.frame(mean=(exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,4]*0.91+MWd2intAll$Sol[,5]*0.91^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,7]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,24])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,41])*0.91^2))-TTmass$mean
TTmass$ash <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,8]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,25])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,42])*0.91^2))-TTmass$mean
TTmass$beech <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,10]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,27])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,44])*0.91^2))-TTmass$mean
TTmass$birch <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,11]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,28])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,45])*0.91^2))-TTmass$mean
TTmass$elm <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,15]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,32])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,49])*0.91^2))-TTmass$mean
TTmass$hazel <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,17]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,34])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,51])*0.91^2))-TTmass$mean
TTmass$oak <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,19]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,36])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,53])*0.91^2))-TTmass$mean
TTmass$rowan <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,20]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,37])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,54])*0.91^2))-TTmass$mean
TTmass$sycamore <- (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,21]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,38])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,55])*0.91^2))-TTmass$mean
TTmass$willow <-   (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,22]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,39])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,56])*0.91^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day160")+
  theme(text = element_text(size=15))

# day 160 compared to each other not the mean
TTmass <- data.frame(mean=(exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,4]*0.91+MWd2intAll$Sol[,5]*0.91^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,7]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,24])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,41])*0.91^2))
TTmass$ash <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,8]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,25])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,42])*0.91^2))
TTmass$beech <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,10]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,27])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,44])*0.91^2))
TTmass$birch <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,11]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,28])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,45])*0.91^2))
TTmass$elm <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,15]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,32])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,49])*0.91^2))
TTmass$hazel <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,17]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,34])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,51])*0.91^2))
TTmass$oak <-      (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,19]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,36])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,53])*0.91^2))
TTmass$rowan <-    (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,20]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,37])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,54])*0.91^2))
TTmass$sycamore <- (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,21]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,38])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,55])*0.91^2))
TTmass$willow <-   (exp(MWd2intAll$Sol[,1]+MWd2intAll$Sol[,22]+(MWd2intAll$Sol[,4]+MWd2intAll$Sol[,39])*0.91+(MWd2intAll$Sol[,5]+MWd2intAll$Sol[,56])*0.91^2))

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  #geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day160")+
  theme(text = element_text(size=15))


#### MWSY : MW model with no date^2 interacion and using SY ####
load("~/Documents/Models/MWSY3.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWSY3$Sol[,1])+mean(MWSY3$Sol[,2])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

TTslope.cropped <- MWSY3$Sol[,14:23] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datescaled.tree.species.","", TTslope$treetaxa)

TTslopecoeff <- ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

TTintercept.cropped <- MWSY3$Sol[,4:13] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- TTslope$treetaxa

TTinterceptcoeff <- ggplot(TTintercept, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))


library(gridExtra)

row1 <- grid.arrange(TTinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(TTslopecoeff, ncol = 1, widths = 1)
TTCoeff <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))

### slope and intercept per tax
TaxaCatGrowth <- data.frame(Taxa=TTslope[,1], Intercept= TTintercept[,2], Slope=TTslope[,2])   

preddayscaled <- seq(0.669,0.983,0.001)
predday <- preddayscaled*175
meanslope <- mean(MWSY3$Sol[,1])+mean(MWSY3$Sol[,2])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
points(predday,meanslope, type="l", lwd=2, col=3)
alder <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,4])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,14])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
ash <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,5])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,15])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
beech <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,6])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,16])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
birch <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,7])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,17])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
elm <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,8])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,18])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
hazel <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,9])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,19])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
oak <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,10])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,20])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
rowan <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,11])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,21])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
sycamore <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,12])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,22])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2
willow <- mean(MWSY3$Sol[,1]+MWSY3$Sol[,13])+mean(MWSY3$Sol[,2]+MWSY3$Sol[,23])*preddayscaled+mean(MWSY3$Sol[,3])*preddayscaled^2

### barchart of lower interval by date
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat$date, exp(cater_habitat$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=1, pch=20)
points(predday,exp(alder), type="l", lwd=2, col=4)
points(predday,exp(ash), type="l", lwd=2, col=4)
points(predday,exp(beech), type="l", lwd=2, col=4)
points(predday,exp(birch), type="l", lwd=2, col=4)
points(predday,exp(elm), type="l", lwd=2, col=4)
points(predday,exp(hazel), type="l", lwd=2, col=4)
points(predday,exp(oak), type="l", lwd=2, col=4)
points(predday,exp(rowan), type="l", lwd=2, col=4)
points(predday,exp(sycamore), type="l", lwd=2, col=4)
points(predday,exp(willow), type="l", lwd=2, col=4)
points(predday,exp(meanslope), type="l", lwd=2, col="red")

#### Plot with data as translucent ####
mycolblack <- rgb(0, 0, 0, max = 250, alpha = 25, names = "blacktrans")
mycolgrey <- rgb(128, 128, 128, max = 250, alpha = 25, names = "greytrans")
#AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass (g)", pch=20, col=mycolgrey, cex=0.3)
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=mycolblack, pch=20, cex=0.3)
points(predday,exp(alder), type="l", lwd=2, col="darkred")
points(predday,exp(ash), type="l", lwd=2, col="firebrick3")
points(predday,exp(beech), type="l", lwd=2, col="chocolate2")
points(predday,exp(birch), type="l", lwd=2, col="goldenrod")
points(predday,exp(elm), type="l", lwd=2, col="olivedrab4")
points(predday,exp(hazel), type="l", lwd=2, col="darkgreen")
points(predday,exp(oak), type="l", lwd=2, col="deepskyblue3")
points(predday,exp(rowan), type="l", lwd=2, col="royalblue4")
points(predday,exp(sycamore), type="l", lwd=2, col="slateblue2")
points(predday,exp(willow), type="l", lwd=2, col="orchid")
#points(predday,exp(meanslope), type="l", lwd=2, col=1)
par(new = T)
hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat$date, exp(cater_habitat$logmpc1), xlab="Date", ylab="Mass (g)", pch=20, col=mycolgrey, cex=0.3)
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=mycolblack, pch=20, cex=0.3)
points(predday,exp(alder), type="l", lwd=2, col="darkred")
points(predday,exp(ash), type="l", lwd=2, col="firebrick3")
points(predday,exp(beech), type="l", lwd=2, col="chocolate2")
points(predday,exp(birch), type="l", lwd=2, col="goldenrod")
points(predday,exp(elm), type="l", lwd=2, col="olivedrab4")
points(predday,exp(hazel), type="l", lwd=2, col="darkgreen")
points(predday,exp(oak), type="l", lwd=2, col="deepskyblue3")
points(predday,exp(rowan), type="l", lwd=2, col="royalblue4")
points(predday,exp(sycamore), type="l", lwd=2, col="slateblue2")
points(predday,exp(willow), type="l", lwd=2, col="orchid")
#points(predday,exp(meanslope), type="l", lwd=2, col=1)


#### Mean curve on data and TT curves on appropraite axis limits ####
plot(intsamples$date, intsamples$mpc1, xlab="Date", ylab="Mass (g)", pch=20, col=mycolblack, cex=0.4, xlim=c(117,172), ylim=c(0,1))
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=mycolblack, pch=20, cex=0.4)
points(predday,exp(meanslope), type="l", lwd=2, col=2)
legend("topleft", legend="Mean curve",
       lty=1, lwd=3, 
       col=2, 
       cex=0.8, seg.len=0.8, bty = "n")
plot(predday,exp(alder), type="l", lwd=1.5, col="darkred", xlab="Date", ylab="Mass (g)", ylim=c(0,0.05))
points(predday,exp(ash), type="l", lwd=1.5, col="firebrick3")
points(predday,exp(beech), type="l", lwd=1.5, col="chocolate2")
points(predday,exp(birch), type="l", lwd=1.5, col="goldenrod")
points(predday,exp(elm), type="l", lwd=1.5, col="olivedrab4")
points(predday,exp(hazel), type="l", lwd=1.5, col="darkgreen")
points(predday,exp(oak), type="l", lwd=1.5, col="deepskyblue3")
points(predday,exp(rowan), type="l", lwd=1.5, col="royalblue4")
points(predday,exp(sycamore), type="l", lwd=1.5, col="slateblue2")
points(predday,exp(willow), type="l", lwd=1.5, col="orchid")
points(predday,exp(meanslope), type="l", lwd=2, col=2, lty="dotted")
points(rep(168, 101), seq(0,0.05, 0.0005), type="l", col=1, lwd=0.5, lty="dashed")
legend("topleft", legend=c("Alder","Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow", "Mean"),
       lty=c(1,1,1,1,1,1,1,1,1,1,3), lwd=3, 
       col=c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid",2), 
       cex=0.8, seg.len=0.8, bty = "n")


# mass on 168- last day on which all tree taxa have had a caterpillar
TTmass <- data.frame(mean=(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$ash <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$beech <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$birch <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$elm <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$hazel <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$oak <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$rowan <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$sycamore <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean
TTmass$willow <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.96+(MWSY3$Sol[,3])*0.96^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day168")+
  theme(text = element_text(size=15))

# day 160
TTmass <- data.frame(mean=(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.91+MWSY3$Sol[,3]*0.91^2)))
colnames(TTmass)[colnames(TTmass)=="var1"] <- "mean"
TTmass$alder <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$ash <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$beech <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$birch <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$elm <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$hazel <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$oak <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$rowan <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$sycamore <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean
TTmass$willow <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.91+(MWSY3$Sol[,3])*0.91^2))-TTmass$mean

TTmass.cropped <- TTmass[,2:11] # crop to just the columns wanted
TTmass2 <- data.frame(treetaxa=c(colnames(TTmass.cropped))) #dataframe with column for yearsite 
TTmass2$coeff <- apply(TTmass.cropped,2, mean) # mean 
for(i in 1:length(TTmass2$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTmass.cropped[,i])
  TTmass2$lowci[i] <- A["var1","lower"] 
  TTmass2$upci[i] <- A["var1","upper"] 
}

ggplot(TTmass2, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  ylab("Difference from mean mass")+
  #ylim(-0.03,0.09)+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Day160")+
  theme(text = element_text(size=15))


#### Day 165 (0.9428571) difference to mean ####
Alder165 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Ash165 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Beech165 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Birch165 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Elm165 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Hazel165 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Oak165 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Rowan165 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Sycamore165 <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))
Willow165 <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.943+MWSY3$Sol[,3]*0.943^2))

Mass165 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                    mean=c(mean(Alder165), mean(Ash165), mean(Beech165), mean(Birch165), mean(Elm165), mean(Hazel165), mean(Oak165), mean(Rowan165), mean(Sycamore165), mean(Willow165)),
                    lowci=c(HPDinterval(Alder165)[1],HPDinterval(Ash165)[1],HPDinterval(Beech165)[1],HPDinterval(Birch165)[1],HPDinterval(Elm165)[1],HPDinterval(Hazel165)[1],HPDinterval(Oak165)[1],HPDinterval(Rowan165)[1],HPDinterval(Sycamore165)[1],HPDinterval(Willow165)[1]),
                    upci=c(HPDinterval(Alder165)[2],HPDinterval(Ash165)[2],HPDinterval(Beech165)[2],HPDinterval(Birch165)[2],HPDinterval(Elm165)[2],HPDinterval(Hazel165)[2],HPDinterval(Oak165)[2],HPDinterval(Rowan165)[2],HPDinterval(Sycamore165)[2],HPDinterval(Willow165)[2]))

Mass1 <- ggplot(Mass165, aes(fct_rev(TT), mean))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to mean (g)")+
  theme_bw()+
  labs(tag = "Day165")+
  theme(text=element_text(size= 20))


#### Difference to oak (on 165: 0.9428571) #### 

AlderOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
AshOD <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
BeechOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
BirchOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
ElmOD <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
HazelOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
RowanOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
SycamoreOD <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))
WillowOD <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.943+(MWSY3$Sol[,3])*0.943^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.943+(MWSY3$Sol[,3])*0.943^2))

OD165 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                   mean=c(mean(AlderOD), mean(AshOD), mean(BeechOD), mean(BirchOD), mean(ElmOD), mean(HazelOD), mean(RowanOD), mean(SycamoreOD), mean(WillowOD)),
                   lowci=c(HPDinterval(AlderOD)[1],HPDinterval(AshOD)[1],HPDinterval(BeechOD)[1],HPDinterval(BirchOD)[1],HPDinterval(ElmOD)[1],HPDinterval(HazelOD)[1],HPDinterval(RowanOD)[1],HPDinterval(SycamoreOD)[1],HPDinterval(WillowOD)[1]),
                   upci=c(HPDinterval(AlderOD)[2],HPDinterval(AshOD)[2],HPDinterval(BeechOD)[2],HPDinterval(BirchOD)[2],HPDinterval(ElmOD)[2],HPDinterval(HazelOD)[2],HPDinterval(RowanOD)[2],HPDinterval(SycamoreOD)[2],HPDinterval(WillowOD)[2]))

Mass2 <- ggplot(OD165, aes(fct_rev(TT), mean))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("")+
  ylab("Mass difference to oak (g)")+
  annotate("text", label="*", x="Beech", y=-0.0018, size=7)+
  annotate("text", label="*", x="Birch", y=-0.0008, size=7)+
  annotate("text", label="*", x="Willow", y=-0.002, size=7)+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 20))

col1 <- grid.arrange(Mass1, nrow = 1, heights = 1)
col2 <- grid.arrange(Mass2, nrow = 1, heights = 1)
MassDif165 <- grid.arrange(col1, col2, ncol = 2, widths = c(1.2,1))

#### Difference to oak (on 168: 0.96) #### 

AlderOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
AshOD <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
BeechOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
BirchOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
ElmOD <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
HazelOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
RowanOD <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
SycamoreOD <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))
WillowOD <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))

OD168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                    mean=c(mean(AlderOD), mean(AshOD), mean(BeechOD), mean(BirchOD), mean(ElmOD), mean(HazelOD), mean(RowanOD), mean(SycamoreOD), mean(WillowOD)),
                    lowci=c(HPDinterval(AlderOD)[1],HPDinterval(AshOD)[1],HPDinterval(BeechOD)[1],HPDinterval(BirchOD)[1],HPDinterval(ElmOD)[1],HPDinterval(HazelOD)[1],HPDinterval(RowanOD)[1],HPDinterval(SycamoreOD)[1],HPDinterval(WillowOD)[1]),
                    upci=c(HPDinterval(AlderOD)[2],HPDinterval(AshOD)[2],HPDinterval(BeechOD)[2],HPDinterval(BirchOD)[2],HPDinterval(ElmOD)[2],HPDinterval(HazelOD)[2],HPDinterval(RowanOD)[2],HPDinterval(SycamoreOD)[2],HPDinterval(WillowOD)[2]))

Mass4 <- ggplot(OD168, aes(fct_rev(TT), mean))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("")+
  ylab("Mass difference to oak (g)")+
  annotate("text", label="*", x="Beech", y=-0.0018, size=7)+
  annotate("text", label="*", x="Birch", y=-0.0006, size=7)+
  annotate("text", label="*", x="Willow", y=-0.0018, size=7)+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 20))

#### Day 168 (0.96) difference to mean ####
Alder168 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,4]+(MWSY3$Sol[,2]+MWSY3$Sol[,14])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Ash168 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,5]+(MWSY3$Sol[,2]+MWSY3$Sol[,15])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Beech168 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,6]+(MWSY3$Sol[,2]+MWSY3$Sol[,16])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Birch168 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,7]+(MWSY3$Sol[,2]+MWSY3$Sol[,17])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Elm168 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,8]+(MWSY3$Sol[,2]+MWSY3$Sol[,18])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Hazel168 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,9]+(MWSY3$Sol[,2]+MWSY3$Sol[,19])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Oak168 <-      (exp(MWSY3$Sol[,1]+MWSY3$Sol[,10]+(MWSY3$Sol[,2]+MWSY3$Sol[,20])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Rowan168 <-    (exp(MWSY3$Sol[,1]+MWSY3$Sol[,11]+(MWSY3$Sol[,2]+MWSY3$Sol[,21])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Sycamore168 <- (exp(MWSY3$Sol[,1]+MWSY3$Sol[,12]+(MWSY3$Sol[,2]+MWSY3$Sol[,22])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))
Willow168 <-   (exp(MWSY3$Sol[,1]+MWSY3$Sol[,13]+(MWSY3$Sol[,2]+MWSY3$Sol[,23])*0.96+(MWSY3$Sol[,3])*0.96^2))-(exp(MWSY3$Sol[,1]+MWSY3$Sol[,2]*0.96+MWSY3$Sol[,3]*0.96^2))

Mass168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                      mean=c(mean(Alder168), mean(Ash168), mean(Beech168), mean(Birch168), mean(Elm168), mean(Hazel168), mean(Oak168), mean(Rowan168), mean(Sycamore168), mean(Willow168)),
                      lowci=c(HPDinterval(Alder168)[1],HPDinterval(Ash168)[1],HPDinterval(Beech168)[1],HPDinterval(Birch168)[1],HPDinterval(Elm168)[1],HPDinterval(Hazel168)[1],HPDinterval(Oak168)[1],HPDinterval(Rowan168)[1],HPDinterval(Sycamore168)[1],HPDinterval(Willow168)[1]),
                      upci=c(HPDinterval(Alder168)[2],HPDinterval(Ash168)[2],HPDinterval(Beech168)[2],HPDinterval(Birch168)[2],HPDinterval(Elm168)[2],HPDinterval(Hazel168)[2],HPDinterval(Oak168)[2],HPDinterval(Rowan168)[2],HPDinterval(Sycamore168)[2],HPDinterval(Willow168)[2]))

Mass3 <- ggplot(Mass168, aes(fct_rev(TT), mean))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to mean (g)")+
  theme_bw()+
  labs(tag = "Day168")+
  theme(text=element_text(size= 20))

col1 <- grid.arrange(Mass3, nrow = 1, heights = 1)
col2 <- grid.arrange(Mass4, nrow = 1, heights = 1)
MassDif168 <- grid.arrange(col1, col2, ncol = 2, widths = c(1.2,1))
#############################
#### Looking at siteyear ####
#############################
load("~/Documents/Models/MWSY.RData")

#### Plotting growth curves ####
preddayscaled <- seq(0.66,1,0.01)
predday <- preddayscaled*max(cater_habitat$date)
meanslope <- mean(MWSY$Sol[,1])+mean(MWSY$Sol[,2])*preddayscaled+mean(MWSY$Sol[,3])*preddayscaled^2

par(mfcol=c(1,2), cex=1.5)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), log="y", xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")
#par(new = T)
#hist(intsamples$date, breaks=100, axes=F, xlab=NA, ylab=NA, ylim=c(0,500), xlim=c(117,175), main="")
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)
points(predday, exp(meanslope), col=2, type="l")

SYslope.cropped <- MWSY$Sol[,156:287] # crop to just the columns wanted
SYslope <- data.frame(siteyear=c(colnames(SYslope.cropped))) #dataframe with column for yearsite 
SYslope$coeff <- apply(SYslope.cropped,2, mean) # mean 
for(i in 1:length(SYslope$siteyear)) {   # loop for CIs
  A <- HPDinterval(SYslope.cropped[,i])
  SYslope$lowci[i] <- A["var1","lower"] 
  SYslope$upci[i] <- A["var1","upper"] 
} 
SYslope$siteyear <- gsub("datescaled.siteyear.","", SYslope$siteyear)

SYslopecoeff <- ggplot(SYslope, aes(siteyear, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site Year")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

SYintercept.cropped <- MWSY$Sol[,24:155] # crop to just the columns wanted
SYintercept <- data.frame(siteyear=c(colnames(SYintercept.cropped))) #dataframe with column for yearsite 
SYintercept$coeff <- apply(SYintercept.cropped,2, mean) # mean 
for(i in 1:length(SYintercept$siteyear)) {   # loop for CIs
  A <- HPDinterval(SYintercept.cropped[,i])
  SYintercept$lowci[i] <- A["var1","lower"] 
  SYintercept$upci[i] <- A["var1","upper"] 
} 
SYintercept$siteyear <- SYslope$siteyear

SYinterceptcoeff <- ggplot(SYintercept, aes(siteyear, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site Year")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Intercept")+
  theme(text = element_text(size=10))


library(gridExtra)

row1 <- grid.arrange(SYinterceptcoeff, ncol = 1, widths = 1)
row2 <- grid.arrange(SYslopecoeff, ncol = 1, widths = 1)
SYCoeff <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))

SYslope$site <- SYslope$siteyear
SYslope$site <- gsub(" 2017","", SYslope$site)
SYslope$site <- gsub(" 2018","", SYslope$site)
SYslope$site <- gsub(" 2019","", SYslope$site)

site <- read.csv("Dropbox/master_data/site/site_details.csv")
colnames(site)[4:6] <- c("latitude", "longitude", "elevation")
store <- pmatch(SYslope$site, site$site, duplicates.ok = TRUE)
SYslope <- cbind(SYslope, site$elevation[store])
colnames(SYslope)[6] <- c("elevation")
SYslope <- SYslope[order(SYslope$elevation),] 

SYslopeElev <- ggplot(SYslope, aes(fct_inorder(siteyear), coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site Year")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
  theme(text = element_text(size=10))

SYintercept$site <- SYintercept$siteyear
SYintercept$site <- gsub(" 2017","", SYintercept$site)
SYintercept$site <- gsub(" 2018","", SYintercept$site)
SYintercept$site <- gsub(" 2019","", SYintercept$site)

store <- pmatch(SYintercept$site, site$site, duplicates.ok = TRUE)
SYintercept <- cbind(SYintercept, site$elevation[store])
colnames(SYintercept)[6] <- c("elevation")
SYintercept <- SYintercept[order(SYintercept$elevation),] 

SYinterceptElev <- ggplot(SYintercept, aes(fct_inorder(siteyear), coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Site Year")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "intercept")+
  theme(text = element_text(size=10))

SYcurves <- data.frame(siteyear=SYintercept$siteyear, site=SYintercept$site, elev=SYintercept$elevation, int=SYintercept$coeff, slope=SYslope$coeff)
dayscal <- seq(0.67,1,0.001)
curve <- mean(MWSY$Sol[,1])+mean(MWSY$Sol[,2])*dayscal+mean(MWSY$Sol[,3])*dayscal^2
days <- dayscal*max(cater_habitat$date)

par(mfcol=c(1,1),mar=c(3.9, 3.8, 1, 1), cex=1.4, las=1)
plot(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc1), xlab="Date", ylab="Mass", pch=20, col="grey")
points(cater_habitat_1719$date, exp(cater_habitat_1719$logmpc2), col=1, pch=20)

for(i in 1:132){
  A <- mean(MWSY$Sol[,1])+SYcurves[i,4]+(mean(MWSY$Sol[,2])+SYcurves[i,5])*dayscal+mean(MWSY$Sol[,3])*dayscal^2
  points(days, exp(A), type="l", col=4, lwd=0.5)
}
points(days,exp(curve), col=2, type="l", lwd=1.5)
