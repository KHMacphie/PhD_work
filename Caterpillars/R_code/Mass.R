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
cater$mpc1 <- revalue(cater$mpc1, c("0.02"="0.002"))
cater$mpc1 <- as.numeric(cater$mpc1)
cater$mpc1 <- cater$mpc1/cater$caterpillars
cater$mpc1 <- ifelse(cater$mpc1 < 0.002,0.002,cater$mpc1)

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


#Mass4<- MCMCglmm(cbind(logmpc1, logmpc2)~ datecent, 
#                 random=~us(1+datecent):tree.species + us(1+datecent):year + us(1+datecent):site + treeID + recorder + siteday, 
#                 family="cengaussian", data=cater.expanded, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)
#save(Mass4, file = "~/Documents/Models/Mass4.RData")
load("~/Documents/Models/Mass4.RData")

TTintercept.cropped <- Mass4$Sol[,3:18] # crop to just the columns wanted
TTintercept <- data.frame(treetaxa=c(colnames(TTintercept.cropped))) #dataframe with column for yearsite 
TTintercept$coeff <- apply(TTintercept.cropped,2, mean) # mean 
for(i in 1:length(TTintercept$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTintercept.cropped[,i])
  TTintercept$lowci[i] <- A["var1","lower"] 
  TTintercept$upci[i] <- A["var1","upper"] 
} 
TTintercept$treetaxa <- gsub("(Intercept).tree.species.","", TTintercept$treetaxa)
TTintercept$treetaxa <- gsub(".tree.species.","", TTintercept$treetaxa)
TTintercept$treetaxa <- gsub("(Intercept)","", TTintercept$treetaxa)
TTintercept$treetaxa <- gsub("()","", TTintercept$treetaxa) # not working

par(mfcol=c(1,1), cex=1.5)
hist(Mass4$VCV[,"(Intercept):(Intercept).tree.species"], breaks=500, xlim=c(0,0.2))
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

TTslope.cropped <- Mass4$Sol[,19:34] # crop to just the columns wanted
TTslope <- data.frame(treetaxa=c(colnames(TTslope.cropped))) #dataframe with column for yearsite 
TTslope$coeff <- apply(TTslope.cropped,2, mean) # mean 
for(i in 1:length(TTslope$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTslope.cropped[,i])
  TTslope$lowci[i] <- A["var1","lower"] 
  TTslope$upci[i] <- A["var1","upper"] 
} 
TTslope$treetaxa <- gsub("datecent.tree.species.","", TTslope$treetaxa)

par(mfcol=c(1,1), cex=1.5)
hist(Mass4$VCV[,"datecent:datecent.tree.species"], breaks=500)
abline(v=0, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
#legend("topright", legend="B", bty="n") 

par(mfcol=c(1,1), cex=0.8)
ggplot(TTslope, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  labs(tag = "Slope")+
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
