#### BOU Conference analysis ####

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
library(DHARMa)
library(ggeffects)
library(LaplacesDemon)
library(schoolmath)

birdphen <- read.csv("~/Dropbox/master_data/blue tits/Bird_Phenology.csv")
nestlings <- read.csv("~/Dropbox/master_data/blue tits/Nestlings.csv")
swap <- read.csv("~/Dropbox/master_data/blue tits/clutchswaps.csv")
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")


nestlings$SYB <- paste(nestlings$site, nestlings$box, nestlings$year)

SYB <- data.frame(SYB=nestlings$SYB, fledged=nestlings$fledged)
SYB <- aggregate(.~SYB, SYB, sum)

store <- pmatch(SYB$SYB,nestlings$SYB)
SYB$site <- nestlings$site[store]
SYB$year <- nestlings$year[store]
SYB$SY <- paste(SYB$site, SYB$year)

fledged <- data.frame(SYB=nestlings$SYB, fledged=nestlings$fledged)

birdphen$SYB <- paste(birdphen$site, birdphen$box, birdphen$year)
store2 <- pmatch(SYB$SYB,birdphen$SYB)
SYB$hd <- as.numeric(birdphen$hd_1.45[store2])

swap$SYB <- paste(swap$site, swap$destination.nest, swap$year)
store7 <- pmatch(SYB$SYB,swap$SYB)
SYB$cs2 <- swap$clutch.size[store7]
SYB$exp <- (2*SYB$cs2)/SYB$cs2
SYB$exp[which(is.na(SYB$exp)==T)] <- 1
SYB$exp <- as.factor(SYB$exp)

SYB$cs1 <- as.numeric(birdphen$cs[store2])
SYB$cs <- ifelse(SYB$exp==2, SYB$cs2, SYB$cs1)
SYB$propfl <- SYB$fledged/SYB$cs

SYB <- SYB[-which(SYB$propfl>1),]
SYB$cs1 <- NULL
SYB$cs2 <- NULL

ggplot(SYB, aes(propfl, fledged, col=exp))+
  geom_jitter()


#####################################
#### Site-year caterpillar peaks ####
#####################################

#### Estimates for each site-year ####
symetrics <- read.csv("Documents/Models/symetrics.csv")
symetrics$SY <- symetrics$siteyear

#Using mcmcglmm estimates
store3<- pmatch(SYB$SY, symetrics$SY, duplicates.ok = TRUE)
SYB$mu <- symetrics$PD[store3]
SYB$logheight <- symetrics$PlH[store3]
SYB$sigma <- symetrics$PS[store3]
#SYB <- na.omit(SYB)
SYB$mu <- as.integer(SYB$mu)
SYB$async <- (SYB$hd+10)-SYB$mu

#### Timing: bird vs caterpillars ####
SYB$birdtime <- (SYB$hd+10)
SYB$site <- as.factor(SYB$site)
SYB$year <- as.factor(SYB$year)

cater$SY <- paste(cater$site, cater$year)
count <- data.frame(SY=cater$SY, cater=cater$caterpillars)
count <- aggregate(.~SY, count, sum)
store9 <- pmatch(SYB$SY, count$SY, duplicates.ok=TRUE)
SYB$cater <- count$cater[store9]
SYB$PSna <- symetrics$PSna[store3]
SYB <- na.omit(SYB)
certcater <- SYB

certcater <- certcater[-which(certcater$PSna>0),] #removing any with NAs in cater peak
certcater <- certcater[-which(certcater$cater<10),] #and any with <10 cater
#n=256

SYBbySY <- data.frame(site=certcater$site, year=certcater$year, sy=certcater$SY, cater=certcater$mu, bird=certcater$birdtime, async=certcater$async)
SYBbySY <- aggregate(.~site+year+sy, SYBbySY, mean)

timemod <- lmer(birdtime~ mu + (1|SY) + (1|site) + (1|year), certcater)#
summary(timemod)

same <- data.frame(x=seq(136,166,1),y=seq(136,166,1))

timepred <- ggpredict(timemod,"mu")
plot(timepred, col=2)+  
  geom_line(data=same, aes(x,y), lty="dashed")+
  geom_point(data=SYBbySY, aes(cater,bird), col="darkgray", pch=20)+
  xlab("Caterpillar timing (ordinal date)")+
  ylab("Bird timing (ordinal date)")+
  ggtitle("")+
  theme(text=element_text(size= 15))

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
                   
timing.mod <- MCMCglmm(birdtime~mu, random=~site+year+SY, data=certcater, nitt=50000, burnin=15000, thin=20, pr=TRUE, prior=prior) 
plot(timing.mod)
summary(timing.mod)

pred <- data.frame(mu=seq(min(certcater$mu), max(certcater$mu), 3))

for(i in 1:length(pred$mu)){
A <- timing.mod$Sol[,1]+timing.mod$Sol[,2]*pred$mu[i]
pred$mean[i] <- mean(A)
pred$lowci[i] <- HPDinterval(mcmc(A))[1]
pred$upci[i] <- HPDinterval(mcmc(A))[2]
}

mycolblack <- rgb(100, 100, 100, max = 250, alpha = 100, names = "blacktrans")

rounded <- data.frame(term=c("intercept", "mu"), mean=c(mean(timing.mod$Sol[,1]),mean(timing.mod$Sol[,2])), lowci=c(HPDinterval(timing.mod$Sol[,1])[1],HPDinterval(timing.mod$Sol[,2])[1]),upci=c(HPDinterval(timing.mod$Sol[,1])[2],HPDinterval(timing.mod$Sol[,2])[2]))
rounded[,2:4] <- round(rounded[,2:4], 2)

text <- paste("Slope: ",rounded[2,2]," CIs: (",rounded[2,3]," - ",rounded[2,4],")")

ggplot(data=pred, aes(mu, mean))+
  xlim(min(certcater$mu), max(certcater$mu))+
  ylim(min(certcater$birdtime), max(certcater$birdtime))+
  geom_point(data=certcater, aes(mu,birdtime), col=mycolblack, pch=20)+
  geom_line(data=same, aes(x,y), lty="dashed", col=2)+
  geom_line(data=pred, aes(mu, mean), lwd=1)+
  geom_ribbon(data=pred, aes(x=mu, ymin=lowci, ymax=upci), alpha=0.15)+
  xlab("Caterpillar timing (ordinal date)")+
  ylab("Bird timing (ordinal date)")+
  geom_label(label=text, x=142,y=172, label.size = 0.5, label.padding = unit(0.55, "lines"))+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"


ggplot(SYB, aes(async, fledged))+geom_point()+geom_smooth(model=lm) #go with fledged, l0ok at with an without numbe of eggs/hatched 
ggplot(SYB, aes(async, propfl))+geom_point()+geom_smooth(model=lm)

SYB$logheight.scal <- scale(SYB$logheight) # cent: -2.438593, scale: 1.15297
SYB$sigma.scal <- scale(SYB$sigma) # cent: 11.18239, scale: 2.860743
SYB$cs.cent <- SYB$cs-as.integer(mean(SYB$cs)) #mean as integer= 8

SYB$failed <- SYB$cs-SYB$fledged

#### Models ####
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

fledged <- MCMCglmm(fledged~async+I(async^2)+logheight.scal+sigma.scal, random=~site + year + SY, data=SYB, nitt=50000, burnin=15000, thin=20, prior=prior)
plot(fledged)
summary(fledged)

fledged2 <- MCMCglmm(fledged~async+I(async^2)+logheight.scal+sigma.scal+cs.cent, random=~site + year + SY, data=SYB, nitt=50000, burnin=15000, thin=20, prior=prior)
plot(fledged2)
summary(fledged2)

propfledg <- MCMCglmm(propfl~async+I(async^2)+logheight.scal+sigma.scal, random=~site + year + SY, data=SYB, nitt=50000, burnin=15000, thin=20, prior=prior)
plot(propfledg)
summary(propfledg)

prior2<-list(R=list(V=1,fix=1),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

binomfledg <- MCMCglmm(cbind(fledged,failed) ~async+I(async^2)+logheight.scal+sigma.scal, random=~site + year + SY, data=SYB, family="multinomial2", nitt=50000, burnin=15000, thin=20, prior=prior2)
plot(binomfledg)
summary(binomfledg)

prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

binomfledg2 <- MCMCglmm(cbind(fledged,failed) ~async+I(async^2)+logheight.scal+sigma.scal, random=~site + year + SY, data=SYB, family="multinomial2", nitt=50000, burnin=15000, thin=20, prior=prior3)
plot(binomfledg2)
summary(binomfledg2)


#### Summarise SYB into async bins of 2 days ####
bins <- data.frame(async=SYB$async, fledged=SYB$fledged, propfl=SYB$propfl)
bins$bins <- ifelse(is.even(bins$async)==T, bins$async-0.5, bins$async+0.5)
bins$async <- NULL
bins <- aggregate(.~bins,bins,mean)

ggplot(data=p.fledge, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$fledged), 10))+
  geom_point(data=bins, aes(bins, fledged), pch=20)+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

ggplot(data=p.fledge, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$fledged), max(SYB$fledged)))+
  geom_point(data=bins, aes(bins, fledged), pch=20, size=3)+
  geom_point(data=SYB, aes(async, fledged), col=mycolblack, pch=20)+
    xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

ggplot(data=p.propfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=bins, aes(bins, propfl), pch=20)+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

ggplot(data=p.propfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=bins, aes(bins, propfl), pch=20, size=3)+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

#### Chicks fledged ####
p.fledge <- data.frame(async=seq(min(SYB$async), max(SYB$async), 2))# 18, 4))

for(i in 1:length(p.fledge$async)){
  
  A <- fledged$Sol[,1]+fledged$Sol[,2]*p.fledge$async[i]+fledged$Sol[,3]*p.fledge$async[i]^2
  
  p.fledge$mean[i] <- mean(A)
  p.fledge$lowci[i] <- HPDinterval(mcmc(A))[1]
  p.fledge$upci[i] <- HPDinterval(mcmc(A))[2]
}


#general mismatch
ggplot(data=p.fledge, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$fledged), max(SYB$fledged)))+
  geom_point(data=SYB, aes(async, fledged), col=mycolblack, pch=20)+
  geom_line(data=p.fledge, aes(async, mean), lwd=1)+
  geom_ribbon(data=p.fledge, aes(x=async, ymin=lowci, ymax=upci), alpha=0.15)+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

#peak height
for(i in 1:length(p.fledge$async)){
  
  A <- fledged$Sol[,1]+fledged$Sol[,2]*p.fledge$async[i]+fledged$Sol[,3]*p.fledge$async[i]^2 + 2*fledged$Sol[,4]
  B <- fledged$Sol[,1]+fledged$Sol[,2]*p.fledge$async[i]+fledged$Sol[,3]*p.fledge$async[i]^2 - 2*fledged$Sol[,4]
  
  p.fledge$himean[i] <- mean(A)
  p.fledge$hilowci[i] <- HPDinterval(mcmc(A))[1]
  p.fledge$hiupci[i] <- HPDinterval(mcmc(A))[2]
  p.fledge$lomean[i] <- mean(B)
  p.fledge$lolowci[i] <- HPDinterval(mcmc(B))[1]
  p.fledge$loupci[i] <- HPDinterval(mcmc(B))[2]
}

ggplot(data=p.fledge, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$fledged), max(SYB$fledged)))+
  geom_point(data=SYB, aes(async, fledged), col=mycolblack, pch=20)+
  #geom_line(data=p.fledge, aes(async, mean), lwd=1)+
  geom_line(data=p.fledge, aes(async, himean), lwd=1, col=2)+
  geom_ribbon(data=p.fledge, aes(x=async, ymin=hilowci, ymax=hiupci), alpha=0.15, fill="red")+
  geom_line(data=p.fledge, aes(async, lomean), lwd=1, col="blue")+
  geom_ribbon(data=p.fledge, aes(x=async, ymin=lolowci, ymax=loupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7

#peak width
for(i in 1:length(p.fledge$async)){
  
  A <- fledged$Sol[,1]+fledged$Sol[,2]*p.fledge$async[i]+fledged$Sol[,3]*p.fledge$async[i]^2 + 1.8*fledged$Sol[,5]
  B <- fledged$Sol[,1]+fledged$Sol[,2]*p.fledge$async[i]+fledged$Sol[,3]*p.fledge$async[i]^2 - 1.8*fledged$Sol[,5]
  
  p.fledge$widemean[i] <- mean(A)
  p.fledge$widelowci[i] <- HPDinterval(mcmc(A))[1]
  p.fledge$wideupci[i] <- HPDinterval(mcmc(A))[2]
  p.fledge$slimmean[i] <- mean(B)
  p.fledge$slimlowci[i] <- HPDinterval(mcmc(B))[1]
  p.fledge$slimupci[i] <- HPDinterval(mcmc(B))[2]
}

ggplot(data=p.fledge, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$fledged), max(SYB$fledged)))+
  geom_point(data=SYB, aes(async, fledged), col=mycolblack, pch=20)+
  #geom_line(data=p.fledge, aes(async, mean), lwd=1)+
  #geom_ribbon(data=p.fledge, aes(x=async, ymin=lowci, ymax=upci), alpha=0.15)+
  geom_line(data=p.fledge, aes(async, widemean), lwd=1, col=2)+
  geom_ribbon(data=p.fledge, aes(x=async, ymin=widelowci, ymax=wideupci), alpha=0.15, fill="red")+
  geom_line(data=p.fledge, aes(async, slimmean), lwd=1, col="blue")+
  geom_ribbon(data=p.fledge, aes(x=async, ymin=slimlowci, ymax=slimupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7

#### Prop fledged ####
p.propfl <- data.frame(async=seq(min(SYB$async), max(SYB$async), 2))# 18, 4))

for(i in 1:length(p.propfl$async)){
  
  A <- propfledg$Sol[,1]+propfledg$Sol[,2]*p.propfl$async[i]+propfledg$Sol[,3]*p.propfl$async[i]^2
  
  p.propfl$mean[i] <- mean(A)
  p.propfl$lowci[i] <- HPDinterval(mcmc(A))[1]
  p.propfl$upci[i] <- HPDinterval(mcmc(A))[2]
}

#general mismatch
ggplot(data=p.propfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  geom_line(data=p.propfl, aes(async, mean), lwd=1)+
  geom_ribbon(data=p.propfl, aes(x=async, ymin=lowci, ymax=upci), alpha=0.15)+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

#peak height
for(i in 1:length(p.propfl$async)){
  
  A <- propfledg$Sol[,1]+propfledg$Sol[,2]*p.propfl$async[i]+propfledg$Sol[,3]*p.propfl$async[i]^2 + 2*propfledg$Sol[,4]
  B <- propfledg$Sol[,1]+propfledg$Sol[,2]*p.propfl$async[i]+propfledg$Sol[,3]*p.propfl$async[i]^2 - 2*propfledg$Sol[,4]
  
  p.propfl$himean[i] <- mean(A)
  p.propfl$hilowci[i] <- HPDinterval(mcmc(A))[1]
  p.propfl$hiupci[i] <- HPDinterval(mcmc(A))[2]
  p.propfl$lomean[i] <- mean(B)
  p.propfl$lolowci[i] <- HPDinterval(mcmc(B))[1]
  p.propfl$loupci[i] <- HPDinterval(mcmc(B))[2]
}

ggplot(data=p.propfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  #geom_line(data=p.propfl, aes(async, mean), lwd=1)+
  geom_line(data=p.propfl, aes(async, himean), lwd=1, col=2)+
  geom_ribbon(data=p.propfl, aes(x=async, ymin=hilowci, ymax=hiupci), alpha=0.15, fill="red")+
  geom_line(data=p.propfl, aes(async, lomean), lwd=1, col="blue")+
  geom_ribbon(data=p.propfl, aes(x=async, ymin=lolowci, ymax=loupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7

#peak width
for(i in 1:length(p.propfl$async)){
  
  A <- propfledg$Sol[,1]+propfledg$Sol[,2]*p.propfl$async[i]+propfledg$Sol[,3]*p.propfl$async[i]^2 + 1.8*propfledg$Sol[,5]
  B <- propfledg$Sol[,1]+propfledg$Sol[,2]*p.propfl$async[i]+propfledg$Sol[,3]*p.propfl$async[i]^2 - 1.8*propfledg$Sol[,5]
  
  p.propfl$widemean[i] <- mean(A)
  p.propfl$widelowci[i] <- HPDinterval(mcmc(A))[1]
  p.propfl$wideupci[i] <- HPDinterval(mcmc(A))[2]
  p.propfl$slimmean[i] <- mean(B)
  p.propfl$slimlowci[i] <- HPDinterval(mcmc(B))[1]
  p.propfl$slimupci[i] <- HPDinterval(mcmc(B))[2]
}

ggplot(data=p.propfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  #geom_line(data=p.propfl, aes(async, mean), lwd=1)+
  geom_line(data=p.propfl, aes(async, widemean), lwd=1, col=2)+
  geom_ribbon(data=p.propfl, aes(x=async, ymin=widelowci, ymax=wideupci), alpha=0.15, fill="red")+
  geom_line(data=p.propfl, aes(async, slimmean), lwd=1, col="blue")+
  geom_ribbon(data=p.propfl, aes(x=async, ymin=slimlowci, ymax=slimupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7

#### Prop fledged- binomial####
p.binomfl <- data.frame(async=seq(min(SYB$async), max(SYB$async), 2))# 18, 4))

for(i in 1:length(p.binomfl$async)){
  
  A <- binomfledg2$Sol[,1]+binomfledg2$Sol[,2]*p.binomfl$async[i]+binomfledg2$Sol[,3]*p.binomfl$async[i]^2
  
  p.binomfl$mean[i] <- invlogit(mean(A))
  p.binomfl$lowci[i] <- invlogit(HPDinterval(mcmc(A))[1])
  p.binomfl$upci[i] <- invlogit(HPDinterval(mcmc(A))[2])
}

#general mismatch
ggplot(data=p.binomfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  geom_line(data=p.binomfl, aes(async, mean), lwd=1)+
  geom_ribbon(data=p.binomfl, aes(x=async, ymin=lowci, ymax=upci), alpha=0.15)+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6"x7"

#peak height
for(i in 1:length(p.binomfl$async)){
  
  A <- binomfledg2$Sol[,1]+binomfledg2$Sol[,2]*p.binomfl$async[i]+binomfledg2$Sol[,3]*p.binomfl$async[i]^2 + 2*binomfledg2$Sol[,4]
  B <- binomfledg2$Sol[,1]+binomfledg2$Sol[,2]*p.binomfl$async[i]+binomfledg2$Sol[,3]*p.binomfl$async[i]^2 - 2*binomfledg2$Sol[,4]
  
  p.binomfl$himean[i] <- invlogit(mean(A))
  p.binomfl$hilowci[i] <- invlogit(HPDinterval(mcmc(A))[1])
  p.binomfl$hiupci[i] <- invlogit(HPDinterval(mcmc(A))[2])
  p.binomfl$lomean[i] <- invlogit(mean(B))
  p.binomfl$lolowci[i] <- invlogit(HPDinterval(mcmc(B))[1])
  p.binomfl$loupci[i] <- invlogit(HPDinterval(mcmc(B))[2])
}

ggplot(data=p.binomfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  #geom_line(data=p.binomfl, aes(async, mean), lwd=1)+
  geom_line(data=p.binomfl, aes(async, himean), lwd=1, col=2)+
  geom_ribbon(data=p.binomfl, aes(x=async, ymin=hilowci, ymax=hiupci), alpha=0.15, fill="red")+
  geom_line(data=p.binomfl, aes(async, lomean), lwd=1, col="blue")+
  geom_ribbon(data=p.binomfl, aes(x=async, ymin=lolowci, ymax=loupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7

#peak width
for(i in 1:length(p.binomfl$async)){
  
  A <- binomfledg2$Sol[,1]+binomfledg2$Sol[,2]*p.binomfl$async[i]+binomfledg2$Sol[,3]*p.binomfl$async[i]^2 + 1.8*binomfledg2$Sol[,5]
  B <- binomfledg2$Sol[,1]+binomfledg2$Sol[,2]*p.binomfl$async[i]+binomfledg2$Sol[,3]*p.binomfl$async[i]^2 - 1.8*binomfledg2$Sol[,5]
  
  p.binomfl$widemean[i] <- invlogit(mean(A))
  p.binomfl$widelowci[i] <- invlogit(HPDinterval(mcmc(A))[1])
  p.binomfl$wideupci[i] <- invlogit(HPDinterval(mcmc(A))[2])
  p.binomfl$slimmean[i] <- invlogit(mean(B))
  p.binomfl$slimlowci[i] <- invlogit(HPDinterval(mcmc(B))[1])
  p.binomfl$slimupci[i] <- invlogit(HPDinterval(mcmc(B))[2])
}

ggplot(data=p.binomfl, aes(async, mean))+
  coord_cartesian(xlim = c(min(SYB$async), max(SYB$async)), ylim = c(min(SYB$propfl), max(SYB$propfl)))+
  geom_point(data=SYB, aes(async, propfl), col=mycolblack, pch=20)+
  #geom_line(data=p.binomfl, aes(async, mean), lwd=1)+
  geom_line(data=p.binomfl, aes(async, widemean), lwd=1, col=2)+
  geom_ribbon(data=p.binomfl, aes(x=async, ymin=widelowci, ymax=wideupci), alpha=0.15, fill="red")+
  geom_line(data=p.binomfl, aes(async, slimmean), lwd=1, col="blue")+
  geom_ribbon(data=p.binomfl, aes(x=async, ymin=slimlowci, ymax=slimupci), alpha=0.15, fill="blue")+
  xlab("Asynchrony (days)")+
  ylab("Prop. of chicks fledged")+
  theme_bw()+
  theme(text=element_text(size= 15)) #save as 6x7


#### Site-year caterpillar peaks ####

load("~/Documents/Models/SiteYear2.RData")

keepend <- function(x, n){  #function to keep last n number of characters in a string
  substr(x, nchar(x)-n+1, nchar(x))
}

#Dataframe with siteyears
modelcols <- data.frame(siteyear=keepend(colnames(SiteYear2$Sol)[669:924],8))
store <- pmatch(modelcols$siteyear,cater$siteyear)
modelcols$site <- cater$site[store]
modelcols$year <- cater$year[store]

syeffects <- SiteYear2$Sol[,1:924] #reduced columns

#Data frames for combined coefficients for each element or site year peaks
C <- data.frame(matrix(ncol = 256, nrow = length(syeffects[,1])))
colnames(C) <- keepend(colnames(syeffects)[669:924],8)
B <- data.frame(matrix(ncol = 256, nrow = length(syeffects[,1])))
colnames(B) <- keepend(colnames(syeffects)[669:924],8)
A <- data.frame(matrix(ncol = 256, nrow = length(syeffects[,1])))
colnames(A) <- keepend(colnames(syeffects)[669:924],8)

#columns showing which column number has each coefficient PD in sy effects
modelcols$C.I <- rep(1,length(modelcols$siteyear))
modelcols$B.I <- rep(2,length(modelcols$siteyear))
modelcols$A.I <- rep(3,length(modelcols$siteyear))

for(i in 1:length(modelcols$siteyear)){
  modelcols$C.S[i] <- which(keepend(colnames(syeffects),3)==modelcols$site[i])[1]
  modelcols$B.S[i] <- which(keepend(colnames(syeffects),3)==modelcols$site[i])[2]
  modelcols$A.S[i] <- which(keepend(colnames(syeffects),3)==modelcols$site[i])[3]
  modelcols$C.Y[i] <- which(keepend(colnames(syeffects),4)==modelcols$year[i] & substr(colnames(syeffects), nchar(colnames(syeffects))-8, nchar(colnames(syeffects))-5)=="year")[1]
  modelcols$B.Y[i] <- which(keepend(colnames(syeffects),4)==modelcols$year[i] & substr(colnames(syeffects), nchar(colnames(syeffects))-8, nchar(colnames(syeffects))-5)=="year")[2]
  modelcols$A.Y[i] <- which(keepend(colnames(syeffects),4)==modelcols$year[i] & substr(colnames(syeffects), nchar(colnames(syeffects))-8, nchar(colnames(syeffects))-5)=="year")[3]
  modelcols$C.SY[i] <- which(keepend(colnames(syeffects),8)==modelcols$siteyear[i])[1]
  modelcols$B.SY[i] <- which(keepend(colnames(syeffects),8)==modelcols$siteyear[i])[2]
  modelcols$A.SY[i] <- which(keepend(colnames(syeffects),8)==modelcols$siteyear[i])[3]
}

#combining relevant columns for each part of peak equation
for(i in 1:length(modelcols$siteyear)){
  C[,i] <- syeffects[,modelcols[i,"C.I"]]+syeffects[,modelcols[i,"C.S"]]+syeffects[,modelcols[i,"C.Y"]]+syeffects[,modelcols[i,"C.SY"]]
  B[,i] <- syeffects[,modelcols[i,"B.I"]]+syeffects[,modelcols[i,"B.S"]]+syeffects[,modelcols[i,"B.Y"]]+syeffects[,modelcols[i,"B.SY"]]
  A[,i] <- syeffects[,modelcols[i,"A.I"]]+syeffects[,modelcols[i,"A.S"]]+syeffects[,modelcols[i,"A.Y"]]+syeffects[,modelcols[i,"A.SY"]]
}

#new dataframe for metrics to loop into
symeans <- data.frame(siteyear=modelcols$siteyear, site=modelcols$site, year=modelcols$year)

for(i in 1:length(symeans$siteyear)){
  symeans$A[i] <- mean(A[,i])
  symeans$B[i] <- mean(B[,i])
  symeans$C[i] <- mean(C[,i])
}

p.peaks <- data.frame(matrix(ncol = 256, nrow = length(seq(-2,2.5,0.05))))
colnames(p.peaks) <- keepend(colnames(syeffects)[669:924],8)
p.peaks$date.scal <- seq(-2,2.5,0.05)

for(i in 1:256){
  p.peaks[,i] <- exp(mean(A[,i])*p.peaks$date.scal^2 + mean(B[,i])*p.peaks$date.scal + mean(C[,i]))
}

plotpeak <- data.frame(s1=p.peaks[,6], s2=p.peaks[,7], s3=p.peaks[,10], s4=p.peaks[,11], s5=p.peaks[,13], s6=p.peaks[,14], s7=p.peaks[,16],s8=p.peaks[,40],s9=p.peaks[,80],s10=p.peaks[,84])
plotpeak$date <- p.peaks$date.scal*14+147
plotpeak <- gather(plotpeak, key="sy", value="abund", 1:10)

ggplot(data=plotpeak, aes(date, abund, colour=sy))+
  geom_line(lwd=1)+
  xlab("Ordinal date")+
  ylab("Caterpillar abundance")+
  scale_colour_brewer(palette="Paired")+
  guides(color = "none", linetype="none")+
  theme_bw()+
  theme(text=element_text(size= 15), axis.text.y=element_blank()) #save as 6x7
