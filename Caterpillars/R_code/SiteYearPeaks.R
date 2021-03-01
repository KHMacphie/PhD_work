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

#### Dataframe ####
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater$treeID <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$datescaled <- scale(cater$date) #unscale(x, center= 146.77, scale=14.04305)   now from -2.11991 to 2.01024
mean(cater$date) # 146.77
sd(cater$date) # 14.04305


#############################################################
#### Model: Caterpillar abundance peaks in each SiteYear ####
#############################################################

#k<-10000
#prior<-list(R=list(V=1,nu=0.002),
#             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#SiteYear<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + us(1+datescaled+I(datescaled^2)):siteyear + siteday + recorder, 
#                         family="poisson", data=cater, prior=prior, nitt=500000, burnin=50000, pr=TRUE, thin=20)
#save(SiteYear, file = "~/Documents/Models/SiteYear.RData")
#load("~/Documents/Models/SiteYear.RData")

#k<-10000
#prior<-list(R=list(V=1,nu=0.002),
#            G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                   G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                   G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#SiteYear2<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + us(1+datescaled+I(datescaled^2)):siteyear + siteday + treeID + recorder, 
#                         family="poisson", data=cater, prior=prior, nitt=750000, burnin=50000, pr=TRUE, thin=50)
#save(SiteYear2, file = "~/Documents/Models/SiteYear2.RData")
load("~/Documents/Models/SiteYear2.RData")

#### Checking model fits the data and converged ####
#plot(SiteYear2) #look at fixed effect and random term trace plots 
#SiteYear2.Sim<-simulate(SiteYear2,nsim=1000) #simulate 1000 times
#par(mfcol=c(1,1))
#hist(apply(SiteYear2.Sim,2,sum), breaks=10000) #histogram of simulation predictions for total abundance
#hist(apply(SiteYear2.Sim,2,sum), breaks=10000, xlim=c(0,100000))
#abline(v=sum(cater$caterpillars),col=2) # red line for observed value in data

#plot(seq(1,length(apply(SiteYear2.Sim,2,sum)),1), apply(SiteYear2.Sim,2,sum), ylim=c(0,1000000))
#abline(h=sum(cater$caterpillars),col=2) 

#propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
#hist(apply(SiteYear2.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
#abline(v=propzero(cater$caterpillars), col="red") # red line for observed proportion in data

# each site year = (intercept + site + year + siteyear)+(date + date:site + date:year + date:siteyear)+(date^2 + date^2:site + date^2:year + date^2:siteyear)

#intercept: SiteYear2$Sol[,1]
#site intercepts: SiteYear2$Sol[,4:47]
#year intercepts: SiteYear2$Sol[,136:142]
#siteyear intercepts: SiteYear2$Sol[,157:412]

#date: SiteYear2$Sol[,2]
#site dates: SiteYear2$Sol[,48:91]
#year dates: SiteYear2$Sol[,143:149]
#siteyear dates: SiteYear2$Sol[,413:668]

#date^2: SiteYear2$Sol[,3]
#site date^2s: SiteYear2$Sol[,92:135]
#year date^2s: SiteYear2$Sol[,150:156]
#siteyear date^2s: SiteYear2$Sol[,669:924]

###################
#### Functions ####
###################
keepend <- function(x, n){  #function to keep last n number of characters in a string
  substr(x, nchar(x)-n+1, nchar(x))
}

forestplot <- function(d, x, y, l, u, xlab, ylab){  #function for forest plot: dataframe then column numbers
  xaxis <- d[,x]
  yaxis <- d[,y]
  lci <- d[,l]
  uci <- d[,u]
  
  ggplot(d, aes(fct_rev(fct_inorder(xaxis)), yaxis))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymax=uci, ymin=lci, width=0.5))+
    theme_bw()+
    xlab(xlab)+
    ylab(ylab)+
    #geom_hline(yintercept=1, linetype="dashed", color = "black")+
    theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
    guides(color = "none")+
    coord_flip()
}
# d=dataframe
# x=xaxis (column number in d)
# y=yaxis (column number in d)
# l=lower CI (column number in d) same as y if no CIs
# u=upper CI (column number in d) same as y if no CIs
# xlab/ylab axis title

###################

#Dataframe with siteyears
modelcols <- data.frame(siteyear=keepend(colnames(SiteYear2$Sol)[669:924],8))
store <- pmatch(modelcols$siteyear,cater$siteyear)
modelcols$site <- cater$site[store]
modelcols$year <- cater$year[store]

syeffects <- SiteYear2$Sol[,1:924] #reduced columns

## testing code
#keepend(colnames(SiteYear2$Sol)[885],8) 
#which(keepend(colnames(syeffects),3)==modelcols$site[5]) 
#substr(colnames(syeffects), nchar(colnames(syeffects))-8, nchar(colnames(syeffects))-5) 

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
symetrics <- data.frame(siteyear=modelcols$siteyear, site=modelcols$site, year=modelcols$year)
 i <- 221                         
#### Peak date per siteyear ####  x= -b/2a
#### Peak height per siteyear #### y and peak x
#### Peak width per siteyear (half each PH) ####    (-b +/- sqrt(b^2 - 4ac))/2a
for(i in 1:length(symetrics$siteyear)){
  PD <- (-B[,i]/(2*A[,i]))
  symetrics$PD[i] <- unscale(median(PD, na.rm=TRUE), scale=14.04305, center= 146.77)$V1
  symetrics$PDvar[i] <- unscale(var(PD, na.rm=TRUE), scale=14.04305)$V1 
  symetrics$PDna[i] <- (length(which(is.na(PD)==TRUE)))/length(PD) #proportion NA
  PH <- (C[,i] +B[,i]*(-B[,i]/(2*A[,i])) + A[,i]*(-B[,i]/(2*A[,i]))^2)
  symetrics$PH[i] <- exp(median(PH, na.rm=TRUE))
  symetrics$PlH[i] <- median(PH, na.rm=TRUE)
  symetrics$PlHvar[i] <- var(PH, na.rm=TRUE)
  symetrics$PHna[i] <- (length(which(is.na(PH)==TRUE)))/length(PH) #proportion NA
  PS <- sqrt(-1/(2*A[,i])) 
  symetrics$PS[i] <- unscale(median(PS, na.rm=TRUE), scale=14.04305)$V1
  symetrics$PSvar[i] <- unscale(var(PS, na.rm=TRUE), scale=14.04305)$V1
  symetrics$PSna[i] <- (length(which(is.na(PS)==TRUE)))/length(PS) #proportion NA
}

hist((-B[,i]/(2*A[,i])), breaks=100)
abline(v=median(-B[,i]/(2*A[,i])), col="darkgreen")
abline(v=median(-B[,i]/(2*A[,i]))+sd(-B[,i]/(2*A[,i])), col="red")
abline(v=median(-B[,i]/(2*A[,i]))-sd(-B[,i]/(2*A[,i])), col="red")

unscalPD <- c()
for(j in 1:length(PD)){
  unscalPD[j] <- unscale((-B[,i]/(2*A[,i]))[j], scale=14.04305, center= 146.77)$V1
}

hist(unscalPD, breaks=100)
abline(v=unscale(median(PD, na.rm=TRUE), scale=14.04305, center= 146.77)$V1, col="darkgreen")
abline(v=unscale(median(PD, na.rm=TRUE), scale=14.04305, center= 146.77)$V1+unscale(sd(PD, na.rm=TRUE), scale=14.04305)$V1 , col="red")
abline(v=unscale(median(PD, na.rm=TRUE), scale=14.04305, center= 146.77)$V1-unscale(sd(PD, na.rm=TRUE), scale=14.04305)$V1 , col="red")

### Plotting each SY peak estimate ###
unscale(2.01024, center= 146.77, scale=14.04305)  # now from -2.11991 to 2.01024
datescal <- seq(-2.12, 2.01, 0.01)
date <- unscale(datescal, center= 146.77, scale=14.04305)$V1
#i <- 190
for (i in 1:length(symetrics$siteyear)){
  X <- subset(cater, siteyear==symetrics$siteyear[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  abund <- exp(median(A[,i])*datescal^2+median(B[,i])*datescal+median(C[,i]))
  
  mypath <- file.path(paste("~/Documents/Models/MCMCglmmSYpeaks/",symetrics$siteyear[i],".pdf"))
  pdf(file = mypath, width=7, height=6)
  
  plot(jitter(X$date,1), jitter(X$cater,0.05), ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund")
  points(date, abund, type="l")
  title(main=symetrics$siteyear[i])
  
  dev.off()
  
}

#### Checking getting gussian coeffs from quadratic model ####
datescal <- seq(-2.12, 2.01, 0.01)
date <- unscale(datescal, center= 146.77, scale=14.04305)$V1

a <- median(A[,204])
b <- median(B[,204])
c <- median(C[,204])

mu <- -b/(2*a)
var <- -1/(2*a)
height <- exp(c+((mu^2)/(2*var)))
height2 <- exp(c-((a*(b^2))/(4*(a^2))))

gausabund <- exp((-((datescal-mu)^2)/(2*var))+log(height))
abund <- exp(a*datescal^2+b*datescal+c)
plot(date, abund, type="l")
points(date,gausabund, type="l", col=2)

######################################################
#### Optim to predict site-year caterpillar peaks ####
######################################################

#cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
#cater$siteyear <- paste(cater$site, cater$year)

gaussianpeak<-function(par,data){
  date <- data[,1]
  cater <- data[,2]
  
  mu <- par[1]
  sigma <- par[2]
  logheight <- par[3]
  
  return(-sum(dpois(cater, exp((-((date-mu)^2)/(2*sigma^2)))*exp(logheight), log=T)))
}

#### Estimates for each site-year ####
SYs <- colnames(A)
SYpeaks <- data.frame(SY=colnames(A))
datescal <- seq(-2.12, 2.01, 0.01)
date <- unscale(datescal, center= 146.77, scale=14.04305)$V1
par(mfrow=c(1,1))
#i <- 221

for(i in 1:length(SYs)){
  
  X <- subset(cater, siteyear==SYs[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  peak.X<-optim(par=c(mean(X$date),10,mean(X$cater)), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks$n[i] <- length(X$date)
  SYpeaks$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks$cater[i] <- sum(X$cater)
  SYpeaks$mu[i] <- peak.X$par[1]
  SYpeaks$mu.se[i] <- sqrt(diag(solve(peak.X$hessian)))[1]
  SYpeaks$sigma[i] <- peak.X$par[2]
  SYpeaks$sigma.se[i] <- sqrt(diag(solve(peak.X$hessian)))[2]
  SYpeaks$logheight[i] <- peak.X$par[3]
  SYpeaks$logheight.se[i] <- sqrt(diag(solve(peak.X$hessian)))[3]
  SYpeaks$loglik[i] <- peak.X$value
  SYpeaks$converge[i] <- peak.X$convergence
 
  optimabund <- exp((-((date-peak.X$par[1])^2)/(2*peak.X$par[2]^2)))*exp(peak.X$par[3])
  mcmcabund <- exp(median(A[,i])*datescal^2+median(B[,i])*datescal+median(C[,i]))
  
  mypath <- file.path(paste("~/Documents/Models/SYpeaks/",symetrics$siteyear[i],".pdf"))
  pdf(file = mypath, width=8, height=6)
  
  plot(jitter(X$date,1), jitter(X$cater,0.05), xlim=c(117,175),ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund")
  points(date, optimabund,type="l", col=2)
  points(date, mcmcabund, type="l")
  text(x=116+(nchar("MCMCglmm:")/2),y=max(X$cater), "MCMCglmm:", cex=0.8)
  text(x=116+(nchar(paste("height =", round(symetrics$PH[i], 2), "(", round(symetrics$PHLCI[i],2), "-", round(symetrics$PHUCI[i],2), ")"))/3),y=max(X$cater)-0.05*max(X$cater), paste("height =", round(symetrics$PH[i], 2), "(", round(symetrics$PHLCI[i],2), "-", round(symetrics$PHUCI[i],2), ")"), cex=0.8)
  text(x=116+(nchar(paste("date =", round(symetrics$PD[i], 2), "(", round(symetrics$PDLCI[i],2), "-", round(symetrics$PDUCI[i],2), ")"))/3),y=max(X$cater)-0.1*max(X$cater), paste("date =", round(symetrics$PD[i], 2), "(", round(symetrics$PDLCI[i],2), "-", round(symetrics$PDUCI[i],2), ")"), cex=0.8)
  text(x=116+(nchar(paste("width =", round(symetrics$PW[i], 2), "(", round(symetrics$PWLCI[i],2), "-", round(symetrics$PWUCI[i],2), ")"))/3),y=max(X$cater)-0.15*max(X$cater), paste("width =", round(symetrics$PW[i], 2), "(", round(symetrics$PWLCI[i],2), "-", round(symetrics$PWUCI[i],2), ")"), cex=0.8)
  text(x=116+(nchar(paste("NAs =", round(symetrics$PWna[i], 2)))/2),y=max(X$cater)-0.2*max(X$cater), paste("NAs =", round(symetrics$PWna[i], 2)), cex=0.8)
  text(x=116+(nchar("Optim:")/2),y=max(X$cater)-0.27*max(X$cater), "Optim:", col=2, cex=0.8)
  text(x=116+(nchar(paste("exp(logheight) =", round(exp(peak.X$par[3]), 2), "exp(se) =", round(exp(sqrt(diag(solve(peak.X$hessian)))[3]),2)))/3),y=max(X$cater)-0.32*max(X$cater), paste("exp(logheight) =", round(exp(peak.X$par[3]), 2), "exp(se) =", round(exp(sqrt(diag(solve(peak.X$hessian)))[3]),2)), col=2, cex=0.8)
  text(x=116+(nchar(paste("mu =", round(peak.X$par[1], 2), "se =", round(sqrt(diag(solve(peak.X$hessian)))[1],2)))/2.5),y=max(X$cater)-0.37*max(X$cater), paste("mu =", round(peak.X$par[1], 2), "se =", round(sqrt(diag(solve(peak.X$hessian)))[1],2)), col=2, cex=0.8)
  text(x=116+(nchar(paste("sigma =", round(peak.X$par[2], 2), "se =", round(sqrt(diag(solve(peak.X$hessian)))[2],2)))/2.5),y=max(X$cater)-0.42*max(X$cater), paste("sigma =", round(peak.X$par[2], 2), "se =", round(sqrt(diag(solve(peak.X$hessian)))[2],2)), col=2, cex=0.8)
  text(x=116+(nchar(paste("converg =", peak.X$convergence))/2),y=max(X$cater)-0.47*max(X$cater), paste("converg =", peak.X$convergence), col=2, cex=0.8)
  
  title(main=symetrics$siteyear[i])
  
  dev.off()
  
  X <- NULL
  peak.X <- NULL
}

#### MCMCglmm predictions####


# PD plots
symetrics$PDwidth <- symetrics$PDUCI-symetrics$PDLCI
symetrics <- symetrics[order(symetrics$PDwidth),] 
hist(symetrics$PD, breaks=200)
forestplot(symetrics, 1, 4, 5, 6, "Site-Year", "Peak Date") 
forestplot(symetrics, x=1, y=14, l=14, u=14, xlab="Site-Year", ylab="PD credible interval width")

# PH plots
symetrics$PHwidth <- symetrics$PHUCI-symetrics$PHLCI
symetrics <- symetrics[order(symetrics$PHwidth),] 
hist(symetrics$PH, breaks=300)
forestplot(symetrics, 1, 7, 8, 9, "Site-Year", "Peak Height") 
forestplot(symetrics, 1, 7, 8, 9, "Site-Year", "Peak Height")+coord_flip(ylim = c(0, 5))
forestplot(symetrics, x=1, y=15, l=15, u=15, xlab="Site-Year", ylab="PH credible interval width")
forestplot(symetrics, x=1, y=15, l=15, u=15, xlab="Site-Year", ylab="PH credible interval width")+coord_flip(ylim = c(0, 5))

hasheight <- subset(symetrics, PHLCI!=0)
forestplot(hasheight, 1, 7, 8, 9, "Site-Year", "Peak Height")+coord_flip(ylim = c(0, 1))


# PW plots
hist(symetrics$PWna, breaks=100)
symetrics$PWwidth <- symetrics$PWUCI-symetrics$PWLCI
symetrics <- symetrics[order(symetrics$PWwidth),] 
hist(symetrics$PW, breaks=300)
forestplot(symetrics, 1, 10, 11, 12, "Site-Year", "Peak Width (days)") 
forestplot(symetrics, x=1, y=16, l=16, u=16, xlab="Site-Year", ylab="PW credible interval width")
forestplot(symetrics, x=1, y=16, l=16, u=16, xlab="Site-Year", ylab="PW credible interval width")+coord_flip(ylim = c(0, 50))








##########################
#### Area under curve ####  
##########################

f<-function(x, i, b1, b2){exp(i+b1*x+b2*x^2)} #calc area under curve
area <- data.frame(matrix(ncol = 256, nrow = length(syeffects[,1]))) #dataframe for post dists of SY areas
colnames(area) <- keepend(colnames(syeffects)[669:924],8) #name columns
Apos <- c()
for(x in 1:length(colnames(A))){
  Apos[x] <- length(which(A[,x]>0)) 
} #number of positive date^2 coefficients in each siteyear posterior

A[A > 0] <- NA #any positive date^2 iterations to NA

for(j in 1:length(symetrics$siteyear)){
  for(i in 1:length(A[,1])){ 
    area[i,j] <- ifelse(
                    is.numeric(try(ifelse(A[i,j]=="NA", 
                                          yes= NA, 
                                          no=integrate(f, lower=-Inf, upper=Inf, i=C[i,j], b1=B[i,j], b2=A[i,j])$value), silent=T))==F,
                    yes= NA,
                    no= integrate(f, lower=-Inf, upper=Inf, i=C[i,j], b1=B[i,j], b2=A[i,j])$value)
}} # calculating areas and putting NA in any that would have an error

areaNA <- c()
for(x in 1:length(colnames(A))){
  areaNA[x] <- length(which(is.na(area[,x]))) 
} # number of NAs in each siteyear area post dist

NAs <- data.frame(Apos =Apos, areaNA=areaNA)
NAs$dif <- NAs[,2]-NAs[,1] #checking for extreme differences.. none

#p.datescal <- seq(-2.12,2.01,0.01)
#p.date <- unscale(p.datescal, center= 146.77, scale=14.04305)$V1
#plot(p.date, exp(mean(A[18081,11])*(p.datescal^2)+mean(B[18081,11])*p.datescal+mean(C[18081,11]))) #very weird extreme broad peak

for(y in 1:length(symetrics$siteyear)){
  symetrics$PA[y] <- median(area[,y], na.rm=TRUE)
  symetrics$PALCI[y] <- HPDinterval(mcmc(area[,y]))[1]
  symetrics$PAUCI[y] <- HPDinterval(mcmc(area[,y]))[2]
  symetrics$PAna[y] <- length(which(is.na(area[,y])))/length(area[,1]) #proportion NA
} # extracting median and CIs

# PA plots
hist(symetrics$PAna, breaks=100)
symetrics$PAwidth <- symetrics$PAUCI-symetrics$PALCI
symetrics <- symetrics[order(symetrics$PAwidth),] 
hist(symetrics$PA, breaks=300)
forestplot(symetrics, 1, 14, 15, 16, "Site-Year", "Peak Area") 
forestplot(symetrics, 1, 14, 15, 16, "Site-Year", "Peak Area")+coord_flip(ylim = c(0, 1))

forestplot(symetrics, x=1, y=18, l=18, u=18, xlab="Site-Year", ylab="PA credible interval width")
forestplot(symetrics, x=1, y=18, l=18, u=18, xlab="Site-Year", ylab="PA credible interval width")+coord_flip(ylim = c(0, 3))

############################
#### Parameters to test ####
############################

## Responses
# nest average mass day 12
# average mass increase day 12-6
# number fledged
# proportion fledged (+ clutchsize)
# female mass
# male mass

# reltiming = hatch+10 - caterpillar peak
# abreltiming = | reltiming |
# height
# width
# peakarea

# response ~ abreltiming*height + abreltimign*width
# response ~ abreltiming*peakarea

## Random
# residual = nest
# site
# year
# siteyear

# popdensity = mean centred proportion of boxes occupied (have eggs?)
# motherID
# fatherID
# number of chicks fledged (factor)
