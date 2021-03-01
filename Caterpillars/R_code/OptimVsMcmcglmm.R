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

###########################
#### Optim vs MCMCglmm ####
###########################

#### MCMCglmm model####
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater$treeID <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$datescaled <- scale(cater$date) #unscale(x, center= 146.77, scale=14.04305)   now from -2.11991 to 2.01024
mean(cater$date) # 146.77
sd(cater$date) # 14.04305

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
symetrics <- data.frame(siteyear=modelcols$siteyear, site=modelcols$site, year=modelcols$year)

for(i in 1:length(symetrics$siteyear)){
  PD <- (-B[,i]/(2*A[,i]))
  symetrics$PD[i] <- unscale(median(PD, na.rm=TRUE), scale=14.04305, center= 146.77)$V1
  symetrics$PDvar[i] <- unscale(var(PD, na.rm=TRUE), scale=14.04305)$V1 
  symetrics$PDna[i] <- (length(which(is.na(PD)==TRUE)))/length(PD) #proportion NA
  PH <- (C[,i] +B[,i]*(-B[,i]/(2*A[,i])) + A[,i]*(-B[,i]/(2*A[,i]))^2)
  #symetrics$PH[i] <- exp(median(PH, na.rm=TRUE))
  symetrics$PlH[i] <- median(PH, na.rm=TRUE)
  symetrics$PlHvar[i] <- var(PH, na.rm=TRUE)
  symetrics$PHna[i] <- (length(which(is.na(PH)==TRUE)))/length(PH) #proportion NA
  PS <- sqrt(-1/(2*A[,i])) 
  symetrics$PS[i] <- unscale(median(PS, na.rm=TRUE), scale=14.04305)$V1
  symetrics$PSvar[i] <- unscale(var(PS, na.rm=TRUE), scale=14.04305)$V1
  symetrics$PSna[i] <- (length(which(is.na(PS)==TRUE)))/length(PS) #proportion NA
}

write.csv(symetrics, "~/Documents/Models/symetrics.csv", row.names = FALSE) 

#### Optim ####

gaussianpeak<-function(par,data){
  date <- data[,1]
  cater <- data[,2]
  
  mu <- par[1]
  sigma <- par[2]
  logheight <- par[3]
  
  return(-sum(dpois(cater, exp((-((date-mu)^2)/(2*sigma^2)))*exp(logheight), log=T)))
}

# Estimates for each site-year choosing highest liklihood option 
SYs <- data.frame(SY=unique(cater$siteyear))
for(j in 1:length(SYs$SY)){
  SYs$cater[j] <-  sum(cater$caterpillars[which(cater$siteyear==SYs$SY[j])])
}
SYs <- SYs$SY[SYs$cater>0]

SYpeaks <- data.frame(siteyear=SYs) # with highest likelihood
#peak.Y= mean,10,mean  
#peak.Z= mean,10,1     
#peak.W= mean,10,0     
#peak.F= mean,5,mean   
#peak.G= mean,5,1      
#peak.H= mean,5,0      
#peak.I= mean,12,mean  
#peak.J= mean,12,1     
#peak.K= mean,12,0     
#peak.L= mean,2,mean   
#peak.M= mean,2,1      
#peak.N= mean,2,0      

for(i in 1:length(SYs)){
  
  X <- subset(cater, siteyear==SYs[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  
  peak.Y<-optim(par=c(mean(X$date),10,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.Z<-optim(par=c(mean(X$date),10,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.W<-optim(par=c(mean(X$date),10,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.F<-optim(par=c(mean(X$date),5,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.G<-optim(par=c(mean(X$date),5,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.H<-optim(par=c(mean(X$date),5,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.I<-optim(par=c(mean(X$date),12,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.J<-optim(par=c(mean(X$date),12,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.K<-optim(par=c(mean(X$date),12,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.L<-optim(par=c(mean(X$date),2,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.M<-optim(par=c(mean(X$date),2,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  peak.N<-optim(par=c(mean(X$date),2,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  df <- data.frame(round=c("peak.Y","peak.Z","peak.W","peak.F","peak.G","peak.H","peak.I","peak.J","peak.K","peak.L","peak.M","peak.N"),
                   loglik=c(peak.Y$value,peak.Z$value,peak.W$value,peak.F$value,peak.G$value,peak.H$value,peak.I$value,peak.J$value,peak.K$value,peak.L$value,peak.M$value,peak.N$value))
  
  for(j in 1:12){
    df$NAorNCon[j] <- as.numeric(tryCatch(ifelse(is.na(sqrt(diag(solve(get(paste(df$round[j]))$hessian)))[1] || sqrt(diag(solve(get(paste(df$round[j]))$hessian)))[2] || sqrt(diag(solve(get(paste(df$round[j]))$hessian)))[3])==T ||
                                                   (get(paste(df$round[j]))$convergence==1)==T, 1, 0), error=function(e){print("1")}))
  }
  
  
  SYpeaks$n[i] <- length(X$date)
  SYpeaks$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks$cater[i] <- sum(X$cater)
  SYpeaks$mu[i] <- get(paste(df$round[which.min(df$loglik)]))$par[1]
  SYpeaks$mu.var[i] <- tryCatch(diag(solve(get(paste(df$round[which.min(df$loglik)]))$hessian))[1], error=function(e){print("Error")})
  SYpeaks$sigma[i] <- get(paste(df$round[which.min(df$loglik)]))$par[2]
  SYpeaks$sigma.var[i] <- tryCatch(diag(solve(get(paste(df$round[which.min(df$loglik)]))$hessian))[2], error=function(e){print("Error")})
  SYpeaks$logheight[i] <- get(paste(df$round[which.min(df$loglik)]))$par[3]
  SYpeaks$logheight.var[i] <- tryCatch(diag(solve(get(paste(df$round[which.min(df$loglik)]))$hessian))[3], error=function(e){print("Error")})
  SYpeaks$loglik[i] <- get(paste(df$round[which.min(df$loglik)]))$value
  SYpeaks$converge[i] <- get(paste(df$round[which.min(df$loglik)]))$convergence  
  SYpeaks$which[i] <- df$round[which.min(df$loglik)] 
  SYpeaks$NAorNConv[i] <- df$NAorNCon[which.min(df$loglik)]==0 ## true= converged/noNAs, false = didnt converge/hasNAs
  SYpeaks$OtherOK[i] <- sum(df$NAorNCon)<12 ## true= at least one round converged/noNAs, false = all had a problem
  
  peak.Y<-NULL
  peak.Z<-NULL
  peak.W<-NULL
  peak.F<-NULL
  peak.G<-NULL
  peak.H<-NULL
  peak.I<-NULL
  peak.J<-NULL
  peak.K<-NULL
  peak.L<-NULL
  peak.M<-NULL
  peak.N<-NULL
  df <- NULL
  
}


#### Combining the two sets of estimates ####

store <- pmatch(SYpeaks$siteyear, symetrics$siteyear)
both <- SYpeaks
both[,16:24] <- symetrics[store,4:12]
plot(both$mu, both$PD)
plot(both$mu, both$PD, xlim=c(100,200))
points(seq(100,200,1),seq(100,200,1),type="l")
plot(both$logheight, both$PlH)
plot(both$logheight, both$PlH, xlim=c(-5,2))
points(seq(-5,2,1),seq(-5,2,1),type="l")
plot(both$sigma, both$PS)
plot(both$sigma, both$PS, xlim=c(5,30))
points(seq(5,30,1),seq(5,30,1),type="l")

# Plotting peak estimates
datescal <- seq(-2.12, 2.01, 0.01)
date <- unscale(datescal, center= 146.77, scale=14.04305)$V1
par(mfrow=c(1,1))

for(i in 1:length(both$siteyear)){

  X <- subset(cater, siteyear==both$siteyear[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  
optimabund <- exp((-((date-both$mu[i])^2)/(2*both$sigma[i]^2)))*exp(both$logheight[i])
mcmcabund <- exp((-((date-both$PD[i])^2)/(2*both$PS[i]^2)))*exp(both$PlH[i])

mypath <- file.path(paste("~/Documents/Models/SYpeaks/",both$siteyear[i],".pdf"))
pdf(file = mypath, width=8, height=6)

plot(jitter(X$date,1), jitter(X$cater,0.05), xlim=c(117,175),ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund")
points(date, optimabund,type="l", col=2)
points(date, mcmcabund, type="l")

text(x=126,y=max(X$cater), "MCMCglmm:", cex=0.8)
text(x=126,y=max(X$cater)-0.04*max(X$cater), paste("mu =", round(both$PD[i], 2), "var =", round(both$PDvar[i], 4)), col="black", cex=0.8)
text(x=126,y=max(X$cater)-0.08*max(X$cater), paste("logheight =", round(both$PlH[i], 3), "var =", round(both$PlHvar[i], 4)), col="black", cex=0.8)
text(x=126,y=max(X$cater)-0.12*max(X$cater), paste("sigma =", round(both$PS[i], 3), "var =", round(both$PSvar[i], 4)), col="black", cex=0.8)
text(x=126,y=max(X$cater)-0.16*max(X$cater), paste("NAs? =", ifelse((both$PSna[i]>0)==TRUE, paste("yes ( prop =", round(both$PSna[i],5), ")"), "no")), col="black", cex=0.8)

text(x=126,y=max(X$cater)-0.24*max(X$cater), "Optim:", cex=0.8, col="red")
text(x=126,y=max(X$cater)-0.28*max(X$cater), paste("mu =", round(both$mu[i], 2), "var =", round(both$mu.var[i], 4)), col="red", cex=0.8)
text(x=126,y=max(X$cater)-0.32*max(X$cater), paste("logheight =", round(both$logheight[i], 3), "var =", round(both$logheight.var[i], 4)), col="red", cex=0.8)
text(x=126,y=max(X$cater)-0.36*max(X$cater), paste("sigma =", round(both$sigma[i], 3), "var =", round(both$sigma.var[i], 4)), col="red", cex=0.8)
text(x=126,y=max(X$cater)-0.4*max(X$cater), paste("NotConv/NAs? =", ifelse((both$NAorNConv[i]>0)==TRUE, "no", "yes")), col="red", cex=0.8)

title(main=both$siteyear[i])

dev.off()

}
