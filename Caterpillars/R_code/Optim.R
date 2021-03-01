rm(list=ls())
setwd('/Users/s1205615/')

######################################################
#### Optim to predict site-year caterpillar peaks ####
######################################################

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$sy <- paste(cater$site, cater$year)

gaussianpeak<-function(par,data){
  date <- data[,1]
  cater <- data[,2]
  
  mu <- par[1]
  sigma <- par[2]
  logheight <- par[3]
   
  return(-sum(dpois(cater, exp((-((date-mu)^2)/(2*sigma^2)))*exp(logheight), log=T)))
}

#### Testing examples ####
pit19 <- subset(cater, sy=="PIT 2019")
pit19 <- data.frame(date=pit19$date, cater=pit19$caterpillars)
pth19 <- subset(cater, sy=="PTH 2019")
pth19 <- data.frame(date=pth19$date, cater=pth19$caterpillars)
plot(pit19$date, pit19$cater, xlim=c(118,170))
plot(pth19$date, pth19$cater, xlim=c(118,170))
date <- seq(118,170,0.1)

peak.pit19<-optim(par=c(mean(pit19$date),10,log(1)), fn=gaussianpeak,data=pit19,method="Nelder-Mead",hessian=T)
peak.pit19
plot(date, exp((-((date-peak.pit19$par[1])^2)/(2*peak.pit19$par[2]^2)))*exp(peak.pit19$par[3]),type="l", ylim=c(0,max(pit19$cater)))
points(jitter(pit19$date,1), jitter(pit19$cater,1)) # jittered so you can see density of points
sqrt(diag(solve(peak.pit19$hessian)))

## gaussian curve as standard quadratic equation
#points(date, exp((-1/(2*peak.pit19$par[2]^2))*date^2 + (peak.pit19$par[1]/peak.pit19$par[2]^2)*date - (peak.pit19$par[1]^2/(2*peak.pit19$par[2]^2)) + peak.pit19$par[3]), type="l", col=2)

peak.pit19<-optim(par=c(mean(pit19$date),10,mean(pit19$cater)), fn=gaussianpeak,data=pit19,method="Nelder-Mead",hessian=T)
peak.pit19
plot(date, exp((-((date-peak.pit19$par[1])^2)/(2*peak.pit19$par[2]^2)))*exp(peak.pit19$par[3]),type="l", ylim=c(0,max(pit19$cater)))
points(jitter(pit19$date,1), jitter(pit19$cater,1)) # jittered so you can see density of points
sqrt(diag(solve(peak.pit19$hessian)))

peak.pth19<-optim(par=c(mean(pth19$date),10,1), fn=gaussianpeak,data=pth19,method="Nelder-Mead",hessian=T)
peak.pth19
plot(date, exp((-((date-peak.pth19$par[1])^2)/(2*peak.pth19$par[2]^2)))*exp(peak.pth19$par[3]),type="l", ylim=c(0,max(pth19$cater)))
points(jitter(pth19$date,1), jitter(pth19$cater,1))
sqrt(diag(solve(peak.pth19$hessian)))

peak.pth19<-optim(par=c(mean(pth19$date),10,mean(pth19$cater)), fn=gaussianpeak,data=pth19,method="Nelder-Mead",hessian=T)
peak.pth19
plot(date, exp((-((date-peak.pth19$par[1])^2)/(2*peak.pth19$par[2]^2)))*exp(peak.pth19$par[3]),type="l", ylim=c(0,max(pth19$cater)))
points(jitter(pth19$date,1), jitter(pth19$cater,1))
sqrt(diag(solve(peak.pth19$hessian)))

date <- seq(118,170,0.1)
i <- 21
X <- subset(cater, sy==SYs[i])
X <- data.frame(date=X$date, cater=X$caterpillars)
peak.X<-optim(par=c(mean(X$date),10,mean(X$cater)), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
peak.X
sqrt(diag(solve(peak.X$hessian)))
plot(date, exp((-((date-peak.X$par[1])^2)/(2*peak.X$par[2]^2)))*exp(peak.X$par[3]),type="l", ylim=c(0,max(X$cater)))
points(X$date, X$cater) # jittered so you can see density of points

#### Estimates for each site-year ####
SYs <- unique(cater$sy)
SYpeaks1 <- data.frame(SY=unique(cater$sy))
date <- seq(118,175,0.1)
par(mfrow=c(1,1))
 
for(i in 1:length(SYs)){
  
    X <- subset(cater, sy==SYs[i])
    X <- data.frame(date=X$date, cater=X$caterpillars)
    peak.X<-optim(par=c(mean(X$date),1,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
    
    SYpeaks1$n[i] <- length(X$date)
    SYpeaks1$propzero[i] <- length(which(X$cater==0))/length(X$cater)
    SYpeaks1$cater[i] <- sum(X$cater)
    SYpeaks1$mu[i] <- peak.X$par[1]
    SYpeaks1$mu.var[i] <- tryCatch(diag(solve(peak.X$hessian))[1], error=function(e){print("Error")})
    SYpeaks1$sigma[i] <- peak.X$par[2]
    SYpeaks1$sigma.var[i] <- tryCatch(diag(solve(peak.X$hessian))[2], error=function(e){print("Error")})
    SYpeaks1$logheight[i] <- peak.X$par[3]
    SYpeaks1$logheight.var[i] <- tryCatch(diag(solve(peak.X$hessian))[3], error=function(e){print("Error")})
    SYpeaks1$loglik[i] <- peak.X$value
    SYpeaks1$converge[i] <- peak.X$convergence
    
    mypath <- file.path(paste("~/Documents/Models/OptimSYpeaks/",SYs[i],".pdf"))
    
    pdf(file = mypath, width=7, height=6)
    
    plot(date, exp((-((date-peak.X$par[1])^2)/(2*peak.X$par[2]^2)))*exp(peak.X$par[3]),type="l", ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund")
    points(jitter(X$date,1), jitter(X$cater,0.05))
    text(x=126,y=max(X$cater), paste("mu =", round(peak.X$par[1], 2), "var =", ifelse(SYpeaks1$mu.var[i]=="Error", SYpeaks1$mu.var[i], round(as.numeric(SYpeaks1$mu.var[i]),2))))
    text(x=126,y=max(X$cater)-0.05*max(X$cater), paste("sigma =", round(peak.X$par[2], 2), "var =", ifelse(SYpeaks1$sigma.var[i]=="Error", SYpeaks1$sigma.var[i], round(as.numeric(SYpeaks1$sigma.var[i]),2))))
    text(x=127.5,y=max(X$cater)-0.1*max(X$cater), paste("logheight =", round(peak.X$par[3], 2), "var =", ifelse(SYpeaks1$logheight.var[i]=="Error", SYpeaks1$logheight.var[i], round(as.numeric(SYpeaks1$logheight.var[i]),2))))
    text(x=129,y=max(X$cater)-0.15*max(X$cater), paste("loglik =", round(peak.X$value, 2), "convergence =", peak.X$convergence))
    title(main=SYs[i])
    
    dev.off()
    
    X <- NULL
    peak.X <- NULL
}

#### Estimates for each site-year with varied start parameters ####
SYs <- data.frame(SY=unique(cater$sy))
for(j in 1:length(SYs$SY)){
 SYs$cater[j] <-  sum(cater$caterpillars[which(cater$sy==SYs$SY[j])])
}
SYs <- SYs$SY[SYs$cater>0]

SYpeaks1 <- data.frame(SY=SYs) #mean,10,mean Y black
SYpeaks2 <- data.frame(SY=SYs) #mean,10,1    Z gray
SYpeaks3 <- data.frame(SY=SYs) #mean,10,0    W red 
SYpeaks4 <- data.frame(SY=SYs) #mean,5,mean  F orange
SYpeaks5 <- data.frame(SY=SYs) #mean,5,1     G yellow
SYpeaks6 <- data.frame(SY=SYs) #mean,5,0     H magenta

date <- seq(118,175,0.1)
par(mfrow=c(1,1))
#i <- 15

for(i in 1:length(SYs)){
  
  X <- subset(cater, sy==SYs[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  
  peak.Y<-optim(par=c(mean(X$date),10,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks1$n[i] <- length(X$date)
  SYpeaks1$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks1$cater[i] <- sum(X$cater)
  SYpeaks1$mu[i] <- peak.Y$par[1]
  SYpeaks1$mu.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[1], error=function(e){print("Error")})
  SYpeaks1$sigma[i] <- peak.Y$par[2]
  SYpeaks1$sigma.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[2], error=function(e){print("Error")})
  SYpeaks1$logheight[i] <- peak.Y$par[3]
  SYpeaks1$logheight.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[3], error=function(e){print("Error")})
  SYpeaks1$loglik[i] <- peak.Y$value
  SYpeaks1$converge[i] <- peak.Y$convergence
  
  peak.Z<-optim(par=c(mean(X$date),10,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks2$n[i] <- length(X$date)
  SYpeaks2$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks2$cater[i] <- sum(X$cater)
  SYpeaks2$mu[i] <- peak.Z$par[1]
  SYpeaks2$mu.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[1], error=function(e){print("Error")})
  SYpeaks2$sigma[i] <- peak.Z$par[2]
  SYpeaks2$sigma.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[2], error=function(e){print("Error")})
  SYpeaks2$logheight[i] <- peak.Z$par[3]
  SYpeaks2$logheight.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[3], error=function(e){print("Error")})
  SYpeaks2$loglik[i] <- peak.Z$value
  SYpeaks2$converge[i] <- peak.Z$convergence
  
  peak.W<-optim(par=c(mean(X$date),10,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks2$n[i] <- length(X$date)
  SYpeaks3$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks3$cater[i] <- sum(X$cater)
  SYpeaks3$mu[i] <- peak.W$par[1]
  SYpeaks3$mu.var[i] <- tryCatch(diag(solve(peak.W$hessian))[1], error=function(e){print("Error")})
  SYpeaks3$sigma[i] <- peak.W$par[2]
  SYpeaks3$sigma.var[i] <- tryCatch(diag(solve(peak.W$hessian))[2], error=function(e){print("Error")})
  SYpeaks3$logheight[i] <- peak.W$par[3]
  SYpeaks3$logheight.var[i] <- tryCatch(diag(solve(peak.W$hessian))[3], error=function(e){print("Error")})
  SYpeaks3$loglik[i] <- peak.W$value
  SYpeaks3$converge[i] <- peak.W$convergence
  
  peak.F<-optim(par=c(mean(X$date),5,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks4$n[i] <- length(X$date)
  SYpeaks4$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks4$cater[i] <- sum(X$cater)
  SYpeaks4$mu[i] <- peak.F$par[1]
  SYpeaks4$mu.var[i] <- tryCatch(diag(solve(peak.F$hessian))[1], error=function(e){print("Error")})
  SYpeaks4$sigma[i] <- peak.F$par[2]
  SYpeaks4$sigma.var[i] <- tryCatch(diag(solve(peak.F$hessian))[2], error=function(e){print("Error")})
  SYpeaks4$logheight[i] <- peak.F$par[3]
  SYpeaks4$logheight.var[i] <- tryCatch(diag(solve(peak.F$hessian))[3], error=function(e){print("Error")})
  SYpeaks4$loglik[i] <- peak.F$value
  SYpeaks4$converge[i] <- peak.F$convergence
  
  peak.G<-optim(par=c(mean(X$date),5,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks5$n[i] <- length(X$date)
  SYpeaks5$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks5$cater[i] <- sum(X$cater)
  SYpeaks5$mu[i] <- peak.G$par[1]
  SYpeaks5$mu.var[i] <- tryCatch(diag(solve(peak.G$hessian))[1], error=function(e){print("Error")})
  SYpeaks5$sigma[i] <- peak.G$par[2]
  SYpeaks5$sigma.var[i] <- tryCatch(diag(solve(peak.G$hessian))[2], error=function(e){print("Error")})
  SYpeaks5$logheight[i] <- peak.G$par[3]
  SYpeaks5$logheight.var[i] <- tryCatch(diag(solve(peak.G$hessian))[3], error=function(e){print("Error")})
  SYpeaks5$loglik[i] <- peak.G$value
  SYpeaks5$converge[i] <- peak.G$convergence
  
  peak.H<-optim(par=c(mean(X$date),5,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks6$n[i] <- length(X$date)
  SYpeaks6$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks6$cater[i] <- sum(X$cater)
  SYpeaks6$mu[i] <- peak.H$par[1]
  SYpeaks6$mu.var[i] <- tryCatch(diag(solve(peak.H$hessian))[1], error=function(e){print("Error")})
  SYpeaks6$sigma[i] <- peak.H$par[2]
  SYpeaks6$sigma.var[i] <- tryCatch(diag(solve(peak.H$hessian))[2], error=function(e){print("Error")})
  SYpeaks6$logheight[i] <- peak.H$par[3]
  SYpeaks6$logheight.var[i] <- tryCatch(diag(solve(peak.H$hessian))[3], error=function(e){print("Error")})
  SYpeaks6$loglik[i] <- peak.H$value
  SYpeaks6$converge[i] <- peak.H$convergence
  
  mypath <- file.path(paste("~/Documents/Models/OptimSYpeaks/",SYs[i],".pdf"))
  
  pdf(file = mypath, width=10, height=20)
  
  plot(date, exp((-((date-peak.Y$par[1])^2)/(2*peak.Y$par[2]^2)))*exp(peak.Y$par[3]),type="l", ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund", col="black")
  points(jitter(X$date,1), jitter(X$cater,0.05))
  text(x=126,y=max(X$cater), paste("mu =", round(peak.Y$par[1], 2), "var =", ifelse(SYpeaks1$mu.var[i]=="Error", SYpeaks1$mu.var[i], round(as.numeric(SYpeaks1$mu.var[i]),2))), col="black")
  text(x=126,y=max(X$cater)-0.01*max(X$cater), paste("sigma =", round(peak.Y$par[2], 2), "var =", ifelse(SYpeaks1$sigma.var[i]=="Error", SYpeaks1$sigma.var[i], round(as.numeric(SYpeaks1$sigma.var[i]),2))), col="black")
  text(x=127.5,y=max(X$cater)-0.02*max(X$cater), paste("logheight =", round(peak.Y$par[3], 2), "var =", ifelse(SYpeaks1$logheight.var[i]=="Error", SYpeaks1$logheight.var[i], round(as.numeric(SYpeaks1$logheight.var[i]),2))), col="black")
  text(x=129,y=max(X$cater)-0.03*max(X$cater), paste("loglik =", round(peak.Y$value, 2), "convergence =", peak.Y$convergence), col="black")
  
  points(date, exp((-((date-peak.Z$par[1])^2)/(2*peak.Z$par[2]^2)))*exp(peak.Z$par[3]),type="l", col="gray")
  text(x=126,y=max(X$cater)-0.05*max(X$cater), paste("mu =", round(peak.Z$par[1], 2), "var =", ifelse(SYpeaks2$mu.var[i]=="Error", SYpeaks2$mu.var[i], round(as.numeric(SYpeaks2$mu.var[i]),2))), col="gray")
  text(x=126,y=max(X$cater)-0.06*max(X$cater), paste("sigma =", round(peak.Z$par[2], 2), "var =", ifelse(SYpeaks2$sigma.var[i]=="Error", SYpeaks2$sigma.var[i], round(as.numeric(SYpeaks2$sigma.var[i]),2))), col="gray")
  text(x=127.5,y=max(X$cater)-0.07*max(X$cater), paste("logheight =", round(peak.Z$par[3], 2), "var =", ifelse(SYpeaks2$logheight.var[i]=="Error", SYpeaks2$logheight.var[i], round(as.numeric(SYpeaks2$logheight.var[i]),2))), col="gray")
  text(x=129,y=max(X$cater)-0.08*max(X$cater), paste("loglik =", round(peak.Z$value, 2), "convergence =", peak.Z$convergence), col="gray")
  
  points(date, exp((-((date-peak.W$par[1])^2)/(2*peak.W$par[2]^2)))*exp(peak.W$par[3]),type="l", col="red")
  text(x=126,y=max(X$cater)-0.10*max(X$cater), paste("mu =", round(peak.W$par[1], 2), "var =", ifelse(SYpeaks3$mu.var[i]=="Error", SYpeaks3$mu.var[i], round(as.numeric(SYpeaks3$mu.var[i]),2))), col="red")
  text(x=126,y=max(X$cater)-0.11*max(X$cater), paste("sigma =", round(peak.W$par[2], 2), "var =", ifelse(SYpeaks3$sigma.var[i]=="Error", SYpeaks3$sigma.var[i], round(as.numeric(SYpeaks3$sigma.var[i]),2))), col="red")
  text(x=127.5,y=max(X$cater)-0.12*max(X$cater), paste("logheight =", round(peak.W$par[3], 2), "var =", ifelse(SYpeaks3$logheight.var[i]=="Error", SYpeaks3$logheight.var[i], round(as.numeric(SYpeaks3$logheight.var[i]),2))), col="red")
  text(x=129,y=max(X$cater)-0.13*max(X$cater), paste("loglik =", round(peak.W$value, 2), "convergence =", peak.W$convergence), col="red")
  
  points(date, exp((-((date-peak.F$par[1])^2)/(2*peak.F$par[2]^2)))*exp(peak.F$par[3]),type="l", col="orange")
  text(x=126,y=max(X$cater)-0.15*max(X$cater), paste("mu =", round(peak.F$par[1], 2), "var =", ifelse(SYpeaks4$mu.var[i]=="Error", SYpeaks4$mu.var[i], round(as.numeric(SYpeaks4$mu.var[i]),2))), col="orange")
  text(x=126,y=max(X$cater)-0.16*max(X$cater), paste("sigma =", round(peak.F$par[2], 2), "var =", ifelse(SYpeaks4$sigma.var[i]=="Error", SYpeaks4$sigma.var[i], round(as.numeric(SYpeaks4$sigma.var[i]),2))), col="orange")
  text(x=127.5,y=max(X$cater)-0.17*max(X$cater), paste("logheight =", round(peak.F$par[3], 2), "var =", ifelse(SYpeaks4$logheight.var[i]=="Error", SYpeaks4$logheight.var[i], round(as.numeric(SYpeaks4$logheight.var[i]),2))), col="orange")
  text(x=129,y=max(X$cater)-0.18*max(X$cater), paste("loglik =", round(peak.F$value, 2), "convergence =", peak.F$convergence), col="orange")
  
  points(date, exp((-((date-peak.G$par[1])^2)/(2*peak.G$par[2]^2)))*exp(peak.G$par[3]),type="l", col="yellow")
  text(x=126,y=max(X$cater)-0.20*max(X$cater), paste("mu =", round(peak.G$par[1], 2), "var =", ifelse(SYpeaks5$mu.var[i]=="Error", SYpeaks5$mu.var[i], round(as.numeric(SYpeaks5$mu.var[i]),2))), col="yellow")
  text(x=126,y=max(X$cater)-0.21*max(X$cater), paste("sigma =", round(peak.G$par[2], 2), "var =", ifelse(SYpeaks5$sigma.var[i]=="Error", SYpeaks5$sigma.var[i], round(as.numeric(SYpeaks5$sigma.var[i]),2))), col="yellow")
  text(x=127.5,y=max(X$cater)-0.22*max(X$cater), paste("logheight =", round(peak.G$par[3], 2), "var =", ifelse(SYpeaks5$logheight.var[i]=="Error", SYpeaks5$logheight.var[i], round(as.numeric(SYpeaks5$logheight.var[i]),2))), col="yellow")
  text(x=129,y=max(X$cater)-0.23*max(X$cater), paste("loglik =", round(peak.G$value, 2), "convergence =", peak.G$convergence), col="yellow")
  
  points(date, exp((-((date-peak.H$par[1])^2)/(2*peak.H$par[2]^2)))*exp(peak.H$par[3]),type="l", col="magenta")
  text(x=126,y=max(X$cater)-0.25*max(X$cater), paste("mu =", round(peak.H$par[1], 2), "var =", ifelse(SYpeaks6$mu.var[i]=="Error", SYpeaks6$mu.var[i], round(as.numeric(SYpeaks6$mu.var[i]),2))), col="magenta")
  text(x=126,y=max(X$cater)-0.26*max(X$cater), paste("sigma =", round(peak.H$par[2], 2), "var =", ifelse(SYpeaks6$sigma.var[i]=="Error", SYpeaks6$sigma.var[i], round(as.numeric(SYpeaks6$sigma.var[i]),2))), col="magenta")
  text(x=127.5,y=max(X$cater)-0.27*max(X$cater), paste("logheight =", round(peak.H$par[3], 2), "var =", ifelse(SYpeaks6$logheight.var[i]=="Error", SYpeaks6$logheight.var[i], round(as.numeric(SYpeaks6$logheight.var[i]),2))), col="magenta")
  text(x=129,y=max(X$cater)-0.28*max(X$cater), paste("loglik =", round(peak.H$value, 2), "convergence =", peak.H$convergence), col="magenta")
  
 title(main=SYs[i])

  dev.off()
  
}

OptimCompare <- data.frame(SYpeaks=c(1:6), 
                           Errors=c(length(which(SYpeaks1=="Error")),
                                    length(which(SYpeaks2=="Error")),
                                    length(which(SYpeaks3=="Error")),
                                    length(which(SYpeaks4=="Error")),
                                    length(which(SYpeaks5=="Error")),
                                    length(which(SYpeaks6=="Error"))),
                           NotConverge=c(sum(SYpeaks1$converge),
                                         sum(SYpeaks2$converge),
                                         sum(SYpeaks3$converge),
                                         sum(SYpeaks4$converge),
                                         sum(SYpeaks5$converge),
                                         sum(SYpeaks6$converge)),
                           VarNAs=c((length(which(is.na(sqrt(as.numeric(SYpeaks1$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks1$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks1$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks1$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks2$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks2$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks2$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks2$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks3$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks3$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks3$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks3$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks4$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks4$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks4$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks4$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks5$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks5$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks5$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks5$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks6$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks6$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks6$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks6$mu.var))==0)))),
                           RowNAs=c(length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$mu.var)))))))),
                           RowNAsORNotConv=c(length(unique(c(which(SYpeaks1$converge==1),which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$mu.var))))))),
                                             length(unique(c(which(SYpeaks2$converge==1),which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$mu.var))))))),
                                             length(unique(c(which(SYpeaks3$converge==1),which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$mu.var))))))),
                                             length(unique(c(which(SYpeaks4$converge==1),which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$mu.var))))))),
                                             length(unique(c(which(SYpeaks5$converge==1),which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$mu.var))))))),
                                             length(unique(c(which(SYpeaks6$converge==1),which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$mu.var)))))))))

#### Estimates for each site-year with varied start parameters ####
SYs <- data.frame(SY=unique(cater$sy))
for(j in 1:length(SYs$SY)){
  SYs$cater[j] <-  sum(cater$caterpillars[which(cater$sy==SYs$SY[j])])
}
SYs <- SYs$SY[SYs$cater>0]

#SYpeaks <- data.frame(SY=SYs) # with highest likelihood
SYpeaks1 <- data.frame(SY=SYs) #mean,10,mean Y black
SYpeaks2 <- data.frame(SY=SYs) #mean,10,1    Z gray
SYpeaks3 <- data.frame(SY=SYs) #mean,10,0    W red 
SYpeaks4 <- data.frame(SY=SYs) #mean,5,mean  F orange
SYpeaks5 <- data.frame(SY=SYs) #mean,5,1     G yellow
SYpeaks6 <- data.frame(SY=SYs) #mean,5,0     H magenta
SYpeaks7 <- data.frame(SY=SYs) #mean,12,mean  I purple
SYpeaks8 <- data.frame(SY=SYs) #mean,12,1     J blue
SYpeaks9 <- data.frame(SY=SYs) #mean,12,0     K darkblue
SYpeaks10 <- data.frame(SY=SYs) #mean,2,mean  L green
SYpeaks11 <- data.frame(SY=SYs) #mean,2,1     M darkgreen
SYpeaks12 <- data.frame(SY=SYs) #mean,2,0     N maroon
date <- seq(118,175,0.1)
par(mfrow=c(1,1))
i <- 218

for(i in 1:length(SYs)){
  
  X <- subset(cater, sy==SYs[i])
  X <- data.frame(date=X$date, cater=X$caterpillars)
  
  peak.Y<-optim(par=c(mean(X$date),10,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks1$n[i] <- length(X$date)
  SYpeaks1$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks1$cater[i] <- sum(X$cater)
  SYpeaks1$mu[i] <- peak.Y$par[1]
  SYpeaks1$mu.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[1], error=function(e){print("Error")})
  SYpeaks1$sigma[i] <- peak.Y$par[2]
  SYpeaks1$sigma.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[2], error=function(e){print("Error")})
  SYpeaks1$logheight[i] <- peak.Y$par[3]
  SYpeaks1$logheight.var[i] <- tryCatch(diag(solve(peak.Y$hessian))[3], error=function(e){print("Error")})
  SYpeaks1$loglik[i] <- peak.Y$value
  SYpeaks1$converge[i] <- peak.Y$convergence
  
  peak.Z<-optim(par=c(mean(X$date),10,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks2$n[i] <- length(X$date)
  SYpeaks2$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks2$cater[i] <- sum(X$cater)
  SYpeaks2$mu[i] <- peak.Z$par[1]
  SYpeaks2$mu.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[1], error=function(e){print("Error")})
  SYpeaks2$sigma[i] <- peak.Z$par[2]
  SYpeaks2$sigma.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[2], error=function(e){print("Error")})
  SYpeaks2$logheight[i] <- peak.Z$par[3]
  SYpeaks2$logheight.var[i] <- tryCatch(diag(solve(peak.Z$hessian))[3], error=function(e){print("Error")})
  SYpeaks2$loglik[i] <- peak.Z$value
  SYpeaks2$converge[i] <- peak.Z$convergence
  
  peak.W<-optim(par=c(mean(X$date),10,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks2$n[i] <- length(X$date)
  SYpeaks3$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks3$cater[i] <- sum(X$cater)
  SYpeaks3$mu[i] <- peak.W$par[1]
  SYpeaks3$mu.var[i] <- tryCatch(diag(solve(peak.W$hessian))[1], error=function(e){print("Error")})
  SYpeaks3$sigma[i] <- peak.W$par[2]
  SYpeaks3$sigma.var[i] <- tryCatch(diag(solve(peak.W$hessian))[2], error=function(e){print("Error")})
  SYpeaks3$logheight[i] <- peak.W$par[3]
  SYpeaks3$logheight.var[i] <- tryCatch(diag(solve(peak.W$hessian))[3], error=function(e){print("Error")})
  SYpeaks3$loglik[i] <- peak.W$value
  SYpeaks3$converge[i] <- peak.W$convergence
  
  peak.F<-optim(par=c(mean(X$date),5,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks4$n[i] <- length(X$date)
  SYpeaks4$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks4$cater[i] <- sum(X$cater)
  SYpeaks4$mu[i] <- peak.F$par[1]
  SYpeaks4$mu.var[i] <- tryCatch(diag(solve(peak.F$hessian))[1], error=function(e){print("Error")})
  SYpeaks4$sigma[i] <- peak.F$par[2]
  SYpeaks4$sigma.var[i] <- tryCatch(diag(solve(peak.F$hessian))[2], error=function(e){print("Error")})
  SYpeaks4$logheight[i] <- peak.F$par[3]
  SYpeaks4$logheight.var[i] <- tryCatch(diag(solve(peak.F$hessian))[3], error=function(e){print("Error")})
  SYpeaks4$loglik[i] <- peak.F$value
  SYpeaks4$converge[i] <- peak.F$convergence
  
  peak.G<-optim(par=c(mean(X$date),5,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks5$n[i] <- length(X$date)
  SYpeaks5$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks5$cater[i] <- sum(X$cater)
  SYpeaks5$mu[i] <- peak.G$par[1]
  SYpeaks5$mu.var[i] <- tryCatch(diag(solve(peak.G$hessian))[1], error=function(e){print("Error")})
  SYpeaks5$sigma[i] <- peak.G$par[2]
  SYpeaks5$sigma.var[i] <- tryCatch(diag(solve(peak.G$hessian))[2], error=function(e){print("Error")})
  SYpeaks5$logheight[i] <- peak.G$par[3]
  SYpeaks5$logheight.var[i] <- tryCatch(diag(solve(peak.G$hessian))[3], error=function(e){print("Error")})
  SYpeaks5$loglik[i] <- peak.G$value
  SYpeaks5$converge[i] <- peak.G$convergence
  
  peak.H<-optim(par=c(mean(X$date),5,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks6$n[i] <- length(X$date)
  SYpeaks6$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks6$cater[i] <- sum(X$cater)
  SYpeaks6$mu[i] <- peak.H$par[1]
  SYpeaks6$mu.var[i] <- tryCatch(diag(solve(peak.H$hessian))[1], error=function(e){print("Error")})
  SYpeaks6$sigma[i] <- peak.H$par[2]
  SYpeaks6$sigma.var[i] <- tryCatch(diag(solve(peak.H$hessian))[2], error=function(e){print("Error")})
  SYpeaks6$logheight[i] <- peak.H$par[3]
  SYpeaks6$logheight.var[i] <- tryCatch(diag(solve(peak.H$hessian))[3], error=function(e){print("Error")})
  SYpeaks6$loglik[i] <- peak.H$value
  SYpeaks6$converge[i] <- peak.H$convergence
  
  peak.I<-optim(par=c(mean(X$date),12,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
    
  SYpeaks7$n[i] <- length(X$date)
  SYpeaks7$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks7$cater[i] <- sum(X$cater)
  SYpeaks7$mu[i] <- peak.I$par[1]
  SYpeaks7$mu.var[i] <- tryCatch(diag(solve(peak.I$hessian))[1], error=function(e){print("Error")})
  SYpeaks7$sigma[i] <- peak.I$par[2]
  SYpeaks7$sigma.var[i] <- tryCatch(diag(solve(peak.I$hessian))[2], error=function(e){print("Error")})
  SYpeaks7$logheight[i] <- peak.I$par[3]
  SYpeaks7$logheight.var[i] <- tryCatch(diag(solve(peak.I$hessian))[3], error=function(e){print("Error")})
  SYpeaks7$loglik[i] <- peak.I$value
  SYpeaks7$converge[i] <- peak.I$convergence
  
  peak.J<-optim(par=c(mean(X$date),12,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks8$n[i] <- length(X$date)
  SYpeaks8$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks8$cater[i] <- sum(X$cater)
  SYpeaks8$mu[i] <- peak.J$par[1]
  SYpeaks8$mu.var[i] <- tryCatch(diag(solve(peak.J$hessian))[1], error=function(e){print("Error")})
  SYpeaks8$sigma[i] <- peak.J$par[2]
  SYpeaks8$sigma.var[i] <- tryCatch(diag(solve(peak.J$hessian))[2], error=function(e){print("Error")})
  SYpeaks8$logheight[i] <- peak.J$par[3]
  SYpeaks8$logheight.var[i] <- tryCatch(diag(solve(peak.J$hessian))[3], error=function(e){print("Error")})
  SYpeaks8$loglik[i] <- peak.J$value
  SYpeaks8$converge[i] <- peak.J$convergence
  
  peak.K<-optim(par=c(mean(X$date),12,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks9$n[i] <- length(X$date)
  SYpeaks9$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks9$cater[i] <- sum(X$cater)
  SYpeaks9$mu[i] <- peak.K$par[1]
  SYpeaks9$mu.var[i] <- tryCatch(diag(solve(peak.K$hessian))[1], error=function(e){print("Error")})
  SYpeaks9$sigma[i] <- peak.K$par[2]
  SYpeaks9$sigma.var[i] <- tryCatch(diag(solve(peak.K$hessian))[2], error=function(e){print("Error")})
  SYpeaks9$logheight[i] <- peak.K$par[3]
  SYpeaks9$logheight.var[i] <- tryCatch(diag(solve(peak.K$hessian))[3], error=function(e){print("Error")})
  SYpeaks9$loglik[i] <- peak.K$value
  SYpeaks9$converge[i] <- peak.K$convergence
  
  peak.L<-optim(par=c(mean(X$date),2,log(mean(X$cater))), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks10$n[i] <- length(X$date)
  SYpeaks10$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks10$cater[i] <- sum(X$cater)
  SYpeaks10$mu[i] <- peak.L$par[1]
  SYpeaks10$mu.var[i] <- tryCatch(diag(solve(peak.L$hessian))[1], error=function(e){print("Error")})
  SYpeaks10$sigma[i] <- peak.L$par[2]
  SYpeaks10$sigma.var[i] <- tryCatch(diag(solve(peak.L$hessian))[2], error=function(e){print("Error")})
  SYpeaks10$logheight[i] <- peak.L$par[3]
  SYpeaks10$logheight.var[i] <- tryCatch(diag(solve(peak.L$hessian))[3], error=function(e){print("Error")})
  SYpeaks10$loglik[i] <- peak.L$value
  SYpeaks10$converge[i] <- peak.L$convergence
  
  peak.M<-optim(par=c(mean(X$date),2,1), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks11$n[i] <- length(X$date)
  SYpeaks11$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks11$cater[i] <- sum(X$cater)
  SYpeaks11$mu[i] <- peak.M$par[1]
  SYpeaks11$mu.var[i] <- tryCatch(diag(solve(peak.M$hessian))[1], error=function(e){print("Error")})
  SYpeaks11$sigma[i] <- peak.M$par[2]
  SYpeaks11$sigma.var[i] <- tryCatch(diag(solve(peak.M$hessian))[2], error=function(e){print("Error")})
  SYpeaks11$logheight[i] <- peak.M$par[3]
  SYpeaks11$logheight.var[i] <- tryCatch(diag(solve(peak.M$hessian))[3], error=function(e){print("Error")})
  SYpeaks11$loglik[i] <- peak.M$value
  SYpeaks11$converge[i] <- peak.M$convergence
  
  peak.N<-optim(par=c(mean(X$date),2,0), fn=gaussianpeak,data=X,method="Nelder-Mead",hessian=T)
  
  SYpeaks12$n[i] <- length(X$date)
  SYpeaks12$propzero[i] <- length(which(X$cater==0))/length(X$cater)
  SYpeaks12$cater[i] <- sum(X$cater)
  SYpeaks12$mu[i] <- peak.N$par[1]
  SYpeaks12$mu.var[i] <- tryCatch(diag(solve(peak.N$hessian))[1], error=function(e){print("Error")})
  SYpeaks12$sigma[i] <- peak.N$par[2]
  SYpeaks12$sigma.var[i] <- tryCatch(diag(solve(peak.N$hessian))[2], error=function(e){print("Error")})
  SYpeaks12$logheight[i] <- peak.N$par[3]
  SYpeaks12$logheight.var[i] <- tryCatch(diag(solve(peak.N$hessian))[3], error=function(e){print("Error")})
  SYpeaks12$loglik[i] <- peak.N$value
  SYpeaks12$converge[i] <- peak.N$convergence
  
  mypath <- file.path(paste("~/Documents/Models/OptimSYpeaks/",SYs[i],".pdf"))
  
  pdf(file = mypath, width=10, height=20)
  
  plot(date, exp((-((date-peak.Y$par[1])^2)/(2*peak.Y$par[2]^2)))*exp(peak.Y$par[3]),type="l", ylim=c(0,max(X$cater)), xlab="Date", ylab="Abund", col="black")
  points(jitter(X$date,1), jitter(X$cater,0.05))
  text(x=126,y=max(X$cater), paste("mu =", round(peak.Y$par[1], 2), "var =", ifelse(SYpeaks1$mu.var[i]=="Error", SYpeaks1$mu.var[i], round(as.numeric(SYpeaks1$mu.var[i]),2))), col="black")
  text(x=126,y=max(X$cater)-0.01*max(X$cater), paste("sigma =", round(peak.Y$par[2], 2), "var =", ifelse(SYpeaks1$sigma.var[i]=="Error", SYpeaks1$sigma.var[i], round(as.numeric(SYpeaks1$sigma.var[i]),2))), col="black")
  text(x=127.5,y=max(X$cater)-0.02*max(X$cater), paste("logheight =", round(peak.Y$par[3], 2), "var =", ifelse(SYpeaks1$logheight.var[i]=="Error", SYpeaks1$logheight.var[i], round(as.numeric(SYpeaks1$logheight.var[i]),2))), col="black")
  text(x=129,y=max(X$cater)-0.03*max(X$cater), paste("loglik =", round(peak.Y$value, 2), "convergence =", peak.Y$convergence), col="black")
  
  points(date, exp((-((date-peak.Z$par[1])^2)/(2*peak.Z$par[2]^2)))*exp(peak.Z$par[3]),type="l", col="gray")
  text(x=126,y=max(X$cater)-0.05*max(X$cater), paste("mu =", round(peak.Z$par[1], 2), "var =", ifelse(SYpeaks2$mu.var[i]=="Error", SYpeaks2$mu.var[i], round(as.numeric(SYpeaks2$mu.var[i]),2))), col="gray")
  text(x=126,y=max(X$cater)-0.06*max(X$cater), paste("sigma =", round(peak.Z$par[2], 2), "var =", ifelse(SYpeaks2$sigma.var[i]=="Error", SYpeaks2$sigma.var[i], round(as.numeric(SYpeaks2$sigma.var[i]),2))), col="gray")
  text(x=127.5,y=max(X$cater)-0.07*max(X$cater), paste("logheight =", round(peak.Z$par[3], 2), "var =", ifelse(SYpeaks2$logheight.var[i]=="Error", SYpeaks2$logheight.var[i], round(as.numeric(SYpeaks2$logheight.var[i]),2))), col="gray")
  text(x=129,y=max(X$cater)-0.08*max(X$cater), paste("loglik =", round(peak.Z$value, 2), "convergence =", peak.Z$convergence), col="gray")
  
  points(date, exp((-((date-peak.W$par[1])^2)/(2*peak.W$par[2]^2)))*exp(peak.W$par[3]),type="l", col="red")
  text(x=126,y=max(X$cater)-0.10*max(X$cater), paste("mu =", round(peak.W$par[1], 2), "var =", ifelse(SYpeaks3$mu.var[i]=="Error", SYpeaks3$mu.var[i], round(as.numeric(SYpeaks3$mu.var[i]),2))), col="red")
  text(x=126,y=max(X$cater)-0.11*max(X$cater), paste("sigma =", round(peak.W$par[2], 2), "var =", ifelse(SYpeaks3$sigma.var[i]=="Error", SYpeaks3$sigma.var[i], round(as.numeric(SYpeaks3$sigma.var[i]),2))), col="red")
  text(x=127.5,y=max(X$cater)-0.12*max(X$cater), paste("logheight =", round(peak.W$par[3], 2), "var =", ifelse(SYpeaks3$logheight.var[i]=="Error", SYpeaks3$logheight.var[i], round(as.numeric(SYpeaks3$logheight.var[i]),2))), col="red")
  text(x=129,y=max(X$cater)-0.13*max(X$cater), paste("loglik =", round(peak.W$value, 2), "convergence =", peak.W$convergence), col="red")
  
  points(date, exp((-((date-peak.F$par[1])^2)/(2*peak.F$par[2]^2)))*exp(peak.F$par[3]),type="l", col="orange")
  text(x=126,y=max(X$cater)-0.15*max(X$cater), paste("mu =", round(peak.F$par[1], 2), "var =", ifelse(SYpeaks4$mu.var[i]=="Error", SYpeaks4$mu.var[i], round(as.numeric(SYpeaks4$mu.var[i]),2))), col="orange")
  text(x=126,y=max(X$cater)-0.16*max(X$cater), paste("sigma =", round(peak.F$par[2], 2), "var =", ifelse(SYpeaks4$sigma.var[i]=="Error", SYpeaks4$sigma.var[i], round(as.numeric(SYpeaks4$sigma.var[i]),2))), col="orange")
  text(x=127.5,y=max(X$cater)-0.17*max(X$cater), paste("logheight =", round(peak.F$par[3], 2), "var =", ifelse(SYpeaks4$logheight.var[i]=="Error", SYpeaks4$logheight.var[i], round(as.numeric(SYpeaks4$logheight.var[i]),2))), col="orange")
  text(x=129,y=max(X$cater)-0.18*max(X$cater), paste("loglik =", round(peak.F$value, 2), "convergence =", peak.F$convergence), col="orange")
  
  points(date, exp((-((date-peak.G$par[1])^2)/(2*peak.G$par[2]^2)))*exp(peak.G$par[3]),type="l", col="chocolate1")
  text(x=126,y=max(X$cater)-0.20*max(X$cater), paste("mu =", round(peak.G$par[1], 2), "var =", ifelse(SYpeaks5$mu.var[i]=="Error", SYpeaks5$mu.var[i], round(as.numeric(SYpeaks5$mu.var[i]),2))), col="chocolate1")
  text(x=126,y=max(X$cater)-0.21*max(X$cater), paste("sigma =", round(peak.G$par[2], 2), "var =", ifelse(SYpeaks5$sigma.var[i]=="Error", SYpeaks5$sigma.var[i], round(as.numeric(SYpeaks5$sigma.var[i]),2))), col="chocolate1")
  text(x=127.5,y=max(X$cater)-0.22*max(X$cater), paste("logheight =", round(peak.G$par[3], 2), "var =", ifelse(SYpeaks5$logheight.var[i]=="Error", SYpeaks5$logheight.var[i], round(as.numeric(SYpeaks5$logheight.var[i]),2))), col="chocolate1")
  text(x=129,y=max(X$cater)-0.23*max(X$cater), paste("loglik =", round(peak.G$value, 2), "convergence =", peak.G$convergence), col="chocolate1")
  
  points(date, exp((-((date-peak.H$par[1])^2)/(2*peak.H$par[2]^2)))*exp(peak.H$par[3]),type="l", col="magenta")
  text(x=126,y=max(X$cater)-0.25*max(X$cater), paste("mu =", round(peak.H$par[1], 2), "var =", ifelse(SYpeaks6$mu.var[i]=="Error", SYpeaks6$mu.var[i], round(as.numeric(SYpeaks6$mu.var[i]),2))), col="magenta")
  text(x=126,y=max(X$cater)-0.26*max(X$cater), paste("sigma =", round(peak.H$par[2], 2), "var =", ifelse(SYpeaks6$sigma.var[i]=="Error", SYpeaks6$sigma.var[i], round(as.numeric(SYpeaks6$sigma.var[i]),2))), col="magenta")
  text(x=127.5,y=max(X$cater)-0.27*max(X$cater), paste("logheight =", round(peak.H$par[3], 2), "var =", ifelse(SYpeaks6$logheight.var[i]=="Error", SYpeaks6$logheight.var[i], round(as.numeric(SYpeaks6$logheight.var[i]),2))), col="magenta")
  text(x=129,y=max(X$cater)-0.28*max(X$cater), paste("loglik =", round(peak.H$value, 2), "convergence =", peak.H$convergence), col="magenta")
  
  points(date, exp((-((date-peak.I$par[1])^2)/(2*peak.I$par[2]^2)))*exp(peak.I$par[3]),type="l", col="purple")
  text(x=126,y=max(X$cater)-0.30*max(X$cater), paste("mu =", round(peak.I$par[1], 2), "var =", ifelse(SYpeaks7$mu.var[i]=="Error", SYpeaks7$mu.var[i], round(as.numeric(SYpeaks7$mu.var[i]),2))), col="purple")
  text(x=126,y=max(X$cater)-0.31*max(X$cater), paste("sigma =", round(peak.I$par[2], 2), "var =", ifelse(SYpeaks7$sigma.var[i]=="Error", SYpeaks7$sigma.var[i], round(as.numeric(SYpeaks7$sigma.var[i]),2))), col="purple")
  text(x=127.5,y=max(X$cater)-0.32*max(X$cater), paste("logheight =", round(peak.I$par[3], 2), "var =", ifelse(SYpeaks7$logheight.var[i]=="Error", SYpeaks7$logheight.var[i], round(as.numeric(SYpeaks7$logheight.var[i]),2))), col="purple")
  text(x=129,y=max(X$cater)-0.33*max(X$cater), paste("loglik =", round(peak.I$value, 2), "convergence =", peak.I$convergence), col="purple")
  
  points(date, exp((-((date-peak.J$par[1])^2)/(2*peak.J$par[2]^2)))*exp(peak.J$par[3]),type="l", col="blue")
  text(x=126,y=max(X$cater)-0.35*max(X$cater), paste("mu =", round(peak.J$par[1], 2), "var =", ifelse(SYpeaks8$mu.var[i]=="Error", SYpeaks8$mu.var[i], round(as.numeric(SYpeaks8$mu.var[i]),2))), col="blue")
  text(x=126,y=max(X$cater)-0.36*max(X$cater), paste("sigma =", round(peak.J$par[2], 2), "var =", ifelse(SYpeaks8$sigma.var[i]=="Error", SYpeaks8$sigma.var[i], round(as.numeric(SYpeaks8$sigma.var[i]),2))), col="blue")
  text(x=127.5,y=max(X$cater)-0.37*max(X$cater), paste("logheight =", round(peak.J$par[3], 2), "var =", ifelse(SYpeaks8$logheight.var[i]=="Error", SYpeaks8$logheight.var[i], round(as.numeric(SYpeaks8$logheight.var[i]),2))), col="blue")
  text(x=129,y=max(X$cater)-0.38*max(X$cater), paste("loglik =", round(peak.J$value, 2), "convergence =", peak.J$convergence), col="blue")
  
  points(date, exp((-((date-peak.K$par[1])^2)/(2*peak.K$par[2]^2)))*exp(peak.K$par[3]),type="l", col="darkblue")
  text(x=126,y=max(X$cater)-0.40*max(X$cater), paste("mu =", round(peak.K$par[1], 2), "var =", ifelse(SYpeaks9$mu.var[i]=="Error", SYpeaks9$mu.var[i], round(as.numeric(SYpeaks9$mu.var[i]),2))), col="darkblue")
  text(x=126,y=max(X$cater)-0.41*max(X$cater), paste("sigma =", round(peak.K$par[2], 2), "var =", ifelse(SYpeaks9$sigma.var[i]=="Error", SYpeaks9$sigma.var[i], round(as.numeric(SYpeaks9$sigma.var[i]),2))), col="darkblue")
  text(x=127.5,y=max(X$cater)-0.42*max(X$cater), paste("logheight =", round(peak.K$par[3], 2), "var =", ifelse(SYpeaks9$logheight.var[i]=="Error", SYpeaks9$logheight.var[i], round(as.numeric(SYpeaks9$logheight.var[i]),2))), col="darkblue")
  text(x=129,y=max(X$cater)-0.43*max(X$cater), paste("loglik =", round(peak.K$value, 2), "convergence =", peak.K$convergence), col="darkblue")
  
  points(date, exp((-((date-peak.L$par[1])^2)/(2*peak.L$par[2]^2)))*exp(peak.L$par[3]),type="l", col="green")
  text(x=126,y=max(X$cater)-0.45*max(X$cater), paste("mu =", round(peak.L$par[1], 2), "var =", ifelse(SYpeaks10$mu.var[i]=="Error", SYpeaks10$mu.var[i], round(as.numeric(SYpeaks10$mu.var[i]),2))), col="green")
  text(x=126,y=max(X$cater)-0.46*max(X$cater), paste("sigma =", round(peak.L$par[2], 2), "var =", ifelse(SYpeaks10$sigma.var[i]=="Error", SYpeaks10$sigma.var[i], round(as.numeric(SYpeaks10$sigma.var[i]),2))), col="green")
  text(x=127.5,y=max(X$cater)-0.47*max(X$cater), paste("logheight =", round(peak.L$par[3], 2), "var =", ifelse(SYpeaks10$logheight.var[i]=="Error", SYpeaks10$logheight.var[i], round(as.numeric(SYpeaks10$logheight.var[i]),2))), col="green")
  text(x=129,y=max(X$cater)-0.48*max(X$cater), paste("loglik =", round(peak.L$value, 2), "convergence =", peak.L$convergence), col="green")
  
  points(date, exp((-((date-peak.M$par[1])^2)/(2*peak.M$par[2]^2)))*exp(peak.M$par[3]),type="l", col="darkgreen")
  text(x=126,y=max(X$cater)-0.50*max(X$cater), paste("mu =", round(peak.M$par[1], 2), "var =", ifelse(SYpeaks11$mu.var[i]=="Error", SYpeaks11$mu.var[i], round(as.numeric(SYpeaks11$mu.var[i]),2))), col="darkgreen")
  text(x=126,y=max(X$cater)-0.51*max(X$cater), paste("sigma =", round(peak.M$par[2], 2), "var =", ifelse(SYpeaks11$sigma.var[i]=="Error", SYpeaks11$sigma.var[i], round(as.numeric(SYpeaks11$sigma.var[i]),2))), col="darkgreen")
  text(x=127.5,y=max(X$cater)-0.52*max(X$cater), paste("logheight =", round(peak.M$par[3], 2), "var =", ifelse(SYpeaks11$logheight.var[i]=="Error", SYpeaks11$logheight.var[i], round(as.numeric(SYpeaks11$logheight.var[i]),2))), col="darkgreen")
  text(x=129,y=max(X$cater)-0.53*max(X$cater), paste("loglik =", round(peak.M$value, 2), "convergence =", peak.M$convergence), col="darkgreen")
  
  points(date, exp((-((date-peak.N$par[1])^2)/(2*peak.N$par[2]^2)))*exp(peak.N$par[3]),type="l", col="maroon")
  text(x=126,y=max(X$cater)-0.55*max(X$cater), paste("mu =", round(peak.N$par[1], 2), "var =", ifelse(SYpeaks12$mu.var[i]=="Error", SYpeaks12$mu.var[i], round(as.numeric(SYpeaks12$mu.var[i]),2))), col="maroon")
  text(x=126,y=max(X$cater)-0.56*max(X$cater), paste("sigma =", round(peak.N$par[2], 2), "var =", ifelse(SYpeaks12$sigma.var[i]=="Error", SYpeaks12$sigma.var[i], round(as.numeric(SYpeaks12$sigma.var[i]),2))), col="maroon")
  text(x=127.5,y=max(X$cater)-0.57*max(X$cater), paste("logheight =", round(peak.N$par[3], 2), "var =", ifelse(SYpeaks12$logheight.var[i]=="Error", SYpeaks12$logheight.var[i], round(as.numeric(SYpeaks12$logheight.var[i]),2))), col="maroon")
  text(x=129,y=max(X$cater)-0.58*max(X$cater), paste("loglik =", round(peak.N$value, 2), "convergence =", peak.N$convergence), col="maroon")
  
  title(main=SYs[i])
  
  dev.off()
  
}

OptimCompare <- data.frame(SYpeaks=c(1:12), 
                           Errors=c(length(which(SYpeaks1=="Error")),
                                    length(which(SYpeaks2=="Error")),
                                    length(which(SYpeaks3=="Error")),
                                    length(which(SYpeaks4=="Error")),
                                    length(which(SYpeaks5=="Error")),
                                    length(which(SYpeaks6=="Error")),
                                    length(which(SYpeaks7=="Error")),
                                    length(which(SYpeaks8=="Error")),
                                    length(which(SYpeaks9=="Error")),
                                    length(which(SYpeaks10=="Error")),
                                    length(which(SYpeaks11=="Error")),
                                    length(which(SYpeaks12=="Error"))),
                           NotConverge=c(sum(SYpeaks1$converge),
                                         sum(SYpeaks2$converge),
                                         sum(SYpeaks3$converge),
                                         sum(SYpeaks4$converge),
                                         sum(SYpeaks5$converge),
                                         sum(SYpeaks6$converge),
                                         sum(SYpeaks7$converge),
                                         sum(SYpeaks8$converge),
                                         sum(SYpeaks9$converge),
                                         sum(SYpeaks10$converge),
                                         sum(SYpeaks11$converge),
                                         sum(SYpeaks12$converge)),
                           VarNAs=c((length(which(is.na(sqrt(as.numeric(SYpeaks1$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks1$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks1$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks1$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks2$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks2$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks2$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks2$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks3$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks3$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks3$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks3$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks4$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks4$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks4$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks4$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks5$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks5$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks5$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks5$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks6$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks6$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks6$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks6$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks7$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks7$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks7$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks7$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks7$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks7$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks8$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks8$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks8$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks8$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks8$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks8$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks9$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks9$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks9$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks9$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks9$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks9$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks10$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks10$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks10$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks10$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks10$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks10$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks11$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks11$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks11$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks11$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks11$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks11$mu.var))==0))),
                                    (length(which(is.na(sqrt(as.numeric(SYpeaks12$mu.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks12$sigma.var)))))+length(which(is.na(sqrt(as.numeric(SYpeaks12$logheight.var)))))+length(which(sqrt(as.numeric(SYpeaks12$logheight.var))==0))+length(which(sqrt(as.numeric(SYpeaks12$sigma.var))==0))+length(which(sqrt(as.numeric(SYpeaks12$mu.var))==0)))),
                           RowNAs=c(length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks7$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks7$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks7$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks8$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks8$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks8$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks9$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks9$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks9$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks10$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks10$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks10$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks11$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks11$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks11$mu.var))))))),
                                    length(unique(c(which(is.na(sqrt(as.numeric(SYpeaks12$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks12$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks12$mu.var)))))))),
                           RowNAsORNotConv=c(length(unique(c(which(SYpeaks1$converge==1),which(is.na(sqrt(as.numeric(SYpeaks1$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks1$mu.var))))))),
                                    length(unique(c(which(SYpeaks2$converge==1),which(is.na(sqrt(as.numeric(SYpeaks2$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks2$mu.var))))))),
                                    length(unique(c(which(SYpeaks3$converge==1),which(is.na(sqrt(as.numeric(SYpeaks3$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks3$mu.var))))))),
                                    length(unique(c(which(SYpeaks4$converge==1),which(is.na(sqrt(as.numeric(SYpeaks4$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks4$mu.var))))))),
                                    length(unique(c(which(SYpeaks5$converge==1),which(is.na(sqrt(as.numeric(SYpeaks5$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks5$mu.var))))))),
                                    length(unique(c(which(SYpeaks6$converge==1),which(is.na(sqrt(as.numeric(SYpeaks6$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks6$mu.var))))))),
                                    length(unique(c(which(SYpeaks7$converge==1),which(is.na(sqrt(as.numeric(SYpeaks7$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks7$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks7$mu.var))))))),
                                    length(unique(c(which(SYpeaks8$converge==1),which(is.na(sqrt(as.numeric(SYpeaks8$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks8$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks8$mu.var))))))),
                                    length(unique(c(which(SYpeaks9$converge==1),which(is.na(sqrt(as.numeric(SYpeaks9$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks9$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks9$mu.var))))))),
                                    length(unique(c(which(SYpeaks10$converge==1),which(is.na(sqrt(as.numeric(SYpeaks10$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks10$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks10$mu.var))))))),
                                    length(unique(c(which(SYpeaks11$converge==1),which(is.na(sqrt(as.numeric(SYpeaks11$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks11$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks11$mu.var))))))),
                                    length(unique(c(which(SYpeaks12$converge==1),which(is.na(sqrt(as.numeric(SYpeaks12$logheight.var)))),which(is.na(sqrt(as.numeric(SYpeaks12$sigma.var)))),which(is.na(sqrt(as.numeric(SYpeaks12$mu.var)))))))))


#### Estimates for each site-year choosing highest liklihood option ####
SYs <- data.frame(SY=unique(cater$sy))
for(j in 1:length(SYs$SY)){
  SYs$cater[j] <-  sum(cater$caterpillars[which(cater$sy==SYs$SY[j])])
}
SYs <- SYs$SY[SYs$cater>0]

SYpeaks <- data.frame(SY=SYs) # with highest likelihood
#peak.Y= mean,10,mean  Y black
#peak.Z= mean,10,1     Z gray
#peak.W= mean,10,0     W red 
#peak.F= mean,5,mean   F orange
#peak.G= mean,5,1      G chocolate1
#peak.H= mean,5,0      H magenta
#peak.I= mean,12,mean  I purple
#peak.J= mean,12,1     J blue
#peak.K= mean,12,0     K darkblue
#peak.L= mean,2,mean   L green
#peak.M= mean,2,1      M darkgreen
#peak.N= mean,2,0      N maroon
#date <- seq(118,175,0.1)

#i <- 17

for(i in 1:length(SYs)){
  
  X <- subset(cater, sy==SYs[i])
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
