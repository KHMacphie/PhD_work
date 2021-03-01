#load Dataframe.R script rm(list=ls()) setwd('/Users/s1205615/')

############################################################
#### Model: Asymmetry in the caterpillar abundance peak ####
############################################################

a<-10000
prior<-list(R=list(V=diag(1), nu=0.002), 
            G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

#CurveShape20<- MCMCglmm(caterpillars~datescaled+I(datescaled^2)+I(datescaled^3), 
#                 random=~us(1+datescaled+I(datescaled^2)):siteyear+recorder+siteday+treeID,
#                 family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000, thin=1000)
#save(CurveShape20, file = "~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/CurveShape20.RData") #sample size 1950
load("~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/CurveShape20.RData")
summary(CurveShape20)

#### Checking model fits the data and converged ####
plot(CurveShape20) #look at fixed effect and random term trace plots 
CurveShape20.Sim<-simulate(CurveShape20,nsim=1000) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(CurveShape20.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(CurveShape20.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data

## Differentiation to calculate peak date
# (-b +/- sqrt(b^2-4*a*c))/(2*a)
cdf <- data.frame(CurveShape20$Sol[,1:4])
cdf$a <- 3*CurveShape20$Sol[,4]
cdf$b <- 2*CurveShape20$Sol[,3]
cdf$c <- CurveShape20$Sol[,2]

cdf$pd <- ((-cdf$b - sqrt(cdf$b^2-4*cdf$a*cdf$c))/(2*cdf$a))

## Calculate peak height
cdf$ph <- CurveShape20$Sol[,1]+CurveShape20$Sol[,2]*cdf$pd+CurveShape20$Sol[,3]*cdf$pd^2+CurveShape20$Sol[,4]*cdf$pd^3

## Calculate width either side at 50% peak height
# manipulate peak height to find dates the peak is at a certain proportion of the peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/2))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.5[i] <- A[1]
  cdf$r2.5[i] <- A[2]
  cdf$r3.5[i] <- A[3]
}

# turn into real numbers
cdf$r1.5 <- Re(cdf$r1.5)[abs(Im(cdf$r1.5)) < 1e-6]
cdf$r2.5 <- Re(cdf$r2.5)[abs(Im(cdf$r2.5)) < 1e-6]
cdf$r3.5 <- Re(cdf$r3.5)[abs(Im(cdf$r3.5)) < 1e-6]

#remove the root that is outside of the data range (cubic curve comes back up)
cdf$r1.5 <- ifelse(cdf$r1.5<min(cater_habitat$datescaled),NA,cdf$r1.5)
cdf$r2.5 <- ifelse(cdf$r2.5<min(cater_habitat$datescaled),NA,cdf$r2.5)
cdf$r3.5 <- ifelse(cdf$r3.5<min(cater_habitat$datescaled),NA,cdf$r3.5)

#get real dates
roots.5 <- data.frame(r1.5=cdf$r1.5,r2.5=cdf$r2.5,r3.5=cdf$r3.5)
cdf$root1.5 <- ggfortify::unscale((apply(roots.5, 1, min, na.rm=TRUE)), center= 146.7727, scale=14.04083)  
cdf$root2.5 <- ggfortify::unscale((apply(roots.5, 1, max, na.rm=TRUE)), center= 146.7727, scale=14.04083) 

#days either side of peak
cdf$left.5 <- ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083)-cdf$root1.5
cdf$left.5 <- cdf$left.5$var1
cdf$right.5 <- cdf$root2.5-(ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083))
cdf$right.5 <- cdf$right.5$V1

#proportion of width
cdf$width.5 <- cdf$left.5+cdf$right.5
cdf$propleft.5 <- cdf$left.5/cdf$width.5
cdf$propright.5 <- cdf$right.5/cdf$width.5

mean(cdf$propleft.5) # 0.5436443 
HPDinterval(mcmc(cdf$propleft.5)) # 0.5318943 0.5551997
mean(cdf$propright.5) # 0.4563557
HPDinterval(mcmc(cdf$propright.5)) # 0.4448003 0.4681057
mean(cdf$width.5) # 24.20633
HPDinterval(mcmc(cdf$width.5)) # 22.70739 25.83797

#### Calculate width either side at 25% peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/4))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.25[i] <- A[1]
  cdf$r2.25[i] <- A[2]
  cdf$r3.25[i] <- A[3]
}

cdf$r1.25 <- Re(cdf$r1.25)[abs(Im(cdf$r1.25)) < 1e-6]
cdf$r2.25 <- Re(cdf$r2.25)[abs(Im(cdf$r2.25)) < 1e-6]
cdf$r3.25 <- Re(cdf$r3.25)[abs(Im(cdf$r3.25)) < 1e-6]

cdf$r1.25 <- ifelse(cdf$r1.25<min(cater_habitat$datescaled),NA,cdf$r1.25)
cdf$r2.25 <- ifelse(cdf$r2.25<min(cater_habitat$datescaled),NA,cdf$r2.25)
cdf$r3.25 <- ifelse(cdf$r3.25<min(cater_habitat$datescaled),NA,cdf$r3.25)

roots.25 <- data.frame(r1.25=cdf$r1.25,r2.25=cdf$r2.25,r3.25=cdf$r3.25)
cdf$root1.25 <- ggfortify::unscale((apply(roots.25, 1, min, na.rm=TRUE)), center= 146.7727, scale=14.04083)   
cdf$root2.25 <- ggfortify::unscale((apply(roots.25, 1, max, na.rm=TRUE)), center= 146.7727, scale=14.04083) 

cdf$left.25 <- ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083)-cdf$root1.25
cdf$left.25 <- cdf$left.25$var1
cdf$right.25 <- cdf$root2.25-ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083)
cdf$right.25 <- cdf$right.25$V1
cdf$width.25 <- cdf$left.25+cdf$right.25
cdf$propleft.25 <- cdf$left.25/cdf$width.25
cdf$propright.25 <- cdf$right.25/cdf$width.25

mean(cdf$propleft.25) # 0.5651141
HPDinterval(mcmc(cdf$propleft.25)) # 0.5457612 0.584862
mean(cdf$propright.25) # 0.4348859
HPDinterval(mcmc(cdf$propright.25)) # 0.415138 0.4542388
mean(cdf$width.25) # 35.04902
HPDinterval(mcmc(cdf$width.25)) # 32.7612 37.32882

#### Calculate width either side at 75% peak height
for(i in 1:length(cdf$a)){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])*0.75))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.75[i] <- A[1]
  cdf$r2.75[i] <- A[2]
  cdf$r3.75[i] <- A[3]
}

cdf$r1.75 <- Re(cdf$r1.75)[abs(Im(cdf$r1.75)) < 1e-6]
cdf$r2.75 <- Re(cdf$r2.75)[abs(Im(cdf$r2.75)) < 1e-6]
cdf$r3.75 <- Re(cdf$r3.75)[abs(Im(cdf$r3.75)) < 1e-6]

cdf$r1.75 <- ifelse(cdf$r1.75<min(cater_habitat$datescaled),NA,cdf$r1.75)
cdf$r2.75 <- ifelse(cdf$r2.75<min(cater_habitat$datescaled),NA,cdf$r2.75)
cdf$r3.75 <- ifelse(cdf$r3.75<min(cater_habitat$datescaled),NA,cdf$r3.75)

roots.75 <- data.frame(r1.75=cdf$r1.75,r2.75=cdf$r2.75,r3.75=cdf$r3.75)
cdf$root1.75 <- ggfortify::unscale((apply(roots.75, 1, min, na.rm=TRUE)), center= 146.7727, scale=14.04083)   
cdf$root2.75 <- ggfortify::unscale((apply(roots.75, 1, max, na.rm=TRUE)), center= 146.7727, scale=14.04083)  

cdf$left.75 <- ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083)-cdf$root1.75
cdf$left.75 <- cdf$left.75$var1
cdf$right.75 <- cdf$root2.75-ggfortify::unscale(cdf$pd, center= 146.7727, scale=14.04083)
cdf$right.75 <- cdf$right.75$V1
cdf$width.75 <- cdf$left.75+cdf$right.75
cdf$propleft.75 <- cdf$left.75/cdf$width.75
cdf$propright.75 <- cdf$right.75/cdf$width.75

mean(cdf$propleft.75) # 0.5273824
HPDinterval(mcmc(cdf$propleft.75)) # 0.520277 0.5341462
mean(cdf$propright.75) # 0.4726176
HPDinterval(mcmc(cdf$propright.75)) # 0.4658538 0.479723
mean(cdf$width.75) # 15.41359
HPDinterval(mcmc(cdf$width.75)) # 14.34555 16.386

#### Plot curves
# colour: 
mycol <- rgb(0, 120, 0, max = 250, alpha = 15, names = "greentrans")

# mean trend
dayscal <- seq(-2.1199,2.0102,0.001)
curve <- mean(CurveShape20$Sol[,1])+mean(CurveShape20$Sol[,2])*dayscal+mean(CurveShape20$Sol[,3])*dayscal^2+mean(CurveShape20$Sol[,4])*dayscal^3
days <- unscale(dayscal, center= 146.7727, scale=14.04083)$V1

# metrics for lines 
quart <- data.frame(qd=seq(mean(cdf$root1.25$V1),(mean(cdf$root2.25$V1)),0.1))
quart$qh <- mean(exp(cdf$ph)/4)
half <- data.frame(hd=seq((mean(cdf$root1.5$V1)+0.2),(mean(cdf$root2.5$V1)+0.2),0.1))
half$hh <- mean(exp(cdf$ph)/2)
tquart <- data.frame(tqd=seq((mean(cdf$root1.75$V1)+0.45),(mean(cdf$root2.75$V1)+0.45),0.1))
tquart$tqh <- mean(exp(cdf$ph)*0.75)

## Plotting posterior distribution curves with percentage of distribution to either side of the peak
par(mfcol=c(1,1),mar=c(3.9, 3.8, 1, 1), cex=1.4, las=1)
plot(days,exp(curve), type="l", ylim=c(0,0.095), xlab="Ordinal Date", ylab="Abundance", yaxs="i")

for(i in 1:1950){ 
  A <- CurveShape20$Sol[i,1]+CurveShape20$Sol[i,2]*dayscal+CurveShape20$Sol[i,3]*dayscal^2+CurveShape20$Sol[i,4]*dayscal^3
  points(days, exp(A), type="l", col=mycol, lwd=0.5)
}

points(quart$qd, quart$qh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(half$hd, half$hh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(tquart$tqd, tquart$tqh, type="l", lty="dashed", lwd=0.7, col="gray66")
abline(v=mean(unscale(cdf$pd, center= 146.7727, scale=14.04083)$var1), lwd=0.8, lty="dashed", col="gray66")
points(days, exp(curve), type="l")

text(150, 0.0135, "56.5%", cex=0.9, col="gray40")
text(158, 0.0135, "43.5%", cex=0.9, col="gray40")
text(150, 0.03, "54.4%", cex=0.9, col="gray40")
text(158, 0.03, "45.6%", cex=0.9, col="gray40")
text(150.3, 0.046, "52.7%", cex=0.9, col="gray40")
text(158, 0.046, "47.3%", cex=0.9, col="gray40")

text(125, mean(quart$qh), "0.25", cex=0.9, col="gray66")
text(125.4, mean(half$hh), "0.5", cex=0.9, col="gray66")
text(125, mean(tquart$tqh), "0.75", cex=0.9, col="gray66")

arrows(x0=127.5, y0=mean(half$hh), x1=136, y1=mean(half$hh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(quart$qh), x1=130.25, y1=mean(quart$qh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(tquart$tqh), x1=140, y1=mean(tquart$tqh), length=0.1, col="gray66") #Saved as 8"x8"

############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(CurveShape20$Sol[,1]),3)," (",
                      round(HPDinterval(CurveShape20$Sol[,1])[1],3)," - ",
                      round(HPDinterval(CurveShape20$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(CurveShape20$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(CurveShape20$Sol[,2]),3)," (",
                          round(HPDinterval(CurveShape20$Sol[,2])[1],3)," - ",
                          round(HPDinterval(CurveShape20$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(CurveShape20$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(CurveShape20$Sol[,3]),3)," (",
                           round(HPDinterval(CurveShape20$Sol[,3])[1],3)," - ",
                           round(HPDinterval(CurveShape20$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(CurveShape20$Sol[,3]))),
  c("Date³ (scaled)",paste(round(mean(CurveShape20$Sol[,4]),3)," (",
                           round(HPDinterval(CurveShape20$Sol[,4])[1],3)," - ",
                           round(HPDinterval(CurveShape20$Sol[,4])[2],3),")",sep=""),
    round(effectiveSize(CurveShape20$Sol[,4]))))

#### random terms
column<-1
siteyear1<-c("SiteYear- Intercept var",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                             round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-2
siteyear2<-c("SiteYear- Intercept:Date slope covar",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                                          round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-3
siteyear3<-c("SiteYear- Intercept:Date² slope covar",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                                           round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-5
siteyear5<-c("SiteYear- Date slope var",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                              round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-6
siteyear6<-c("SiteYear- Date slope:Date² slope covar",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                                            round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-9
siteyear9<-c("SiteYear- Date² slope var",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                                               round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(CurveShape20$VCV[, column])))

column<-10
recorder<-c("Recorder",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape20$VCV[, column])))


column<-11
siteday<-c("Site Day",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                            round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                            round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(CurveShape20$VCV[, column])))


column<-12
treeID<-c("Tree ID",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                          round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                          round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(CurveShape20$VCV[, column])))

column<-13
residual<-c("Residual",paste(round(posterior.mode(CurveShape20$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape20$VCV[, column])))

random<-rbind(siteyear1,siteyear2,siteyear3,siteyear5,siteyear6,siteyear9, recorder, siteday, treeID, residual)

#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/Inc2020/TableCurveShape20.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)