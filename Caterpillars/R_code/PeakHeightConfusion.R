rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)
library(forcats)
library(gridExtra)

###############################
#### Peak height confusion ####
###############################

load("/Users/s1205615/Documents/PhD (stopped using 12:2:19)/GitHub/R/Caterpillar analysis/results/TempPoisson.RData")


pred_day<-seq(120,175,0.5)
curve_2 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(2*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(2*TempPoisson$Sol[,"Apr"]) #Apr
curve_4 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(4*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(4*TempPoisson$Sol[,"Apr"]) #Apr
curve_4.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(4.5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(4.5*TempPoisson$Sol[,"Apr"]) #Apr
curve_5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(5*TempPoisson$Sol[,"Apr"]) #Apr
curve_5.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(5.5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(5.5*TempPoisson$Sol[,"Apr"]) #Apr
curve_6 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6*TempPoisson$Sol[,"Apr"]) #Apr
curve_6.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(6.5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(6.5*TempPoisson$Sol[,"Apr"]) #Apr
curve_7 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7*TempPoisson$Sol[,"Apr"]) #Apr
curve_7.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(7.5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(7.5*TempPoisson$Sol[,"Apr"]) #Apr
curve_8 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8*TempPoisson$Sol[,"Apr"]) #Apr
curve_8.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.5*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.5*TempPoisson$Sol[,"Apr"]) #Apr

plot(pred_day, exp(curve_4), type="l")
points(pred_day, exp(curve_4.5), type="l")
points(pred_day, exp(curve_5), type="l")
points(pred_day, exp(curve_5.5), type="l")
points(pred_day, exp(curve_6), type="l")
points(pred_day, exp(curve_6.5), type="l")
points(pred_day, exp(curve_7), type="l")
points(pred_day, exp(curve_7.5), type="l")
points(pred_day, exp(curve_8), type="l")
points(pred_day, exp(curve_8.5), type="l")


peak <- data.frame(temp=c(4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15))

for(i in 1:length(peak$temp)){
  peak$date[i] <- mean(-((TempPoisson$Sol[,"date"])+(peak$temp[i]*TempPoisson$Sol[,"date:Apr"]))/
                         (2*TempPoisson$Sol[,"I(date^2)"]))
}

for(i in 1:length(peak$temp)){
  peak$height[i] <- exp(mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                           mean(TempPoisson$Sol[,"date"])*peak$date[i] +  #date+yearchange
                           mean(TempPoisson$Sol[,"I(date^2)"])*peak$date[i]^2 +  #I(date^2)
                           mean(TempPoisson$Sol[,"date:Apr"])*peak$temp[i]*peak$date[i] +  #date*Apr
                           mean(TempPoisson$Sol[,"Apr"])*peak$temp[i]) #Apr
}

for(i in 1:length(peak$temp)){
  peak$logheight[i] <- (mean(TempPoisson$Sol[,"(Intercept)"]) +  #Intercept+yearchange
                          mean(TempPoisson$Sol[,"date"])*peak$date[i] +  #date+yearchange
                          mean(TempPoisson$Sol[,"I(date^2)"])*peak$date[i]^2 +  #I(date^2)
                          mean(TempPoisson$Sol[,"date:Apr"])*peak$temp[i]*peak$date[i] +  #date*Apr
                          mean(TempPoisson$Sol[,"Apr"])*peak$temp[i]) #Apr
}

plot(peak$temp, peak$date, type="l")
plot(peak$date, peak$height, type="l")
plot(peak$temp, peak$height, type="l")
plot(peak$date, peak$logheight, type="l")
plot(peak$temp, peak$logheight, type="l")

b1 <- mean(TempPoisson$Sol[,"date"])
b2 <- mean(TempPoisson$Sol[,"I(date^2)"])
b3 <- mean(TempPoisson$Sol[,"Apr"])
b4 <- mean(TempPoisson$Sol[,"date:Apr"])
int <- mean(TempPoisson$Sol[,"(Intercept)"])
A <- -b1/(2*b2)
B <- -b4/(2*b2)

for(i in 1:length(peak$temp)){
  peak$eqheight[i] <- (b1*A + b2*(A^2) + int + (b1*B + b2*2*A*B + b3 + b4*A)*peak$temp[i] + (b2*(B^2) + b4*B)*(peak$temp[i]^2))
}

plot(peak$temp, peak$logheight)
points(peak$temp, peak$eqheight, type="l")

plot(peak$temp, peak$height)
points(peak$temp, exp(peak$eqheight), type="l")

#### change in date but not height- can it deal with it? ####

curves <- data.frame(x=seq(-7,3,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y5  <- -curves$x^2 + -4*curves$x -2.75
curves$y4  <- -curves$x^2 + -5*curves$x -5
curves$y7  <- -curves$x^2 + -2*curves$x +0.25
curves$y6  <- -curves$x^2 + -3*curves$x -1
curves$y3  <- -curves$x^2 + -6*curves$x -7.75
curves$y2  <- -curves$x^2 + -7*curves$x -11
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y8  <- -curves$x^2 + -1*curves$x +1
curves$y9  <- -curves$x^2 + 0*curves$x +1.25
curves$y10 <- -curves$x^2 + 1*curves$x +1

# figure checking heights are uniform
peak <- 4/-2
height <- exp(-peak^2 + -4*peak -2.75)
curves$height <- height
plot(curves$x, exp(curves$y1), lwd=3,type="l",xlim=c(-7,3), ylim=c(0,5), col=1)
points(curves$x, curves$height, type="l", lty="dashed", col=2, lwd=3)
points(curves$x, exp(curves$y2), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y3), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y4), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y5), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y6), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y7), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y8), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y9), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y10), type="l", col=2, lwd=3)

#finalising data set
curves$date <- curves$x+7
curves$x=NULL

curveslong <- gather(data=curves, key="contvar", value="abundance", select=1:10)
curveslong$contvar <- gsub("y","", curveslong$contvar)
curveslong$expabund <- exp(curveslong$abundance)
curveslong$expabund <- round(curveslong$expabund, digits=2)
curveslong$contvar <- as.numeric(curveslong$contvar)
curveslong$countabund <- as.integer(curveslong$expabund*100)

#model 

k<-1000
prior<-list(R=list(V=1,nu=0.002))

HeightCheck1<- MCMCglmm(expabund~contvar*date+I(date^2)+I(contvar^2), family="gaussian", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck1, file = "~/Documents/Models/HeightCheck1.RData")
load("~/Documents/Models/HeightCheck1.RData")
summary(HeightCheck1)

-(mean(HeightCheck1$Sol[,6])^2)/mean(HeightCheck1$Sol[,4])-(mean(HeightCheck1$Sol[,6])^2)/mean(HeightCheck1$Sol[,4])

date <- seq(0,10,0.1)
contvar1 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*1) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*1^2) +
  mean(HeightCheck1$Sol[,6]*1)*date

contvar2 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*2) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*2^2) +
  mean(HeightCheck1$Sol[,6]*2)*date

contvar3 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*3) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*3^2) +
  mean(HeightCheck1$Sol[,6]*3)*date

contvar4 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*4) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*4^2) +
  mean(HeightCheck1$Sol[,6]*4)*date

contvar5 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*5) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*5^2) +
  mean(HeightCheck1$Sol[,6]*5)*date

contvar6 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*6) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*6^2) +
  mean(HeightCheck1$Sol[,6]*6)*date

contvar7 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*7) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*7^2) +
  mean(HeightCheck1$Sol[,6]*7)*date

contvar8 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*8) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*8^2) +
  mean(HeightCheck1$Sol[,6]*8)*date

contvar9 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*9) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*9^2) +
  mean(HeightCheck1$Sol[,6]*9)*date

contvar10 <- mean(HeightCheck1$Sol[,1]) +
  mean(HeightCheck1$Sol[,2]*10) +
  mean(HeightCheck1$Sol[,3])*date +
  mean(HeightCheck1$Sol[,4])*date^2 +
  mean(HeightCheck1$Sol[,5]*10^2) +
  mean(HeightCheck1$Sol[,6]*10)*date

plot(date, contvar1, type="l")
points(date, contvar2, type="l")
points(date, contvar3, type="l")
points(date, contvar3, type="l")
points(date, contvar4, type="l")
points(date, contvar5, type="l")
points(date, contvar6, type="l")
points(date, contvar7, type="l")
points(date, contvar8, type="l")
points(date, contvar9, type="l")
points(date, contvar10, type="l")

HeightCheck2<- MCMCglmm(expabund~contvar*date+I(date^2), family="gaussian", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck2, file = "~/Documents/Models/HeightCheck2.RData")
load("~/Documents/Models/HeightCheck2.RData")
summary(HeightCheck2)

date <- seq(0,10,0.1)
contvar1.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*1) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*1)*date

contvar2.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*2) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*2)*date

contvar3.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*3) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*3)*date

contvar4.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*4) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*4)*date

contvar5.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*5) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*5)*date

contvar6.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*6) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*6)*date

contvar7.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*7) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*7)*date

contvar8.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*8) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*8)*date

contvar9.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*9) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*9)*date

contvar10.2 <- mean(HeightCheck2$Sol[,1]) +
  mean(HeightCheck2$Sol[,2]*10) +
  mean(HeightCheck2$Sol[,3])*date +
  mean(HeightCheck2$Sol[,4])*date^2 +
  mean(HeightCheck2$Sol[,5]*10)*date

plot(date, contvar1, type="l", ylim=c(0,4))
points(date, contvar2.2, type="l")
points(date, contvar3.2, type="l")
points(date, contvar4.2, type="l")
points(date, contvar5.2, type="l")
points(date, contvar6.2, type="l")
points(date, contvar7.2, type="l")
points(date, contvar8.2, type="l")
points(date, contvar9.2, type="l")
points(date, contvar10.2, type="l")
points(date, contvar1, type="l", col=2)
points(date, contvar2, type="l", col=2)
points(date, contvar3, type="l", col=2)
points(date, contvar4, type="l", col=2)
points(date, contvar5, type="l", col=2)
points(date, contvar6, type="l", col=2)
points(date, contvar7, type="l", col=2)
points(date, contvar8, type="l", col=2)
points(date, contvar9, type="l", col=2)
points(date, contvar10, type="l", col=2)
points(curves$date, exp(curves$y1), type="l", col=3)
points(curves$date, exp(curves$y2), type="l", col=3)
points(curves$date, exp(curves$y3), type="l", col=3)
points(curves$date, exp(curves$y4), type="l", col=3)
points(curves$date, exp(curves$y5), type="l", col=3)
points(curves$date, exp(curves$y6), type="l", col=3)
points(curves$date, exp(curves$y7), type="l", col=3)
points(curves$date, exp(curves$y8), type="l", col=3)
points(curves$date, exp(curves$y9), type="l", col=3)
points(curves$date, exp(curves$y10), type="l", col=3)

intHC1 <- mean(HeightCheck1$Sol[,1])
b1HC1 <- mean(HeightCheck1$Sol[,3])
b2HC1 <- mean(HeightCheck1$Sol[,4])
b3HC1 <- mean(HeightCheck1$Sol[,2])
b4HC1 <- mean(HeightCheck1$Sol[,6])
b5HC1 <- mean(HeightCheck1$Sol[,5])
AHC1 <- mean(-HeightCheck1$Sol[,3]/(2*HeightCheck1$Sol[,4]))
BHC1 <- mean(-HeightCheck1$Sol[,6]/(2*HeightCheck1$Sol[,4]))

intHC2 <- mean(HeightCheck2$Sol[,1])
b1HC2 <- mean(HeightCheck2$Sol[,3])
b2HC2 <- mean(HeightCheck2$Sol[,4])
b3HC2 <- mean(HeightCheck2$Sol[,2])
b4HC2 <- mean(HeightCheck2$Sol[,5])
AHC2 <- mean(-HeightCheck2$Sol[,3]/(2*HeightCheck2$Sol[,4]))
BHC2 <- mean(-HeightCheck2$Sol[,5]/(2*HeightCheck2$Sol[,4]))

x <- seq(-200,205,0.1)
HC1 <- b1HC1*AHC1 + b2HC1*(AHC1^2) + intHC1 + (b1HC1*BHC1 + b2HC1*2*AHC1*BHC1 + b3HC1 + b4HC1*AHC1)*x + (b2HC1*(BHC1^2) + b4HC1*BHC1 + b5HC1)*(x^2)
HC2 <- b1HC2*AHC2 + b2HC2*(AHC2^2) + intHC2 + (b1HC2*BHC2 + b2HC2*2*AHC2*BHC2 + b3HC2 + b4HC2*AHC2)*x + (b2HC2*(BHC2^2) + b4HC2*BHC2)*(x^2)
plot(x, HC1, type="l")
points(x, HC2, type="l", col=2)


#### spread heights out further ####

curves <- data.frame(x=seq(-7,13,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y2  <- -curves$x^2 + -5*curves$x -5
curves$y3  <- -curves$x^2 + -2*curves$x +0.25
curves$y4  <- -curves$x^2 + 1*curves$x +1
curves$y5  <- -curves$x^2 + 4*curves$x -2.75
curves$y6  <- -curves$x^2 + 7*curves$x -11
curves$y7  <- -curves$x^2 + 10*curves$x -23.75
curves$y8  <- -curves$x^2 + 13*curves$x -41
curves$y9  <- -curves$x^2 + 16*curves$x -62.75
curves$y10 <- -curves$x^2 + 19*curves$x -89

# figure checking heights are uniform
peak <- 4/-2
height <- exp(-peak^2 + -4*peak -2.75)
curves$height <- height
plot(curves$x, exp(curves$y1), lwd=3,type="l",xlim=c(-7,13), ylim=c(0,5), col=1)
points(curves$x, curves$height, type="l", lty="dashed", col=2, lwd=3)
points(curves$x, exp(curves$y2), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y3), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y4), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y5), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y6), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y7), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y8), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y9), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y10), type="l", col=1, lwd=3)

curves$date <- curves$x+7
curves$x=NULL

curveslong <- gather(data=curves, key="contvar", value="abundance", select=1:10)
curveslong$contvar <- gsub("y","", curveslong$contvar)
curveslong$expabund <- exp(curveslong$abundance)
curveslong$contvar <- as.numeric(curveslong$contvar)

#model 

k<-1000
prior<-list(R=list(V=1,nu=0.002))

HeightCheck3<- MCMCglmm(expabund~contvar*date+I(date^2)+I(contvar^2), family="gaussian", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck3, file = "~/Documents/Models/HeightCheck3.RData")
load("~/Documents/Models/HeightCheck3.RData")
summary(HeightCheck3)

HeightCheck4<- MCMCglmm(expabund~contvar*date+I(date^2), family="gaussian", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck4, file = "~/Documents/Models/HeightCheck4.RData")
load("~/Documents/Models/HeightCheck4.RData")
summary(HeightCheck4)


date <- seq(0,20,0.1)
contvar1 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*1) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*1^2) +
  mean(HeightCheck3$Sol[,6]*1)*date

contvar2 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*2) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*2^2) +
  mean(HeightCheck3$Sol[,6]*2)*date

contvar3 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*3) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*3^2) +
  mean(HeightCheck3$Sol[,6]*3)*date

contvar4 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*4) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*4^2) +
  mean(HeightCheck3$Sol[,6]*4)*date

contvar5 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*5) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*5^2) +
  mean(HeightCheck3$Sol[,6]*5)*date

contvar6 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*6) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*6^2) +
  mean(HeightCheck3$Sol[,6]*6)*date

contvar7 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*7) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*7^2) +
  mean(HeightCheck3$Sol[,6]*7)*date

contvar8 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*8) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*8^2) +
  mean(HeightCheck3$Sol[,6]*8)*date

contvar9 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*9) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*9^2) +
  mean(HeightCheck3$Sol[,6]*9)*date

contvar10 <- mean(HeightCheck3$Sol[,1]) +
  mean(HeightCheck3$Sol[,2]*10) +
  mean(HeightCheck3$Sol[,3])*date +
  mean(HeightCheck3$Sol[,4])*date^2 +
  mean(HeightCheck3$Sol[,5]*10^2) +
  mean(HeightCheck3$Sol[,6]*10)*date

contvar1.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*1) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*1)*date

contvar2.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*2) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*2)*date

contvar3.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*3) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*3)*date

contvar4.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*4) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*4)*date

contvar5.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*5) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*5)*date

contvar6.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*6) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*6)*date

contvar7.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*7) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*7)*date

contvar8.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*8) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*8)*date

contvar9.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*9) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*9)*date

contvar10.2 <- mean(HeightCheck4$Sol[,1]) +
  mean(HeightCheck4$Sol[,2]*10) +
  mean(HeightCheck4$Sol[,3])*date +
  mean(HeightCheck4$Sol[,4])*date^2 +
  mean(HeightCheck4$Sol[,5]*10)*date

plot(date, contvar1.2, type="l", ylim=c(0,4))
points(date, contvar2.2, type="l")
points(date, contvar3.2, type="l")
points(date, contvar4.2, type="l")
points(date, contvar5.2, type="l")
points(date, contvar6.2, type="l")
points(date, contvar7.2, type="l")
points(date, contvar8.2, type="l")
points(date, contvar9.2, type="l")
points(date, contvar10.2, type="l")
points(date, contvar1, type="l", col=2)
points(date, contvar2, type="l", col=2)
points(date, contvar3, type="l", col=2)
points(date, contvar4, type="l", col=2)
points(date, contvar5, type="l", col=2)
points(date, contvar6, type="l", col=2)
points(date, contvar7, type="l", col=2)
points(date, contvar8, type="l", col=2)
points(date, contvar9, type="l", col=2)
points(date, contvar10, type="l", col=2)
points(curves$date, exp(curves$y1), type="l", col=3)
points(curves$date, exp(curves$y2), type="l", col=3)
points(curves$date, exp(curves$y3), type="l", col=3)
points(curves$date, exp(curves$y4), type="l", col=3)
points(curves$date, exp(curves$y5), type="l", col=3)
points(curves$date, exp(curves$y6), type="l", col=3)
points(curves$date, exp(curves$y7), type="l", col=3)
points(curves$date, exp(curves$y8), type="l", col=3)
points(curves$date, exp(curves$y9), type="l", col=3)
points(curves$date, exp(curves$y10), type="l", col=3)

intHC3 <- mean(HeightCheck3$Sol[,1])
b1HC3 <- mean(HeightCheck3$Sol[,3])
b2HC3 <- mean(HeightCheck3$Sol[,4])
b3HC3 <- mean(HeightCheck3$Sol[,2])
b4HC3 <- mean(HeightCheck3$Sol[,6])
b5HC3 <- mean(HeightCheck3$Sol[,5])
AHC3 <- mean(-HeightCheck3$Sol[,3]/(2*HeightCheck3$Sol[,4]))
BHC3 <- mean(-HeightCheck3$Sol[,6]/(2*HeightCheck3$Sol[,4]))

intHC4 <- mean(HeightCheck4$Sol[,1])
b1HC4 <- mean(HeightCheck4$Sol[,3])
b2HC4 <- mean(HeightCheck4$Sol[,4])
b3HC4 <- mean(HeightCheck4$Sol[,2])
b4HC4 <- mean(HeightCheck4$Sol[,5])
AHC4 <- mean(-HeightCheck4$Sol[,3]/(2*HeightCheck4$Sol[,4]))
BHC4 <- mean(-HeightCheck4$Sol[,5]/(2*HeightCheck4$Sol[,4]))

x <- seq(-200,205,0.1)
HC3 <- b1HC3*AHC3 + b2HC3*(AHC3^2) + intHC3 + (b1HC3*BHC3 + b2HC3*2*AHC3*BHC3 + b3HC3 + b4HC3*AHC3)*x + (b2HC3*(BHC3^2) + b4HC3*BHC3 + b5HC3)*(x^2)
HC4 <- b1HC4*AHC4 + b2HC4*(AHC4^2) + intHC4 + (b1HC4*BHC4 + b2HC4*2*AHC4*BHC4 + b3HC4 + b4HC4*AHC4)*x + (b2HC4*(BHC4^2) + b4HC4*BHC4)*(x^2)
plot(x, HC3, type="l")#, xlim=c(0,20), ylim=c(0,4))
points(x, HC4, type="l", col=2)
points(x, HC1, type="l", col=3)


### Poisson model ####

prior<-list(R=list(V=1,nu=0.002))

HeightCheck5<- MCMCglmm(countabund~contvar*date+I(date^2)+I(contvar^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck5, file = "~/Documents/Models/HeightCheck5.RData")
load("~/Documents/Models/HeightCheck5.RData")
summary(HeightCheck5)

HeightCheck6<- MCMCglmm(countabund~contvar*date+I(date^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck6, file = "~/Documents/Models/HeightCheck6.RData")
load("~/Documents/Models/HeightCheck6.RData")
summary(HeightCheck6)

date <- seq(0,10,0.1)
contvar1 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*1) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*1^2) +
  mean(HeightCheck5$Sol[,6]*1)*date

contvar2 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*2) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*2^2) +
  mean(HeightCheck5$Sol[,6]*2)*date

contvar3 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*3) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*3^2) +
  mean(HeightCheck5$Sol[,6]*3)*date

contvar4 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*4) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*4^2) +
  mean(HeightCheck5$Sol[,6]*4)*date

contvar5 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*5) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*5^2) +
  mean(HeightCheck5$Sol[,6]*5)*date

contvar6 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*6) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*6^2) +
  mean(HeightCheck5$Sol[,6]*6)*date

contvar7 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*7) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*7^2) +
  mean(HeightCheck5$Sol[,6]*7)*date

contvar8 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*8) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*8^2) +
  mean(HeightCheck5$Sol[,6]*8)*date

contvar9 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*9) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*9^2) +
  mean(HeightCheck5$Sol[,6]*9)*date

contvar10 <- mean(HeightCheck5$Sol[,1]) +
  mean(HeightCheck5$Sol[,2]*10) +
  mean(HeightCheck5$Sol[,3])*date +
  mean(HeightCheck5$Sol[,4])*date^2 +
  mean(HeightCheck5$Sol[,5]*10^2) +
  mean(HeightCheck5$Sol[,6]*10)*date

plot(date, contvar1, type="l")
points(date, contvar2, type="l")
points(date, contvar3, type="l")
points(date, contvar3, type="l")
points(date, contvar4, type="l")
points(date, contvar5, type="l")
points(date, contvar6, type="l")
points(date, contvar7, type="l")
points(date, contvar8, type="l")
points(date, contvar9, type="l")
points(date, contvar10, type="l")


contvar1.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*1) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*1)*date

contvar2.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*2) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*2)*date

contvar3.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*3) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*3)*date

contvar4.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*4) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*4)*date

contvar5.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*5) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*5)*date

contvar6.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*6) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*6)*date

contvar7.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*7) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*7)*date

contvar8.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*8) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*8)*date

contvar9.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*9) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*9)*date

contvar10.2 <- mean(HeightCheck6$Sol[,1]) +
  mean(HeightCheck6$Sol[,2]*10) +
  mean(HeightCheck6$Sol[,3])*date +
  mean(HeightCheck6$Sol[,4])*date^2 +
  mean(HeightCheck6$Sol[,5]*10)*date

plot(date, exp(contvar1.2), type="l", ylim=c(0,400))
points(date, exp(contvar2.2), type="l")
points(date, exp(contvar3.2), type="l")
points(date, exp(contvar4.2), type="l")
points(date, exp(contvar5.2), type="l")
points(date, exp(contvar6.2), type="l")
points(date, exp(contvar7.2), type="l")
points(date, exp(contvar8.2), type="l")
points(date, exp(contvar9.2), type="l")
points(date, exp(contvar10.2), type="l")
points(date, exp(contvar1), type="l", col=2)
points(date, exp(contvar2), type="l", col=2)
points(date, exp(contvar3), type="l", col=2)
points(date, exp(contvar4), type="l", col=2)
points(date, exp(contvar5), type="l", col=2)
points(date, exp(contvar6), type="l", col=2)
points(date, exp(contvar7), type="l", col=2)
points(date, exp(contvar8), type="l", col=2)
points(date, exp(contvar9), type="l", col=2)
points(date, exp(contvar10), type="l", col=2)
points(curves$date, exp(curves$y1)*100, type="l", col=3)
points(curves$date, exp(curves$y2)*100, type="l", col=3)
points(curves$date, exp(curves$y3)*100, type="l", col=3)
points(curves$date, exp(curves$y4)*100, type="l", col=3)
points(curves$date, exp(curves$y5)*100, type="l", col=3)
points(curves$date, exp(curves$y6)*100, type="l", col=3)
points(curves$date, exp(curves$y7)*100, type="l", col=3)
points(curves$date, exp(curves$y8)*100, type="l", col=3)
points(curves$date, exp(curves$y9)*100, type="l", col=3)
points(curves$date, exp(curves$y10)*100, type="l", col=3)
points(x, exp(HC5), type="l")

intHC5 <- mean(HeightCheck5$Sol[,1])
b1HC5 <- mean(HeightCheck5$Sol[,3])
b2HC5 <- mean(HeightCheck5$Sol[,4])
b3HC5 <- mean(HeightCheck5$Sol[,2])
b4HC5 <- mean(HeightCheck5$Sol[,6])
b5HC5 <- mean(HeightCheck5$Sol[,5])
AHC5 <- mean(-HeightCheck5$Sol[,3]/(2*HeightCheck5$Sol[,4]))
BHC5 <- mean(-HeightCheck5$Sol[,6]/(2*HeightCheck5$Sol[,4]))

intHC6 <- mean(HeightCheck6$Sol[,1])
b1HC6 <- mean(HeightCheck6$Sol[,3])
b2HC6 <- mean(HeightCheck6$Sol[,4])
b3HC6 <- mean(HeightCheck6$Sol[,2])
b4HC6 <- mean(HeightCheck6$Sol[,5])
AHC6 <- mean(-HeightCheck6$Sol[,3]/(2*HeightCheck6$Sol[,4]))
BHC6 <- mean(-HeightCheck6$Sol[,5]/(2*HeightCheck6$Sol[,4]))

x <- seq(-200,205,0.1)
HC5 <- b1HC5*AHC5 + b2HC5*(AHC5^2) + intHC5 + (b1HC5*BHC5 + b2HC5*2*AHC5*BHC5 + b3HC5 + b4HC5*AHC5)*x + (b2HC5*(BHC5^2) + b4HC5*BHC5 + b5HC5)*(x^2)
HC6 <- b1HC6*AHC6 + b2HC6*(AHC6^2) + intHC6 + (b1HC6*BHC6 + b2HC6*2*AHC6*BHC6 + b3HC6 + b4HC6*AHC6)*x + (b2HC6*(BHC6^2) + b4HC6*BHC6)*(x^2)
plot(x, exp(HC5), type="l", xlim=c(-50,50))
plot(x, exp(HC6), type="l", col=2)

#################################
### Peak height linear slope ####

curves <- data.frame(x=seq(-7,3,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y2  <- -curves$x^2 + -7*curves$x -10.87
curves$y3  <- -curves$x^2 + -6*curves$x -7.5
curves$y4  <- -curves$x^2 + -5*curves$x -4.65
curves$y5  <- -curves$x^2 + -4*curves$x -2.3
curves$y6  <- -curves$x^2 + -3*curves$x -0.46
curves$y7  <- -curves$x^2 + -2*curves$x +0.86
curves$y8  <- -curves$x^2 + -1*curves$x +1.69
curves$y9  <- -curves$x^2 + 0*curves$x +2.01
curves$y10 <- -curves$x^2 + 1*curves$x +1.83

# figure checking heights are uniform
peak <- 4/-2
height <- exp(-peak^2 + -4*peak -2.75)
curves$height <- curves$x+7.5
plot(curves$x, exp(curves$y1), lwd=3,type="l",xlim=c(-7,3), ylim=c(0,15), col=1)
points(curves$x, curves$height, type="l", col=2, lwd=3)
points(curves$x, exp(curves$y2), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y3), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y4), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y5), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y6), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y7), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y8), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y9), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y10), type="l", col=1, lwd=3)

#finalising data set
curves$date <- curves$x+7
curves$x=NULL

curveslong <- gather(data=curves, key="contvar", value="abundance", select=1:10)
curveslong$contvar <- gsub("y","", curveslong$contvar)
curveslong$expabund <- exp(curveslong$abundance)
curveslong$expabund <- round(curveslong$expabund, digits=2)
curveslong$contvar <- as.numeric(curveslong$contvar)
curveslong$countabund <- as.integer(curveslong$expabund*100)

### Poisson model ####

prior<-list(R=list(V=1,nu=0.002))

HeightCheck7<- MCMCglmm(countabund~contvar*date+I(date^2)+I(contvar^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck7, file = "~/Documents/Models/HeightCheck7.RData")
load("~/Documents/Models/HeightCheck7.RData")
summary(HeightCheck7)

HeightCheck8<- MCMCglmm(countabund~contvar*date+I(date^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck8, file = "~/Documents/Models/HeightCheck8.RData")
load("~/Documents/Models/HeightCheck8.RData")
summary(HeightCheck8)

date <- seq(0,10,0.1)
contvar1 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*1) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*1^2) +
  mean(HeightCheck7$Sol[,6]*1)*date

contvar2 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*2) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*2^2) +
  mean(HeightCheck7$Sol[,6]*2)*date

contvar3 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*3) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*3^2) +
  mean(HeightCheck7$Sol[,6]*3)*date

contvar4 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*4) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*4^2) +
  mean(HeightCheck7$Sol[,6]*4)*date

contvar5 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*5) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*5^2) +
  mean(HeightCheck7$Sol[,6]*5)*date

contvar6 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*6) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*6^2) +
  mean(HeightCheck7$Sol[,6]*6)*date

contvar7 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*7) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*7^2) +
  mean(HeightCheck7$Sol[,6]*7)*date

contvar8 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*8) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*8^2) +
  mean(HeightCheck7$Sol[,6]*8)*date

contvar9 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*9) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*9^2) +
  mean(HeightCheck7$Sol[,6]*9)*date

contvar10 <- mean(HeightCheck7$Sol[,1]) +
  mean(HeightCheck7$Sol[,2]*10) +
  mean(HeightCheck7$Sol[,3])*date +
  mean(HeightCheck7$Sol[,4])*date^2 +
  mean(HeightCheck7$Sol[,5]*10^2) +
  mean(HeightCheck7$Sol[,6]*10)*date

plot(date, contvar1, type="l")
points(date, contvar2, type="l")
points(date, contvar3, type="l")
points(date, contvar3, type="l")
points(date, contvar4, type="l")
points(date, contvar5, type="l")
points(date, contvar6, type="l")
points(date, contvar7, type="l")
points(date, contvar8, type="l")
points(date, contvar9, type="l")
points(date, contvar10, type="l")


contvar1.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*1) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*1)*date

contvar2.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*2) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*2)*date

contvar3.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*3) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*3)*date

contvar4.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*4) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*4)*date

contvar5.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*5) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*5)*date

contvar6.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*6) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*6)*date

contvar7.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*7) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*7)*date

contvar8.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*8) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*8)*date

contvar9.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*9) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*9)*date

contvar10.2 <- mean(HeightCheck8$Sol[,1]) +
  mean(HeightCheck8$Sol[,2]*10) +
  mean(HeightCheck8$Sol[,3])*date +
  mean(HeightCheck8$Sol[,4])*date^2 +
  mean(HeightCheck8$Sol[,5]*10)*date

plot(date, exp(contvar1.2), type="l", ylim=c(0,900))
points(date, exp(contvar2.2), type="l")
points(date, exp(contvar3.2), type="l")
points(date, exp(contvar4.2), type="l")
points(date, exp(contvar5.2), type="l")
points(date, exp(contvar6.2), type="l")
points(date, exp(contvar7.2), type="l")
points(date, exp(contvar8.2), type="l")
points(date, exp(contvar9.2), type="l")
points(date, exp(contvar10.2), type="l")
points(date, exp(contvar1), type="l", col=2)
points(date, exp(contvar2), type="l", col=2)
points(date, exp(contvar3), type="l", col=2)
points(date, exp(contvar4), type="l", col=2)
points(date, exp(contvar5), type="l", col=2)
points(date, exp(contvar6), type="l", col=2)
points(date, exp(contvar7), type="l", col=2)
points(date, exp(contvar8), type="l", col=2)
points(date, exp(contvar9), type="l", col=2)
points(date, exp(contvar10), type="l", col=2)
points(curves$date, exp(curves$y1)*100, type="l", col=3)
points(curves$date, exp(curves$y2)*100, type="l", col=3)
points(curves$date, exp(curves$y3)*100, type="l", col=3)
points(curves$date, exp(curves$y4)*100, type="l", col=3)
points(curves$date, exp(curves$y5)*100, type="l", col=3)
points(curves$date, exp(curves$y6)*100, type="l", col=3)
points(curves$date, exp(curves$y7)*100, type="l", col=3)
points(curves$date, exp(curves$y8)*100, type="l", col=3)
points(curves$date, exp(curves$y9)*100, type="l", col=3)
points(curves$date, exp(curves$y10)*100, type="l", col=3)
points(x, exp(HC7), type="l", lwd=2)

intHC7 <- mean(HeightCheck7$Sol[,1])
b1HC7 <- mean(HeightCheck7$Sol[,3])
b2HC7 <- mean(HeightCheck7$Sol[,4])
b3HC7 <- mean(HeightCheck7$Sol[,2])
b4HC7 <- mean(HeightCheck7$Sol[,6])
b5HC7 <- mean(HeightCheck7$Sol[,5])
AHC7 <- mean(-HeightCheck7$Sol[,3]/(2*HeightCheck7$Sol[,4]))
BHC7 <- mean(-HeightCheck7$Sol[,6]/(2*HeightCheck7$Sol[,4]))

intHC8 <- mean(HeightCheck8$Sol[,1])
b1HC8 <- mean(HeightCheck8$Sol[,3])
b2HC8 <- mean(HeightCheck8$Sol[,4])
b3HC8 <- mean(HeightCheck8$Sol[,2])
b4HC8 <- mean(HeightCheck8$Sol[,5])
AHC8 <- mean(-HeightCheck8$Sol[,3]/(2*HeightCheck8$Sol[,4]))
BHC8 <- mean(-HeightCheck8$Sol[,5]/(2*HeightCheck8$Sol[,4]))

x <- seq(-200,205,0.1)
HC7 <- b1HC7*AHC7 + b2HC7*(AHC7^2) + intHC7 + (b1HC7*BHC7 + b2HC7*2*AHC7*BHC7 + b3HC7 + b4HC7*AHC7)*x + (b2HC7*(BHC7^2) + b4HC7*BHC7 + b5HC7)*(x^2)
HC8 <- b1HC8*AHC8 + b2HC8*(AHC8^2) + intHC8 + (b1HC8*BHC8 + b2HC8*2*AHC8*BHC8 + b3HC8 + b4HC8*AHC8)*x + (b2HC8*(BHC8^2) + b4HC8*BHC8)*(x^2)
plot(x, HC7, type="l", lwd=2, xlim=c(0,50), ylim=c(0,300))
plot(x, HC8, type="l", col=2)


#### spreading out the data points ####
curves <- data.frame(x=seq(-7,13,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y2  <- -curves$x^2 + -5*curves$x -4.65
curves$y3  <- -curves$x^2 + -2*curves$x +0.865
curves$y4  <- -curves$x^2 + 1*curves$x +1.83
curves$y5  <- -curves$x^2 + 4*curves$x -1.75
curves$y6  <- -curves$x^2 + 7*curves$x -9.85
curves$y7  <- -curves$x^2 + 10*curves$x -22.475
curves$y8  <- -curves$x^2 + 13*curves$x -39.61
curves$y9  <- -curves$x^2 + 16*curves$x -61.26
curves$y10 <- -curves$x^2 + 19*curves$x -87.415


# figure checking heights are uniform
peak <- 4/-2
height <- exp(-peak^2 + -4*peak -2.75)
curves$height <- curves$x+7.5
plot(curves$x, exp(curves$y1), lwd=3,type="l",xlim=c(-7,13), ylim=c(0,20), col=1)
points(curves$x, curves$height, type="l", col=2, lwd=3)
points(curves$x, exp(curves$y2), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y3), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y4), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y5), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y6), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y7), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y8), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y9), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y10), type="l", col=1, lwd=3)

#finalising data set
curves$date <- curves$x+7
curves$x=NULL

curveslong <- gather(data=curves, key="contvar", value="abundance", select=1:10)
curveslong$contvar <- gsub("y","", curveslong$contvar)
curveslong$expabund <- exp(curveslong$abundance)
curveslong$expabund <- round(curveslong$expabund, digits=2)
curveslong$contvar <- as.numeric(curveslong$contvar)
curveslong$countabund <- as.integer(curveslong$expabund*100)

prior<-list(R=list(V=1,nu=0.002))

HeightCheck9<- MCMCglmm(countabund~contvar*date+I(date^2)+I(contvar^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck9, file = "~/Documents/Models/HeightCheck9.RData")
load("~/Documents/Models/HeightCheck9.RData")
summary(HeightCheck9)

HeightCheck10<- MCMCglmm(countabund~contvar*date+I(date^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck10, file = "~/Documents/Models/HeightCheck10.RData")
load("~/Documents/Models/HeightCheck10.RData")
summary(HeightCheck10)

date <- seq(0,20,0.1)
contvar1 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*1) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*1^2) +
  mean(HeightCheck9$Sol[,6]*1)*date

contvar2 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*2) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*2^2) +
  mean(HeightCheck9$Sol[,6]*2)*date

contvar3 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*3) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*3^2) +
  mean(HeightCheck9$Sol[,6]*3)*date

contvar4 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*4) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*4^2) +
  mean(HeightCheck9$Sol[,6]*4)*date

contvar5 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*5) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*5^2) +
  mean(HeightCheck9$Sol[,6]*5)*date

contvar6 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*6) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*6^2) +
  mean(HeightCheck9$Sol[,6]*6)*date

contvar7 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*7) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*7^2) +
  mean(HeightCheck9$Sol[,6]*7)*date

contvar8 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*8) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*8^2) +
  mean(HeightCheck9$Sol[,6]*8)*date

contvar9 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*9) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*9^2) +
  mean(HeightCheck9$Sol[,6]*9)*date

contvar10 <- mean(HeightCheck9$Sol[,1]) +
  mean(HeightCheck9$Sol[,2]*10) +
  mean(HeightCheck9$Sol[,3])*date +
  mean(HeightCheck9$Sol[,4])*date^2 +
  mean(HeightCheck9$Sol[,5]*10^2) +
  mean(HeightCheck9$Sol[,6]*10)*date

plot(date, contvar1, type="l")
points(date, contvar2, type="l")
points(date, contvar3, type="l")
points(date, contvar3, type="l")
points(date, contvar4, type="l")
points(date, contvar5, type="l")
points(date, contvar6, type="l")
points(date, contvar7, type="l")
points(date, contvar8, type="l")
points(date, contvar9, type="l")
points(date, contvar10, type="l")


contvar1.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*1) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*1)*date

contvar2.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*2) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*2)*date

contvar3.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*3) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*3)*date

contvar4.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*4) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*4)*date

contvar5.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*5) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*5)*date

contvar6.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*6) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*6)*date

contvar7.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*7) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*7)*date

contvar8.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*8) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*8)*date

contvar9.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*9) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*9)*date

contvar10.2 <- mean(HeightCheck10$Sol[,1]) +
  mean(HeightCheck10$Sol[,2]*10) +
  mean(HeightCheck10$Sol[,3])*date +
  mean(HeightCheck10$Sol[,4])*date^2 +
  mean(HeightCheck10$Sol[,5]*10)*date

plot(date, exp(contvar1.2), type="l", ylim=c(0,1500))
points(date, exp(contvar2.2), type="l")
points(date, exp(contvar3.2), type="l")
points(date, exp(contvar4.2), type="l")
points(date, exp(contvar5.2), type="l")
points(date, exp(contvar6.2), type="l")
points(date, exp(contvar7.2), type="l")
points(date, exp(contvar8.2), type="l")
points(date, exp(contvar9.2), type="l")
points(date, exp(contvar10.2), type="l")
points(date, exp(contvar1), type="l", col=2)
points(date, exp(contvar2), type="l", col=2)
points(date, exp(contvar3), type="l", col=2)
points(date, exp(contvar4), type="l", col=2)
points(date, exp(contvar5), type="l", col=2)
points(date, exp(contvar6), type="l", col=2)
points(date, exp(contvar7), type="l", col=2)
points(date, exp(contvar8), type="l", col=2)
points(date, exp(contvar9), type="l", col=2)
points(date, exp(contvar10), type="l", col=2)
points(curves$date, exp(curves$y1)*100, type="l", col=3)
points(curves$date, exp(curves$y2)*100, type="l", col=3)
points(curves$date, exp(curves$y3)*100, type="l", col=3)
points(curves$date, exp(curves$y4)*100, type="l", col=3)
points(curves$date, exp(curves$y5)*100, type="l", col=3)
points(curves$date, exp(curves$y6)*100, type="l", col=3)
points(curves$date, exp(curves$y7)*100, type="l", col=3)
points(curves$date, exp(curves$y8)*100, type="l", col=3)
points(curves$date, exp(curves$y9)*100, type="l", col=3)
points(curves$date, exp(curves$y10)*100, type="l", col=3)
#points(x, exp(HC9), type="l", lwd=2)
#points(x, exp(HC10), type="l", lwd=2) these arent right- x axis is date not variable

intHC9 <- mean(HeightCheck9$Sol[,1])
b1HC9 <- mean(HeightCheck9$Sol[,3])
b2HC9 <- mean(HeightCheck9$Sol[,4])
b3HC9 <- mean(HeightCheck9$Sol[,2])
b4HC9 <- mean(HeightCheck9$Sol[,6])
b5HC9 <- mean(HeightCheck9$Sol[,5])
AHC9 <- mean(-HeightCheck9$Sol[,3]/(2*HeightCheck9$Sol[,4]))
BHC9 <- mean(-HeightCheck9$Sol[,6]/(2*HeightCheck9$Sol[,4]))

intHC10 <- mean(HeightCheck10$Sol[,1])
b1HC10 <- mean(HeightCheck10$Sol[,3])
b2HC10 <- mean(HeightCheck10$Sol[,4])
b3HC10 <- mean(HeightCheck10$Sol[,2])
b4HC10 <- mean(HeightCheck10$Sol[,5])
AHC10 <- mean(-HeightCheck10$Sol[,3]/(2*HeightCheck10$Sol[,4]))
BHC10 <- mean(-HeightCheck10$Sol[,5]/(2*HeightCheck10$Sol[,4]))

x <- seq(-200,205,0.1)
HC9 <- b1HC9*AHC9 + b2HC9*(AHC9^2) + intHC9 + (b1HC9*BHC9 + b2HC9*2*AHC9*BHC9 + b3HC9 + b4HC9*AHC9)*x + (b2HC9*(BHC9^2) + b4HC9*BHC9 + b5HC9)*(x^2)
HC10 <- b1HC10*AHC10 + b2HC10*(AHC10^2) + intHC10 + (b1HC10*BHC10 + b2HC10*2*AHC10*BHC10 + b3HC10 + b4HC10*AHC10)*x + (b2HC10*(BHC10^2) + b4HC10*BHC10)*(x^2)
plot(x, exp(HC9), type="l", lwd=2, xlim=c(-30,40), ylim=c(0,900))
plot(x, HC10, type="l", col=2)fv

#############################################
#### Height consistent, increasing width ####
curves <- data.frame(x=seq(-7,17,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y2  <- -0.9*curves$x^2 +-4.5*curves$x - 4.375# -5*curves$x -5
curves$y3  <- -0.8*curves$x^2 +-1.6*curves$x + 0.45# -2*curves$x +0.25
curves$y4  <- -0.7*curves$x^2 +0.75*curves$x + 1.05# 1*curves$x +1
curves$y5  <- -0.6*curves$x^2 +2.357*curves$x-1.065# 4*curves$x -2.75
curves$y6  <- -0.5*curves$x^2 +3.5*curves$x -4.875# 7*curves$x -11
curves$y7  <- -0.4*curves$x^2 +4*curves$x -8.75# 10*curves$x -23.75
curves$y8  <- -0.3*curves$x^2 +3.9*curves$x -11.425# 13*curves$x -41
curves$y9  <- -0.2*curves$x^2 +3.2*curves$x-11.55# 16*curves$x -62.75
curves$y10 <- -0.1*curves$x^2 +1.9*curves$x -7.775# 19*curves$x -89

# figure checking heights are uniform
peak <- -4/-2
height <- exp(-peak^2 + -4*peak -2.75)
curves$height <- height
plot(curves$x, curves$height, lwd=3,type="l",lty="dashed", xlim=c(-7,15), ylim=c(0,5), col=2)


points(curves$x, exp(curves$y1), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y2), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y3), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y4), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y5), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y6), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y7), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y8), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y9), type="l", col=1, lwd=3)
points(curves$x, exp(curves$y10), type="l", col=1, lwd=3)

curves$date <- curves$x+7
curves$x=NULL

curveslong <- gather(data=curves, key="contvar", value="abundance", select=1:10)
curveslong$contvar <- gsub("y","", curveslong$contvar)
curveslong$expabund <- exp(curveslong$abundance)
curveslong$expabund <- round(curveslong$expabund, digits=2)
curveslong$contvar <- as.numeric(curveslong$contvar)
curveslong$countabund <- as.integer(curveslong$expabund*100)

prior<-list(R=list(V=1,nu=0.002))

HeightCheck11<- MCMCglmm(countabund~contvar*date+I(date^2)*contvar+I(contvar^2)+I(contvar^3), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck11, file = "~/Documents/Models/HeightCheck11.RData")
load("~/Documents/Models/HeightCheck11.RData")
summary(HeightCheck11)

HeightCheck12<- MCMCglmm(countabund~contvar*date+I(date^2)*contvar+I(contvar^2), family="poisson", data=curveslong, prior=prior, nitt=200000, burnin=20000)
save(HeightCheck12, file = "~/Documents/Models/HeightCheck12.RData")
load("~/Documents/Models/HeightCheck12.RData")
summary(HeightCheck12)


date <- seq(0,24,0.1)
contvar1 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*1) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*1^2) +
  mean(HeightCheck11$Sol[,6]*1^3) +
  mean(HeightCheck11$Sol[,7]*1)*date +
  mean(HeightCheck11$Sol[,8]*1)*date^2

contvar2 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*2) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*2^2) +
  mean(HeightCheck11$Sol[,6]*2^3) +
  mean(HeightCheck11$Sol[,7]*2)*date +
  mean(HeightCheck11$Sol[,8]*2)*date^2

contvar3 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*3) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*3^2) +
  mean(HeightCheck11$Sol[,6]*3^3) +
  mean(HeightCheck11$Sol[,7]*3)*date +
  mean(HeightCheck11$Sol[,8]*3)*date^2

contvar4 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*4) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*4^2) +
  mean(HeightCheck11$Sol[,6]*4^3) +
  mean(HeightCheck11$Sol[,7]*4)*date +
  mean(HeightCheck11$Sol[,8]*4)*date^2

contvar5 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*5) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*5^2) +
  mean(HeightCheck11$Sol[,6]*5^3) +
  mean(HeightCheck11$Sol[,7]*5)*date  +
  mean(HeightCheck11$Sol[,8]*5)*date^2

contvar6 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*6) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*6^2) +
  mean(HeightCheck11$Sol[,6]*6^3) +
  mean(HeightCheck11$Sol[,7]*6)*date +
  mean(HeightCheck11$Sol[,8]*6)*date^2

contvar7 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*7) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*7^2) +
  mean(HeightCheck11$Sol[,6]*7^3) +
  mean(HeightCheck11$Sol[,7]*7)*date +
  mean(HeightCheck11$Sol[,8]*7)*date^2

contvar8 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*8) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*8^2) +
  mean(HeightCheck11$Sol[,6]*8^3) +
  mean(HeightCheck11$Sol[,7]*8)*date +
  mean(HeightCheck11$Sol[,8]*8)*date^2

contvar9 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*9) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*9^2) +
  mean(HeightCheck11$Sol[,6]*9^3) +
  mean(HeightCheck11$Sol[,7]*9)*date +
  mean(HeightCheck11$Sol[,8]*9)*date^2

contvar10 <- mean(HeightCheck11$Sol[,1]) +
  mean(HeightCheck11$Sol[,2]*10) +
  mean(HeightCheck11$Sol[,3])*date +
  mean(HeightCheck11$Sol[,4])*date^2 +
  mean(HeightCheck11$Sol[,5]*10^2) +
  mean(HeightCheck11$Sol[,6]*10^3) +
  mean(HeightCheck11$Sol[,7]*10)*date +
  mean(HeightCheck11$Sol[,8]*10)*date^2

plot(date, contvar1, type="l")
points(date, contvar2, type="l")
points(date, contvar3, type="l")
points(date, contvar3, type="l")
points(date, contvar4, type="l")
points(date, contvar5, type="l")
points(date, contvar6, type="l")
points(date, contvar7, type="l")
points(date, contvar8, type="l")
points(date, contvar9, type="l")
points(date, contvar10, type="l")


contvar1.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*1) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*1^2) +
  mean(HeightCheck12$Sol[,6]*1)*date +
  mean(HeightCheck12$Sol[,7]*1)*date^2 

contvar2.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*2) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*2^2) +
  mean(HeightCheck12$Sol[,6]*2)*date +
  mean(HeightCheck12$Sol[,7]*2)*date^2 

contvar3.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*3) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*3^2) +
  mean(HeightCheck12$Sol[,6]*3)*date +
  mean(HeightCheck12$Sol[,7]*3)*date^2

contvar4.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*4) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*4^2) +
  mean(HeightCheck12$Sol[,6]*4)*date +
  mean(HeightCheck12$Sol[,7]*4)*date^2 

contvar5.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*5) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*5^2) +
  mean(HeightCheck12$Sol[,6]*5)*date +
  mean(HeightCheck12$Sol[,7]*5)*date^2 

contvar6.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*6) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*6^2) +
  mean(HeightCheck12$Sol[,6]*6)*date +
  mean(HeightCheck12$Sol[,7]*6)*date^2 

contvar7.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*7) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*7^2) +
  mean(HeightCheck12$Sol[,6]*7)*date +
  mean(HeightCheck12$Sol[,7]*7)*date^2 

contvar8.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*8) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*8^2) +
  mean(HeightCheck12$Sol[,6]*8)*date +
  mean(HeightCheck12$Sol[,7]*8)*date^2 

contvar9.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*9) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*9^2) +
  mean(HeightCheck12$Sol[,6]*9)*date +
  mean(HeightCheck12$Sol[,7]*9)*date^2 

contvar10.2 <- mean(HeightCheck12$Sol[,1]) +
  mean(HeightCheck12$Sol[,2]*10) +
  mean(HeightCheck12$Sol[,3])*date +
  mean(HeightCheck12$Sol[,4])*date^2 +
  mean(HeightCheck12$Sol[,5]*10^2) +
  mean(HeightCheck12$Sol[,6]*10)*date +
  mean(HeightCheck12$Sol[,7]*10)*date^2

plot(date, exp(contvar1.2), type="l", ylim=c(0,400))
points(date, exp(contvar2.2), type="l")
points(date, exp(contvar3.2), type="l")
points(date, exp(contvar4.2), type="l")
points(date, exp(contvar5.2), type="l")
points(date, exp(contvar6.2), type="l")
points(date, exp(contvar7.2), type="l")
points(date, exp(contvar8.2), type="l")
points(date, exp(contvar9.2), type="l")
points(date, exp(contvar10.2), type="l")
points(date, exp(contvar1), type="l", col=2)
points(date, exp(contvar2), type="l", col=2)
points(date, exp(contvar3), type="l", col=2)
points(date, exp(contvar4), type="l", col=2)
points(date, exp(contvar5), type="l", col=2)
points(date, exp(contvar6), type="l", col=2)
points(date, exp(contvar7), type="l", col=2)
points(date, exp(contvar8), type="l", col=2)
points(date, exp(contvar9), type="l", col=2)
points(date, exp(contvar10), type="l", col=2)
points(curves$date, exp(curves$y1)*100, type="l", col=3)
points(curves$date, exp(curves$y2)*100, type="l", col=3)
points(curves$date, exp(curves$y3)*100, type="l", col=3)
points(curves$date, exp(curves$y4)*100, type="l", col=3)
points(curves$date, exp(curves$y5)*100, type="l", col=3)
points(curves$date, exp(curves$y6)*100, type="l", col=3)
points(curves$date, exp(curves$y7)*100, type="l", col=3)
points(curves$date, exp(curves$y8)*100, type="l", col=3)
points(curves$date, exp(curves$y9)*100, type="l", col=3)
points(curves$date, exp(curves$y10)*100, type="l", col=3)
points(x, exp(HC11), type="l")
abline(v=(8/-2)+7, col=2)
abline(v=(5/-2)+7, col=2)
abline(v=(2/-2)+7, col=2)
abline(v=(-1/-2)+7, col=2)
abline(v=(-4/-2)+7, col=2)
abline(v=(-7/-2)+7, col=2)
abline(v=(-10/-2)+7, col=2)
abline(v=(-13/-2)+7, col=2)
abline(v=(-16/-2)+7, col=2)
abline(v=(-19/-2)+7, col=2)

# Peak date and height are incorrect but is width correct?

#Peak date
HC <- data.frame(contvar=c(seq(-100,100,0.1)))#1,2,3,4,5,6,7,8,9,10))

mean(-(HeightCheck12$Sol[,"date"] + HeightCheck12$Sol[,"contvar:date"]*1)/(2*(HeightCheck12$Sol[,"I(date^2)"] + HeightCheck12$Sol[,"contvar:I(date^2)"]*1)))

for(i in 1:length(HC$contvar)){
  HC$peakdate[i] <- mean(-(HeightCheck12$Sol[,"date"] + HeightCheck12$Sol[,"contvar:date"]*HC$contvar[i])/(2*(HeightCheck12$Sol[,"I(date^2)"] + HeightCheck12$Sol[,"contvar:I(date^2)"]*HC$contvar[i])))
}

plot(HC$peakdate, HC$contvar)

  for(i in 1:length(PHCI$Lat)){
  PHCI$height[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                               Biogeog$Sol[,"datecentred"]*PHCI$Peak[i]+
                               Biogeog$Sol[,"latitude"]*PHCI$Lat[i]+
                               Biogeog$Sol[,"elevation"]*PHCI$Elev[i]+
                               Biogeog$Sol[,"I(datecentred^2)"]*(PHCI$Peak[i]^2)+ 
                               Biogeog$Sol[,"datecentred:latitude"]*PHCI$Lat*PHCI$Peak[i]+
                               Biogeog$Sol[,"datecentred:elevation"]*PHCI$Elev*PHCI$Peak[i]))
}

# understanding width
curves <- data.frame(x=seq(-7,17,0.1))  #### Curves, peak becoming later with increasing "y.." peak height consistent
curves$y1  <- -curves$x^2 + -8*curves$x -14.75
curves$y2  <- -0.9*curves$x^2 +-4.5*curves$x - 4.375# -5*curves$x -5
curves$y3  <- -0.8*curves$x^2 +-1.6*curves$x + 0.45# -2*curves$x +0.25
curves$y4  <- -0.7*curves$x^2 +0.75*curves$x + 1.05# 1*curves$x +1
curves$y5  <- -0.6*curves$x^2 +2.357*curves$x-1.065# 4*curves$x -2.75
curves$y6  <- -0.5*curves$x^2 +3.5*curves$x -4.875# 7*curves$x -11
curves$y7  <- -0.4*curves$x^2 +4*curves$x -8.75# 10*curves$x -23.75
curves$y8  <- -0.3*curves$x^2 +3.9*curves$x -11.425# 13*curves$x -41
curves$y9  <- -0.2*curves$x^2 +3.2*curves$x-11.55# 16*curves$x -62.75
curves$y10 <- -0.1*curves$x^2 +1.9*curves$x -7.775# 19*curves$x -89

plot(curves$x, curves$height, lwd=3,type="l",lty="dashed", xlim=c(-7,15), ylim=c(0,5), col=2)
plot(curves$x, curves$y1, type="l", col=1, lwd=3)
points(curves$x, curves$y2, type="l", col=1, lwd=3)
points(curves$x, curves$y3, type="l", col=1, lwd=3)
points(curves$x, curves$y4, type="l", col=1, lwd=3)
points(curves$x, curves$y5, type="l", col=1, lwd=3)
points(curves$x, curves$y6, type="l", col=1, lwd=3)
points(curves$x, curves$y7, type="l", col=1, lwd=3)
points(curves$x, curves$y8, type="l", col=1, lwd=3)
points(curves$x, curves$y9, type="l", col=1, lwd=3)
points(curves$x, curves$y10, type="l", col=1, lwd=3)

a <- -1
b <- -8
c <- -14.75
date <- -b/(2*a)
height <- a*date^2 + b*date + c
ch <- c-(height/2)
w1 <- (-b+sqrt(b^2-4*a*c)/(2*a))
w2 <- (-b-sqrt(b^2-4*a*c)/(2*a))
