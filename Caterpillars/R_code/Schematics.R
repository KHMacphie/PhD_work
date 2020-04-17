####################
#### Schematics ####
####################

#### Effect of temp on Phenology ####
# line graph: temp sensitivity from thackeray 2016
phentemp <- data.frame(temp=seq(4,10,0.5))
phentemp$prod <- 120+(-4.1*phentemp$temp)
phentemp$prim <- 130+(-3.7*phentemp$temp)
phentemp$seco <- 135+(-1.9*phentemp$temp)

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(phentemp$temp, phentemp$seco, type="l", lwd=3, col=2, ylim=c(75, 135), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(phentemp$temp, phentemp$prim, type="l", lwd=3, col=1)
points(phentemp$temp, phentemp$prod, type="l", lwd=3, col=3)
title(xlab="Temperature (°C)",  outer=TRUE, line=0)
title(ylab="Phenology", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 
text(5.5,132, "A", col="red")
text(8.5,132, "B", col="red")
text(4.6, 82.5, "producer", col=3, cex=0.7)
text(5.15, 79, "primary consumer", col=1, cex=0.7)
text(5.35, 75.5, "secondary consumer", col=2, cex=0.7) # 9'X7.5"

# curves for caterpillar and bird at 2 temps
curves <- data.frame(x=seq(-30,30,0.1))
curves$y1 <- -curves$x^2 + -4*curves$x -2.75
curves$y2 <- -curves$x^2 + -4.5*curves$x -3.814
curves$y3 <- -curves$x^2 + -2*curves$x +0.25

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y1), lwd=3,type="l",xlim=c(-4,1), ylim=c(0,5), col=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y3), lwd=3,type="l", col=2, lty="dashed")
points(curves$x, exp(curves$y2), lwd=3,type="l")
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 
text(-2,4, "A", col="red")
text(-1,4, "B", col="red")

#### Buffering hypothesis ####

curves <- data.frame(x=seq(-30,30,0.1))
curves$y4 <- -curves$x^2 + curves$x +1
curves$y5 <- -curves$x^2 + curves$x +1.5
curves$y6 <- -(0.4)*(curves$x)^2 + 0.4*curves$x +1.152
curves$line <- 2.7

# base curve
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y4), type="l", ylim=c(0,6), xlim=c(-2.5,3.5), lwd=3, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, curves$line, type="l", lty="dotted", col=2, lwd=2)
title(ylab="Abundance", outer=TRUE, line=0)
#title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 
text(3.3,5.5, "A", cex=0.8)

# plus higher peak
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y4), type="l", ylim=c(0,6), xlim=c(-2.5,3.5), lwd=3, lty="dashed", xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y5), type="l", lwd=3, col=2)
points(curves$x, curves$line, type="l", lty="dotted", col=2, lwd=2)
title(ylab="Abundance", outer=TRUE, line=0)
#title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 
text(3.3,5.5, "B", cex=0.8)

# plus broader peak
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y4), type="l", ylim=c(0,6), xlim=c(-2.5,3.5), lwd=3, lty="dashed", xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y6), type="l", lwd=3, col=2)
points(curves$x, curves$line, type="l", lty="dotted", col=2, lwd=2)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 
text(3.3,5.5, "C", cex=0.8)


#### Schematic for core model ####
# basic curve
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 

#intercept +/-
curves <- data.frame(x=seq(-30,30,0.1))
curves$y7 <- -curves$x^2 + curves$x +1
curves$y8 <- -curves$x^2 + curves$x +1.2
curves$y9 <- -curves$x^2 + curves$x +0.8

plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y8), type="l", col="red2", lwd=2, lty=5)
points(curves$x, exp(curves$y9), type="l", col="blue2", lwd=2, lty=5)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 

#quadratic +/-
curves <- data.frame(x=seq(-30,30,0.1))
curves$y7 <- -curves$x^2 + curves$x +1
curves$y10 <- -(1.2)*(curves$x)^2 + curves$x +1
curves$y11<- -(0.8)*(curves$x)^2 + curves$x +1

plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y10), type="l", col="red2", lwd=2, lty=5)
points(curves$x, exp(curves$y11), type="l", col="blue2", lwd=2, lty=5)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 

#date +/-
curves <- data.frame(x=seq(-30,30,0.1))
curves$y7  <- -curves$x^2 + curves$x +1
curves$y12 <- -curves$x^2 + 1.2*curves$x +1
curves$y13 <- -curves$x^2 + 0.8*curves$x +1

plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y12), type="l", col="red2", lwd=2, lty=5)
points(curves$x, exp(curves$y13), type="l", col="blue2", lwd=2, lty=5)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 

#date:year +/-
curves <- data.frame(x=seq(-30,30,0.1))
curves$y7  <- -curves$x^2 + curves$x +1
curves$y14 <- -curves$x^2 + 4*curves$x -3
curves$y15 <- -curves$x^2 + -curves$x +1.5

plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y14), type="l", col="blueviolet", lwd=2, lty=5)
points(curves$x, exp(curves$y15), type="l", col="darkgreen", lwd=2, lty=5)  
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 

#contvar:date
curves <- data.frame(x=seq(-30,30,0.1))
curves$y7  <- -curves$x^2 + curves$x +1 
curves$y16 <- -curves$x^2 + 1.6*curves$x +1.06
curves$y17 <- -curves$x^2 + 1.4*curves$x +1.04
curves$y18 <- -curves$x^2 + 1.2*curves$x +1.02
curves$y19 <- -curves$x^2 + 0.8*curves$x +0.98
curves$y20 <- -curves$x^2 + 0.6*curves$x +0.96
curves$y21 <- -curves$x^2 + 0.4*curves$x +0.94

plot(curves$x, exp(curves$y7), type="l", ylim=c(0,6), xlim=c(-2,3), lwd=2, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y16), type="l", col="red3", lwd=2, lty=5)
points(curves$x, exp(curves$y17), type="l", col="red", lwd=2, lty=5)
points(curves$x, exp(curves$y18), type="l", col="orangered", lwd=2, lty=5)
points(curves$x, exp(curves$y19), type="l", col="royalblue", lwd=2, lty=5)
points(curves$x, exp(curves$y20), type="l", col="blue2", lwd=2, lty=5)
points(curves$x, exp(curves$y21), type="l", col="blue4", lwd=2, lty=5)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F) 
axis(side=2,labels=F) 



#### different tree taxa peaks ####

x <- seq(-30,30,0.01)
y1  <- -x^2 + x +1
y2 <- -x^2 + 3*x -1.5
y3 <- -4*x^2 + -x +1.5
y4 <- -2*x^2 + 10*x -11.5
sum1 <- exp(y1)+exp(y2)+exp(y3)+exp(y4)
  
par(mfcol=c(2,1), mar=c(1,1,1,1), cex=1.5)#, oma=c(2,2,1,1),cex=2)
plot(x, exp(y1), type="l", col="blue", ylim=c(0,8), xlim=c(-3,5.5), lwd=2, xlab="", ylab="", xaxt='n', yaxt='n')
points(x, exp(y2), type="l", col="green3", lwd=2)
points(x, exp(y3), type="l", col="darkgreen", lwd=2)  
points(x, exp(y4), type="l", col="purple", lwd=2) 
points(x,sum1, type="l", lwd=2, lty=2)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
#plot(x, exp(y5), type="l", col="white", ylim=c(0,8), xlim=c(-3,5.5), lwd=2, xlab="", ylab="", xaxt='n', yaxt='n')

#axis(side=1,labels=F) 
#axis(side=2,labels=F) 

y5 <- -0.5*x^2 + 3.5*x -5.5
y6 <- -0.3*x^2 + x 
y7 <- -0.5*x^2 + -0.2*x +0.8
y8 <- -4*x^2 + 14*x -11
sum2 <- exp(y5)+exp(y6)+exp(y7)+exp(y8)
sum3 <- exp(y5)+exp(y6)+exp(y7)+exp(y8)+exp(y1)+exp(y2)+exp(y3)+exp(y4)

#par(mfcol=c(1,2), mar=c(1,1,1,1), cex=1.5)#, oma=c(2,2,1,1),cex=2)
plot(x, exp(y5), type="l", col="darkred", ylim=c(0,8), xlim=c(-3,5.5), lwd=2, xlab="", ylab="", xaxt='n', yaxt='n')
points(x, exp(y6), type="l", col="red2", lwd=2)
points(x, exp(y7), type="l", col="orange", lwd=2)  
points(x, exp(y8), type="l", col="gold2", lwd=2) 
points(x,sum2, type="l", lwd=2, lty=2)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)  #saving as 7"x5"

#par(mfcol=c(1,1), mar=c(1,1,1,1), cex=1.5)#, oma=c(2,2,1,1),cex=2)
plot(x, exp(y5), type="l", col="darkred", ylim=c(0,9), xlim=c(-3,5.5), lwd=2, xlab="", ylab="", xaxt='n', yaxt='n')
points(x, exp(y6), type="l", col="red2", lwd=2)
points(x, exp(y7), type="l", col="orange", lwd=2)  
points(x, exp(y8), type="l", col="gold2", lwd=2) 
points(x, exp(y1), type="l", col="blue", lwd=2)
points(x, exp(y2), type="l", col="green3", lwd=2)
points(x, exp(y3), type="l", col="darkgreen", lwd=2)  
points(x, exp(y4), type="l", col="purple", lwd=2) 
points(x,sum3, type="l", lwd=2, lty=2)
