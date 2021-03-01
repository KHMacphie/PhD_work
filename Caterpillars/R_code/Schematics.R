####################
#### Schematics ####
####################

#### Effect of temp on Phenology ####
# line graph: temp sensitivity from thackeray 2016
phentemp <- data.frame(temp=seq(4,10,0.5))
phentemp$prod <- 120+(-4.1*phentemp$temp)
phentemp$prim <- 130+(-3.7*phentemp$temp)
phentemp$seco <- 135+(-1.9*phentemp$temp)

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.5)
plot(phentemp$temp, phentemp$seco, type="l", lwd=3, col=2, ylim=c(75, 135), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(phentemp$temp, phentemp$prim, type="l", lwd=3, col=1)
points(phentemp$temp, phentemp$prod, type="l", lwd=3, col=3)
title(xlab="Temperature (°C)",  outer=TRUE, line=0)
title(ylab="Phenology", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
text(5.5,132, "A", col="red")
text(8.5,132, "B", col="red")
text(4.6, 82.5, "producer", col=3, cex=0.7)
text(5.15, 79, "primary consumer", col=1, cex=0.7)
text(5.35, 75.5, "secondary consumer", col=2, cex=0.7) # 9'X7.5"

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(phentemp$temp, phentemp$prod, type="l", lwd=3, col=3, ylim=c(75, 135), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(xlab="Temperature (°C)",  outer=TRUE, line=0)
title(ylab="Phenology (timing)", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(phentemp$temp, phentemp$prim, type="l", lwd=3, col=1)
points(phentemp$temp, phentemp$seco, type="l", lwd=3, col=2) #6x6
text(4.59, 82.5, "producer", col=3, cex=0.7)
text(5.2, 79, "primary consumer", col=1, cex=0.7)
text(5.39, 75.5, "secondary consumer", col=2, cex=0.7)
text(5.5,132, "A", col="red")
text(8.5,132, "B", col="red")

# curves for caterpillar and bird at 2 temps
curves <- data.frame(x=seq(-30,30,0.1))
curves$y1 <- -curves$x^2 + -4*curves$x -2.75
curves$y2 <- -curves$x^2 + -4.5*curves$x -3.814
curves$y3 <- -curves$x^2 + -2*curves$x +0.25

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y2), lwd=3,type="l",xlim=c(-4,1), ylim=c(0,5), col=1, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(curves$x, exp(curves$y1), lwd=3,type="l", col=2)
points(curves$x, exp(curves$y3), lwd=3,type="l", col=2, lty="dashed")
text(-2,4, "A", col="red")
text(-1,4, "B", col="red") #5.5x6

curves$y2 <- -0.5*curves$x^2 + -2.2*curves$x -0.9
points(curves$x, exp(curves$y2), type="l", lwd=3, lty="dashed")
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y2), lwd=3,type="l",xlim=c(-4,1), ylim=c(0,5), col=1, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(curves$x, exp(curves$y1), lwd=3,type="l", col=2)
points(curves$x, exp(curves$y3), lwd=3,type="l", col=2, lty="dashed")

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
#title(ylab="Abundance", outer=TRUE, line=0)
#title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F,tck=-0.02) 
axis(side=2,labels=F,tck=-0.02) #3.5x4"
#text(3.3,5.5, "A", cex=0.8)

# plus higher peak
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y4), type="l", ylim=c(0,6), xlim=c(-2.5,3.5), lwd=3, lty="dashed", xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y5), type="l", lwd=3, col=2)
points(curves$x, curves$line, type="l", lty="dotted", col=2, lwd=2)
#title(ylab="Abundance", outer=TRUE, line=0)
#title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F,tck=-0.02) 
axis(side=2,labels=F,tck=-0.02) 
#text(3.3,5.5, "B", cex=0.8)

# plus broader peak
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=2)
plot(curves$x, exp(curves$y4), type="l", ylim=c(0,6), xlim=c(-2.5,3.5), lwd=3, lty="dashed", xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y6), type="l", lwd=3, col=2)
points(curves$x, curves$line, type="l", lty="dotted", col=2, lwd=2)
#title(ylab="Abundance", outer=TRUE, line=0)
#title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F,tck=-0.02) 
axis(side=2,labels=F,tck=-0.02) 
#text(3.3,5.5, "C", cex=0.8)


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



####  Mismatch buffering model options ####

x <- seq(-10,10,0.5)
y <- -5*x^2+x
y1 <- -5*x^2+x+30
y2 <- -5*x^2+x+100
y3 <- -5*x^2+x-30
y4 <- -5*x^2+x-100
par(mfcol=c(1,1),cex=1.5)# mar=c(3,3,1,1), cex=1.5, oma=c(3,3,1,1))
plot(x,y, type="l", xlab="Mistiming (days)", ylab="Fitness trait e.g. mean chick mass",  xlim=c(-10,10), ylim=c(-620,110), yaxt='n', lwd=2)
points(x, y1, type="l", col=2, lty="dashed", lwd=2)
points(x, y3, type="l", col="blue", lty="dashed", lwd=2)
points(x, y2, type="l", col=2, lty="dotted", lwd=2)
points(x, y4, type="l", col="blue", lty="dotted", lwd=2) #6"x6"

plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c("Average peak shape","Height", "Width", "Increase","Decrease"), 
       lty=c("solid","dotted", "dashed","solid","solid"), 
       col=c(1,1,1,2,"blue"), 
       lwd=3, seg.len=2, bty='n')

x <- seq(0,10,0.5)
y <- -x
y1 <- -5*x^2+x+30
y2 <- -5*x^2+x+100
y3 <- -5*x^2+x-30
y4 <- -5*x^2+x-100
par(mfcol=c(1,1),cex=1.5)# mar=c(3,3,1,1), cex=1.5, oma=c(3,3,1,1))
plot(x,y, type="l", xlab="Mistiming (days)", ylab="Fitness trait e.g. mean chick mass",  xlim=c(-10,10), ylim=c(-620,110), yaxt='n', lwd=2)
points(x, y1, type="l", col=2, lty="dashed", lwd=2)
points(x, y3, type="l", col="blue", lty="dashed", lwd=2)
points(x, y2, type="l", col=2, lty="dotted", lwd=2)
points(x, y4, type="l", col="blue", lty="dotted", lwd=2) #6"x6"


### OR ####
# fitness ~ absolute relative timing * area * sign

x <- seq(0,10,0.1)
negrtim <- 3*-x+2
posrtim <- 2*-x+2

par(mfcol=c(1,1),cex=1.5)
plot(x, posrtim, type="l", col="red", yaxt='n', lwd=2, ylab="Fitness trait", xlab="Absolute relative timing", ylim=c(-40,10))  
points(x, negrtim, type="l", col="blue", lwd=2)


#### Multimembership ####
# line graph: crossing at intercept
line <- data.frame(temp=seq(-10,10,0.5))
line$prod <- 1+(0.5*line$temp)
line$prim <- 1+(-0.3*line$temp)
line$seco <- 1+(1*line$temp)


par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(line$temp, line$prod, type="l", lwd=3, col=1, ylim=c(-15,30),xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(xlab="Foliage score (small trees)",  outer=TRUE, line=0)
title(ylab="Abundance per branch", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(line$temp, line$prim, type="l", lwd=3, col="blue", lty="dashed")
points(line$temp, line$seco, type="l", lwd=3, col=2, lty="dashed") #5.5x6


# Example caterpillar peaks
curves <- data.frame(x=seq(-30,30,0.01))
curves$y1 <- -0.8*curves$x^2 + -0.5*curves$x
curves$y2 <- -curves$x^2 + -4.5*curves$x -3.814
curves$y3 <- -1.9*curves$x^2 + -4*curves$x -0.6
curves$y4 <- -curves$x^2 + -2.5*curves$x -0.6
points(curves$x, exp(curves$y4), lwd=3,type="l")

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y4), lwd=3,type="l",xlim=c(-4,1), ylim=c(0,5), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="y", outer=TRUE, line=0)
title(xlab="x", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(curves$x, exp(curves$y2), lwd=3,type="l", col="orange", lty="dashed")
points(curves$x, exp(curves$y1), lwd=3,type="l", col=2, lty="dashed")
points(curves$x, exp(curves$y3), lwd=3,type="l", col="blue", lty="dashed")#5.5x6

curves$y1 <- -0.8*curves$x^2 + -0.5*curves$x +0.5
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y3), lwd=3,type="l",xlim=c(-4,1), col="blue", ylim=c(0,5), xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(curves$x, exp(curves$y2), lwd=3,type="l", col="orange")
points(curves$x, exp(curves$y1), lwd=3,type="l", col=2) #4x6



curves$y2 <- -0.5*curves$x^2 + -2.2*curves$x -0.9
points(curves$x, exp(curves$y2), type="l", lwd=3)
par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y2), lwd=3,type="l",xlim=c(-4,1), ylim=c(0,5), col=1, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) 
points(curves$x, exp(curves$y1), lwd=3,type="l", col=2)
points(curves$x, exp(curves$y3), lwd=3,type="l", col=2, lty="dashed")


####Problem with GLMM approach

curves <- data.frame(x=seq(-30,30,0.05))
curves$y7  <- -curves$x^2 + curves$x +1 
curves$y16 <- -curves$x^2 + 1.6*curves$x +1.06
curves$y17 <- -curves$x^2 + 1.4*curves$x +1.04
curves$y18 <- -curves$x^2 + 1.2*curves$x +1.02
curves$y19 <- -curves$x^2 + 0.8*curves$x +0.98
curves$y20 <- -curves$x^2 + 0.6*curves$x +0.96
curves$y21 <- -curves$x^2 + 0.4*curves$x +0.94

curves$y7  <- -curves$x^2 + curves$x +1 
curves$y16 <- -curves$x^2 + 1.6*curves$x +0.94
curves$y17 <- -curves$x^2 + 1.4*curves$x +0.96
curves$y18 <- -curves$x^2 + 1.2*curves$x +0.98
curves$y19 <- -curves$x^2 + 0.8*curves$x +1.02
curves$y20 <- -curves$x^2 + 0.6*curves$x +1.04
curves$y21 <- -curves$x^2 + 0.4*curves$x +1.06
curves$y22 <- -curves$x^2 + 0.2*curves$x +1.08
curves$y23 <- -curves$x^2 + 0*curves$x +1.1
curves$y24 <- -curves$x^2 + -0.2*curves$x +1.12
curves$y25 <- -curves$x^2 + -0.4*curves$x +1.14
curves$y26 <- -curves$x^2 + -0.6*curves$x +1.16
curves$y27 <- -curves$x^2 + -0.8*curves$x +1.18

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(curves$x, exp(curves$y16), type="l",col="#27408B", ylim=c(0,6), xlim=c(-2.3,2.7), lwd=1.5, xlab=NA, ylab=NA, xaxt='n', yaxt='n')
points(curves$x, exp(curves$y17), type="l", col="#3A5FCD", lwd=1.5)
points(curves$x, exp(curves$y18), type="l", col="#436EEE", lwd=1.5)
points(curves$x, exp(curves$y7), type="l", col="#4876FF", lwd=1.5)
points(curves$x, exp(curves$y19), type="l", col="#7A67EE", lwd=1.5)
points(curves$x, exp(curves$y20), type="l", col="#8968CD", lwd=1.5)
points(curves$x, exp(curves$y21), type="l", col="#9F79EE", lwd=1.5)
points(curves$x, exp(curves$y22), type="l", col="#AB82FF", lwd=1.5)
points(curves$x, exp(curves$y23), type="l", col="#D15FEE", lwd=1.5)
points(curves$x, exp(curves$y24), type="l", col="#E066FF", lwd=1.5)
points(curves$x, exp(curves$y25), type="l", col="#EE7AE9", lwd=1.5)
points(curves$x, exp(curves$y26), type="l", col="#FF83FA", lwd=1.5)
points(curves$x, exp(curves$y27), type="l", col="#F29EEE", lwd=1.5)
title(ylab="Abundance", outer=TRUE, line=0)
title(xlab="Date", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02) #5.5x7

df <- data.frame(a=rep(-1, 13), b=seq(1.6,-0.8,-0.2), c=seq(0.94,1.18,0.02))
df$pd <- -df$b/(2*df$a)
df$ph <- df$a*df$pd^2 + df$b*df$pd + df$c

par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(df$c,df$pd, xlab=NA, ylab=NA, xaxt='n', yaxt='n', col=c("#27408B", "#3A5FCD", "#436EEE", "#4876FF", "#7A67EE", "#8968CD", "#9F79EE", "#AB82FF", "#D15FEE", "#E066FF", "#EE7AE9", "#FF83FA", "#F29EEE"), pch=20)
title(ylab="Peak date", outer=TRUE, line=0)
title(xlab="Temperature", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02)


par(mfcol=c(1,1), mar=c(1,1,1,1), oma=c(2,2,1,1),cex=1.8)
plot(df$c,df$ph, xlab=NA, ylab=NA, xaxt='n', yaxt='n', col=c("#27408B", "#3A5FCD", "#436EEE", "#4876FF", "#7A67EE", "#8968CD", "#9F79EE", "#AB82FF", "#D15FEE", "#E066FF", "#EE7AE9", "#FF83FA", "#F29EEE"), pch=20)
title(ylab="Peak height", outer=TRUE, line=0)
title(xlab="Temperature", outer=TRUE, line=0)
axis(side=1,labels=F, tck=-0.02) 
axis(side=2,labels=F, tck=-0.02)

c("#27408B", "#3A5FCD", "#436EEE", "#4876FF", "#7A67EE", "#8968CD", "#9F79EE", "#AB82FF", "#D15FEE", "#E066FF", "#EE7AE9", "#FF83FA", "#F29EEE")