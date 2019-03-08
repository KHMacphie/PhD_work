###################################
#### Extra Presentation Graphs ####
###################################

#### Abundance/peak height variation
pred_day <- seq(100,200,1) # different to all other graphs!
H.TempPoisson.2018 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2018.2 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.3) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])  #Apr

H.TempPoisson.2018.3 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.6) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.4 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]+0.9) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.5 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]-0.3) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

H.TempPoisson.2018.6 <-  
  mean(TempPoisson$Sol[,"(Intercept)"]+TempPoisson$Sol[,"year2018"]-0.6) +  #Intercept+yearchange
  mean(TempPoisson$Sol[,"date"]+TempPoisson$Sol[,"date:year2018"])*pred_day +  #date+yearchange
  mean(TempPoisson$Sol[,"I(date^2)"])*pred_day^2 +  #I(date^2)
  mean(8.442*TempPoisson$Sol[,"date:Apr"])*pred_day +  #date*Apr
  mean(8.442*TempPoisson$Sol[,"Apr"])

par(mfcol=c(1,1),mar=c(2,1,1,1), oma=c(1,3,0,0),cex=1.8)
plot(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=8, xlab="Date", ylab="Abundance", ylim=c(0,0.2)) #black
points(pred_day, exp(H.TempPoisson.2018.2), type="l",lwd=3, col=4) 
points(pred_day, exp(H.TempPoisson.2018.3), type="l",lwd=3,col=5) 
points(pred_day, exp(H.TempPoisson.2018.4), type="l", lwd=3, col=3)
points(pred_day, exp(H.TempPoisson.2018.5), type="l", lwd=3, col=6)
points(pred_day, exp(H.TempPoisson.2018.6), type="l", lwd=3, col=2)
title(ylab="Abundance", outer=TRUE, line = 1)
title( xlab="Date", outer=TRUE, line = 0)

plot(pred_day, exp(H.TempPoisson.2018), type="l", lwd=3, col=2, xlab="Date", ylab="Abundance", ylim=c(0,0.12), xlim=c(100,200)) 
