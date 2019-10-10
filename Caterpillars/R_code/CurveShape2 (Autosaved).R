rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)

cater$siteday <- paste(cater$site, cater$day, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$datecent <- cater$date-mean(cater$date)
cater$datescaled <- cater$date/max(cater$date)

a<-1000
prior<-list(R=list(V=diag(1), nu=0.002), 
               G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
               		  G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

DateCubed<- MCMCglmm(caterpillars~datecent*year+I(datecent^2)+I(datecent^3), random=~recorder+siteday+treeID+site, family="poisson", data=cater, prior=prior, nitt=250000, burnin=25000)
save(DateCubed, file = "~/Documents/Models/DateCubed.RData")
load("~/Documents/Models/DateCubed.RData") # DIC 20664.23 

DateSquared<- MCMCglmm(caterpillars~datecent*year+I(datecent^2), 
 random=~recorder+siteday+treeID+site, family="poisson", data=cater, prior=prior, nitt=250000, burnin=25000)
save(DateSquared, file = "~/Documents/Models/DateSquared.RData")
load("~/Documents/Models/DateSquared.RData") # DIC 20689.56

predday <- seq(-30,30, 0.1)
cubic2014 <-mean(DateCubed$Sol[,1])+
			mean(DateCubed$Sol[,2])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3
			
cubic2015 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,3])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,10])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2016 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,4])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,11])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2017 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,5])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,12])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2018 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,6])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,13])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

cubic2019 <-mean(DateCubed$Sol[,1]+DateCubed$Sol[,7])+
			mean(DateCubed$Sol[,2]+DateCubed$Sol[,14])*predday+
			mean(DateCubed$Sol[,8])*predday^2+
			mean(DateCubed$Sol[,9])*predday^3

			
sqrd2014 <- mean(DateSquared$Sol[,1])+
			mean(DateSquared$Sol[,2])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2015 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,3])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,9])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2016 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,4])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,10])*predday+
			mean(DateSquared$Sol[,8])*predday^2
			
sqrd2017 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,5])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,11])*predday+
			mean(DateSquared$Sol[,8])*predday^2

sqrd2018 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,6])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,12])*predday+
			mean(DateSquared$Sol[,8])*predday^2

sqrd2019 <- mean(DateSquared$Sol[,1]+DateSquared$Sol[,7])+
			mean(DateSquared$Sol[,2]+DateSquared$Sol[,13])*predday+
			mean(DateSquared$Sol[,8])*predday^2
	
	
par(mfcol=c(1,2))			
plot(predday,exp(sqrd2014), type="l", col=2, ylim=c(0,0.12))
points(predday,exp(sqrd2015), type="l", col=2)
points(predday,exp(sqrd2016), type="l", col=2)
points(predday,exp(sqrd2017), type="l", col=2)
points(predday,exp(sqrd2018), type="l", col=2)
points(predday,exp(sqrd2019), type="l", col=2)

plot(predday, exp(cubic2014), type="l", ylim=c(0,0.12))
points(predday,exp(cubic2015), type="l", col=1)
points(predday,exp(cubic2016), type="l", col=1)
points(predday,exp(cubic2017), type="l", col=1)
points(predday,exp(cubic2018), type="l", col=1)
points(predday,exp(cubic2019), type="l", col=1)

#### Peak dates for cubic
a <- 3*mean(DateCubed$Sol[,9])
b <- 2*mean(DateCubed$Sol[,8])
c <- mean(DateCubed$Sol[,2])

abline(v=((-b - sqrt(b^2-4*a*c))/(2*a)), col=2)
peak14c <- 

############################################
#### Cubic model with year interactions ####

#DateCubedYear<- MCMCglmm(caterpillars~datescaled*year+I(datescaled^2)*year+I(datescaled^3)*year, random=~recorder+siteday+treeID+site, family="poisson", data=cater, prior=prior, nitt=250000, burnin=25000)
#save(DateCubedYear, file = "~/Documents/Models/DateCubedYear.RData")
load("~/Documents/Models/DateCubedYear.RData")


#### Plot yearly curves

dayscal <- seq(0.67,1,0.001)
curve14 <- mean(DateCubedYear$Sol[,1])+mean(DateCubedYear$Sol[,2])*dayscal+mean(DateCubedYear$Sol[,8])*dayscal^2+mean(DateCubedYear$Sol[,9])*dayscal^3
curve15 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,3])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,10])*dayscal+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,15])*dayscal^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,20])*dayscal^3
curve16 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,4])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,11])*dayscal+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,16])*dayscal^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,21])*dayscal^3
curve17 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,5])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,12])*dayscal+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,17])*dayscal^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,22])*dayscal^3
curve18 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,6])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,13])*dayscal+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,18])*dayscal^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,23])*dayscal^3
curve19 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,7])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,14])*dayscal+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,19])*dayscal^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,24])*dayscal^3
days <- dayscal*max(cater$date)

plot(days,exp(curve14), type="l")
points(days,exp(curve15),type="l", col=2)
points(days,exp(curve16),type="l", col=3)
points(days,exp(curve17),type="l", col=4)
points(days,exp(curve18),type="l", col=5)
points(days,exp(curve19),type="l", col=6)

#### Calculate yearly peak date
# (-b +/- sqrt(b^2-4*a*c))/(2*a)
y14 <- data.frame(intercept=DateCubedYear$Sol[,1], date=DateCubedYear$Sol[,2], date2=DateCubedYear$Sol[,8], date3=DateCubedYear$Sol[,9])
y14$a <- 3*DateCubedYear$Sol[,9]
y14$b <- 2*DateCubedYear$Sol[,8]
y14$c <- DateCubedYear$Sol[,2]

y15 <- data.frame(intercept=(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,3]), date=(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,10]), date2=(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,15]), date3=(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,20]))
y15$a <- 3*(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,20])
y15$b <- 2*(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,15])
y15$c <- (DateCubedYear$Sol[,2]+DateCubedYear$Sol[,10])

y16 <- data.frame(intercept=(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,4]), date=(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,11]), date2=(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,16]), date3=(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,21]))
y16$a <- 3*(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,21])
y16$b <- 2*(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,16])
y16$c <- (DateCubedYear$Sol[,2]+DateCubedYear$Sol[,11])

y17 <- data.frame(intercept=(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,5]), date=(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,12]), date2=(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,17]), date3=(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,22]))
y17$a <- 3*(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,22])
y17$b <- 2*(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,17])
y17$c <- (DateCubedYear$Sol[,2]+DateCubedYear$Sol[,12])

y18 <- data.frame(intercept=(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,6]), date=(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,13]), date2=(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,18]), date3=(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,23]))
y18$a <- 3*(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,23])
y18$b <- 2*(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,18])
y18$c <- (DateCubedYear$Sol[,2]+DateCubedYear$Sol[,13])

y19 <- data.frame(intercept=(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,7]), date=(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,14]), date2=(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,19]), date3=(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,24]))
y19$a <- 3*(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,24])
y19$b <- 2*(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,19])
y19$c <- (DateCubedYear$Sol[,2]+DateCubedYear$Sol[,14])


# Not reverted to day of the year for height calcs
y14$pd <- ((-y14$b - sqrt(y14$b^2-4*y14$a*y14$c))/(2*y14$a))
y15$pd <- ((-y15$b - sqrt(y15$b^2-4*y15$a*y15$c))/(2*y15$a))
y16$pd <- ((-y16$b - sqrt(y16$b^2-4*y16$a*y16$c))/(2*y16$a))
y17$pd <- ((-y17$b - sqrt(y17$b^2-4*y17$a*y17$c))/(2*y17$a))
y18$pd <- ((-y18$b - sqrt(y18$b^2-4*y18$a*y18$c))/(2*y18$a))
y19$pd <- ((-y19$b - sqrt(y19$b^2-4*y19$a*y19$c))/(2*y19$a))

# NaNs in 14-16 not sure what to do  --- all iterations dont have a peak
y14$sqrtbit <- y14$b^2-4*y14$a*y14$c
y14crop <- subset(y14, sqrtbit<0)
seq <- seq(-100,100,0.1)
plot(seq,(y14crop$var1[150]+y14crop$var1.1[150]*seq+y14crop$var1.2[150]*seq^2+y14$var1.3[150]*seq^3), type="l")
#### Calculate yearly peak height

ph14 <-  mean(DateCubedYear$Sol[,1])+mean(DateCubedYear$Sol[,2])*pd14+mean(DateCubedYear$Sol[,8])*pd14^2+mean(DateCubedYear$Sol[,9])*pd14^3
ph15 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,3])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,10])*pd15+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,15])*pd15^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,20])*pd15^3
ph16 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,4])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,11])*pd16+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,16])*pd16^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,21])*pd16^3
ph17 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,5])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,12])*pd17+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,17])*pd17^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,22])*pd17^3
ph18 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,6])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,13])*pd18+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,18])*pd18^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,23])*pd18^3
ph19 <- mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,7])+mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,14])*pd19+mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,19])*pd19^2+mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,24])*pd19^3

abline(h=exp(ph14), col=1, lty="dashed")
abline(h=exp(ph15), col=2, lty="dashed")
abline(h=exp(ph16), col=3, lty="dashed")
abline(h=exp(ph17), col=4, lty="dashed")
abline(h=exp(ph18), col=5, lty="dashed")
abline(h=exp(ph19), col=6, lty="dashed")

#### Calculate width either side at 50% peak height

roots14 <- polyroot(c((mean(DateCubedYear$Sol[,1])-(log(exp(ph14)/2))),mean(DateCubedYear$Sol[,2]),mean(DateCubedYear$Sol[,8]),mean(DateCubedYear$Sol[,9])))
roots15 <- polyroot(c((mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,3])-(log(exp(ph15)/2))),mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,10]),mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,15]),mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,20])))
roots16 <- polyroot(c((mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,4])-(log(exp(ph16)/2))),mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,11]),mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,16]),mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,21])))
roots17 <- polyroot(c((mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,5])-(log(exp(ph17)/2))),mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,12]),mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,17]),mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,22])))
roots18 <- polyroot(c((mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,6])-(log(exp(ph18)/2))),mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,13]),mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,18]),mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,23])))
roots19 <- polyroot(c((mean(DateCubedYear$Sol[,1]+DateCubedYear$Sol[,7])-(log(exp(ph19)/2))),mean(DateCubedYear$Sol[,2]+DateCubedYear$Sol[,14]),mean(DateCubedYear$Sol[,8]+DateCubedYear$Sol[,19]),mean(DateCubedYear$Sol[,9]+DateCubedYear$Sol[,24])))
roots14 <- Re(roots14)[abs(Im(roots14)) < 1e-6]
roots15 <- Re(roots15)[abs(Im(roots15)) < 1e-6]
roots16 <- Re(roots16)[abs(Im(roots16)) < 1e-6]
roots17 <- Re(roots17)[abs(Im(roots17)) < 1e-6]
roots18 <- Re(roots18)[abs(Im(roots18)) < 1e-6]
roots19 <- Re(roots19)[abs(Im(roots19)) < 1e-6]

# each has one real root and 2 complex so:
r14 <- Re(roots14)[abs(Im(roots14)) < 1e-6]
r15 <- Re(roots15)[abs(Im(roots15)) < 1e-6]
r16 <- Re(roots16)[abs(Im(roots16)) < 1e-6]
r17 <- Re(roots17)[abs(Im(roots17)) < 1e-6]
r18 <- Re(roots18)[abs(Im(roots18)) < 1e-6]
r19 <- Re(roots19)[abs(Im(roots19)) < 1e-6]

abline(h=exp(ph14)/2, col=1, lty="dashed")
abline(v=(r14[2]*max(cater$date)), col=1, lty="dashed")
abline(v=(r14[3]*max(cater$date)), col=1, lty="dashed")
abline(h=exp(ph15)/2, col=2, lty="dashed")
abline(v=(r15[2]*max(cater$date)), col=2, lty="dashed")
abline(v=(r15[3]*max(cater$date)), col=2, lty="dashed")
abline(h=exp(ph16)/2, col=3, lty="dashed")
abline(v=(r16[2]*max(cater$date)), col=3, lty="dashed")
abline(v=(r16[3]*max(cater$date)), col=3, lty="dashed")
abline(h=exp(ph17)/2, col=4, lty="dashed")
abline(v=(r17[2]*max(cater$date)), col=4, lty="dashed")
abline(v=(r17[3]*max(cater$date)), col=4, lty="dashed")
abline(h=exp(ph18)/2, col=5, lty="dashed")
abline(v=(r18[2]*max(cater$date)), col=5, lty="dashed")
abline(v=(r18[3]*max(cater$date)), col=5, lty="dashed")
abline(h=exp(ph19)/2, col=6, lty="dashed")
abline(v=(r19[1]*max(cater$date)), col=6, lty="dashed")
abline(v=(r19[3]*max(cater$date)), col=6, lty="dashed")

r1.14 <- r14[2]*max(cater$date)
r2.14 <- r14[3]*max(cater$date)
r1.15 <- r15[2]*max(cater$date)
r2.15 <- r15[3]*max(cater$date)
r1.16 <- r16[2]*max(cater$date)
r2.16 <- r16[3]*max(cater$date)
r1.17 <- r17[2]*max(cater$date)
r2.17 <- r17[3]*max(cater$date)
r1.18 <- r18[2]*max(cater$date)
r2.18 <- r18[3]*max(cater$date)
r1.19 <- r19[1]*max(cater$date)
r2.19 <- r19[3]*max(cater$date)

cubicpeaks <- data.frame(year=c("2014","2015","2016","2017","2018","2019"), 
                         pd=c(pd14*max(cater$date),pd15*max(cater$date),pd16*max(cater$date),pd17*max(cater$date),pd18*max(cater$date),pd19*max(cater$date)),
                         r1=c(r1.14,r1.15,r1.16,r1.17,r1.18,r1.19),
                         r2=c(r2.14,r2.15,r2.16,r2.17,r2.18,r2.19))
cubicpeaks$left <- cubicpeaks$pd-cubicpeaks$r1
cubicpeaks$right <- cubicpeaks$r2-cubicpeaks$pd
cubicpeaks$width <- cubicpeaks$r2-cubicpeaks$r1
cubicpeaks$propleft <- cubicpeaks$left/cubicpeaks$width
cubicpeaks$propright <- cubicpeaks$right/cubicpeaks$width

#### SiteYear cubic ####
cater$siteyear <- paste(cater$site,cater$year)
a<-10000
prior<-list(R=list(V=diag(1), nu=0.002), 
            G=list(G1=list(V=diag(4), nu=4, alpha.mu=c(0,0,0,0), alpha.V=diag(4)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

DateCubedSY2<- MCMCglmm(caterpillars~datescaled+I(datescaled^2)+I(datescaled^3), random=~us(1+datescaled+I(datescaled^2)+I(datescaled^3)):siteyear+recorder+siteday+treeID, family="poisson", data=cater, prior=prior, nitt=1000000, burnin=75000, pr=TRUE, thin=100)
save(DateCubedSY2, file = "~/Documents/Models/DateCubedSY2.RData")
load("~/Documents/Models/DateCubedSY2.RData")

#### Plot curves

dayscal <- seq(0.67,1,0.001)
curve <- mean(DateCubedSY$Sol[,1])+mean(DateCubedSY$Sol[,2])*dayscal+mean(DateCubedSY$Sol[,3])*dayscal^2+mean(DateCubedSY$Sol[,4])*dayscal^3
days <- dayscal*max(cater$date)

plot(days,exp(curve), type="l")


#### Calculate peak date
# (-b +/- sqrt(b^2-4*a*c))/(2*a)
cdf <- data.frame(DateCubedSY$Sol[,1:4])
cdf$a <- 3*DateCubedSY$Sol[,4]
cdf$b <- 2*DateCubedSY$Sol[,3]
cdf$c <- DateCubedSY$Sol[,2]
#abline(v=((-b - sqrt(b^2-4*a*c))/(2*a))*max(cater$date), col=1, lty="dashed")


# Not reverted to day of the year for height calcs
cdf$pd <- ((-cdf$b - sqrt(cdf$b^2-4*cdf$a*cdf$c))/(2*cdf$a))

#### Calculate peak height

cdf$ph <- DateCubedSY$Sol[,1]+DateCubedSY$Sol[,2]*cdf$pd+DateCubedSY$Sol[,3]*cdf$pd^2+DateCubedSY$Sol[,4]*cdf$pd^3


#### Calculate width either side at 50% peak height
for(i in 1:22500){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/2))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1[i] <- A[1]
  cdf$r2[i] <- A[2]
  cdf$r3[i] <- A[3]
  }

cdf$r1 <- Re(cdf$r1)[abs(Im(cdf$r1)) < 1e-6]
cdf$r2 <- Re(cdf$r2)[abs(Im(cdf$r1)) < 1e-6]
cdf$r3 <- Re(cdf$r3)[abs(Im(cdf$r1)) < 1e-6]

cdf$r1 <- ifelse(cdf$r1<min(cater$datescaled),NA,cdf$r1)
cdf$r2 <- ifelse(cdf$r2<min(cater$datescaled),NA,cdf$r2)
cdf$r3 <- ifelse(cdf$r3<min(cater$datescaled),NA,cdf$r3)


roots <- data.frame(r1=cdf$r1,r2=cdf$r2,r3=cdf$r3)
cdf$root1 <- (apply(roots, 1, min, na.rm=TRUE))*max(cater$date)   
cdf$root2 <- (apply(roots, 1, max, na.rm=TRUE))*max(cater$date)  

cdf$left <- (cdf$pd*max(cater$date))-cdf$root1
cdf$right <- cdf$root2-(cdf$pd*max(cater$date))
cdf$width <- cdf$left+cdf$right
cdf$propleft <- cdf$left/cdf$width
cdf$propright <- cdf$right/cdf$width

mean(cdf$propleft)
HPDinterval(cdf$propleft)
mean(cdf$propright)
HPDinterval(cdf$propright)
mean(cdf$width)
HPDinterval(cdf$width)

#### DateCubedSY3- SY not with date^3 to see if converges better
a<-10000
prior<-list(R=list(V=diag(1), nu=0.002), 
            G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

DateCubedSY3<- MCMCglmm(caterpillars~datescaled+I(datescaled^2)+I(datescaled^3), random=~us(1+datescaled+I(datescaled^2)):siteyear+recorder+siteday+treeID, family="poisson", data=cater, prior=prior, nitt=300000, burnin=30000, pr=TRUE, thin=50)
save(DateCubedSY3, file = "~/Documents/Models/DateCubedSY3.RData")
load("~/Documents/Models/DateCubedSY3.RData")
