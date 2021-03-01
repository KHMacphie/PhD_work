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

birdphen <- read.csv("~/Dropbox/master_data/blue tits/Bird_Phenology.csv")
nestlings <- read.csv("~/Dropbox/master_data/blue tits/Nestlings.csv")

colnames(birdphen)
colnames(nestlings)

birdphen$v1v2 <- birdphen$v2date-birdphen$v1date
nestlings$v1v2 <- nestlings$v2date-nestlings$v1date
table(birdphen$v1v2)
table(nestlings$v1v2)
nestlings$SYB <- paste(nestlings$site, nestlings$box, nestlings$year)
nestlings$v2mass <- gsub("N/A",NA, nestlings$v2mass)
SYB <- data.frame(SYB=nestlings$SYB, v2mass=nestlings$v2mass)

SYB$v2mass <- as.numeric(SYB$v2mass)
SYB <- aggregate(.~SYB, SYB, mean)
store <- pmatch(SYB$SYB,nestlings$SYB)
SYB$site <- nestlings$site[store]
SYB$year <- nestlings$year[store]
SYB$SY <- paste(SYB$site, SYB$year)

fledged <- data.frame(SYB=nestlings$SYB, fledged=nestlings$fledged)
fledged <- aggregate(.~SYB, fledged, sum)

store1 <- pmatch(SYB$SYB,fledged$SYB)
SYB$fledged <- fledged$fledged[store1] 

birdphen$SYB <- paste(birdphen$site, birdphen$box, birdphen$year)
store2 <- pmatch(SYB$SYB,birdphen$SYB)
SYB$hd <- as.numeric(birdphen$hd_1.45[store2])
SYB$v2date <- birdphen$v2date[store2]
SYB$v2age <- SYB$v2date-SYB$hd

biomassinc <- data.frame(SYB=nestlings$SYB, v1mass=nestlings$v1mass, v2mass=as.numeric(nestlings$v2mass))
biomassinc <- aggregate(.~SYB, biomassinc, sum)
biomassinc$biomassinc <- biomassinc$v2mass - biomassinc$v1mass
store4 <- pmatch(SYB$SYB,biomassinc$SYB)
SYB$biomassinc <- biomassinc$biomassinc[store4]
SYB$v2biomass <- biomassinc$v2mass[store4]

######################################################
#### Site-year caterpillar peaks ####
######################################################

#### Estimates for each site-year ####
both <- read.csv("Documents/Models/both.csv")
both$SY <- both$siteyear

#Using mcmcglmm estimates
store3<- pmatch(SYB$SY, both$SY, duplicates.ok = TRUE)
SYB$mu <- both$PD[store3]
SYB$logheight <- both$PlH[store3]
SYB$sigma <- both$PS[store3]
SYB <- na.omit(SYB)
SYB$mu <- as.integer(SYB$mu)
SYB$async <- (SYB$hd+10)-SYB$mu

ggplot(SYB, aes(async, v2mass))+geom_point()+geom_smooth(model=lm)
ggplot(SYB, aes(async, fledged))+geom_point()+geom_smooth(model=lm)
ggplot(SYB, aes(async, v2biomass))+geom_point()+geom_smooth(model=lm)
ggplot(SYB, aes(async, biomassinc))+geom_point()+geom_smooth(model=lm)
SYB$v2age <- as.factor(SYB$v2age)
SYB$v2age <- relevel(SYB$v2age, ref="12")
SYB$logheight.cent <- SYB$logheight-mean(SYB$logheight)
SYB$logheight.scal <- scale(SYB$logheight) # cent: -2.348996, scale: 1.123537
SYB$sigma.scal <- scale(SYB$sigma) # cent: 11.11513, scale: 2.796189


#### preliminary analysis
mcmcmod1 <- MCMCglmm(v2biomass~async+I(async^2)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=30000, thin=20)
summary(mcmcmod1)
plot(mcmcmod1)

mod1 <- lmer(v2biomass~async+I(async^2)+logheight.scal+sigma.scal+v2age+(1|site)+(1|year)+(1|SY), data=SYB)
summary(mod1)
plot(mod1)

cube.mod1 <- lmer(v2biomass~async+I(async^2)+I(async^3)+logheight.scal+sigma.scal+v2age+(1|site)+(1|year)+(1|SY), data=SYB)
summary(cube.mod1)

mcmcmod2 <- MCMCglmm(v2mass~async+I(async^2)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=30000, thin=20)
summary(mcmcmod2)
plot(mcmcmod2)

mod2 <- lmer(v2mass~async+I(async^2)+logheight.scal+sigma.scal+v2age+(1|site)+(1|year)+(1|SY), data=SYB)
summary(mod2)
plot(mod2)

mcmcmod3 <- MCMCglmm(fledged~async+I(async^2)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=30000, thin=20)
summary(mcmcmod3)
plot(mcmcmod3)

mod3 <- lmer(fledged~async+I(async^2)+logheight.scal+sigma.scal+v2age+(1|site)+(1|year)+(1|SY), data=SYB)
summary(mod3)
plot(mod3)

mycolblack <- rgb(100, 100, 100, max = 250, alpha = 100, names = "blacktrans") # 0,0,0=black higher =lighter
async <- seq(-30,30,0.2)

biomass.avpeak <- mean(mcmcmod1$Sol[,"logheight.scal"])*0 + mean(mcmcmod1$Sol[,"sigma.scal"])*0 + mean(mcmcmod1$Sol[,"(Intercept)"]) + mean(mcmcmod1$Sol[,"async"])*async + mean(mcmcmod1$Sol[,"I(async^2)"])*(async^2) 
biomass.highpeak <- mean(mcmcmod1$Sol[,"logheight.scal"])*2 + mean(mcmcmod1$Sol[,"sigma.scal"])*0 + mean(mcmcmod1$Sol[,"(Intercept)"]) + mean(mcmcmod1$Sol[,"async"])*async + mean(mcmcmod1$Sol[,"I(async^2)"])*(async^2)
biomass.lowpeak <- mean(mcmcmod1$Sol[,"logheight.scal"])*-2 + mean(mcmcmod1$Sol[,"sigma.scal"])*0 + mean(mcmcmod1$Sol[,"(Intercept)"]) + mean(mcmcmod1$Sol[,"async"])*async + mean(mcmcmod1$Sol[,"I(async^2)"])*(async^2)
biomass.widepeak <- mean(mcmcmod1$Sol[,"logheight.scal"])*0 + mean(mcmcmod1$Sol[,"sigma.scal"])*2 + mean(mcmcmod1$Sol[,"(Intercept)"]) + mean(mcmcmod1$Sol[,"async"])*async + mean(mcmcmod1$Sol[,"I(async^2)"])*(async^2)
biomass.narrowpeak <- mean(mcmcmod1$Sol[,"logheight.scal"])*0 + mean(mcmcmod1$Sol[,"sigma.scal"])*-2 + mean(mcmcmod1$Sol[,"(Intercept)"]) + mean(mcmcmod1$Sol[,"async"])*async + mean(mcmcmod1$Sol[,"I(async^2)"])*(async^2)
plot(SYB$async, SYB$v2biomass, pch=20, col=mycolblack)
points(async, biomass.avpeak, type="l", col="black")
points(async, biomass.lowpeak, type="l", col="blue", lty="dotdash")
points(async, biomass.highpeak, type="l", col="red", lty="dotdash")
points(async, biomass.narrowpeak, type="l", col="blue", lty="dashed")
points(async, biomass.widepeak, type="l", col="red", lty="dashed")

mass.avpeak <- mean(mcmcmod2$Sol[,"logheight.scal"])*0 + mean(mcmcmod2$Sol[,"sigma.scal"])*0 + mean(mcmcmod2$Sol[,"(Intercept)"]) + mean(mcmcmod2$Sol[,"async"])*async + mean(mcmcmod2$Sol[,"I(async^2)"])*(async^2) 
mass.highpeak <- mean(mcmcmod2$Sol[,"logheight.scal"])*2 + mean(mcmcmod2$Sol[,"sigma.scal"])*0 + mean(mcmcmod2$Sol[,"(Intercept)"]) + mean(mcmcmod2$Sol[,"async"])*async + mean(mcmcmod2$Sol[,"I(async^2)"])*(async^2)
mass.lowpeak <- mean(mcmcmod2$Sol[,"logheight.scal"])*-2 + mean(mcmcmod2$Sol[,"sigma.scal"])*0 + mean(mcmcmod2$Sol[,"(Intercept)"]) + mean(mcmcmod2$Sol[,"async"])*async + mean(mcmcmod2$Sol[,"I(async^2)"])*(async^2)
mass.widepeak <- mean(mcmcmod2$Sol[,"logheight.scal"])*0 + mean(mcmcmod2$Sol[,"sigma.scal"])*2 + mean(mcmcmod2$Sol[,"(Intercept)"]) + mean(mcmcmod2$Sol[,"async"])*async + mean(mcmcmod2$Sol[,"I(async^2)"])*(async^2)
mass.narrowpeak <- mean(mcmcmod2$Sol[,"logheight.scal"])*0 + mean(mcmcmod2$Sol[,"sigma.scal"])*-2 + mean(mcmcmod2$Sol[,"(Intercept)"]) + mean(mcmcmod2$Sol[,"async"])*async + mean(mcmcmod2$Sol[,"I(async^2)"])*(async^2)
plot(SYB$async, SYB$v2mass, pch=20, col=mycolblack)
points(async, mass.avpeak, type="l", col="black")
points(async, mass.lowpeak, type="l", col="blue", lty="dotdash")
points(async, mass.highpeak, type="l", col="red", lty="dotdash")
points(async, mass.narrowpeak, type="l", col="blue", lty="dashed")
points(async, mass.widepeak, type="l", col="red", lty="dashed")

fledged.avpeak <- mean(mcmcmod3$Sol[,"logheight.scal"])*0 + mean(mcmcmod3$Sol[,"sigma.scal"])*0 + mean(mcmcmod3$Sol[,"(Intercept)"]) + mean(mcmcmod3$Sol[,"async"])*async + mean(mcmcmod3$Sol[,"I(async^2)"])*(async^2) 
fledged.highpeak <- mean(mcmcmod3$Sol[,"logheight.scal"])*2 + mean(mcmcmod3$Sol[,"sigma.scal"])*0 + mean(mcmcmod3$Sol[,"(Intercept)"]) + mean(mcmcmod3$Sol[,"async"])*async + mean(mcmcmod3$Sol[,"I(async^2)"])*(async^2)
fledged.lowpeak <- mean(mcmcmod3$Sol[,"logheight.scal"])*-2 + mean(mcmcmod3$Sol[,"sigma.scal"])*0 + mean(mcmcmod3$Sol[,"(Intercept)"]) + mean(mcmcmod3$Sol[,"async"])*async + mean(mcmcmod3$Sol[,"I(async^2)"])*(async^2)
fledged.widepeak <- mean(mcmcmod3$Sol[,"logheight.scal"])*0 + mean(mcmcmod3$Sol[,"sigma.scal"])*2 + mean(mcmcmod3$Sol[,"(Intercept)"]) + mean(mcmcmod3$Sol[,"async"])*async + mean(mcmcmod3$Sol[,"I(async^2)"])*(async^2)
fledged.narrowpeak <- mean(mcmcmod3$Sol[,"logheight.scal"])*0 + mean(mcmcmod3$Sol[,"sigma.scal"])*-2 + mean(mcmcmod3$Sol[,"(Intercept)"]) + mean(mcmcmod3$Sol[,"async"])*async + mean(mcmcmod3$Sol[,"I(async^2)"])*(async^2)
plot(SYB$async, SYB$fledged, pch=20, col=mycolblack)
points(async, fledged.avpeak, type="l", col="black")
points(async, fledged.lowpeak, type="l", col="blue", lty="dotdash")
points(async, fledged.highpeak, type="l", col="red", lty="dotdash")
points(async, fledged.narrowpeak, type="l", col="blue", lty="dashed")
points(async, fledged.widepeak, type="l", col="red", lty="dashed")

### Model1 (biomass) fit ###
# Plot fitted values and residuals 
mu.mod1 <- fitted(mod1) #Get fitted values (covariate effects + random effects)
resid.mod1  <- resid(mod1)  #Get residuals

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu.mod1, 
     y = resid.mod1,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Homegeneity")
abline(h = 0, lty = 2)     


# Plot residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("year", "site")
SYB$resid.1 <- resid.mod1
source("~/Documents/R-INLA Course/AllData/HighstatLibV13.R")
MyMultipanel.ggp2(Z = SYB, 
                  varx = MyVar, 
                  vary = "resid.1", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = TRUE,
                  addHorizontalLine = TRUE)

par(mfrow = c(2,2))
plot(SYB$async, 
     SYB$resid.1,
     xlab = "async",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$mu, 
     SYB$resid.1,
     xlab = "mu",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$logheight.scal, 
     SYB$resid.1,
     xlab = "logheight.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2) 
plot(SYB$sigma.scal, 
     SYB$resid.1,
     xlab = "sigma.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2)


# Plot the residuals versus year and site.
par(mfrow=c(1,2))
boxplot(resid.1 ~ year, data = SYB)
boxplot(resid.1 ~ site, SYB)


# Visualise the residuals on a map
SYB$Sign.1 <- ifelse(SYB$resid.1 >=0, "positive", "negative") #red for positive residuals and
#black for negative residuals
SYB$Size.1 <-  round(5 * sqrt(abs(SYB$resid.1) /max(SYB$resid.1))) + 1 #Size of residual
site <- read.csv("Dropbox/master_data/site/site_details.csv")
store5 <- pmatch(SYB$site,site$site,duplicates.ok = TRUE)
SYB$Longitude <- site$Mean.Long[store5]
SYB$Latitude <- site$Mean.Lat[store5]

p2 <- ggplot()
p2 <- p2 + theme(text = element_text(size=13)) 
p2 <- p2 + geom_jitter(aes(x = Longitude, 
                           y = Latitude,
                           size = Size.1,
                           col = Sign.1),
                       shape = 1,
                       data = SYB) 
p2 <- p2 + xlab("Longitude") + ylab("Latitude")  
p2



### Model2 (mass) fit ###
# Plot fitted values and residuals 
mu.mod2 <- fitted(mod2) #Get fitted values (covariate effects + random effects)
resid.mod2  <- resid(mod2)  #Get residuals

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu.mod2, 
     y = resid.mod2,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Homegeneity")
abline(h = 0, lty = 2)     


# Plot residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("year", "site")
SYB$resid.2 <- resid.mod2
source("~/Documents/R-INLA Course/AllData/HighstatLibV13.R")
MyMultipanel.ggp2(Z = SYB, 
                  varx = MyVar, 
                  vary = "resid.2", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = TRUE,
                  addHorizontalLine = TRUE)

par(mfrow = c(2,2))
plot(SYB$async, 
     SYB$resid.2,
     xlab = "async",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$mu, 
     SYB$resid.2,
     xlab = "mu",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$logheight.scal, 
     SYB$resid.2,
     xlab = "logheight.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2) 
plot(SYB$sigma.scal, 
     SYB$resid.2,
     xlab = "sigma.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2)


# Plot the residuals versus year and site.
par(mfrow=c(1,2))
boxplot(resid.2 ~ year, data = SYB)
boxplot(resid.2 ~ site, SYB)


# Visualise the residuals on a map
SYB$Sign.2 <- ifelse(SYB$resid.2 >=0, "positive", "negative") #red for positive residuals and
#black for negative residuals
SYB$Size.2 <-  round(5 * sqrt(abs(SYB$resid.2) /max(SYB$resid.2))) + 1 #Size of residual


p2 <- ggplot()
p2 <- p2 + theme(text = element_text(size=13)) 
p2 <- p2 + geom_jitter(aes(x = Longitude, 
                           y = Latitude,
                           size = Size.2,
                           col = Sign.2),
                       shape = 1,
                       data = SYB) 
p2 <- p2 + xlab("Longitude") + ylab("Latitude")  
p2



### Model2 (mass) fit ###
# Plot fitted values and residuals 
mu.mod3 <- fitted(mod3) #Get fitted values (covariate effects + random effects)
resid.mod3  <- resid(mod3)  #Get residuals

par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.5)
plot(x = mu.mod3, 
     y = resid.mod3,
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Homegeneity")
abline(h = 0, lty = 2)     


# Plot residuals vs each covariate (in and not in the model).
# Check for non-linear patterns in the residuals.
MyVar <- c("year", "site")
SYB$resid.3 <- resid.mod3
source("~/Documents/R-INLA Course/AllData/HighstatLibV13.R")
MyMultipanel.ggp2(Z = SYB, 
                  varx = MyVar, 
                  vary = "resid.3", 
                  ylab = "Residuals",
                  addSmoother = TRUE,
                  addRegressionLine = TRUE,
                  addHorizontalLine = TRUE)

par(mfrow = c(2,2))
plot(SYB$async, 
     SYB$resid.3,
     xlab = "async",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$mu, 
     SYB$resid.3,
     xlab = "mu",
     ylab = "Residuals")
abline(h = 0, lty = 2)
plot(SYB$logheight.scal, 
     SYB$resid.3,
     xlab = "logheight.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2) 
plot(SYB$sigma.scal, 
     SYB$resid.3,
     xlab = "sigma.scal",
     ylab = "Residuals")
abline(h = 0, lty = 2)


# Plot the residuals versus year and site.
par(mfrow=c(1,2))
boxplot(resid.3 ~ year, data = SYB)
boxplot(resid.3 ~ site, SYB)


# Visualise the residuals on a map
SYB$Sign.3 <- ifelse(SYB$resid.3 >=0, "positive", "negative") #red for positive residuals and
#black for negative residuals
SYB$Size.3 <-  round(5 * sqrt(abs(SYB$resid.3) /max(SYB$resid.3))) + 1 #Size of residual


p2 <- ggplot()
p2 <- p2 + theme(text = element_text(size=13)) 
p2 <- p2 + geom_jitter(aes(x = Longitude, 
                           y = Latitude,
                           size = Size.3,
                           col = Sign.3),
                       shape = 1,
                       data = SYB) 
p2 <- p2 + xlab("Longitude") + ylab("Latitude")  
p2


### cubic models

mod1.cubic <- MCMCglmm(v2biomass~async+I(async^2)+I(async^3)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=40000, thin=20)
summary(mod1.cubic)
plot(mod1.cubic)

mod2.cubic <- MCMCglmm(v2mass~async+I(async^2)+I(async^3)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=40000, thin=20)
summary(mod2.cubic)
plot(mod2.cubic)

mod3.cubic <- MCMCglmm(fledged~async+I(async^2)+I(async^3)+logheight.scal+sigma.scal+v2age, random=~site + year + SY, data=SYB, burnin=10000, nitt=40000, thin=20)
summary(mod3.cubic)
plot(mod3.cubic)



async <- seq(-30,30,0.2)
biomass.avpeak <- mean(mod1.cubic$Sol[,"logheight.scal"])*0 + mean(mod1.cubic$Sol[,"sigma.scal"])*0 + mean(mod1.cubic$Sol[,"(Intercept)"]) + mean(mod1.cubic$Sol[,"async"])*async + mean(mod1.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod1.cubic$Sol[,"I(async^3)"])*(async^3) 
biomass.highpeak <- mean(mod1.cubic$Sol[,"logheight.scal"])*2 + mean(mod1.cubic$Sol[,"sigma.scal"])*0 + mean(mod1.cubic$Sol[,"(Intercept)"]) + mean(mod1.cubic$Sol[,"async"])*async + mean(mod1.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod1.cubic$Sol[,"I(async^3)"])*(async^3)
biomass.lowpeak <- mean(mod1.cubic$Sol[,"logheight.scal"])*-2 + mean(mod1.cubic$Sol[,"sigma.scal"])*0 + mean(mod1.cubic$Sol[,"(Intercept)"]) + mean(mod1.cubic$Sol[,"async"])*async + mean(mod1.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod1.cubic$Sol[,"I(async^3)"])*(async^3)
biomass.widepeak <- mean(mod1.cubic$Sol[,"logheight.scal"])*0 + mean(mod1.cubic$Sol[,"sigma.scal"])*2 + mean(mod1.cubic$Sol[,"(Intercept)"]) + mean(mod1.cubic$Sol[,"async"])*async + mean(mod1.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod1.cubic$Sol[,"I(async^3)"])*(async^3)
biomass.narrowpeak <- mean(mod1.cubic$Sol[,"logheight.scal"])*0 + mean(mod1.cubic$Sol[,"sigma.scal"])*-2 + mean(mod1.cubic$Sol[,"(Intercept)"]) + mean(mod1.cubic$Sol[,"async"])*async + mean(mod1.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod1.cubic$Sol[,"I(async^3)"])*(async^3)
plot(SYB$async, SYB$v2biomass, pch=20, col=mycolblack)
points(async, biomass.avpeak, type="l", col="black")
points(async, biomass.lowpeak, type="l", col="blue", lty="dotdash")
points(async, biomass.highpeak, type="l", col="red", lty="dotdash")
points(async, biomass.narrowpeak, type="l", col="blue", lty="dashed")
points(async, biomass.widepeak, type="l", col="red", lty="dashed")

mass.avpeak <- mean(mod2.cubic$Sol[,"logheight.scal"])*0 + mean(mod2.cubic$Sol[,"sigma.scal"])*0 + mean(mod2.cubic$Sol[,"(Intercept)"]) + mean(mod2.cubic$Sol[,"async"])*async + mean(mod2.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod2.cubic$Sol[,"I(async^3)"])*(async^3)
mass.highpeak <- mean(mod2.cubic$Sol[,"logheight.scal"])*2 + mean(mod2.cubic$Sol[,"sigma.scal"])*0 + mean(mod2.cubic$Sol[,"(Intercept)"]) + mean(mod2.cubic$Sol[,"async"])*async + mean(mod2.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod2.cubic$Sol[,"I(async^3)"])*(async^3)
mass.lowpeak <- mean(mod2.cubic$Sol[,"logheight.scal"])*-2 + mean(mod2.cubic$Sol[,"sigma.scal"])*0 + mean(mod2.cubic$Sol[,"(Intercept)"]) + mean(mod2.cubic$Sol[,"async"])*async + mean(mod2.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod2.cubic$Sol[,"I(async^3)"])*(async^3)
mass.widepeak <- mean(mod2.cubic$Sol[,"logheight.scal"])*0 + mean(mod2.cubic$Sol[,"sigma.scal"])*2 + mean(mod2.cubic$Sol[,"(Intercept)"]) + mean(mod2.cubic$Sol[,"async"])*async + mean(mod2.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod2.cubic$Sol[,"I(async^3)"])*(async^3)
mass.narrowpeak <- mean(mod2.cubic$Sol[,"logheight.scal"])*0 + mean(mod2.cubic$Sol[,"sigma.scal"])*-2 + mean(mod2.cubic$Sol[,"(Intercept)"]) + mean(mod2.cubic$Sol[,"async"])*async + mean(mod2.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod2.cubic$Sol[,"I(async^3)"])*(async^3)
plot(SYB$async, SYB$v2mass, pch=20, col=mycolblack)
points(async, mass.avpeak, type="l", col="black")
points(async, mass.lowpeak, type="l", col="blue", lty="dotdash")
points(async, mass.highpeak, type="l", col="red", lty="dotdash")
points(async, mass.narrowpeak, type="l", col="blue", lty="dashed")
points(async, mass.widepeak, type="l", col="red", lty="dashed")

fledged.avpeak <- mean(mod3.cubic$Sol[,"logheight.scal"])*0 + mean(mod3.cubic$Sol[,"sigma.scal"])*0 + mean(mod3.cubic$Sol[,"(Intercept)"]) + mean(mod3.cubic$Sol[,"async"])*async + mean(mod3.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod3.cubic$Sol[,"I(async^3)"])*(async^3) 
fledged.highpeak <- mean(mod3.cubic$Sol[,"logheight.scal"])*2 + mean(mod3.cubic$Sol[,"sigma.scal"])*0 + mean(mod3.cubic$Sol[,"(Intercept)"]) + mean(mod3.cubic$Sol[,"async"])*async + mean(mod3.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod3.cubic$Sol[,"I(async^3)"])*(async^3) 
fledged.lowpeak <- mean(mod3.cubic$Sol[,"logheight.scal"])*-2 + mean(mod3.cubic$Sol[,"sigma.scal"])*0 + mean(mod3.cubic$Sol[,"(Intercept)"]) + mean(mod3.cubic$Sol[,"async"])*async + mean(mod3.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod3.cubic$Sol[,"I(async^3)"])*(async^3) 
fledged.widepeak <- mean(mod3.cubic$Sol[,"logheight.scal"])*0 + mean(mod3.cubic$Sol[,"sigma.scal"])*2 + mean(mod3.cubic$Sol[,"(Intercept)"]) + mean(mod3.cubic$Sol[,"async"])*async + mean(mod3.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod3.cubic$Sol[,"I(async^3)"])*(async^3) 
fledged.narrowpeak <- mean(mod3.cubic$Sol[,"logheight.scal"])*0 + mean(mod3.cubic$Sol[,"sigma.scal"])*-2 + mean(mod3.cubic$Sol[,"(Intercept)"]) + mean(mod3.cubic$Sol[,"async"])*async + mean(mod3.cubic$Sol[,"I(async^2)"])*(async^2) + mean(mod3.cubic$Sol[,"I(async^3)"])*(async^3) 
plot(SYB$async, SYB$fledged, pch=20, col=mycolblack)
points(async, fledged.avpeak, type="l", col="black")
points(async, fledged.lowpeak, type="l", col="blue", lty="dotdash")
points(async, fledged.highpeak, type="l", col="red", lty="dotdash")
points(async, fledged.narrowpeak, type="l", col="blue", lty="dashed")
points(async, fledged.widepeak, type="l", col="red", lty="dashed")
