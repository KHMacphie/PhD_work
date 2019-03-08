######################
#### Biogeography ####
######################

rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe #####

site <- read.csv("Dropbox/master_data/site/site_details.csv")
colnames(site)[4:6] <- c("latitude", "longitude", "elevation")

cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
cater$year <- as.factor(cater$year)

pmatch(cater$site,site$site,duplicates.ok=TRUE)
all_data<- merge(cater, site, by="site", duplicates.ok=TRUE)

all_data$treespecies <- all_data$`tree species`
all_data$yearsite<- paste(all_data$site, all_data$year)
all_data$sitetree <- paste(all_data$site, all_data$tree)
all_data$siteday <- paste(all_data$site, all_data$date, all_data$year)
all_data$obs<-as.factor(seq(1,length(all_data[,1])))
all_data$sitetreesp <- paste(all_data$site, all_data$treespecies)
all_data$datecentred <- all_data$date-mean(all_data$date)
######################################
#### Model with date*biogeography ####
######################################

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

Biogeog<- MCMCglmm(caterpillars~datecentred*year+latitude*datecentred+elevation*datecentred+I(datecentred^2), random=~site+sitetree+siteday+recorder, family="poisson", data=all_data, prior=prior2, nitt=250000, burnin=25000)
save(Biogeog, file = "~/Documents/Models/Biogeog.RData")
load("~/Documents/Models/Biogeog.RData") 

summary(Biogeog)

plot(Biogeog$VCV) #random effects
plot(Biogeog$Sol) #fixedeffects
autocorr(Biogeog$Sol) #looking to level of auto correlation in fixed variables

#check if model generates sensible results
Biogeog.Sim<-simulate(Biogeog,nsim=100)
sum(all_data$caterpillars)
par(mfcol=c(1,1))
hist(apply(Biogeog.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(Biogeog.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")

### Plotting curves across elevation and latitude for 2014 at mean elv/lat
# latitude: 55.98-57.89, mean=56.93
# elevation 10-433, mean=152.64
centpredday <- seq(-23,23,0.5)
biogeog <- data.frame(xcentred=centpredday)
biogeog$L56 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*152)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*152)*centpredday

biogeog$L56.5 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56.5)+
  mean(Biogeog$Sol[,"elevation"]*152)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+
  mean(Biogeog$Sol[,"datecentred:latitude"]*56.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*152)*centpredday

biogeog$L57 <- mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*152)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*152)*centpredday

biogeog$L57.5 <- mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*152)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*152)*centpredday

biogeog$L57.9 <- mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.9)+
  mean(Biogeog$Sol[,"elevation"]*152)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.9)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*152)*centpredday

plot(biogeog$xcentred, biogeog$L56, type="l")
points(biogeog$xcentred, biogeog$L56.5, type="l")
points(biogeog$xcentred, biogeog$L57, type="l")
points(biogeog$xcentred, biogeog$L57.5, type="l")
points(biogeog$xcentred, biogeog$L57.9, type="l")

plot(biogeog$xcentred, exp(biogeog$L56), type="l")
points(biogeog$xcentred, exp(biogeog$L56.5), type="l")
points(biogeog$xcentred, exp(biogeog$L57), type="l")
points(biogeog$xcentred, exp(biogeog$L57.5), type="l")
points(biogeog$xcentred, exp(biogeog$L57.9), type="l")
#10 to 430
biogeog$E10 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"])*56.9+
  mean(Biogeog$Sol[,"elevation"])*10+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"])*56.9*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"])*10*centpredday

biogeog$E170 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"])*56.9+
  mean(Biogeog$Sol[,"elevation"])*170+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"])*56.9*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"])*170*centpredday

biogeog$E90 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"])*56.9+
  mean(Biogeog$Sol[,"elevation"])*90+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"])*56.9*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"])*90*centpredday

biogeog$E250 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"])*56.9+
  mean(Biogeog$Sol[,"elevation"])*250+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"])*56.9*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"])*250*centpredday

biogeog$E330 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56.9)+
  mean(Biogeog$Sol[,"elevation"]*330)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56.9)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*330)*centpredday

biogeog$E410 <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"])*56.9+
  mean(Biogeog$Sol[,"elevation"])*410+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"])*56.9*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"])*410*centpredday


plot(biogeog$xcentred, exp(biogeog$E410), type="l")
points(biogeog$xcentred, exp(biogeog$E10), type="l")
points(biogeog$xcentred, exp(biogeog$E90), type="l")
points(biogeog$xcentred, exp(biogeog$E170), type="l")
points(biogeog$xcentred, exp(biogeog$E250), type="l")
points(biogeog$xcentred, exp(biogeog$E330), type="l")
points(biogeog$xcentred, exp(biogeog$E410), type="l")

plot(biogeog$xcentred, biogeog$E410, type="l")
points(biogeog$xcentred, biogeog$E10, type="l")
points(biogeog$xcentred, biogeog$E90, type="l")
points(biogeog$xcentred, biogeog$E170, type="l")
points(biogeog$xcentred, biogeog$E250, type="l")
points(biogeog$xcentred, biogeog$E330, type="l")

##### plot peak date by lat/elev and then height by lat/elev

PHCI <- data.frame(Lat=c(56, 56.5, 57, 57.5, 57.9, 56.9,56.9,56.9,56.9,56.9), 
                   Elev=c(152, 152, 152, 152, 152, 40, 130, 220, 310, 400))

for(i in 1:length(PHCI$Lat)){
PHCI$Peak[i] <- mean(-((Biogeog$Sol[,"datecentred"])+(PHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(PHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                      (2*Biogeog$Sol[,"I(datecentred^2)"]))
}

for(i in 1:length(PHCI$Lat)){
PHCI$Height[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                         Biogeog$Sol[,"datecentred"]*PHCI$Peak[i]+
                         Biogeog$Sol[,"latitude"]*PHCI$Lat[i]+
                         Biogeog$Sol[,"elevation"]*PHCI$Elev[i]+
                         Biogeog$Sol[,"I(datecentred^2)"]*(PHCI$Peak[i]^2)+ 
                         Biogeog$Sol[,"datecentred:latitude"]*PHCI$Lat*PHCI$Peak[i]+
                         Biogeog$Sol[,"datecentred:elevation"]*PHCI$Elev*PHCI$Peak[i]))
}

for(i in 1:length(PHCI$Lat)){
  PHCI$logHeight[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                               Biogeog$Sol[,"datecentred"]*PHCI$Peak[i]+
                               Biogeog$Sol[,"latitude"]*PHCI$Lat[i]+
                               Biogeog$Sol[,"elevation"]*PHCI$Elev[i]+
                               Biogeog$Sol[,"I(datecentred^2)"]*(PHCI$Peak[i]^2)+ 
                               Biogeog$Sol[,"datecentred:latitude"]*PHCI$Lat*PHCI$Peak[i]+
                               Biogeog$Sol[,"datecentred:elevation"]*PHCI$Elev*PHCI$Peak[i])
}

for(i in 1:length(PHCI$Lat)){
  PHCI$Height3[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                              Biogeog$Sol[,"datecentred"]*7+
                              Biogeog$Sol[,"latitude"]*PHCI$Lat[i]+
                              Biogeog$Sol[,"elevation"]*PHCI$Elev[i]+
                              Biogeog$Sol[,"I(datecentred^2)"]*(7^2)+ 
                              Biogeog$Sol[,"datecentred:latitude"]*PHCI$Lat*7+
                              Biogeog$Sol[,"datecentred:elevation"]*PHCI$Elev*7))
}

### CIs

for(i in 1:length(PHCI$Lat)) {   # loop for CIs
  A <- HPDinterval((-((Biogeog$Sol[,"datecentred"])+(PHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(PHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                      (2*Biogeog$Sol[,"I(datecentred^2)"])))
  PHCI$Plowci[i] <- A["var1","lower"] 
  PHCI$Pupci[i] <- A["var1","upper"] 
}

for(i in 1:length(PHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*PHCI$Peak[i]+
                     Biogeog$Sol[,"latitude"]*PHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*PHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(PHCI$Peak[i]^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*PHCI$Lat*PHCI$Peak[i]+
                     Biogeog$Sol[,"datecentred:elevation"]*PHCI$Elev*PHCI$Peak[i])
  PHCI$Hlowci[i] <- exp(A["var1","lower"]) 
  PHCI$Hupci[i] <- exp(A["var1","upper"]) 
}

#### Each seperately
# Latitude
LatPHCI <- data.frame(Lat=c(56, 56.5, 57, 57.5, 57.9), 
                   Elev=c(152, 152, 152, 152, 152))

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$Peak[i] <- mean(-((Biogeog$Sol[,"datecentred"])+(LatPHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(LatPHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                         (2*Biogeog$Sol[,"I(datecentred^2)"]))
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$Height[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                               Biogeog$Sol[,"datecentred"]*LatPHCI$Peak[i]+
                               Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                               Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                               Biogeog$Sol[,"I(datecentred^2)"]*(LatPHCI$Peak[i]^2)+ 
                               Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*LatPHCI$Peak[i]+
                               Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*LatPHCI$Peak[i]))
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$logHeight[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                              Biogeog$Sol[,"datecentred"]*LatPHCI$Peak[i]+
                              Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                              Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                              Biogeog$Sol[,"I(datecentred^2)"]*(LatPHCI$Peak[i]^2)+ 
                              Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*LatPHCI$Peak[i]+
                              Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*LatPHCI$Peak[i])
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$Height3[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                                Biogeog$Sol[,"datecentred"]*2+
                                Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                                Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                                Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                                Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*2+
                                Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*2))
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$logHeight2[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                   Biogeog$Sol[,"datecentred"]*2+
                                   Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                                   Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                                   Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                                   Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*2+
                                   Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*2)
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$logHeightneg2[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                  Biogeog$Sol[,"datecentred"]*-2+
                                  Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                                  Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                                  Biogeog$Sol[,"I(datecentred^2)"]*(-2^2)+ 
                                  Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*-2+
                                  Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*-2)
}

for(i in 1:length(LatPHCI$Lat)){
  LatPHCI$logHeight6[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                  Biogeog$Sol[,"datecentred"]*6+
                                  Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                                  Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                                  Biogeog$Sol[,"I(datecentred^2)"]*(6^2)+ 
                                  Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*6+
                                  Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*6)
}


### CIs

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval((-((Biogeog$Sol[,"datecentred"])+(LatPHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(LatPHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                      (2*Biogeog$Sol[,"I(datecentred^2)"])))
  LatPHCI$Plowci[i] <- A["var1","lower"] 
  LatPHCI$Pupci[i] <- A["var1","upper"] 
}

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*LatPHCI$Peak[i]+
                     Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(LatPHCI$Peak[i]^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*LatPHCI$Peak[i]+
                     Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*LatPHCI$Peak[i])
  LatPHCI$Hlowci[i] <- exp(A["var1","lower"]) 
  LatPHCI$Hupci[i] <- exp(A["var1","upper"]) 
}

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*2+
                     Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*2+
                     Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*2)
  LatPHCI$H3lowci[i] <- exp(A["var1","lower"]) 
  LatPHCI$H3upci[i] <- exp(A["var1","upper"]) 
}

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*2+
                     Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*2+
                     Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*2)
  LatPHCI$logH2lowci[i] <- A["var1","lower"]
  LatPHCI$logH2upci[i] <- A["var1","upper"]
}

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*-2+
                     Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(-2^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*-2+
                     Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*-2)
  LatPHCI$logHneg2lowci[i] <- A["var1","lower"]
  LatPHCI$logHneg2upci[i] <- A["var1","upper"]
}

for(i in 1:length(LatPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*6+
                     Biogeog$Sol[,"latitude"]*LatPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*LatPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(6^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*LatPHCI$Lat[i]*6+
                     Biogeog$Sol[,"datecentred:elevation"]*LatPHCI$Elev[i]*6)
  LatPHCI$logH6lowci[i] <- A["var1","lower"]
  LatPHCI$logH6upci[i] <- A["var1","upper"]
}

# Elevation
ElevPHCI <- data.frame(Lat=c(56.5,56.5,56.5,56.5,56.5), 
                   Elev=c(40, 130, 220, 310, 400))

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$Peak[i] <- mean(-((Biogeog$Sol[,"datecentred"])+(ElevPHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(ElevPHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                            (2*Biogeog$Sol[,"I(datecentred^2)"]))
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$Height[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                                  Biogeog$Sol[,"datecentred"]*ElevPHCI$Peak[i]+
                                  Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                  Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                  Biogeog$Sol[,"I(datecentred^2)"]*(ElevPHCI$Peak[i]^2)+ 
                                  Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*ElevPHCI$Peak[i]+
                                  Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*ElevPHCI$Peak[i]))
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$logHeight[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                 Biogeog$Sol[,"datecentred"]*ElevPHCI$Peak[i]+
                                 Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                 Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                 Biogeog$Sol[,"I(datecentred^2)"]*(ElevPHCI$Peak[i]^2)+ 
                                 Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*ElevPHCI$Peak[i]+
                                 Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*ElevPHCI$Peak[i])
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$Height3[i] <- exp(mean(Biogeog$Sol[,"(Intercept)"]+
                                   Biogeog$Sol[,"datecentred"]*2+
                                   Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                   Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                   Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                                   Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*2+
                                   Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*2))
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$logHeight2[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                    Biogeog$Sol[,"datecentred"]*2+
                                    Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                    Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                    Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                                    Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*2+
                                    Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*2)
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$logHeightneg2[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                   Biogeog$Sol[,"datecentred"]*-2+
                                   Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                   Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                   Biogeog$Sol[,"I(datecentred^2)"]*(-2^2)+ 
                                   Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*-2+
                                   Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*-2)
}

for(i in 1:length(ElevPHCI$Lat)){
  ElevPHCI$logHeight6[i] <- mean(Biogeog$Sol[,"(Intercept)"]+
                                   Biogeog$Sol[,"datecentred"]*6+
                                   Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                                   Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                                   Biogeog$Sol[,"I(datecentred^2)"]*(6^2)+ 
                                   Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*6+
                                   Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*6)
}

### CIs

for(i in 1:length(ElevPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval((-((Biogeog$Sol[,"datecentred"])+(ElevPHCI$Lat[i]*Biogeog$Sol[,"datecentred:latitude"])+(ElevPHCI$Elev[i]*Biogeog$Sol[,"datecentred:elevation"]))/
                      (2*Biogeog$Sol[,"I(datecentred^2)"])))
  ElevPHCI$Plowci[i] <- A["var1","lower"] 
  ElevPHCI$Pupci[i] <- A["var1","upper"] 
}

for(i in 1:length(ElevPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*ElevPHCI$Peak[i]+
                     Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(ElevPHCI$Peak[i]^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*ElevPHCI$Peak[i]+
                     Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*ElevPHCI$Peak[i])
  ElevPHCI$Hlowci[i] <- exp(A["var1","lower"]) 
  ElevPHCI$Hupci[i] <- exp(A["var1","upper"]) 
}

for(i in 1:length(ElevPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*2+
                     Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*2+
                     Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*2)
  ElevPHCI$H3lowci[i] <- exp(A["var1","lower"]) 
  ElevPHCI$H3upci[i] <- exp(A["var1","upper"]) 
}

for(i in 1:length(ElevPHCI$Lat)) {   # loop for CIs
  A <- HPDinterval(Biogeog$Sol[,"(Intercept)"]+
                     Biogeog$Sol[,"datecentred"]*2+
                     Biogeog$Sol[,"latitude"]*ElevPHCI$Lat[i]+
                     Biogeog$Sol[,"elevation"]*ElevPHCI$Elev[i]+
                     Biogeog$Sol[,"I(datecentred^2)"]*(2^2)+ 
                     Biogeog$Sol[,"datecentred:latitude"]*ElevPHCI$Lat[i]*2+
                     Biogeog$Sol[,"datecentred:elevation"]*ElevPHCI$Elev[i]*2)
  ElevPHCI$logH3lowci[i] <- A["var1","lower"]
  ElevPHCI$logH3upci[i] <- A["var1","upper"]
}

## reverse date centering
# all_data$datecentred <- all_data$date-mean(all_data$date)
LatPHCI$PeakDate <- LatPHCI$Peak+mean(all_data$date)
LatPHCI$PDlowci <- LatPHCI$Plowci+mean(all_data$date)
LatPHCI$PDupci <- LatPHCI$Pupci+mean(all_data$date)
ElevPHCI$PeakDate <- ElevPHCI$Peak+mean(all_data$date)
ElevPHCI$PDlowci <- ElevPHCI$Plowci+mean(all_data$date)
ElevPHCI$PDupci <- ElevPHCI$Pupci+mean(all_data$date)

LatTiming <- ggplot(LatPHCI, aes(x=Lat, y=PeakDate))+
  geom_point()+
  geom_smooth(stat="identity")+
  geom_errorbar(aes(ymax=PDupci, ymin=PDlowci, width=0))+
  theme_bw()+
  xlab("Latitude (°N)")+
  ylab("Peak Date")+
  ylim(145, 162)+
  theme(text = element_text(size=30)) 
LatTiming

ElevTiming <- ggplot(ElevPHCI, aes(x=Elev, y=PeakDate))+
  geom_point()+
  geom_smooth(stat="identity")+
  geom_errorbar(aes(ymax=PDupci, ymin=PDlowci, width=0))+
  theme_bw()+
  xlab("Elevation (m)")+
  ylab("Peak Date")+
  ylim(145, 162)+
  theme(text = element_text(size=30))
ElevTiming

# Latitude Peak gradient
(LatPHCI[5,3]-LatPHCI[1,3])/(LatPHCI[5,1]-LatPHCI[1,1]) # 3.100646
# Latitude B equivalent
-mean(Biogeog$Sol[,"datecentred:latitude"]/
    (2*Biogeog$Sol[,"I(datecentred^2)"])) # 3.100646
-HPDinterval(Biogeog$Sol[,"datecentred:latitude"]/
        (2*Biogeog$Sol[,"I(datecentred^2)"])) # 4.398177  1.793462

#Elevation Peak gradient
(ElevPHCI[4,3]-ElevPHCI[1,3])/(ElevPHCI[4,2]-ElevPHCI[1,2]) # 0.02942217
# Latitude B equivalent
-mean(Biogeog$Sol[,"datecentred:elevation"]/
        (2*Biogeog$Sol[,"I(datecentred^2)"])) # 0.02942217  
-HPDinterval(Biogeog$Sol[,"datecentred:elevation"]/
        (2*Biogeog$Sol[,"I(datecentred^2)"])) # 0.0359856   0.02282919  

# Latitude height gradient
(LatPHCI[5,19]-LatPHCI[1,19])/(LatPHCI[5,1]-LatPHCI[1,1]) #-0.2732204   day 2
(LatPHCI[5,20]-LatPHCI[1,20])/(LatPHCI[5,1]-LatPHCI[1,1]) #-0.3809531   day -2
(LatPHCI[5,21]-LatPHCI[1,21])/(LatPHCI[5,1]-LatPHCI[1,1]) #-0.1654877   day 6
# Elevation height gradient
(ElevPHCI[5,19]-ElevPHCI[1,19])/(ElevPHCI[5,2]-ElevPHCI[1,2]) #-1.696034e-05 day 2
(ElevPHCI[5,20]-ElevPHCI[1,20])/(ElevPHCI[5,2]-ElevPHCI[1,2]) #-0.001039616  day -2
(ElevPHCI[5,21]-ElevPHCI[1,21])/(ElevPHCI[5,2]-ElevPHCI[1,2]) #0.001005696   day 6


### the relationship between biogeography and abundance is dependent on which day you look at

#mean(Biogeog$Sol[,"datecentred"])  negative slope
#[1] -1.53604
#> mean(Biogeog$Sol[,"elevation"]) seems insignificant to the relationship
#[1] -0.0005282883
#> mean(Biogeog$Sol[,"datecentred:elevation"]) becomes more positive with increasing elevation so beneath a certain elevation its a negative relationship but then shifts becomign more positive 
#[1] 0.000255664

## plot for curves in each year at high and low lat/elev
# Lat low: 56.5, Lat high: 57.5   elev mean: 150
# Elev low: 50, Elev high: 350    lat mean: 57

Lat14.low <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Lat14.high <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Elev14.low <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*50)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*50)*centpredday

Elev14.high <-  mean(Biogeog$Sol[,"(Intercept)"])+
  mean(Biogeog$Sol[,"datecentred"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*350)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*350)*centpredday

Lat15.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2015"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2015"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Lat15.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2015"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2015"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Elev15.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2015"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2015"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*50)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*50)*centpredday

Elev15.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2015"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2015"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*350)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*350)*centpredday

Lat16.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2016"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2016"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Lat16.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2016"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2016"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Elev16.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2016"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2016"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*50)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*50)*centpredday

Elev16.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2016"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2016"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*350)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*350)*centpredday

Lat17.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2017"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2017"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Lat17.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2017"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2017"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Elev17.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2017"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2017"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*50)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*50)*centpredday

Elev17.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2017"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2017"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*350)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*350)*centpredday

Lat18.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2018"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2018"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*56)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*56)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Lat18.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2018"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2018"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57.5)+
  mean(Biogeog$Sol[,"elevation"]*150)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57.5)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*150)*centpredday

Elev18.low <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2018"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2018"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*50)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*50)*centpredday

Elev18.high <-  mean(Biogeog$Sol[,"(Intercept)"]+Biogeog$Sol[,"year2018"])+
  mean(Biogeog$Sol[,"datecentred"]+Biogeog$Sol[,"datecentred:year2018"])*centpredday+
  mean(Biogeog$Sol[,"latitude"]*57)+
  mean(Biogeog$Sol[,"elevation"]*350)+
  mean(Biogeog$Sol[,"I(datecentred^2)"])*centpredday^2+ 
  mean(Biogeog$Sol[,"datecentred:latitude"]*57)*centpredday+
  mean(Biogeog$Sol[,"datecentred:elevation"]*350)*centpredday

predday <- seq(125,171,0.5)
# Latitude figures
par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(predday, exp(Lat14.low), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.17), axes=F) #black
points(predday, exp(Lat15.low), type="l",lwd=3, col=2) #red
points(predday, exp(Lat16.low), type="l",lwd=3,col=3) #green
points(predday, exp(Lat17.low), type="l", lwd=3, col=4)
points(predday, exp(Lat18.low), type="l", lwd=3, col=5)
legend("topright", legend="56°", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(predday, exp(Lat14.high), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.17), axes=F) #black
points(predday, exp(Lat15.high), type="l",lwd=3, col=2) #red
points(predday, exp(Lat16.high), type="l",lwd=3,col=3) #green
points(predday, exp(Lat17.high), type="l", lwd=3, col=4)
points(predday, exp(Lat18.high), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="57.5°", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
title(ylab="Caterpillar Abundance", outer=TRUE, line = 2)
title( xlab="Date", outer=TRUE, line = 0)


# Elevation figures
par(mfcol=c(2,1), mar=c(2,1,1,1), oma=c(2,4,0,0), cex=1.8)
plot(predday, exp(Elev14.low), type="l", lwd=3, xlab=NA, ylab=NA, ylim=c(0,0.17), axes=F) #black
points(predday, exp(Elev15.low), type="l",lwd=3, col=2) #red
points(predday, exp(Elev16.low), type="l",lwd=3,col=3) #green
points(predday, exp(Elev17.low), type="l", lwd=3, col=4)
points(predday, exp(Elev18.low), type="l", lwd=3, col=5)
legend("topright", legend="50m", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)

plot(predday, exp(Elev14.high), type="l", lwd=3,  ylab=NA, xlab=NA, ylim=c(0,0.17), axes=F) #black
points(predday, exp(Elev15.high), type="l",lwd=3, col=2) #red
points(predday, exp(Elev16.high), type="l",lwd=3,col=3) #green
points(predday, exp(Elev17.high), type="l", lwd=3, col=4)
points(predday, exp(Elev18.high), type="l", lwd=3, col=5)
legend("topleft", legend=c("2014","2015", "2016", "2017", "2018"), lty=c(1,1,1,1,1), lwd=3, col=c(1,2,3,4,5), cex=0.85, seg.len=0.8)
legend("topright", legend="350m", bty="n")
box()
axis(side = 1, tck = -.015, labels = NA)
axis(side = 2, tck = -.015, labels = NA)
axis(side = 1, lwd = 0, line = -.4)
axis(side = 2, lwd = 0, line = -.4, las = 1)
title(ylab="Caterpillar Abundance", outer=TRUE, line = 2)
title( xlab="Date", outer=TRUE, line = 0)
