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

####################################################
#### Organising Habitat and TreeTaxa Categories ####
####################################################

# Categories: Alder, Ash, Birch, Beech, Elm, Hazel, Oak, Rowan, Sycamore, Willow
# Work out amount of foliage missed without other deciduous and conifers

#### Habitat data set ####
habitat <- read.csv("~/Dropbox/2018/Habitats_2018_all.csv")
table(habitat$site) #checking all there

#### Stands, thickets and small trees combined as small trees and other deciduous combined
#converting 6 stands into small trees and combining oth.decid and conifers (combining pine, yew, juniper and conifer)

# combining other deciduous
habitat$X6_oth.decid <- rowSums(habitat[,c("X6_cherry", "X6_elder", "X6_holly", "X6_rose")], na.rm=TRUE)
habitat$X21_oth.decid <- rowSums(habitat[,c("X21_blackthorn", "X21_holly")], na.rm=TRUE)
habitat$s_oth.decid <- rowSums(habitat[,c("s_aspen", "s_cherry", "s_chestnut", "s_elder", "s_hawthorn", "s_holly", "s_lime", "s_whitebeam", "s_other")], na.rm=TRUE)
habitat$m_oth.decid <- rowSums(habitat[,c("m_aspen", "m_cherry", "m_chestnut", "m_holly", "m_lime", "m_whitebeam", "m_horsechestnut")], na.rm=TRUE)
habitat$l_oth.decid <- habitat$l_cherry
habitat$z_oth.decid <- rowSums(habitat[,c("z_blackthorn", "z_cherry", "z_other")], na.rm=TRUE)

# combining conifers
habitat$X6_conifer_all <- habitat$X6_juniper
habitat$s_conifer_all <- rowSums(habitat[,c("s_conifer", "s_pine", "s_yew", "s_juniper")], na.rm=TRUE)
habitat$m_conifer_all <- rowSums(habitat[,c("m_conifer", "m_pine", "m_yew")], na.rm=TRUE)
habitat$l_conifer_all <- rowSums(habitat[,c("l_conifer", "l_pine")], na.rm=TRUE)

# converting 6 stands into small trees
habitat$alder6smt <- habitat$X6_alder*0.5
habitat$ash6smt <- habitat$X6_ash*0.5
habitat$birch6smt <- habitat$X6_birch*0.5
habitat$beech6smt <- habitat$X6_beech*0.5
habitat$hazel6smt <- habitat$X6_hazel*0.5
habitat$rowan6smt <- habitat$X6_rowan*0.5
habitat$syc6smt <- habitat$X6_sycamore*0.5
habitat$willow6smt <- habitat$X6_willow*0.5
habitat$oth.decid6smt <- habitat$X6_oth.decid*0.5
habitat$conifer6smt <- habitat$X6_conifer_all*0.5
# thickets already as small trees

# combining stands, thickets and small trees
habitat$Alder_S <- rowSums(habitat[,c("s_alder", "X21_alder", "alder6smt")], na.rm=TRUE)
habitat$Ash_S <- rowSums(habitat[,c("s_ash", "ash6smt")], na.rm=TRUE)
habitat$Beech_S <- rowSums(habitat[,c("s_beech", "beech6smt")], na.rm=TRUE)
habitat$Birch_S <- rowSums(habitat[,c("s_birch", "X21_birch", "birch6smt")], na.rm=TRUE)
habitat$Hazel_S <- rowSums(habitat[,c("s_hazel", "X21_hazel", "hazel6smt")], na.rm=TRUE)
habitat$Rowan_S <- rowSums(habitat[,c("s_rowan", "rowan6smt")], na.rm=TRUE)
habitat$Sycamore_S <- rowSums(habitat[,c("s_sycamore", "syc6smt")], na.rm=TRUE)
habitat$Willow_S <- rowSums(habitat[,c("s_willow", "willow6smt", "X21_willow", "z_willow")], na.rm=TRUE)
habitat$Conifer_S <- rowSums(habitat[,c("X6_conifer_all", "s_conifer_all")], na.rm=TRUE)
habitat$OthDecid_S <- rowSums(habitat[,c("X6_oth.decid", "X21_oth.decid", "s_oth.decid", "z_oth.decid")], na.rm=TRUE)

##### Use the weighting function to calculate foliage per tree type and size at each NB

## making a function for the weighting equation
weighting <- function(x){return(pi*(x/(2*pi))^2)}  ##weighted by cross secitonal area of the trunk in cm^2 (min possible size for each category)
Sm <- weighting(40)
Med <- weighting(100)
Lar <- weighting(250)
S=1
M=Med/Sm
L=Lar/Sm

Habitat_byNB <- data.frame(Site=habitat$site, 
                           NB=habitat$nestbox, 
                           Alder_S=(S*habitat$Alder_S),
                           Alder_M=(M*habitat$m_alder),
                           Ash_L=(L*habitat$l_ash),
                           Ash_M=(M*habitat$m_ash),
                           Ash_S=(S*habitat$Ash_S),
                           Beech_S=(S*habitat$Beech_S),
                           Beech_M=(M*habitat$m_beech),
                           Beech_L=(L*habitat$l_beech),
                           Birch_S=(S*habitat$Birch_S),
                           Birch_M=(M*habitat$m_birch),
                           Birch_L=(L*habitat$l_birch),
                           Elm_S=(S*habitat$s_elm),
                           Elm_M=(M*habitat$m_elm),
                           Elm_L=(L*habitat$l_elm),
                           Hazel_S=(S*habitat$Hazel_S),
                           Oak_S=(S*habitat$s_oak),
                           Oak_M=(M*habitat$m_oak),
                           Oak_L=(L*habitat$l_oak), 
                           Rowan_S=(S*habitat$Rowan_S),
                           Rowan_M=(M*habitat$m_rowan),
                           Sycamore_S=(S*habitat$Sycamore_S),
                           Sycamore_M=(M*habitat$m_sycamore),
                           Sycamore_L=(L*habitat$l_sycamore),
                           Willow_S=(S*habitat$Willow_S),
                           Willow_M=(M*habitat$m_willow),
                           Willow_L=(L*habitat$l_willow),
                           Conifer_S=(S*habitat$Conifer_S),
                           Conifer_M=(M*habitat$m_conifer_all),
                           Conifer_L=(L*habitat$l_conifer_all),
                           OthDecid_S=(S*habitat$OthDecid_S),  
                           OthDecid_M=(M*habitat$m_oth.decid),
                           OthDecid_L=(L*habitat$l_oth.decid))
Habitat_byNB[is.na(Habitat_byNB)] <- 0

#Foliage score for each tree type at each NB
Habitat_byNB$Alder_FS <- rowSums(Habitat_byNB[,c("Alder_S", "Alder_M")], na.rm=TRUE)
Habitat_byNB$Ash_FS <- rowSums(Habitat_byNB[,c("Ash_S", "Ash_M", "Ash_L")], na.rm=TRUE)
Habitat_byNB$Beech_FS <- rowSums(Habitat_byNB[,c("Beech_S", "Beech_M", "Beech_L")], na.rm=TRUE)
Habitat_byNB$Birch_FS <- rowSums(Habitat_byNB[,c("Birch_S", "Birch_M", "Birch_L")], na.rm=TRUE)
Habitat_byNB$Elm_FS <- rowSums(Habitat_byNB[,c("Elm_S", "Elm_M", "Elm_L")], na.rm=TRUE)
Habitat_byNB$Hazel_FS <- Habitat_byNB$Hazel_S
Habitat_byNB$Oak_FS <- rowSums(Habitat_byNB[,c("Oak_S", "Oak_M", "Oak_L")], na.rm=TRUE)
Habitat_byNB$Rowan_FS <- rowSums(Habitat_byNB[,c("Rowan_S", "Rowan_M")], na.rm=TRUE)
Habitat_byNB$Sycamore_FS <- rowSums(Habitat_byNB[,c("Sycamore_S", "Sycamore_M", "Sycamore_L")], na.rm=TRUE)
Habitat_byNB$Willow_FS <- rowSums(Habitat_byNB[,c("Willow_S", "Willow_M", "Willow_L")], na.rm=TRUE)
Habitat_byNB$Conifer_FS <- rowSums(Habitat_byNB[,c("Conifer_S", "Conifer_M", "Conifer_L")], na.rm=TRUE)
Habitat_byNB$OthDecid_FS <- rowSums(Habitat_byNB[,c("OthDecid_S", "OthDecid_M", "OthDecid_L")], na.rm=TRUE)

## combining NB within each site- mean to account for different number of NBs at sites
site <- read.csv("Dropbox/master_data/site/site_details.csv")
Habitat_Site <- Habitat_byNB[,35:46]
Habitat_Site$Site <- Habitat_byNB$Site
Habitat_Site <- aggregate(.~Site, Habitat_Site, mean)

# getting proportions of each tree category
Habitat_Site$Total <- rowSums(Habitat_Site[2:13])
Habitat_Site$Alder_prop <- Habitat_Site$Alder_FS/Habitat_Site$Total
Habitat_Site$Ash_prop <- Habitat_Site$Ash_FS/Habitat_Site$Total
Habitat_Site$Beech_prop <- Habitat_Site$Beech_FS/Habitat_Site$Total
Habitat_Site$Birch_prop <- Habitat_Site$Birch_FS/Habitat_Site$Total
Habitat_Site$Elm_prop <- Habitat_Site$Elm_FS/Habitat_Site$Total
Habitat_Site$Hazel_prop <- Habitat_Site$Hazel_FS/Habitat_Site$Total
Habitat_Site$Oak_prop <- Habitat_Site$Oak_FS/Habitat_Site$Total
Habitat_Site$Rowan_prop <- Habitat_Site$Rowan_FS/Habitat_Site$Total
Habitat_Site$Sycamore_prop <- Habitat_Site$Sycamore_FS/Habitat_Site$Total
Habitat_Site$Willow_prop <- Habitat_Site$Willow_FS/Habitat_Site$Total
Habitat_Site$Conifer_prop <- Habitat_Site$Conifer_FS/Habitat_Site$Total
Habitat_Site$OthDecid_prop <- Habitat_Site$OthDecid_FS/Habitat_Site$Total

# FS scaled
Habitat_Site[,27:38] <- (Habitat_Site[,2:13]/max(Habitat_Site[,2:13])) # Dividing all FS's by max FS
colnames(Habitat_Site)[27:38] <- c("Alder_Scld","Ash_Scld","Beech_Scld","Birch_Scld","Elm_Scld","Hazel_Scld","Oak_Scld","Rowan_Scld","Sycamore_Scld","Willow_Scld","Conifer_Scld","OthDecid_Scld")

##############################
##### TreeTaxa categories ####

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)
# Creating other deciduous
cater$tree.species<- revalue(cater$tree.species, c("?"="OthDecid", "Cherry"="OthDecid", "Aspen"="OthDecid", "Chestnut"="OthDecid", "Lime"="OthDecid", "Field maple"="OthDecid", "Damson"="OthDecid", "Whitebeam"="OthDecid"))

#####################
#### Adding mass ####

cater$caterpillar.mass <- revalue(cater$caterpillar.mass, c("0"=""))
cater$caterpillar.mass2 <- revalue(cater$caterpillar.mass, c("<0.01"="0.02", "0.01"="0.02"))
cater$caterpillar.mass2 <- as.numeric(as.character(cater$caterpillar.mass2))

#mass per caterpillar
cater$mpc2 <- cater$caterpillar.mass2/cater$caterpillars
cater$mpc1 <- as.character(cater$caterpillar.mass2)
cater$mpc1 <- revalue(cater$mpc1, c("0.02"="0.001"))
cater$mpc1 <- as.numeric(cater$mpc1)
cater$mpc1 <- cater$mpc1/cater$caterpillars
cater$mpc1 <- ifelse(cater$mpc1 < 0.001,0.001,cater$mpc1)

# Sort rows with mass.uncertain
which(cater$mass.uncertain=="1") # 11694 21835 22385 25555 25743 25876 26665 27967 28137
cater$notes[11694]
cater$mpc1[11694] <- NA
cater$mpc2[11694] <- NA
cater[21835,]
cater$mpc2[21835] <- cater$caterpillar.mass2[21835]
cater[22385,]
cater$mpc1[22385] <- cater$caterpillar.mass2[22385]/4
cater$mpc2[22385] <- cater$caterpillar.mass2[22385]/4
cater[25555,]
cater$mpc2[25555] <- cater$caterpillar.mass2[25555]
cater[25743,]
cater$mpc2[25743] <- cater$caterpillar.mass2[25743]
cater[25876,]
cater$mpc1[25876] <- cater$caterpillar.mass2[25876]/3
cater$mpc2[25876] <- cater$caterpillar.mass2[25876]/3
cater[26665,]
cater$mpc1[26665] <- cater$caterpillar.mass2[26665]/7
cater$mpc2[26665] <- cater$caterpillar.mass2[26665]/7
cater[27967,]
cater$mpc1[27967] <- cater$caterpillar.mass2[27967]
cater$mpc2[27967] <- cater$caterpillar.mass2[27967]
cater[28137,]
cater$mpc1[28137] <- cater$caterpillar.mass2[28137]/3
cater$mpc2[28137] <- cater$caterpillar.mass2[28137]/3

#interval censor and log
cater$logmpc1 <- log(cater$mpc1)
cater$logmpc2 <- log(cater$mpc2)

########################
#### Full dataframe ####
########################

Habitat_Site$site <- Habitat_Site$Site
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE)
cater_habitat$treeID <- paste(cater_habitat$tree, cater_habitat$site)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$date, cater_habitat$year)
cater_habitat$siteyear <- paste(cater_habitat$site, cater_habitat$year)
#cater_habitat$datescaled <- cater_habitat$date/max(cater_habitat$date) ## WRONG
cater_habitat$datescaled <- scale(cater_habitat$date) #unscale(x, center= 146.4095, scale=14.19835)   now from -2.071330 to 2.013653
mean(cater_habitat$date) # 146.4095
sd(cater_habitat$date) # 14.19835

#!!!!!!!!!!! if removing OthDecid !!!!!!!!!!!!
cater_habitat<- subset(cater_habitat, tree.species!="OthDecid")

#############################
#### Abundance TT curves ####
#############################

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#AbundTTCurves<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                                 random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):siteyear + treeID + siteday + recorder, 
#                                 family="poisson", data=cater_habitat, prior=prior3, nitt=300000, burnin=30000, pr=TRUE)
#save(AbundTTCurves, file = "~/Documents/Models/AbundTTCurves.RData")
#load("~/Documents/Models/AbundTTCurves.RData")

#AbundTTCurves1<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):siteyear + treeID + siteday + recorder, 
#                         family="poisson", data=cater_habitat, prior=prior3, nitt=750000, burnin=50000, pr=TRUE, thin=50)
#save(AbundTTCurves1, file = "~/Documents/Models/AbundTTCurves1.RData")
#load("~/Documents/Models/AbundTTCurves1.RData")

# was 1million, needs a bit longer
#AbundTTCurves3<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):siteyear + treeID + siteday + recorder, 
#                         family="poisson", data=cater_habitat, prior=prior3, nitt=1200000, burnin=50000, pr=TRUE, thin=30)
#save(AbundTTCurves3, file = "~/Documents/Models/AbundTTCurves3.RData")
#rm(list=ls())
#load("~/Documents/Models/AbundTTCurves3.RData")

#### !!!!!!! Using proper scaling !!!!!!! ####  and changed random regression

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#ATTCscaled<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):siteyear + year + treeID + siteday + recorder, 
#                         family="poisson", data=cater_habitat, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=30)
#save(ATTCscaled, file = "~/Documents/Models/ATTCscaled.RData")
#load("~/Documents/Models/ATTCscaled.RData")

#ATTCyearint<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + siteyear + treeID + siteday + recorder, 
#                         family="poisson", data=cater_habitat, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=30)
#save(ATTCyearint, file = "~/Documents/Models/ATTCyearint.RData")
load("~/Documents/Models/ATTCyearint.RData")

AlderCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,4],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,14], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,24])
AshCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,5],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,15], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,25])
BeechCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,6],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,16], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,26])
BirchCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,7],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,17], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,27])
ElmCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,8],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,18], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,28])
HazelCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,9],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,19], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,29])
OakCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,10],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,20], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,30])
RowanCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,11],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,21], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,31])
SycamoreCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,12],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,22], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,32])
WillowCurve <- data.frame(int=AbundTTCurves3$Sol[,1], TTint=AbundTTCurves3$Sol[,13],date=AbundTTCurves3$Sol[,2], TTdate=AbundTTCurves3$Sol[,23], date2=AbundTTCurves3$Sol[,3], TTdate2=AbundTTCurves3$Sol[,33])



dayscal <- seq(0.67,1,0.001)
curve <- mean(AbundTTCurves3$Sol[,1])+mean(AbundTTCurves3$Sol[,2])*dayscal+mean(AbundTTCurves3$Sol[,3])*dayscal^2
days <- dayscal*max(cater_habitat$date)
mycol <- rgb(0, 150, 0, max = 250, alpha = 10, names = NULL)
mycol2 <- rgb(150, 0, 0, max = 250, alpha = 10, names = NULL)
mycol3 <- rgb(0, 0, 150, max = 250, alpha = 10, names = NULL)
mycol4 <- rgb(75, 75, 0, max = 250, alpha = 10, names = NULL)
mycol5 <- rgb(0, 75, 75, max = 250, alpha = 10, names = NULL)
mycol6 <- rgb(75, 0, 75, max = 250, alpha = 10, names = NULL)
mycol7 <- rgb(50, 50, 50, max = 250, alpha = 10, names = NULL)
mycol8 <- rgb(100, 25, 25, max = 250, alpha = 10, names = NULL)
mycol9 <- rgb(25, 100, 25, max = 250, alpha = 10, names = NULL)
mycol10 <- rgb(25, 25, 100, max = 250, alpha = 10, names = NULL)

par(mfcol=c(1,1),mar=c(3.9, 3.8, 1, 1), cex=1.4, las=1)
plot(days,exp(curve), type="l", ylim=c(0,0.15), xlab="Date", ylab="Abundance", yaxs="i")

for(i in 1:5400){
  A <- AlderCurve[i,1]+AlderCurve[i,2] +(AlderCurve[i,3]+AlderCurve[i,4])*dayscal+(AlderCurve[i,5]+AlderCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol, lwd=0.5)
}
for(i in 1:5400){
  A <- AshCurve[i,1]+AshCurve[i,2] +(AshCurve[i,3]+AshCurve[i,4])*dayscal+(AshCurve[i,5]+AshCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol2, lwd=0.5)
}
for(i in 1:5400){
  A <- BeechCurve[i,1]+BeechCurve[i,2] +(BeechCurve[i,3]+BeechCurve[i,4])*dayscal+(BeechCurve[i,5]+BeechCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol3, lwd=0.5)
}
for(i in 1:5400){
  A <- BirchCurve[i,1]+BirchCurve[i,2] +(BirchCurve[i,3]+BirchCurve[i,4])*dayscal+(BirchCurve[i,5]+BirchCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol4, lwd=0.5)
}
for(i in 1:5400){
  A <- ElmCurve[i,1]+ElmCurve[i,2] +(ElmCurve[i,3]+ElmCurve[i,4])*dayscal+(ElmCurve[i,5]+ElmCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol5, lwd=0.5)
}
for(i in 1:5400){
  A <- HazelCurve[i,1]+HazelCurve[i,2] +(HazelCurve[i,3]+HazelCurve[i,4])*dayscal+(HazelCurve[i,5]+HazelCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol6, lwd=0.5)
}
for(i in 1:5400){
  A <- OakCurve[i,1]+OakCurve[i,2] +(OakCurve[i,3]+OakCurve[i,4])*dayscal+(OakCurve[i,5]+OakCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol7, lwd=0.5)
}
for(i in 1:5400){
  A <- RowanCurve[i,1]+RowanCurve[i,2] +(RowanCurve[i,3]+RowanCurve[i,4])*dayscal+(RowanCurve[i,5]+RowanCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol8, lwd=0.5)
}
for(i in 1:5400){
  A <- SycamoreCurve[i,1]+SycamoreCurve[i,2] +(SycamoreCurve[i,3]+SycamoreCurve[i,4])*dayscal+(SycamoreCurve[i,5]+SycamoreCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol9, lwd=0.5)
}
for(i in 1:5400){
  A <- WillowCurve[i,1]+WillowCurve[i,2] +(WillowCurve[i,3]+WillowCurve[i,4])*dayscal+(WillowCurve[i,5]+WillowCurve[i,6])*dayscal^2
  points(days, exp(A), type="l", col=mycol10, lwd=0.5)
}

AlderLine <- mean(AlderCurve[,1]+AlderCurve[,2]) +mean(AlderCurve[,3]+AlderCurve[,4])*dayscal+mean(AlderCurve[,5]+AlderCurve[,6])*dayscal^2
AshLine <- mean(AshCurve[,1]+AshCurve[,2]) +mean(AshCurve[,3]+AshCurve[,4])*dayscal+mean(AshCurve[,5]+AshCurve[,6])*dayscal^2
BeechLine <- mean(BeechCurve[,1]+BeechCurve[,2]) +mean(BeechCurve[,3]+BeechCurve[,4])*dayscal+mean(BeechCurve[,5]+BeechCurve[,6])*dayscal^2
BirchLine <- mean(BirchCurve[,1]+BirchCurve[,2]) +mean(BirchCurve[,3]+BirchCurve[,4])*dayscal+mean(BirchCurve[,5]+BirchCurve[,6])*dayscal^2
ElmLine <- mean(ElmCurve[,1]+ElmCurve[,2]) +mean(ElmCurve[,3]+ElmCurve[,4])*dayscal+mean(ElmCurve[,5]+ElmCurve[,6])*dayscal^2
HazelLine <- mean(HazelCurve[,1]+HazelCurve[,2]) +mean(HazelCurve[,3]+HazelCurve[,4])*dayscal+mean(HazelCurve[,5]+HazelCurve[,6])*dayscal^2
OakLine <- mean(OakCurve[,1]+OakCurve[,2]) +mean(OakCurve[,3]+OakCurve[,4])*dayscal+mean(OakCurve[,5]+OakCurve[,6])*dayscal^2
RowanLine <- mean(RowanCurve[,1]+RowanCurve[,2]) +mean(RowanCurve[,3]+RowanCurve[,4])*dayscal+mean(RowanCurve[,5]+RowanCurve[,6])*dayscal^2
SycamoreLine <- mean(SycamoreCurve[,1]+SycamoreCurve[,2]) +mean(SycamoreCurve[,3]+SycamoreCurve[,4])*dayscal+mean(SycamoreCurve[,5]+SycamoreCurve[,6])*dayscal^2
WillowLine <- mean(WillowCurve[,1]+WillowCurve[,2]) +mean(WillowCurve[,3]+WillowCurve[,4])*dayscal+mean(WillowCurve[,5]+WillowCurve[,6])*dayscal^2

plot(days, exp(OakLine), type="l", col="#36C2C2", lwd=2, xlab="Date", ylab="Abundance")
points(days, exp(AlderLine), type="l", col="#C70808", lwd=2)
points(days, exp(AshLine), type="l", col="#FA3737", lwd=2)
points(days, exp(BeechLine), type="l", col="#FF792B", lwd=2)
points(days, exp(BirchLine), type="l", col="#FFAD14", lwd=2)
points(days, exp(ElmLine), type="l", col="#85C23A", lwd=2)
points(days, exp(HazelLine), type="l", col="#207A39", lwd=2)
points(days, exp(RowanLine), type="l", col="#2C69EB", lwd=2)
points(days, exp(SycamoreLine), type="l", col="#9A3BFF", lwd=2)
points(days, exp(WillowLine), type="l", col="#FF57FF", lwd=2)
#c("#C70808", "#FA3737", "#FF792B", "#FFAD14", "#85C23A", "#207A39", "#36C2C2", "#2C69EB", "#9A3BFF", "#FF57FF")

###################
#### Peak date ####    -b/2a
###################

AlderCurve$a <- AlderCurve[,5]+AlderCurve[,6]
AshCurve$a <- AshCurve[,5]+AshCurve[,6]
BeechCurve$a <- BeechCurve[,5]+BeechCurve[,6]
BirchCurve$a <- BirchCurve[,5]+BirchCurve[,6]
ElmCurve$a <- ElmCurve[,5]+ElmCurve[,6]
HazelCurve$a <- HazelCurve[,5]+HazelCurve[,6]
OakCurve$a <- OakCurve[,5]+OakCurve[,6]
RowanCurve$a <- RowanCurve[,5]+RowanCurve[,6]
SycamoreCurve$a <- SycamoreCurve[,5]+SycamoreCurve[,6]
WillowCurve$a <- WillowCurve[,5]+WillowCurve[,6]
AlderCurve$b <- AlderCurve[,3]+AlderCurve[,4]
AshCurve$b <- AshCurve[,3]+AshCurve[,4]
BeechCurve$b <- BeechCurve[,3]+BeechCurve[,4]
BirchCurve$b <- BirchCurve[,3]+BirchCurve[,4]
ElmCurve$b <- ElmCurve[,3]+ElmCurve[,4]
HazelCurve$b <- HazelCurve[,3]+HazelCurve[,4]
OakCurve$b <- OakCurve[,3]+OakCurve[,4]
RowanCurve$b <- RowanCurve[,3]+RowanCurve[,4]
SycamoreCurve$b <- SycamoreCurve[,3]+SycamoreCurve[,4]
WillowCurve$b <- WillowCurve[,3]+WillowCurve[,4]
AlderCurve$c <- AlderCurve[,1]+AlderCurve[,2]
AshCurve$c <- AshCurve[,1]+AshCurve[,2]
BeechCurve$c <- BeechCurve[,1]+BeechCurve[,2]
BirchCurve$c <- BirchCurve[,1]+BirchCurve[,2]
ElmCurve$c <- ElmCurve[,1]+ElmCurve[,2]
HazelCurve$c <- HazelCurve[,1]+HazelCurve[,2]
OakCurve$c <- OakCurve[,1]+OakCurve[,2]
RowanCurve$c <- RowanCurve[,1]+RowanCurve[,2]
SycamoreCurve$c <- SycamoreCurve[,1]+SycamoreCurve[,2]
WillowCurve$c <- WillowCurve[,1]+WillowCurve[,2]

AlderCurve$pd <- -AlderCurve$b/(2*AlderCurve$a)
AshCurve$pd <- -AshCurve$b/(2*AshCurve$a)
BeechCurve$pd <- -BeechCurve$b/(2*BeechCurve$a)
BirchCurve$pd <- -BirchCurve$b/(2*BirchCurve$a)
ElmCurve$pd <- -ElmCurve$b/(2*ElmCurve$a)
HazelCurve$pd <- -HazelCurve$b/(2*HazelCurve$a)
OakCurve$pd <- -OakCurve$b/(2*OakCurve$a)
RowanCurve$pd <- -RowanCurve$b/(2*RowanCurve$a)
SycamoreCurve$pd <- -SycamoreCurve$b/(2*SycamoreCurve$a)
WillowCurve$pd <- -WillowCurve$b/(2*WillowCurve$a)

AlderPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])))
AshPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])))
BeechPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])))
BirchPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])))
ElmPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])))
HazelPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])))
OakPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])))
RowanPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])))
SycamorePDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])))
WillowPDCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])))

meanPD <- -(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))
AlderCurve$pddiff <- (-AlderCurve$b/(2*AlderCurve$a))-meanPD
AshCurve$pddiff <- (-AshCurve$b/(2*AshCurve$a))-meanPD
BeechCurve$pddiff <- (-BeechCurve$b/(2*BeechCurve$a))-meanPD
BirchCurve$pddiff <- (-BirchCurve$b/(2*BirchCurve$a))-meanPD
ElmCurve$pddiff <- (-ElmCurve$b/(2*ElmCurve$a))-meanPD
HazelCurve$pddiff <- (-HazelCurve$b/(2*HazelCurve$a))-meanPD
OakCurve$pddiff <- (-OakCurve$b/(2*OakCurve$a))-meanPD
RowanCurve$pddiff <- (-RowanCurve$b/(2*RowanCurve$a))-meanPD
SycamoreCurve$pddiff <- (-SycamoreCurve$b/(2*SycamoreCurve$a))-meanPD
WillowCurve$pddiff <- (-WillowCurve$b/(2*WillowCurve$a))-meanPD

AlderPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
AshPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
BeechPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
BirchPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
ElmPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
HazelPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
OakPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
RowanPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
SycamorePDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))
WillowPDdifCI <- HPDinterval(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))-(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3]))))


TTcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                       PD=c(mean(AlderCurve$pd),mean(AshCurve$pd),mean(BeechCurve$pd),mean(BirchCurve$pd),mean(ElmCurve$pd),mean(HazelCurve$pd),mean(OakCurve$pd),mean(RowanCurve$pd),mean(SycamoreCurve$pd),mean(WillowCurve$pd)),
                       PDLCI=c(AlderPDCI[1], AshPDCI[1], BeechPDCI[1], BirchPDCI[1], ElmPDCI[1], HazelPDCI[1], OakPDCI[1], RowanPDCI[1], SycamorePDCI[1], WillowPDCI[1]),
                       PDUCI=c(AlderPDCI[2], AshPDCI[2], BeechPDCI[2], BirchPDCI[2], ElmPDCI[2], HazelPDCI[2], OakPDCI[2], RowanPDCI[2], SycamorePDCI[2], WillowPDCI[2]),
                       PDdiff=c(mean(AlderCurve$pddiff),mean(AshCurve$pddiff),mean(BeechCurve$pddiff),mean(BirchCurve$pddiff),mean(ElmCurve$pddiff),mean(HazelCurve$pddiff),mean(OakCurve$pddiff),mean(RowanCurve$pddiff),mean(SycamoreCurve$pddiff),mean(WillowCurve$pddiff)),
                       PDdiffLCI=c(AlderPDdifCI[1], AshPDdifCI[1], BeechPDdifCI[1], BirchPDdifCI[1], ElmPDdifCI[1], HazelPDdifCI[1], OakPDdifCI[1], RowanPDdifCI[1], SycamorePDdifCI[1], WillowPDdifCI[1]),
                       PDdiffUCI=c(AlderPDdifCI[2], AshPDdifCI[2], BeechPDdifCI[2], BirchPDdifCI[2], ElmPDdifCI[2], HazelPDdifCI[2], OakPDdifCI[2], RowanPDdifCI[2], SycamorePDdifCI[2], WillowPDdifCI[2]))

#####################
#### Peak Height ####
#####################

meanPH <- AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*meanPD+AbundTTCurves3$Sol[,3]*meanPD^2
AlderCurve$ph <- AlderCurve$a*AlderCurve$pd^2+AlderCurve$b*AlderCurve$pd+AlderCurve$c
AshCurve$ph <- AshCurve$a*AshCurve$pd^2+AshCurve$b*AshCurve$pd+AshCurve$c
BeechCurve$ph <- BeechCurve$a*BeechCurve$pd^2+BeechCurve$b*BeechCurve$pd+BeechCurve$c
BirchCurve$ph <- BirchCurve$a*BirchCurve$pd^2+BirchCurve$b*BirchCurve$pd+BirchCurve$c
ElmCurve$ph <- ElmCurve$a*ElmCurve$pd^2+ElmCurve$b*ElmCurve$pd+ElmCurve$c
HazelCurve$ph <- HazelCurve$a*HazelCurve$pd^2+HazelCurve$b*HazelCurve$pd+HazelCurve$c
OakCurve$ph <- OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c
RowanCurve$ph <- RowanCurve$a*RowanCurve$pd^2+RowanCurve$b*RowanCurve$pd+RowanCurve$c
SycamoreCurve$ph <- SycamoreCurve$a*SycamoreCurve$pd^2+SycamoreCurve$b*SycamoreCurve$pd+SycamoreCurve$c
WillowCurve$ph <- WillowCurve$a*WillowCurve$pd^2+WillowCurve$b*WillowCurve$pd+WillowCurve$c

AlderCurve$phdif <- exp(AlderCurve$a*AlderCurve$pd^2+AlderCurve$b*AlderCurve$pd+AlderCurve$c)-exp(meanPH)
AshCurve$phdif <- exp(AshCurve$a*AshCurve$pd^2+AshCurve$b*AshCurve$pd+AshCurve$c)-exp(meanPH)
BeechCurve$phdif <- exp(BeechCurve$a*BeechCurve$pd^2+BeechCurve$b*BeechCurve$pd+BeechCurve$c)-exp(meanPH)
BirchCurve$phdif <- exp(BirchCurve$a*BirchCurve$pd^2+BirchCurve$b*BirchCurve$pd+BirchCurve$c)-exp(meanPH)
ElmCurve$phdif <- exp(ElmCurve$a*ElmCurve$pd^2+ElmCurve$b*ElmCurve$pd+ElmCurve$c)-exp(meanPH)
HazelCurve$phdif <- exp(HazelCurve$a*HazelCurve$pd^2+HazelCurve$b*HazelCurve$pd+HazelCurve$c)-exp(meanPH)
OakCurve$phdif <- exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)-exp(meanPH)
RowanCurve$phdif <- exp(RowanCurve$a*RowanCurve$pd^2+RowanCurve$b*RowanCurve$pd+RowanCurve$c)-exp(meanPH)
SycamoreCurve$phdif <- exp(SycamoreCurve$a*SycamoreCurve$pd^2+SycamoreCurve$b*SycamoreCurve$pd+SycamoreCurve$c)-exp(meanPH)
WillowCurve$phdif <- exp(WillowCurve$a*WillowCurve$pd^2+WillowCurve$b*WillowCurve$pd+WillowCurve$c)-exp(meanPH)

AlderPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])))^2)
AshPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])))^2)
BeechPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])))^2)
BirchPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])))^2)
ElmPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])))^2)
HazelPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])))^2)
OakPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])))^2)
RowanPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])))^2)
SycamorePHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])))^2)
WillowPHCI <- HPDinterval(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])))^2)

AlderPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
AshPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
BeechPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
BirchPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
ElmPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
HazelPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
OakPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
RowanPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
SycamorePHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))
WillowPHdifCI <- HPDinterval(exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13]+(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])))+(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*(-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])))^2)-exp(AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,2]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))+AbundTTCurves3$Sol[,3]*(-(AbundTTCurves3$Sol[,2])/(2*(AbundTTCurves3$Sol[,3])))^2))

TTcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                       PD=c(mean(AlderCurve$pd),mean(AshCurve$pd),mean(BeechCurve$pd),mean(BirchCurve$pd),mean(ElmCurve$pd),mean(HazelCurve$pd),mean(OakCurve$pd),mean(RowanCurve$pd),mean(SycamoreCurve$pd),mean(WillowCurve$pd)),
                       PDLCI=c(AlderPDCI[1], AshPDCI[1], BeechPDCI[1], BirchPDCI[1], ElmPDCI[1], HazelPDCI[1], OakPDCI[1], RowanPDCI[1], SycamorePDCI[1], WillowPDCI[1]),
                       PDUCI=c(AlderPDCI[2], AshPDCI[2], BeechPDCI[2], BirchPDCI[2], ElmPDCI[2], HazelPDCI[2], OakPDCI[2], RowanPDCI[2], SycamorePDCI[2], WillowPDCI[2]),
                       PDdiff=c(mean(AlderCurve$pddiff),mean(AshCurve$pddiff),mean(BeechCurve$pddiff),mean(BirchCurve$pddiff),mean(ElmCurve$pddiff),mean(HazelCurve$pddiff),mean(OakCurve$pddiff),mean(RowanCurve$pddiff),mean(SycamoreCurve$pddiff),mean(WillowCurve$pddiff)),
                       PDdiffLCI=c(AlderPDdifCI[1], AshPDdifCI[1], BeechPDdifCI[1], BirchPDdifCI[1], ElmPDdifCI[1], HazelPDdifCI[1], OakPDdifCI[1], RowanPDdifCI[1], SycamorePDdifCI[1], WillowPDdifCI[1]),
                       PDdiffUCI=c(AlderPDdifCI[2], AshPDdifCI[2], BeechPDdifCI[2], BirchPDdifCI[2], ElmPDdifCI[2], HazelPDdifCI[2], OakPDdifCI[2], RowanPDdifCI[2], SycamorePDdifCI[2], WillowPDdifCI[2]),
                       PH=c(mean(AlderCurve$ph),mean(AshCurve$ph),mean(BeechCurve$ph),mean(BirchCurve$ph),mean(ElmCurve$ph),mean(HazelCurve$ph),mean(OakCurve$ph),mean(RowanCurve$ph),mean(SycamoreCurve$ph),mean(WillowCurve$ph)),
                       PHLCI=c(AlderPHCI[1], AshPHCI[1], BeechPHCI[1], BirchPHCI[1], ElmPHCI[1], HazelPHCI[1], OakPHCI[1], RowanPHCI[1], SycamorePHCI[1], WillowPHCI[1]),
                       PHUCI=c(AlderPHCI[2], AshPHCI[2], BeechPHCI[2], BirchPHCI[2], ElmPHCI[2], HazelPHCI[2], OakPHCI[2], RowanPHCI[2], SycamorePHCI[2], WillowPHCI[2]),
                       PHdiff=c(mean(AlderCurve$phdif),mean(AshCurve$phdif),mean(BeechCurve$phdif),mean(BirchCurve$phdif),mean(ElmCurve$phdif),mean(HazelCurve$phdif),mean(OakCurve$phdif),mean(RowanCurve$phdif),mean(SycamoreCurve$phdif),mean(WillowCurve$phdif)),
                       PHdiffLCI=c(AlderPHdifCI[1], AshPHdifCI[1], BeechPHdifCI[1], BirchPHdifCI[1], ElmPHdifCI[1], HazelPHdifCI[1], OakPHdifCI[1], RowanPHdifCI[1], SycamorePHdifCI[1], WillowPHdifCI[1]),
                       PHdiffUCI=c(AlderPHdifCI[2], AshPHdifCI[2], BeechPHdifCI[2], BirchPHdifCI[2], ElmPHdifCI[2], HazelPHdifCI[2], OakPHdifCI[2], RowanPHdifCI[2], SycamorePHdifCI[2], WillowPHdifCI[2]))


######################################
#### Peak Width at 0.01 abundance ####  x= (-b +/- sqrt(b^2 - 4ac))/2a
######################################

#Alder.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])
#Ash.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])
#Beech.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])
#Birch.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])
#Elm.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])
#Hazel.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])
#Oak.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])
#Rowan.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])
#Sycamore.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])
#Willow.a <- (AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])

#Alder.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])
#Ash.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])
#Beech.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])
#Birch.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])
#Elm.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])
#Hazel.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])
#Oak.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])
#Rowan.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])
#Sycamore.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])
#Willow.b <- (AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])

#Alder.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))
#Ash.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))
#Beech.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))
#Birch.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))
#Elm.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))
#Hazel.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))
#Oak.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))
#Rowan.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))
#Sycamore.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))
#Willow.c <- ((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))

AlderPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AlderPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AshPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
AshPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
BeechPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BeechPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BirchPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
BirchPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
ElmPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
ElmPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
HazelPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
HazelPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
OakPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
OakPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
RowanPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
RowanPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
SycamorePW1 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
SycamorePW2 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
WillowPW1   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))
WillowPW2   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))

MeanPW1    <- (-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
MeanPW2    <- (-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))

AlderPW1dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
AlderPW2dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
AshPW1dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
AshPW2dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
BeechPW1dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
BeechPW2dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
BirchPW1dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
BirchPW2dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
ElmPW1dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
ElmPW2dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
HazelPW1dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
HazelPW2dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
OakPW1dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
OakPW2dif      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
RowanPW1dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
RowanPW2dif    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
SycamorePW1dif <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
SycamorePW2dif <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
WillowPW1dif   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))-(-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))
WillowPW2dif   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))-(-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]))

#### NaNs in Ash and Mean !

#PW1 lower, PW2 higher

#For left and right dates
#TTcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
#                       PD=c(mean(AlderCurve$pd),mean(AshCurve$pd),mean(BeechCurve$pd),mean(BirchCurve$pd),mean(ElmCurve$pd),mean(HazelCurve$pd),mean(OakCurve$pd),mean(RowanCurve$pd),mean(SycamoreCurve$pd),mean(WillowCurve$pd)),
#                       PDLCI=c(AlderPDCI[1], AshPDCI[1], BeechPDCI[1], BirchPDCI[1], ElmPDCI[1], HazelPDCI[1], OakPDCI[1], RowanPDCI[1], SycamorePDCI[1], WillowPDCI[1]),
#                       PDUCI=c(AlderPDCI[2], AshPDCI[2], BeechPDCI[2], BirchPDCI[2], ElmPDCI[2], HazelPDCI[2], OakPDCI[2], RowanPDCI[2], SycamorePDCI[2], WillowPDCI[2]),
#                       PDdiff=c(mean(AlderCurve$pddiff),mean(AshCurve$pddiff),mean(BeechCurve$pddiff),mean(BirchCurve$pddiff),mean(ElmCurve$pddiff),mean(HazelCurve$pddiff),mean(OakCurve$pddiff),mean(RowanCurve$pddiff),mean(SycamoreCurve$pddiff),mean(WillowCurve$pddiff)),
#                       PDdiffLCI=c(AlderPDdifCI[1], AshPDdifCI[1], BeechPDdifCI[1], BirchPDdifCI[1], ElmPDdifCI[1], HazelPDdifCI[1], OakPDdifCI[1], RowanPDdifCI[1], SycamorePDdifCI[1], WillowPDdifCI[1]),
#                       PDdiffUCI=c(AlderPDdifCI[2], AshPDdifCI[2], BeechPDdifCI[2], BirchPDdifCI[2], ElmPDdifCI[2], HazelPDdifCI[2], OakPDdifCI[2], RowanPDdifCI[2], SycamorePDdifCI[2], WillowPDdifCI[2]),
#                       PH=c(mean(AlderCurve$ph),mean(AshCurve$ph),mean(BeechCurve$ph),mean(BirchCurve$ph),mean(ElmCurve$ph),mean(HazelCurve$ph),mean(OakCurve$ph),mean(RowanCurve$ph),mean(SycamoreCurve$ph),mean(WillowCurve$ph)),
#                       PHLCI=c(AlderPHCI[1], AshPHCI[1], BeechPHCI[1], BirchPHCI[1], ElmPHCI[1], HazelPHCI[1], OakPHCI[1], RowanPHCI[1], SycamorePHCI[1], WillowPHCI[1]),
#                       PHUCI=c(AlderPHCI[2], AshPHCI[2], BeechPHCI[2], BirchPHCI[2], ElmPHCI[2], HazelPHCI[2], OakPHCI[2], RowanPHCI[2], SycamorePHCI[2], WillowPHCI[2]),
#                       PHdiffLCI=c(AlderPHdifCI[1], AshPHdifCI[1], BeechPHdifCI[1], BirchPHdifCI[1], ElmPHdifCI[1], HazelPHdifCI[1], OakPHdifCI[1], RowanPHdifCI[1], SycamorePHdifCI[1], WillowPHdifCI[1]),
#                       PHdiff=c(mean(AlderCurve$phdif),mean(AshCurve$phdif),mean(BeechCurve$phdif),mean(BirchCurve$phdif),mean(ElmCurve$phdif),mean(HazelCurve$phdif),mean(OakCurve$phdif),mean(RowanCurve$phdif),mean(SycamoreCurve$phdif),mean(WillowCurve$phdif)),
#                       PWL=c(mean(AlderPW1), mean(AshPW1,na.rm=TRUE), mean(BeechPW1), mean(BirchPW1), mean(ElmPW1), mean(HazelPW1), mean(OakPW1), mean(RowanPW1), mean(SycamorePW1), mean(WillowPW1)),
#                       PHdiffUCI=c(AlderPHdifCI[2], AshPHdifCI[2], BeechPHdifCI[2], BirchPHdifCI[2], ElmPHdifCI[2], HazelPHdifCI[2], OakPHdifCI[2], RowanPHdifCI[2], SycamorePHdifCI[2], WillowPHdifCI[2]),
#                       PWR=c(mean(AlderPW2), mean(AshPW2,na.rm=TRUE), mean(BeechPW2), mean(BirchPW2), mean(ElmPW2), mean(HazelPW2), mean(OakPW2), mean(RowanPW2), mean(SycamorePW2), mean(WillowPW2)),
#                       PWLLCI=c(HPDinterval(AlderPW1)[1], HPDinterval(AshPW1)[1], HPDinterval(BeechPW1)[1], HPDinterval(BirchPW1)[1], HPDinterval(ElmPW1)[1], HPDinterval(HazelPW1)[1], HPDinterval(OakPW1)[1], HPDinterval(RowanPW1)[1], HPDinterval(SycamorePW1)[1], HPDinterval(WillowPW1)[1]),
#                       PWRLCI=c(HPDinterval(AlderPW2)[1], HPDinterval(AshPW2)[1], HPDinterval(BeechPW2)[1], HPDinterval(BirchPW2)[1], HPDinterval(ElmPW2)[1], HPDinterval(HazelPW2)[1], HPDinterval(OakPW2)[1], HPDinterval(RowanPW2)[1], HPDinterval(SycamorePW2)[1], HPDinterval(WillowPW2)[1]),
#                       PWLUCI=c(HPDinterval(AlderPW1)[2], HPDinterval(AshPW1)[2], HPDinterval(BeechPW1)[2], HPDinterval(BirchPW1)[2], HPDinterval(ElmPW1)[2], HPDinterval(HazelPW1)[2], HPDinterval(OakPW1)[2], HPDinterval(RowanPW1)[2], HPDinterval(SycamorePW1)[2], HPDinterval(WillowPW1)[2]),
#                       PWRUCI=c(HPDinterval(AlderPW2)[2], HPDinterval(AshPW2)[2], HPDinterval(BeechPW2)[2], HPDinterval(BirchPW2)[2], HPDinterval(ElmPW2)[2], HPDinterval(HazelPW2)[2], HPDinterval(OakPW2)[2], HPDinterval(RowanPW2)[2], HPDinterval(SycamorePW2)[2], HPDinterval(WillowPW2)[2]))

# Width calc
TTcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                       PD=c(mean(AlderCurve$pd),mean(AshCurve$pd),mean(BeechCurve$pd),mean(BirchCurve$pd),mean(ElmCurve$pd),mean(HazelCurve$pd),mean(OakCurve$pd),mean(RowanCurve$pd),mean(SycamoreCurve$pd),mean(WillowCurve$pd)),
                       PDLCI=c(AlderPDCI[1], AshPDCI[1], BeechPDCI[1], BirchPDCI[1], ElmPDCI[1], HazelPDCI[1], OakPDCI[1], RowanPDCI[1], SycamorePDCI[1], WillowPDCI[1]),
                       PDUCI=c(AlderPDCI[2], AshPDCI[2], BeechPDCI[2], BirchPDCI[2], ElmPDCI[2], HazelPDCI[2], OakPDCI[2], RowanPDCI[2], SycamorePDCI[2], WillowPDCI[2]),
                       PDdiff=c(mean(AlderCurve$pddiff),mean(AshCurve$pddiff),mean(BeechCurve$pddiff),mean(BirchCurve$pddiff),mean(ElmCurve$pddiff),mean(HazelCurve$pddiff),mean(OakCurve$pddiff),mean(RowanCurve$pddiff),mean(SycamoreCurve$pddiff),mean(WillowCurve$pddiff)),
                       PDdiffLCI=c(AlderPDdifCI[1], AshPDdifCI[1], BeechPDdifCI[1], BirchPDdifCI[1], ElmPDdifCI[1], HazelPDdifCI[1], OakPDdifCI[1], RowanPDdifCI[1], SycamorePDdifCI[1], WillowPDdifCI[1]),
                       PDdiffUCI=c(AlderPDdifCI[2], AshPDdifCI[2], BeechPDdifCI[2], BirchPDdifCI[2], ElmPDdifCI[2], HazelPDdifCI[2], OakPDdifCI[2], RowanPDdifCI[2], SycamorePDdifCI[2], WillowPDdifCI[2]),
                       PH=c(mean(AlderCurve$ph),mean(AshCurve$ph),mean(BeechCurve$ph),mean(BirchCurve$ph),mean(ElmCurve$ph),mean(HazelCurve$ph),mean(OakCurve$ph),mean(RowanCurve$ph),mean(SycamoreCurve$ph),mean(WillowCurve$ph)),
                       PHLCI=c(AlderPHCI[1], AshPHCI[1], BeechPHCI[1], BirchPHCI[1], ElmPHCI[1], HazelPHCI[1], OakPHCI[1], RowanPHCI[1], SycamorePHCI[1], WillowPHCI[1]),
                       PHUCI=c(AlderPHCI[2], AshPHCI[2], BeechPHCI[2], BirchPHCI[2], ElmPHCI[2], HazelPHCI[2], OakPHCI[2], RowanPHCI[2], SycamorePHCI[2], WillowPHCI[2]),
                       PHdiff=c(mean(AlderCurve$phdif),mean(AshCurve$phdif),mean(BeechCurve$phdif),mean(BirchCurve$phdif),mean(ElmCurve$phdif),mean(HazelCurve$phdif),mean(OakCurve$phdif),mean(RowanCurve$phdif),mean(SycamoreCurve$phdif),mean(WillowCurve$phdif)),
                       PHdiffLCI=c(AlderPHdifCI[1], AshPHdifCI[1], BeechPHdifCI[1], BirchPHdifCI[1], ElmPHdifCI[1], HazelPHdifCI[1], OakPHdifCI[1], RowanPHdifCI[1], SycamorePHdifCI[1], WillowPHdifCI[1]),
                       PHdiffUCI=c(AlderPHdifCI[2], AshPHdifCI[2], BeechPHdifCI[2], BirchPHdifCI[2], ElmPHdifCI[2], HazelPHdifCI[2], OakPHdifCI[2], RowanPHdifCI[2], SycamorePHdifCI[2], WillowPHdifCI[2]),
                       PW=c(mean(AlderPW2-AlderPW1), mean(AshPW2-AshPW1,na.rm=TRUE), mean(BeechPW2-BeechPW1), mean(BirchPW2-BirchPW1), mean(ElmPW2-ElmPW1), mean(HazelPW2-HazelPW1), mean(OakPW2-OakPW1), mean(RowanPW2-RowanPW1), mean(SycamorePW2-SycamorePW1), mean(WillowPW2-WillowPW1)),
                       PWLCI=c(HPDinterval(AlderPW2-AlderPW1)[1], HPDinterval(AshPW2-AshPW1)[1], HPDinterval(BeechPW2-BeechPW1)[1], HPDinterval(BirchPW2-BirchPW1)[1], HPDinterval(ElmPW2-ElmPW1)[1], HPDinterval(HazelPW2-HazelPW1)[1], HPDinterval(OakPW2-OakPW1)[1], HPDinterval(RowanPW2-RowanPW1)[1], HPDinterval(SycamorePW2-SycamorePW1)[1], HPDinterval(WillowPW2-WillowPW1)[1]),
                       PWUCI=c(HPDinterval(AlderPW2-AlderPW1)[2], HPDinterval(AshPW2-AshPW1)[2], HPDinterval(BeechPW2-BeechPW1)[2], HPDinterval(BirchPW2-BirchPW1)[2], HPDinterval(ElmPW2-ElmPW1)[2], HPDinterval(HazelPW2-HazelPW1)[2], HPDinterval(OakPW2-OakPW1)[2], HPDinterval(RowanPW2-RowanPW1)[2], HPDinterval(SycamorePW2-SycamorePW1)[2], HPDinterval(WillowPW2-WillowPW1)[2]),
                       PWdiff=c(mean(AlderPW2dif-AlderPW1dif,na.rm=TRUE), mean(AshPW2dif-AshPW1dif,na.rm=TRUE), mean(BeechPW2dif-BeechPW1dif,na.rm=TRUE), mean(BirchPW2dif-BirchPW1dif,na.rm=TRUE), mean(ElmPW2dif-ElmPW1dif,na.rm=TRUE), mean(HazelPW2dif-HazelPW1dif,na.rm=TRUE), mean(OakPW2dif-OakPW1dif,na.rm=TRUE), mean(RowanPW2dif-RowanPW1dif,na.rm=TRUE), mean(SycamorePW2dif-SycamorePW1dif,na.rm=TRUE), mean(WillowPW2dif-WillowPW1dif,na.rm=TRUE)),
                       PWdiffLCI=c(HPDinterval(AlderPW2dif-AlderPW1dif)[1], HPDinterval(AshPW2dif-AshPW1dif)[1], HPDinterval(BeechPW2dif-BeechPW1dif)[1], HPDinterval(BirchPW2dif-BirchPW1dif)[1], HPDinterval(ElmPW2dif-ElmPW1dif)[1], HPDinterval(HazelPW2dif-HazelPW1dif)[1], HPDinterval(OakPW2dif-OakPW1dif)[1], HPDinterval(RowanPW2dif-RowanPW1dif)[1], HPDinterval(SycamorePW2dif-SycamorePW1dif)[1], HPDinterval(WillowPW2dif-WillowPW1dif)[1]),
                       PWdiffUCI=c(HPDinterval(AlderPW2dif-AlderPW1dif)[2], HPDinterval(AshPW2dif-AshPW1dif)[2], HPDinterval(BeechPW2dif-BeechPW1dif)[2], HPDinterval(BirchPW2dif-BirchPW1dif)[2], HPDinterval(ElmPW2dif-ElmPW1dif)[2], HPDinterval(HazelPW2dif-HazelPW1dif)[2], HPDinterval(OakPW2dif-OakPW1dif)[2], HPDinterval(RowanPW2dif-RowanPW1dif)[2], HPDinterval(SycamorePW2dif-SycamorePW1dif)[2], HPDinterval(WillowPW2dif-WillowPW1dif)[2]),
                       PWdiff2=c(mean((AlderPW2-AlderPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((AshPW2-AshPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((BeechPW2-BeechPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((BirchPW2-BirchPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((ElmPW2-ElmPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((HazelPW2-HazelPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((OakPW2-OakPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((RowanPW2-RowanPW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((SycamorePW2-SycamorePW1)-(MeanPW2-MeanPW1),na.rm=TRUE), mean((WillowPW2-WillowPW1)-(MeanPW2-MeanPW1),na.rm=TRUE)))


#### Peak Width at 0.5 height ####

AlderPW1.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(exp(AlderCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AlderPW2.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(exp(AlderCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AshPW1.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(exp(AshCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
AshPW2.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(exp(AshCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
BeechPW1.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(exp(BeechCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BeechPW2.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(exp(BeechCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BirchPW1.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(exp(BirchCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
BirchPW2.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(exp(BirchCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
ElmPW1.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(exp(ElmCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
ElmPW2.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(exp(ElmCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
HazelPW1.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(exp(HazelCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
HazelPW2.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(exp(HazelCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
OakPW1.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(exp(OakCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
OakPW2.5      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(exp(OakCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
RowanPW1.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(exp(RowanCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
RowanPW2.5    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(exp(RowanCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
SycamorePW1.5 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(exp(SycamoreCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
SycamorePW2.5 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(exp(SycamoreCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
WillowPW1.5   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(exp(WillowCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))
WillowPW2.5   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(exp(WillowCurve$ph)/2))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))

MeanPW1.5    <- (-(AbundTTCurves3$Sol[,2])+sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(exp(meanPH)/2))))/(2*(AbundTTCurves3$Sol[,3]))
MeanPW2.5    <- (-(AbundTTCurves3$Sol[,2])-sqrt((AbundTTCurves3$Sol[,2])^2-4*(AbundTTCurves3$Sol[,3])*((AbundTTCurves3$Sol[,1])-log(exp(meanPH)/2))))/(2*(AbundTTCurves3$Sol[,3]))

TTcurves2 <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                       PD=c(mean(AlderCurve$pd),mean(AshCurve$pd),mean(BeechCurve$pd),mean(BirchCurve$pd),mean(ElmCurve$pd),mean(HazelCurve$pd),mean(OakCurve$pd),mean(RowanCurve$pd),mean(SycamoreCurve$pd),mean(WillowCurve$pd)),
                       PDLCI=c(AlderPDCI[1], AshPDCI[1], BeechPDCI[1], BirchPDCI[1], ElmPDCI[1], HazelPDCI[1], OakPDCI[1], RowanPDCI[1], SycamorePDCI[1], WillowPDCI[1]),
                       PDUCI=c(AlderPDCI[2], AshPDCI[2], BeechPDCI[2], BirchPDCI[2], ElmPDCI[2], HazelPDCI[2], OakPDCI[2], RowanPDCI[2], SycamorePDCI[2], WillowPDCI[2]),
                       PDdiff=c(mean(AlderCurve$pddiff),mean(AshCurve$pddiff),mean(BeechCurve$pddiff),mean(BirchCurve$pddiff),mean(ElmCurve$pddiff),mean(HazelCurve$pddiff),mean(OakCurve$pddiff),mean(RowanCurve$pddiff),mean(SycamoreCurve$pddiff),mean(WillowCurve$pddiff)),
                       PDdiffLCI=c(AlderPDdifCI[1], AshPDdifCI[1], BeechPDdifCI[1], BirchPDdifCI[1], ElmPDdifCI[1], HazelPDdifCI[1], OakPDdifCI[1], RowanPDdifCI[1], SycamorePDdifCI[1], WillowPDdifCI[1]),
                       PDdiffUCI=c(AlderPDdifCI[2], AshPDdifCI[2], BeechPDdifCI[2], BirchPDdifCI[2], ElmPDdifCI[2], HazelPDdifCI[2], OakPDdifCI[2], RowanPDdifCI[2], SycamorePDdifCI[2], WillowPDdifCI[2]),
                       PH=c(mean(AlderCurve$ph),mean(AshCurve$ph),mean(BeechCurve$ph),mean(BirchCurve$ph),mean(ElmCurve$ph),mean(HazelCurve$ph),mean(OakCurve$ph),mean(RowanCurve$ph),mean(SycamoreCurve$ph),mean(WillowCurve$ph)),
                       PHLCI=c(AlderPHCI[1], AshPHCI[1], BeechPHCI[1], BirchPHCI[1], ElmPHCI[1], HazelPHCI[1], OakPHCI[1], RowanPHCI[1], SycamorePHCI[1], WillowPHCI[1]),
                       PHUCI=c(AlderPHCI[2], AshPHCI[2], BeechPHCI[2], BirchPHCI[2], ElmPHCI[2], HazelPHCI[2], OakPHCI[2], RowanPHCI[2], SycamorePHCI[2], WillowPHCI[2]),
                       PHdiff=c(mean(AlderCurve$phdif),mean(AshCurve$phdif),mean(BeechCurve$phdif),mean(BirchCurve$phdif),mean(ElmCurve$phdif),mean(HazelCurve$phdif),mean(OakCurve$phdif),mean(RowanCurve$phdif),mean(SycamoreCurve$phdif),mean(WillowCurve$phdif)),
                       PHdiffLCI=c(AlderPHdifCI[1], AshPHdifCI[1], BeechPHdifCI[1], BirchPHdifCI[1], ElmPHdifCI[1], HazelPHdifCI[1], OakPHdifCI[1], RowanPHdifCI[1], SycamorePHdifCI[1], WillowPHdifCI[1]),
                       PHdiffUCI=c(AlderPHdifCI[2], AshPHdifCI[2], BeechPHdifCI[2], BirchPHdifCI[2], ElmPHdifCI[2], HazelPHdifCI[2], OakPHdifCI[2], RowanPHdifCI[2], SycamorePHdifCI[2], WillowPHdifCI[2]),
                       PW=c(mean(AlderPW2-AlderPW1), mean(AshPW2-AshPW1,na.rm=TRUE), mean(BeechPW2-BeechPW1), mean(BirchPW2-BirchPW1), mean(ElmPW2-ElmPW1), mean(HazelPW2-HazelPW1), mean(OakPW2-OakPW1), mean(RowanPW2-RowanPW1), mean(SycamorePW2-SycamorePW1), mean(WillowPW2-WillowPW1)),
                       PWLCI=c(HPDinterval(AlderPW2-AlderPW1)[1], HPDinterval(AshPW2-AshPW1)[1], HPDinterval(BeechPW2-BeechPW1)[1], HPDinterval(BirchPW2-BirchPW1)[1], HPDinterval(ElmPW2-ElmPW1)[1], HPDinterval(HazelPW2-HazelPW1)[1], HPDinterval(OakPW2-OakPW1)[1], HPDinterval(RowanPW2-RowanPW1)[1], HPDinterval(SycamorePW2-SycamorePW1)[1], HPDinterval(WillowPW2-WillowPW1)[1]),
                       PWUCI=c(HPDinterval(AlderPW2-AlderPW1)[2], HPDinterval(AshPW2-AshPW1)[2], HPDinterval(BeechPW2-BeechPW1)[2], HPDinterval(BirchPW2-BirchPW1)[2], HPDinterval(ElmPW2-ElmPW1)[2], HPDinterval(HazelPW2-HazelPW1)[2], HPDinterval(OakPW2-OakPW1)[2], HPDinterval(RowanPW2-RowanPW1)[2], HPDinterval(SycamorePW2-SycamorePW1)[2], HPDinterval(WillowPW2-WillowPW1)[2]),
                       PWdiff=c(mean(AlderPW2dif-AlderPW1dif,na.rm=TRUE), mean(AshPW2dif-AshPW1dif,na.rm=TRUE), mean(BeechPW2dif-BeechPW1dif,na.rm=TRUE), mean(BirchPW2dif-BirchPW1dif,na.rm=TRUE), mean(ElmPW2dif-ElmPW1dif,na.rm=TRUE), mean(HazelPW2dif-HazelPW1dif,na.rm=TRUE), mean(OakPW2dif-OakPW1dif,na.rm=TRUE), mean(RowanPW2dif-RowanPW1dif,na.rm=TRUE), mean(SycamorePW2dif-SycamorePW1dif,na.rm=TRUE), mean(WillowPW2dif-WillowPW1dif,na.rm=TRUE)),
                       PWdiffLCI=c(HPDinterval(AlderPW2dif-AlderPW1dif)[1], HPDinterval(AshPW2dif-AshPW1dif)[1], HPDinterval(BeechPW2dif-BeechPW1dif)[1], HPDinterval(BirchPW2dif-BirchPW1dif)[1], HPDinterval(ElmPW2dif-ElmPW1dif)[1], HPDinterval(HazelPW2dif-HazelPW1dif)[1], HPDinterval(OakPW2dif-OakPW1dif)[1], HPDinterval(RowanPW2dif-RowanPW1dif)[1], HPDinterval(SycamorePW2dif-SycamorePW1dif)[1], HPDinterval(WillowPW2dif-WillowPW1dif)[1]),
                       PWdiffUCI=c(HPDinterval(AlderPW2dif-AlderPW1dif)[2], HPDinterval(AshPW2dif-AshPW1dif)[2], HPDinterval(BeechPW2dif-BeechPW1dif)[2], HPDinterval(BirchPW2dif-BirchPW1dif)[2], HPDinterval(ElmPW2dif-ElmPW1dif)[2], HPDinterval(HazelPW2dif-HazelPW1dif)[2], HPDinterval(OakPW2dif-OakPW1dif)[2], HPDinterval(RowanPW2dif-RowanPW1dif)[2], HPDinterval(SycamorePW2dif-SycamorePW1dif)[2], HPDinterval(WillowPW2dif-WillowPW1dif)[2]),
                       PW.5=c(mean(AlderPW2.5-AlderPW1.5), mean(AshPW2.5-AshPW1.5,na.rm=TRUE), mean(BeechPW2.5-BeechPW1.5), mean(BirchPW2.5-BirchPW1.5), mean(ElmPW2.5-ElmPW1.5), mean(HazelPW2.5-HazelPW1.5), mean(OakPW2.5-OakPW1.5), mean(RowanPW2.5-RowanPW1.5), mean(SycamorePW2.5-SycamorePW1.5), mean(WillowPW2.5-WillowPW1.5)),
                       PWLCI.5=c(HPDinterval(AlderPW2.5-AlderPW1.5)[1], HPDinterval(AshPW2.5-AshPW1.5)[1], HPDinterval(BeechPW2.5-BeechPW1.5)[1], HPDinterval(BirchPW2.5-BirchPW1.5)[1], HPDinterval(ElmPW2.5-ElmPW1.5)[1], HPDinterval(HazelPW2.5-HazelPW1.5)[1], HPDinterval(OakPW2.5-OakPW1.5)[1], HPDinterval(RowanPW2.5-RowanPW1.5)[1], HPDinterval(SycamorePW2.5-SycamorePW1.5)[1], HPDinterval(WillowPW2.5-WillowPW1.5)[1]),
                       PWUCI.5=c(HPDinterval(AlderPW2.5-AlderPW1.5)[2], HPDinterval(AshPW2.5-AshPW1.5)[2], HPDinterval(BeechPW2.5-BeechPW1.5)[2], HPDinterval(BirchPW2.5-BirchPW1.5)[2], HPDinterval(ElmPW2.5-ElmPW1.5)[2], HPDinterval(HazelPW2.5-HazelPW1.5)[2], HPDinterval(OakPW2.5-OakPW1.5)[2], HPDinterval(RowanPW2.5-RowanPW1.5)[2], HPDinterval(SycamorePW2.5-SycamorePW1.5)[2], HPDinterval(WillowPW2.5-WillowPW1.5)[2]),
                       PWdiff.5=c(mean((AlderPW2.5-AlderPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((AshPW2.5-AshPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((BeechPW2.5-BeechPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((BirchPW2.5-BirchPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((ElmPW2.5-ElmPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((HazelPW2.5-HazelPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((OakPW2.5-OakPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((RowanPW2.5-RowanPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((SycamorePW2.5-SycamorePW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE), mean((WillowPW2.5-WillowPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)),
                       PWdiffLCI.5=c(HPDinterval((AlderPW2.5-AlderPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((AshPW2.5-AshPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((BeechPW2.5-BeechPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((BirchPW2.5-BirchPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((ElmPW2.5-ElmPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((HazelPW2.5-HazelPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((OakPW2.5-OakPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((RowanPW2.5-RowanPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((SycamorePW2.5-SycamorePW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1], HPDinterval((WillowPW2.5-WillowPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[1]),
                       PWdiffUCI.5=c(HPDinterval((AlderPW2.5-AlderPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((AshPW2.5-AshPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((BeechPW2.5-BeechPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((BirchPW2.5-BirchPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((ElmPW2.5-ElmPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((HazelPW2.5-HazelPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((OakPW2.5-OakPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((RowanPW2.5-RowanPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((SycamorePW2.5-SycamorePW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2], HPDinterval((WillowPW2.5-WillowPW1.5)-(MeanPW2.5-MeanPW1.5),na.rm=TRUE)[2]))

##########################
#### Curves in ggplot ####
##########################
dayscal <- seq(0.67,1,0.001)
curve <- mean(AbundTTCurves3$Sol[,1])+mean(AbundTTCurves3$Sol[,2])*dayscal+mean(AbundTTCurves3$Sol[,3])*dayscal^2
days <- dayscal*max(cater_habitat$date)
TTlines <- data.frame(dayscal=seq(0.67,1,0.001))
TTlines$days <- TTlines$dayscal*max(cater_habitat$date)
TTlines$Alder <- mean(AlderCurve[,1]+AlderCurve[,2]) +mean(AlderCurve[,3]+AlderCurve[,4])*TTlines$dayscal+mean(AlderCurve[,5]+AlderCurve[,6])*TTlines$dayscal^2
TTlines$Ash <- mean(AshCurve[,1]+AshCurve[,2]) +mean(AshCurve[,3]+AshCurve[,4])*TTlines$dayscal+mean(AshCurve[,5]+AshCurve[,6])*TTlines$dayscal^2
TTlines$Beech <- mean(BeechCurve[,1]+BeechCurve[,2]) +mean(BeechCurve[,3]+BeechCurve[,4])*TTlines$dayscal+mean(BeechCurve[,5]+BeechCurve[,6])*TTlines$dayscal^2
TTlines$Birch <- mean(BirchCurve[,1]+BirchCurve[,2]) +mean(BirchCurve[,3]+BirchCurve[,4])*TTlines$dayscal+mean(BirchCurve[,5]+BirchCurve[,6])*TTlines$dayscal^2
TTlines$Elm <- mean(ElmCurve[,1]+ElmCurve[,2]) +mean(ElmCurve[,3]+ElmCurve[,4])*TTlines$dayscal+mean(ElmCurve[,5]+ElmCurve[,6])*TTlines$dayscal^2
TTlines$Hazel <- mean(HazelCurve[,1]+HazelCurve[,2]) +mean(HazelCurve[,3]+HazelCurve[,4])*TTlines$dayscal+mean(HazelCurve[,5]+HazelCurve[,6])*TTlines$dayscal^2
TTlines$Oak <- mean(OakCurve[,1]+OakCurve[,2]) +mean(OakCurve[,3]+OakCurve[,4])*TTlines$dayscal+mean(OakCurve[,5]+OakCurve[,6])*TTlines$dayscal^2
TTlines$Rowan <- mean(RowanCurve[,1]+RowanCurve[,2]) +mean(RowanCurve[,3]+RowanCurve[,4])*TTlines$dayscal+mean(RowanCurve[,5]+RowanCurve[,6])*TTlines$dayscal^2
TTlines$Sycamore <- mean(SycamoreCurve[,1]+SycamoreCurve[,2]) +mean(SycamoreCurve[,3]+SycamoreCurve[,4])*TTlines$dayscal+mean(SycamoreCurve[,5]+SycamoreCurve[,6])*TTlines$dayscal^2
TTlines$Willow <- mean(WillowCurve[,1]+WillowCurve[,2]) +mean(WillowCurve[,3]+WillowCurve[,4])*TTlines$dayscal+mean(WillowCurve[,5]+WillowCurve[,6])*TTlines$dayscal^2

TTlineslong <- gather(TTlines, key="TreeTaxa", value="logAbund", select= 3:12)
TTlineslong$Abund <- exp(TTlineslong$logAbund)


#################
#### Figures ####
#################

# Colour palette
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")
NoOakCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "royalblue4", "slateblue2", "orchid")


# Tree taxa curves plotted
Curvesplot <- ggplot(TTlineslong, aes(days, Abund, col=TreeTaxa))+
  geom_line(size=0.7)+
  theme_bw()+
  xlab("Date")+
  ylab("Abundance")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Date
PDplot <- ggplot(TTcurves, aes(TT, (PD*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PDUCI*175), ymin=(PDLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak Date")+
  coord_flip()+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Date difference from mean
PDdifplot <- ggplot(TTcurves, aes(TT, (PDdiff*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PDdiffUCI*175), ymin=(PDdiffLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak date difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Height
PHplot <- ggplot(TTcurves, aes(TT, exp(PH), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=exp(PHUCI), ymin=exp(PHLCI), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak Height")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  coord_flip()+
  scale_colour_manual(values=AllTaxaCols)

# Peak Height difference form mean
PHdifplot <- ggplot(TTcurves, aes(TT, PHdiff, colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=PHdiffUCI, ymin=PHdiffLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height difference")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  coord_flip()+
  scale_colour_manual(values=AllTaxaCols)

# Peak date and Peak height
PHPDplot <- ggplot(TTcurves, aes((PD*175), exp(PH), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=exp(PHUCI), ymin=exp(PHLCI), width=0.5))+
  geom_errorbarh(aes(xmax=(PDUCI*175), xmin=(PDLCI*175), width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  xlab("Peak Date")+
  ylab("Peak Height")+
  theme(text = element_text(size=20))+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width 
PWplot <- ggplot(TTcurves, aes(TT, (PW*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWUCI*175), ymin=(PWLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak Width")+
  coord_flip()+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width difference from mean
PWdifplot <- ggplot(TTcurves, aes(TT, (PWdiff*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWdiffUCI*175), ymin=(PWdiffLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak width difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width at half height
PW.5plot <- ggplot(TTcurves2, aes(TT, (PW.5*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWUCI.5*175), ymin=(PWLCI.5*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak Width (0.5 height)")+
  coord_flip()+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width at half height difference from mean 
PW.5difplot <- ggplot(TTcurves2, aes(TT, (PWdiff.5*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWdiffUCI.5*175), ymin=(PWdiffLCI.5*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak width difference (0.5 height)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

library(gridExtra)

col1 <- grid.arrange(Curvesplot, nrow = 1, widths = 1)
col2 <- grid.arrange(PDplot, PDdifplot, nrow = 2, heights = c(1,1))
col3 <- grid.arrange(PHplot, PHdifplot, nrow = 2, heights = c(1,1))
col4 <- grid.arrange(PWplot, PWdifplot, nrow = 2, heights = c(1,1))
col5 <- grid.arrange(PW.5plot, PW.5difplot, nrow = 2, heights = c(1,1))
AbundCurvesFig <- grid.arrange(col1, col2, col3, col4, col5, ncol = 5, widths = c(2,1,1,1,1))


#############################
#### Difference from oak ####
#############################

#### Peak Date - Oak ####

AlderCurve$a <- AlderCurve[,5]+AlderCurve[,6]
AshCurve$a <- AshCurve[,5]+AshCurve[,6]
BeechCurve$a <- BeechCurve[,5]+BeechCurve[,6]
BirchCurve$a <- BirchCurve[,5]+BirchCurve[,6]
ElmCurve$a <- ElmCurve[,5]+ElmCurve[,6]
HazelCurve$a <- HazelCurve[,5]+HazelCurve[,6]
OakCurve$a <- OakCurve[,5]+OakCurve[,6]
RowanCurve$a <- RowanCurve[,5]+RowanCurve[,6]
SycamoreCurve$a <- SycamoreCurve[,5]+SycamoreCurve[,6]
WillowCurve$a <- WillowCurve[,5]+WillowCurve[,6]
AlderCurve$b <- AlderCurve[,3]+AlderCurve[,4]
AshCurve$b <- AshCurve[,3]+AshCurve[,4]
BeechCurve$b <- BeechCurve[,3]+BeechCurve[,4]
BirchCurve$b <- BirchCurve[,3]+BirchCurve[,4]
ElmCurve$b <- ElmCurve[,3]+ElmCurve[,4]
HazelCurve$b <- HazelCurve[,3]+HazelCurve[,4]
OakCurve$b <- OakCurve[,3]+OakCurve[,4]
RowanCurve$b <- RowanCurve[,3]+RowanCurve[,4]
SycamoreCurve$b <- SycamoreCurve[,3]+SycamoreCurve[,4]
WillowCurve$b <- WillowCurve[,3]+WillowCurve[,4]
AlderCurve$c <- AlderCurve[,1]+AlderCurve[,2]
AshCurve$c <- AshCurve[,1]+AshCurve[,2]
BeechCurve$c <- BeechCurve[,1]+BeechCurve[,2]
BirchCurve$c <- BirchCurve[,1]+BirchCurve[,2]
ElmCurve$c <- ElmCurve[,1]+ElmCurve[,2]
HazelCurve$c <- HazelCurve[,1]+HazelCurve[,2]
OakCurve$c <- OakCurve[,1]+OakCurve[,2]
RowanCurve$c <- RowanCurve[,1]+RowanCurve[,2]
SycamoreCurve$c <- SycamoreCurve[,1]+SycamoreCurve[,2]
WillowCurve$c <- WillowCurve[,1]+WillowCurve[,2]

AlderpdOD <- (-AlderCurve$b/(2*AlderCurve$a))-(-OakCurve$b/(2*OakCurve$a))
AshpdOD <- (-AshCurve$b/(2*AshCurve$a))-(-OakCurve$b/(2*OakCurve$a))
BeechpdOD <- (-BeechCurve$b/(2*BeechCurve$a))-(-OakCurve$b/(2*OakCurve$a))
BirchpdOD <- (-BirchCurve$b/(2*BirchCurve$a))-(-OakCurve$b/(2*OakCurve$a))
ElmpdOD <- (-ElmCurve$b/(2*ElmCurve$a))-(-OakCurve$b/(2*OakCurve$a))
HazelpdOD <- (-HazelCurve$b/(2*HazelCurve$a))-(-OakCurve$b/(2*OakCurve$a))
RowanpdOD <- (-RowanCurve$b/(2*RowanCurve$a))-(-OakCurve$b/(2*OakCurve$a))
SycamorepdOD <- (-SycamoreCurve$b/(2*SycamoreCurve$a))-(-OakCurve$b/(2*OakCurve$a))
WillowpdOD <- (-WillowCurve$b/(2*WillowCurve$a))-(-OakCurve$b/(2*OakCurve$a))

AlderphOD <- exp(AlderCurve$a*AlderCurve$pd^2+AlderCurve$b*AlderCurve$pd+AlderCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
AshphOD <- exp(AshCurve$a*AshCurve$pd^2+AshCurve$b*AshCurve$pd+AshCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
BeechphOD <- exp(BeechCurve$a*BeechCurve$pd^2+BeechCurve$b*BeechCurve$pd+BeechCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
BirchphOD <- exp(BirchCurve$a*BirchCurve$pd^2+BirchCurve$b*BirchCurve$pd+BirchCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
ElmphOD <- exp(ElmCurve$a*ElmCurve$pd^2+ElmCurve$b*ElmCurve$pd+ElmCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
HazelphOD <- exp(HazelCurve$a*HazelCurve$pd^2+HazelCurve$b*HazelCurve$pd+HazelCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
RowanphOD <- exp(RowanCurve$a*RowanCurve$pd^2+RowanCurve$b*RowanCurve$pd+RowanCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
SycamorephOD <- exp(SycamoreCurve$a*SycamoreCurve$pd^2+SycamoreCurve$b*SycamoreCurve$pd+SycamoreCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)
WillowphOD <- exp(WillowCurve$a*WillowCurve$pd^2+WillowCurve$b*WillowCurve$pd+WillowCurve$c)-exp(OakCurve$a*OakCurve$pd^2+OakCurve$b*OakCurve$pd+OakCurve$c)

AlderPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AlderPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,14])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,4])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,24]))
AshPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
AshPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,15])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,5])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,25]))
BeechPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BeechPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,16])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,6])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,26]))
BirchPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
BirchPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,17])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,7])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,27]))
ElmPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
ElmPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,18])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,8])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,28]))
HazelPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
HazelPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,19])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,9])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,29]))
OakPW1      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
OakPW2      <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,20])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,10])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,30]))
RowanPW1    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
RowanPW2    <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,21])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,11])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,31]))
SycamorePW1 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
SycamorePW2 <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,22])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,12])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,32]))
WillowPW1   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])+sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))
WillowPW2   <- (-(AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])-sqrt((AbundTTCurves3$Sol[,2]+AbundTTCurves3$Sol[,23])^2-4*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33])*((AbundTTCurves3$Sol[,1]+AbundTTCurves3$Sol[,13])-log(0.01))))/(2*(AbundTTCurves3$Sol[,3]+AbundTTCurves3$Sol[,33]))

AlderpwOD <- (AlderPW2-AlderPW1)-(OakPW2-OakPW1)
AshpwOD <- (AshPW2-AshPW1)-(OakPW2-OakPW1)
BeechpwOD <- (BeechPW2-BeechPW1)-(OakPW2-OakPW1)
BirchpwOD <- (BirchPW2-BirchPW1)-(OakPW2-OakPW1)
ElmpwOD <- (ElmPW2-ElmPW1)-(OakPW2-OakPW1)
HazelpwOD <- (HazelPW2-HazelPW1)-(OakPW2-OakPW1)
RowanpwOD <- (RowanPW2-RowanPW1)-(OakPW2-OakPW1)
SycamorepwOD <- (SycamorePW2-SycamorePW1)-(OakPW2-OakPW1)
WillowpwOD <- (WillowPW2-WillowPW1)-(OakPW2-OakPW1)


TTcurvesOD <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                       PD=c(mean(AlderpdOD),mean(AshpdOD),mean(BeechpdOD),mean(BirchpdOD),mean(ElmpdOD),mean(HazelpdOD),mean(RowanpdOD),mean(SycamorepdOD),mean(WillowpdOD)),
                       PDLCI=c(HPDinterval(mcmc(AlderpdOD))[1], HPDinterval(mcmc(AshpdOD))[1], HPDinterval(mcmc(BeechpdOD))[1], HPDinterval(mcmc(BirchpdOD))[1], HPDinterval(mcmc(ElmpdOD))[1], HPDinterval(mcmc(HazelpdOD))[1],  HPDinterval(mcmc(RowanpdOD))[1], HPDinterval(mcmc(SycamorepdOD))[1], HPDinterval(mcmc(WillowpdOD))[1]),
                       PDUCI=c(HPDinterval(mcmc(AlderpdOD))[2], HPDinterval(mcmc(AshpdOD))[2], HPDinterval(mcmc(BeechpdOD))[2], HPDinterval(mcmc(BirchpdOD))[2], HPDinterval(mcmc(ElmpdOD))[2], HPDinterval(mcmc(HazelpdOD))[2],  HPDinterval(mcmc(RowanpdOD))[2], HPDinterval(mcmc(SycamorepdOD))[2], HPDinterval(mcmc(WillowpdOD))[2]),
                       PH=c(mean(AlderphOD),mean(AshphOD),mean(BeechphOD),mean(BirchphOD),mean(ElmphOD),mean(HazelphOD),mean(RowanphOD),mean(SycamorephOD),mean(WillowphOD)),
                       PHLCI=c(HPDinterval(mcmc(AlderphOD))[1], HPDinterval(mcmc(AshphOD))[1], HPDinterval(mcmc(BeechphOD))[1], HPDinterval(mcmc(BirchphOD))[1], HPDinterval(mcmc(ElmphOD))[1], HPDinterval(mcmc(HazelphOD))[1],  HPDinterval(mcmc(RowanphOD))[1], HPDinterval(mcmc(SycamorephOD))[1], HPDinterval(mcmc(WillowphOD))[1]),
                       PHUCI=c(HPDinterval(mcmc(AlderphOD))[2], HPDinterval(mcmc(AshphOD))[2], HPDinterval(mcmc(BeechphOD))[2], HPDinterval(mcmc(BirchphOD))[2], HPDinterval(mcmc(ElmphOD))[2], HPDinterval(mcmc(HazelphOD))[2],  HPDinterval(mcmc(RowanphOD))[2], HPDinterval(mcmc(SycamorephOD))[2], HPDinterval(mcmc(WillowphOD))[2]),
                       PW=c(mean(AlderpwOD),mean(AshpwOD, na.rm=TRUE),mean(BeechpwOD),mean(BirchpwOD),mean(ElmpwOD),mean(HazelpwOD),mean(RowanpwOD),mean(SycamorepwOD),mean(WillowpwOD)),
                       PWLCI=c(HPDinterval(mcmc(AlderpwOD))[1], HPDinterval(mcmc(AshpwOD))[1], HPDinterval(mcmc(BeechpwOD))[1], HPDinterval(mcmc(BirchpwOD))[1], HPDinterval(mcmc(ElmpwOD))[1], HPDinterval(mcmc(HazelpwOD))[1],  HPDinterval(mcmc(RowanpwOD))[1], HPDinterval(mcmc(SycamorepwOD))[1], HPDinterval(mcmc(WillowpwOD))[1]),
                       PWUCI=c(HPDinterval(mcmc(AlderpwOD))[2], HPDinterval(mcmc(AshpwOD))[2], HPDinterval(mcmc(BeechpwOD))[2], HPDinterval(mcmc(BirchpwOD))[2], HPDinterval(mcmc(ElmpwOD))[2], HPDinterval(mcmc(HazelpwOD))[2],  HPDinterval(mcmc(RowanpwOD))[2], HPDinterval(mcmc(SycamorepwOD))[2], HPDinterval(mcmc(WillowpwOD))[2]))
                       
AlderpwOD.5 <- (AlderPW2.5-AlderPW1.5)-(OakPW2.5-OakPW1.5)
AshpwOD.5 <- (AshPW2.5-AshPW1.5)-(OakPW2.5-OakPW1.5)
BeechpwOD.5 <- (BeechPW2.5-BeechPW1.5)-(OakPW2.5-OakPW1.5)
BirchpwOD.5 <- (BirchPW2.5-BirchPW1.5)-(OakPW2.5-OakPW1.5)
ElmpwOD.5 <- (ElmPW2.5-ElmPW1.5)-(OakPW2.5-OakPW1.5)
HazelpwOD.5 <- (HazelPW2.5-HazelPW1.5)-(OakPW2.5-OakPW1.5)
RowanpwOD.5 <- (RowanPW2.5-RowanPW1.5)-(OakPW2.5-OakPW1.5)
SycamorepwOD.5 <- (SycamorePW2.5-SycamorePW1.5)-(OakPW2.5-OakPW1.5)
WillowpwOD.5 <- (WillowPW2.5-WillowPW1.5)-(OakPW2.5-OakPW1.5)


TTcurvesOD.5 <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                         PD=c(mean(AlderpdOD),mean(AshpdOD),mean(BeechpdOD),mean(BirchpdOD),mean(ElmpdOD),mean(HazelpdOD),mean(RowanpdOD),mean(SycamorepdOD),mean(WillowpdOD)),
                         PDLCI=c(HPDinterval(mcmc(AlderpdOD))[1], HPDinterval(mcmc(AshpdOD))[1], HPDinterval(mcmc(BeechpdOD))[1], HPDinterval(mcmc(BirchpdOD))[1], HPDinterval(mcmc(ElmpdOD))[1], HPDinterval(mcmc(HazelpdOD))[1],  HPDinterval(mcmc(RowanpdOD))[1], HPDinterval(mcmc(SycamorepdOD))[1], HPDinterval(mcmc(WillowpdOD))[1]),
                         PDUCI=c(HPDinterval(mcmc(AlderpdOD))[2], HPDinterval(mcmc(AshpdOD))[2], HPDinterval(mcmc(BeechpdOD))[2], HPDinterval(mcmc(BirchpdOD))[2], HPDinterval(mcmc(ElmpdOD))[2], HPDinterval(mcmc(HazelpdOD))[2],  HPDinterval(mcmc(RowanpdOD))[2], HPDinterval(mcmc(SycamorepdOD))[2], HPDinterval(mcmc(WillowpdOD))[2]),
                         PH=c(mean(AlderphOD),mean(AshphOD),mean(BeechphOD),mean(BirchphOD),mean(ElmphOD),mean(HazelphOD),mean(RowanphOD),mean(SycamorephOD),mean(WillowphOD)),
                         PHLCI=c(HPDinterval(mcmc(AlderphOD))[1], HPDinterval(mcmc(AshphOD))[1], HPDinterval(mcmc(BeechphOD))[1], HPDinterval(mcmc(BirchphOD))[1], HPDinterval(mcmc(ElmphOD))[1], HPDinterval(mcmc(HazelphOD))[1],  HPDinterval(mcmc(RowanphOD))[1], HPDinterval(mcmc(SycamorephOD))[1], HPDinterval(mcmc(WillowphOD))[1]),
                         PHUCI=c(HPDinterval(mcmc(AlderphOD))[2], HPDinterval(mcmc(AshphOD))[2], HPDinterval(mcmc(BeechphOD))[2], HPDinterval(mcmc(BirchphOD))[2], HPDinterval(mcmc(ElmphOD))[2], HPDinterval(mcmc(HazelphOD))[2],  HPDinterval(mcmc(RowanphOD))[2], HPDinterval(mcmc(SycamorephOD))[2], HPDinterval(mcmc(WillowphOD))[2]),
                         PW.5=c(mean(AlderpwOD.5),mean(AshpwOD.5, na.rm=TRUE),mean(BeechpwOD.5),mean(BirchpwOD.5),mean(ElmpwOD.5),mean(HazelpwOD.5),mean(RowanpwOD.5),mean(SycamorepwOD.5),mean(WillowpwOD.5)),
                         PWLCI.5=c(HPDinterval(mcmc(AlderpwOD.5))[1], HPDinterval(mcmc(AshpwOD.5))[1], HPDinterval(mcmc(BeechpwOD.5))[1], HPDinterval(mcmc(BirchpwOD.5))[1], HPDinterval(mcmc(ElmpwOD.5))[1], HPDinterval(mcmc(HazelpwOD.5))[1],  HPDinterval(mcmc(RowanpwOD.5))[1], HPDinterval(mcmc(SycamorepwOD.5))[1], HPDinterval(mcmc(WillowpwOD.5))[1]),
                         PWUCI.5=c(HPDinterval(mcmc(AlderpwOD.5))[2], HPDinterval(mcmc(AshpwOD.5))[2], HPDinterval(mcmc(BeechpwOD.5))[2], HPDinterval(mcmc(BirchpwOD.5))[2], HPDinterval(mcmc(ElmpwOD.5))[2], HPDinterval(mcmc(HazelpwOD.5))[2],  HPDinterval(mcmc(RowanpwOD.5))[2], HPDinterval(mcmc(SycamorepwOD.5))[2], HPDinterval(mcmc(WillowpwOD.5))[2]))


PDODplot <- ggplot(TTcurvesOD, aes(TT, (PD*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PDUCI*175), ymin=(PDLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak date difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=NoOakCols)

PHODplot <- ggplot(TTcurvesOD, aes(TT, PH, colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=PHUCI, ymin=PHLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak height difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=NoOakCols)

PWODplot <- ggplot(TTcurvesOD, aes(TT, (PW*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWUCI*175), ymin=(PWLCI*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak width difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=NoOakCols)

PWOD.5plot <- ggplot(TTcurvesOD.5, aes(TT, (PW.5*175), colour=TT))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=(PWUCI.5*175), ymin=(PWLCI.5*175), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak width (0.5 height) difference")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "red")+
  theme(text = element_text(size=20))+
  guides(color = "none")+
  scale_colour_manual(values=NoOakCols)

col6 <- grid.arrange(PDODplot, nrow = 1, heights = 1)
col7 <- grid.arrange(PHODplot, nrow = 1, heights = 1)
col8 <- grid.arrange(PWODplot, nrow = 1, heights = 1)
col9 <- grid.arrange(PWOD.5plot, nrow = 1, heights = 1)
AbundCurvesODFig <- grid.arrange(col6, col7, col8, col9, ncol = 4, widths = c(1,1,1,1))


############################
#### Model output table ####   !!!!! NEEDS NEW MODEL NAME FOR EFFECTIVE SAMPLE SIZE
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(ATTCyearint$Sol[,1]),3)," (",
                      round(HPDinterval(ATTCyearint$Sol[,1])[1],3)," - ",
                      round(HPDinterval(ATTCyearint$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(ATTCyearint$Sol[,2]),3)," (",
                          round(HPDinterval(ATTCyearint$Sol[,2])[1],3)," - ",
                          round(HPDinterval(ATTCyearint$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(ATTCyearint$Sol[,3]),3)," (",
                           round(HPDinterval(ATTCyearint$Sol[,3])[1],3)," - ",
                           round(HPDinterval(ATTCyearint$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,3]))))

####random 
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-3
treetaxa3<-c("TreeTaxa- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-5
treetaxa5<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-6
treetaxa6<-c("TreeTaxa- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-9
treetaxa9<-c("TreeTaxa- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-10
site10<-c("Site- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-11
site11<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-12
site12<-c("Site- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-14
site14<-c("Site- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-15
site15<-c("Site- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-18
site18<-c("Site- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(ATTCyearint$VCV[, column])))

column<-19
year10<-c("Year- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-20
year11<-c("Year- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-21
year12<-c("Year- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-23
year14<-c("Year- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-24
year15<-c("Year- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-27
year18<-c("Year- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-31
recorder<-c("Recorder",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))


column<-30
siteday<-c("Site-Day",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(ATTCyearint$VCV[, column])))

column<-28
siteyear<-c("Site-Year",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))

column<-29
treeID<-c("Tree ID",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                          round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                          round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-32
residual<-c("Residual",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))




random<-rbind(treetaxa1,treetaxa2,treetaxa3,treetaxa5,treetaxa6,treetaxa9,site10,site11,site12,site14,site15,site18,year10,year11,year12,year14,year15,year18, siteyear, recorder, siteday, treeID, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableATTCyearint.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
