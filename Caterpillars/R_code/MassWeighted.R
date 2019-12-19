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
cater_habitat$datescaled <- cater_habitat$date/max(cater_habitat$date)

#!!!!!!!!!!! if removing OthDecid !!!!!!!!!!!!
cater_habitat<- subset(cater_habitat, tree.species!="OthDecid")

###################
#### Weighting ####
###################

# changing caterpillars to 1 so no infinity values in weighting equation: us(sqrt(1/caterpillars)):units
cater_habitat$weight <- as.numeric(revalue(as.character(cater_habitat$caterpillars), c("0"="1")))

# reordering so year 2017 is first
cater_habitat$year <- relevel(cater_habitat$year, ref="2017")

# Model with full data set
k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

system.time(MassWeightedFull<- MCMCglmm(cbind(logmpc1, logmpc2)~ year + datescaled + I(datescaled^2), 
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
                 family="cengaussian", data=cater_habitat, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=500))
save(MassWeightedFull, file = "~/Documents/Models/MassWeightedFull.RData")
load("~/Documents/Models/MassWeightedFull.RData")

# Model with just 2017-19
cater_habitat_1719 <- subset(cater_habitat, year!="2014")
cater_habitat_1719 <- subset(cater_habitat_1719, year!="2015")
cater_habitat_1719 <- subset(cater_habitat_1719, year!="2016")

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

system.time(MassWeighted1719<- MCMCglmm(cbind(logmpc1, logmpc2)~ year + datescaled + I(datescaled^2), 
                 random=~us(1+datescaled):tree.species + us(1+datescaled):site + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
                 family="cengaussian", data=cater_habitat_1719, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=500))
save(MassWeighted1719, file = "~/Documents/Models/MassWeighted1719.RData")
load("~/Documents/Models/MassWeighted1719.RData")

# With interaction with date squared
k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

system.time(MWd2int<- MCMCglmm(cbind(logmpc1, logmpc2)~ year + datescaled + I(datescaled^2), 
                                        random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
                                        family="cengaussian", data=cater_habitat_1719, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=250))
save(MWd2int, file = "~/Documents/Models/MWd2int.RData")
load("~/Documents/Models/MWd2int.RData")

#MWd2intAll : MWd2int with all tree taxa categories
#MWd2intOth : MWd2int with other deciduous kept in
load("~/Documents/Models/MWd2intAll.RData")
load("~/Documents/Models/MWd2intOth.RData")

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

system.time(MWd2intSY<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
                               random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
                               family="cengaussian", data=cater_habitat_1719, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=250))
save(MWd2intSY, file = "~/Documents/Models/MWd2intSY.RData")
load("~/Documents/Models/MWd2intSY.RData")

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#system.time(MWSY<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
#                            random=~us(1+datescaled):tree.species + us(1+datescaled):siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
#                            family="cengaussian", data=cater_habitat_1719, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=250))
#save(MWSY, file = "~/Documents/Models/MWSY.RData")
load("~/Documents/Models/MWSY.RData")

system.time(MWSY1<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
                            random=~us(1+datescaled):tree.species + us(1+datescaled):siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
                            family="cengaussian", data=cater_habitat_1719, prior=prior3, nitt=1750000, burnin=50000, pr=TRUE, thin=100))
save(MWSY1, file = "~/Documents/Models/MWSY1.RData")
rm(list=ls())
load("~/Documents/Models/MWSY1.RData")

 ############################
#### Model output table ####   !!!!! NEEDS NEW MODEL NAME FOR EFFECTIVE SAMPLE SIZE
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(MWSY$Sol[,1]),3)," (",
                      round(HPDinterval(MWSY$Sol[,1])[1],3)," - ",
                      round(HPDinterval(MWSY$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(MWSY$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(MWSY$Sol[,2]),3)," (",
                          round(HPDinterval(MWSY$Sol[,2])[1],3)," - ",
                          round(HPDinterval(MWSY$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(MWSY$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(MWSY$Sol[,3]),3)," (",
                           round(HPDinterval(MWSY$Sol[,3])[1],3)," - ",
                           round(HPDinterval(MWSY$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(MWSY$Sol[,3]))))

####random 
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                             round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(MWSY$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                                          round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(MWSY$VCV[, column])))


column<-4
treetaxa4<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                              round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(MWSY$VCV[, column])))


column<-5
siteyear5<-c("SiteYear- Intercept var",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                              round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(MWSY$VCV[, column])))

column<-6
siteyear6<-c("SiteYear- Intercept:Date slope covar",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                                           round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(MWSY$VCV[, column])))


column<-8
siteyear8<-c("SiteYear- Date slope var",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                                               round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
              round(effectiveSize(MWSY$VCV[, column])))


column<-11
recorder<-c("Recorder",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                             round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                             round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(MWSY$VCV[, column])))


column<-10
siteday<-c("Site Day",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                            round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                            round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(MWSY$VCV[, column])))


column<-9
treeID<-c("Tree ID",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                          round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                          round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(MWSY$VCV[, column])))

column<-12
weight<-c("Weighting",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                          round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                          round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(MWSY$VCV[, column])))

column<-13
residual<-c("Residual",paste(round(posterior.mode(MWSY$VCV[, column]),3)," (",
                             round(HPDinterval(MWSY$VCV[, column])[1],3)," - ",
                             round(HPDinterval(MWSY$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(MWSY$VCV[, column])))




random<-rbind(treetaxa1,treetaxa2,treetaxa4,siteyear5,siteyear6,siteyear8, recorder, siteday, treeID, weight, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableMWSY.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
