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

#################################################
#### Plotting trends over time for each year ####
#################################################

#### Organising Habitat and TreeTaxa Categories ####

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



cater14 <- subset(cater_habitat, year=="2014", select = c(site, date, caterpillars))
cater15 <- subset(cater_habitat, year=="2015", select = c(site, date, caterpillars))
cater16 <- subset(cater_habitat, year=="2016", select = c(site, date, caterpillars))
cater17 <- subset(cater_habitat, year=="2017", select = c(site, date, caterpillars))
cater18 <- subset(cater_habitat, year=="2018", select = c(site, date, caterpillars))
cater19 <- subset(cater_habitat, year=="2019", select = c(site, date, caterpillars))

means19 <- data.frame(site=cater19$site, date=cater19$date, cater=cater19$caterpillars)
means19$siteday <- paste(means19$site, means19$date)
means19$date <- as.factor(means19$date)
sitemeans19 <- aggregate(means19$cater~means19$siteday, FUN=mean)
sitemeans19$siteday <- sitemeans19$`means19$siteday`
sitemeans19$cater <- sitemeans19$`means19$cater`
store <- pmatch(sitemeans19$siteday, means19$siteday)
sitemeans19 <- cbind(sitemeans19, means19$date[store])
sitemeans19$date <- sitemeans19$`means19$date[store]`
sitemeans19 <- aggregate(sitemeans19$cater~sitemeans19$date, FUN=mean)

means18 <- data.frame(site=cater18$site, date=cater18$date, cater=cater18$caterpillars)
means18$siteday <- paste(means18$site, means18$date)
means18$date <- as.factor(means18$date)
sitemeans18 <- aggregate(means18$cater~means18$siteday, FUN=mean)
sitemeans18$siteday <- sitemeans18$`means18$siteday`
sitemeans18$cater <- sitemeans18$`means18$cater`
store <- pmatch(sitemeans18$siteday, means18$siteday)
sitemeans18 <- cbind(sitemeans18, means18$date[store])
sitemeans18$date <- sitemeans18$`means18$date[store]`
sitemeans18 <- aggregate(sitemeans18$cater~sitemeans18$date, FUN=mean)

means17 <- data.frame(site=cater17$site, date=cater17$date, cater=cater17$caterpillars)
means17$siteday <- paste(means17$site, means17$date)
means17$date <- as.factor(means17$date)
sitemeans17 <- aggregate(means17$cater~means17$siteday, FUN=mean)
sitemeans17$siteday <- sitemeans17$`means17$siteday`
sitemeans17$cater <- sitemeans17$`means17$cater`
store <- pmatch(sitemeans17$siteday, means17$siteday)
sitemeans17 <- cbind(sitemeans17, means17$date[store])
sitemeans17$date <- sitemeans17$`means17$date[store]`
sitemeans17 <- aggregate(sitemeans17$cater~sitemeans17$date, FUN=mean)

means16 <- data.frame(site=cater16$site, date=cater16$date, cater=cater16$caterpillars)
means16$siteday <- paste(means16$site, means16$date)
means16$date <- as.factor(means16$date)
sitemeans16 <- aggregate(means16$cater~means16$siteday, FUN=mean)
sitemeans16$siteday <- sitemeans16$`means16$siteday`
sitemeans16$cater <- sitemeans16$`means16$cater`
store <- pmatch(sitemeans16$siteday, means16$siteday)
sitemeans16 <- cbind(sitemeans16, means16$date[store])
sitemeans16$date <- sitemeans16$`means16$date[store]`
sitemeans16 <- aggregate(sitemeans16$cater~sitemeans16$date, FUN=mean)

means15 <- data.frame(site=cater15$site, date=cater15$date, cater=cater15$caterpillars)
means15$siteday <- paste(means15$site, means15$date)
means15$date <- as.factor(means15$date)
sitemeans15 <- aggregate(means15$cater~means15$siteday, FUN=mean)
sitemeans15$siteday <- sitemeans15$`means15$siteday`
sitemeans15$cater <- sitemeans15$`means15$cater`
store <- pmatch(sitemeans15$siteday, means15$siteday)
sitemeans15 <- cbind(sitemeans15, means15$date[store])
sitemeans15$date <- sitemeans15$`means15$date[store]`
sitemeans15 <- aggregate(sitemeans15$cater~sitemeans15$date, FUN=mean)

means14 <- data.frame(site=cater14$site, date=cater14$date, cater=cater14$caterpillars)
means14$siteday <- paste(means14$site, means14$date)
means14$date <- as.factor(means14$date)
sitemeans14 <- aggregate(means14$cater~means14$siteday, FUN=mean)
sitemeans14$siteday <- sitemeans14$`means14$siteday`
sitemeans14$cater <- sitemeans14$`means14$cater`
store <- pmatch(sitemeans14$siteday, means14$siteday)
sitemeans14 <- cbind(sitemeans14, means14$date[store])
sitemeans14$date <- sitemeans14$`means14$date[store]`
sitemeans14 <- aggregate(sitemeans14$cater~sitemeans14$date, FUN=mean)

par(mfrow=c(1,1))
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="Ordinal Date", ylab="Mean caterpillar abundance", xlim=c(117,176))
points(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, col=2)
points(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, col=3)
points(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, col=4)
points(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, col=5)
points(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, col=6)

par(mfrow=c(3,2),mar=c(4.1, 4.0, 2.5, 1.5), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2014", ylim=c(0,0.65), yaxt='n')
axis(side = 2, at = c(0,.15, .3,.45, .6), labels = c( '0.0','', '0.3','', '0.6'))
legend("topleft", legend="n = 1106", bty = "n")
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,0.22), yaxt='n')
axis(side = 2, at = c(0,.05, .1,.15, .2), labels = c( '0.0','', '0.1','', '0.2'))
legend("topleft", legend="n = 2632", bty = "n")
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2016", ylim=c(0,0.42), yaxt='n')
axis(side = 2, at = c(0,.1, .2,.3, .4), labels = c( '0.0','', '0.2','', '0.4'))
legend("topleft", legend="n = 2351", bty = "n")
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,0.23), yaxt='n')
axis(side = 2, at = c(0,.05, .1,.15, .2), labels = c( '0.0','', '0.1','', '0.2'))
legend("topleft", legend="n = 6877", bty = "n")
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2018", ylim=c(0,0.82), yaxt='n')
axis(side = 2, at = c(0,.2, .4,.6, .8), labels = c( '0.0','', '0.4','', '0.8'))
legend("topleft", legend="n = 6791", bty = "n")
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2019", ylim=c(0,3), yaxt='n')
axis(side = 2, at = c(0,0.75,1.5,2.25, 3), labels = c( '0.0','', '1.5','', '3.0'))
legend("topleft", legend="n = 8310", bty = "n") # saved as 8"x8"

par(mfrow=c(3,2),mar=c(4.1, 4.0, 2.5, 1.5), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2014", ylim=c(0,2.7), xaxt="n")
legend("topleft", legend="n = 1106", bty = "n")
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,2.7))
legend("topleft", legend="n = 2632", bty = "n")
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2016", ylim=c(0,2.7))
legend("topleft", legend="n = 2351", bty = "n")
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,2.7))
legend("topleft", legend="n = 6877", bty = "n")
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2018", ylim=c(0,2.7))
legend("topleft", legend="n = 6791", bty = "n")
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2019", ylim=c(0,2.7))
legend("topleft", legend="n = 8310", bty = "n") # saved as 8"x8"

par(mfrow=c(1,6),mar=c(4.1, 4.0, 2.5, 1), cex=0.9, cex.lab=1.25, cex.axis=1.2)
plot(as.numeric(as.character(sitemeans14$`sitemeans14$date`)),sitemeans14$`sitemeans14$cater`, pch=20, xlab="", ylab="Mean caterpillar abundance", xlim=c(117,176), main="2014", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n=1106")
plot(as.numeric(as.character(sitemeans15$`sitemeans15$date`)),sitemeans15$`sitemeans15$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2015", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 2632")
plot(as.numeric(as.character(sitemeans16$`sitemeans16$date`)),sitemeans16$`sitemeans16$cater`, pch=20, xlab="Ordinal Date", ylab="", xlim=c(117,176), main="2016", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 2351")
plot(as.numeric(as.character(sitemeans17$`sitemeans17$date`)),sitemeans17$`sitemeans17$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2017", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 6877")
plot(as.numeric(as.character(sitemeans18$`sitemeans18$date`)),sitemeans18$`sitemeans18$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2018", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 6791")
plot(as.numeric(as.character(sitemeans19$`sitemeans19$date`)),sitemeans19$`sitemeans19$cater`, pch=20, xlab="", ylab="", xlim=c(117,176), main="2019", ylim=c(0,2.7), xaxt="n")
axis(side = 1, at = c(120,145,170), labels = c( '120','','170'))
text(145,2.7, "n = 8310") # saved as 5"x12"
