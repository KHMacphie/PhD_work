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

###################################################################
#### Model for testing asymmetry in caterpillar abundance peak ####
###################################################################

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
cater$tree.species<- revalue(cater$tree.species, c("?"="OthDecid", "Cherry"="OthDecid", "Aspen"="OthDecid", "Chestnut"="OthDecid", "Lime"="OthDecid", "Field maple"="OthDecid", "Damson"="OthDecid", "Whitebeam"="OthDecid"))

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
cater_habitat$weight <- as.numeric(revalue(as.character(cater_habitat$caterpillars), c("0"="1")))

#!!!!!!!!!!! if removing OthDecid !!!!!!!!!!!!
cater_habitat<- subset(cater_habitat, tree.species!="OthDecid")

##################################################################
#### Model: cubic curve, siteyear date and date^2 interaction ####
##################################################################
a<-10000
prior<-list(R=list(V=diag(1), nu=0.002), 
            G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))

#CurveShape<- MCMCglmm(caterpillars~datescaled+I(datescaled^2)+I(datescaled^3), random=~us(1+datescaled+I(datescaled^2)):siteyear+recorder+siteday+treeID, family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000, thin=50)
#save(CurveShape, file = "~/Documents/Models/CurveShape.RData")
load("~/Documents/Models/CurveShape.RData")


#### Calculate peak date
# (-b +/- sqrt(b^2-4*a*c))/(2*a)
cdf <- data.frame(CurveShape$Sol[,1:4])
cdf$a <- 3*CurveShape$Sol[,4]
cdf$b <- 2*CurveShape$Sol[,3]
cdf$c <- CurveShape$Sol[,2]
#abline(v=((-b - sqrt(b^2-4*a*c))/(2*a))*max(cater$date), col=1, lty="dashed")


# Not reverted to day of the year for height calcs
cdf$pd <- ((-cdf$b - sqrt(cdf$b^2-4*cdf$a*cdf$c))/(2*cdf$a))

#### Calculate peak height

cdf$ph <- CurveShape$Sol[,1]+CurveShape$Sol[,2]*cdf$pd+CurveShape$Sol[,3]*cdf$pd^2+CurveShape$Sol[,4]*cdf$pd^3


#### Calculate width either side at 50% peak height
for(i in 1:5400){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/2))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.5[i] <- A[1]
  cdf$r2.5[i] <- A[2]
  cdf$r3.5[i] <- A[3]
}

cdf$r1.5 <- Re(cdf$r1.5)[abs(Im(cdf$r1.5)) < 1e-6]
cdf$r2.5 <- Re(cdf$r2.5)[abs(Im(cdf$r2.5)) < 1e-6]
cdf$r3.5 <- Re(cdf$r3.5)[abs(Im(cdf$r3.5)) < 1e-6]

cdf$r1.5 <- ifelse(cdf$r1.5<min(cater_habitat$datescaled),NA,cdf$r1.5)
cdf$r2.5 <- ifelse(cdf$r2.5<min(cater_habitat$datescaled),NA,cdf$r2.5)
cdf$r3.5 <- ifelse(cdf$r3.5<min(cater_habitat$datescaled),NA,cdf$r3.5)


roots.5 <- data.frame(r1.5=cdf$r1.5,r2.5=cdf$r2.5,r3.5=cdf$r3.5)
cdf$root1.5 <- (apply(roots.5, 1, min, na.rm=TRUE))*max(cater_habitat$date)   
cdf$root2.5 <- (apply(roots.5, 1, max, na.rm=TRUE))*max(cater_habitat$date)  

cdf$left.5 <- (cdf$pd*max(cater_habitat$date))-cdf$root1.5
cdf$right.5 <- cdf$root2.5-(cdf$pd*max(cater_habitat$date))
cdf$width.5 <- cdf$left.5+cdf$right.5
cdf$propleft.5 <- cdf$left.5/cdf$width.5
cdf$propright.5 <- cdf$right.5/cdf$width.5

mean(cdf$propleft.5) # 0.5136632
HPDinterval(cdf$propleft.5) # 0.4591579 0.5601217
mean(cdf$propright.5) # 0.4863368
HPDinterval(cdf$propright.5) # 0.0.4398783 0.5408421
mean(cdf$width.5) # 24.47275
HPDinterval(cdf$width.5) # 23.62189 25.18829

#### Calculate width either side at 25% peak height
for(i in 1:5400){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])/4))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.25[i] <- A[1]
  cdf$r2.25[i] <- A[2]
  cdf$r3.25[i] <- A[3]
}

cdf$r1.25 <- Re(cdf$r1.25)[abs(Im(cdf$r1.25)) < 1e-6]
cdf$r2.25 <- Re(cdf$r2.25)[abs(Im(cdf$r2.25)) < 1e-6]
cdf$r3.25 <- Re(cdf$r3.25)[abs(Im(cdf$r3.25)) < 1e-6]

cdf$r1.25 <- ifelse(cdf$r1.25<min(cater_habitat$datescaled),NA,cdf$r1.25)
cdf$r2.25 <- ifelse(cdf$r2.25<min(cater_habitat$datescaled),NA,cdf$r2.25)
cdf$r3.25 <- ifelse(cdf$r3.25<min(cater_habitat$datescaled),NA,cdf$r3.25)


roots.25 <- data.frame(r1.25=cdf$r1.25,r2.25=cdf$r2.25,r3.25=cdf$r3.25)
cdf$root1.25 <- (apply(roots.25, 1, min, na.rm=TRUE))*max(cater_habitat$date)   
cdf$root2.25 <- (apply(roots.25, 1, max, na.rm=TRUE))*max(cater_habitat$date)  

cdf$left.25 <- (cdf$pd*max(cater_habitat$date))-cdf$root1.25
cdf$right.25 <- cdf$root2.25-(cdf$pd*max(cater_habitat$date))
cdf$width.25 <- cdf$left.25+cdf$right.25
cdf$propleft.25 <- cdf$left.25/cdf$width.25
cdf$propright.25 <- cdf$right.25/cdf$width.25

mean(cdf$propleft.25) # 0.5520194
HPDinterval(cdf$propleft.25) # 0.5168431 0.5863098
mean(cdf$propright.25) # 0.4479806
HPDinterval(cdf$propright.25) # 0.4136902 0.4831569
mean(cdf$width.25) # 35.79095
HPDinterval(cdf$width.25) # 34.18659 36.48023

#### Calculate width either side at 75% peak height
for(i in 1:5400){
  A <- polyroot(c((cdf$X.Intercept.[i]-(log(exp(cdf$ph[i])*0.75))),cdf$datescaled[i],cdf$I.datescaled.2.[i],cdf$I.datescaled.3.[i]))
  cdf$r1.75[i] <- A[1]
  cdf$r2.75[i] <- A[2]
  cdf$r3.75[i] <- A[3]
}

cdf$r1.75 <- Re(cdf$r1.75)[abs(Im(cdf$r1.75)) < 1e-6]
cdf$r2.75 <- Re(cdf$r2.75)[abs(Im(cdf$r2.75)) < 1e-6]
cdf$r3.75 <- Re(cdf$r3.75)[abs(Im(cdf$r3.75)) < 1e-6]

cdf$r1.75 <- ifelse(cdf$r1.75<min(cater_habitat$datescaled),NA,cdf$r1.75)
cdf$r2.75 <- ifelse(cdf$r2.75<min(cater_habitat$datescaled),NA,cdf$r2.75)
cdf$r3.75 <- ifelse(cdf$r3.75<min(cater_habitat$datescaled),NA,cdf$r3.75)


roots.75 <- data.frame(r1.75=cdf$r1.75,r2.75=cdf$r2.75,r3.75=cdf$r3.75)
cdf$root1.75 <- (apply(roots.75, 1, min, na.rm=TRUE))*max(cater_habitat$date)   
cdf$root2.75 <- (apply(roots.75, 1, max, na.rm=TRUE))*max(cater_habitat$date)  

cdf$left.75 <- (cdf$pd*max(cater_habitat$date))-cdf$root1.75
cdf$right.75 <- cdf$root2.75-(cdf$pd*max(cater_habitat$date))
cdf$width.75 <- cdf$left.75+cdf$right.75
cdf$propleft.75 <- cdf$left.75/cdf$width.75
cdf$propright.75 <- cdf$right.75/cdf$width.75

mean(cdf$propleft.75) # 0.4729583
HPDinterval(cdf$propleft.75) # 0.3856603 0.5449088
mean(cdf$propright.75) # 0.5270417
HPDinterval(cdf$propright.75) # 0.4550912 0.6143397
mean(cdf$width.75) # 15.51931
HPDinterval(cdf$width.75) # 15.06974 16.08495
#### Plot curves
# colour: 
mycol <- rgb(0, 153, 0, max = 250, alpha = 10, names = "greentrans")

dayscal <- seq(0.67,1,0.001)
curve <- mean(CurveShape$Sol[,1])+mean(CurveShape$Sol[,2])*dayscal+mean(CurveShape$Sol[,3])*dayscal^2+mean(CurveShape$Sol[,4])*dayscal^3
days <- dayscal*max(cater_habitat$date)
quart <- data.frame(qd=seq(mean(cdf$root1.25),(mean(cdf$root2.25)-0.5),0.1))
quart$qh <- mean(exp(cdf$ph)/4)
half <- data.frame(hd=seq((mean(cdf$root1.5)-0.5),(mean(cdf$root2.5)-0.75),0.1))
half$hh <- mean(exp(cdf$ph)/2)
tquart <- data.frame(tqd=seq((mean(cdf$root1.75)-0.5),(mean(cdf$root2.75)-0.75),0.1))
tquart$tqh <- mean(exp(cdf$ph)*0.75)

par(mfcol=c(1,1),mar=c(3.9, 3.8, 1, 1), cex=1.4, las=1)
plot(days,exp(curve), type="l", ylim=c(0,0.095), xlab="Date", ylab="Abundance", yaxs="i")

for(i in 1:5400){
  A <- CurveShape$Sol[i,1]+CurveShape$Sol[i,2]*dayscal+CurveShape$Sol[i,3]*dayscal^2+CurveShape$Sol[i,4]*dayscal^3
  points(days, exp(A), type="l", col=mycol, lwd=0.5)
}

points(quart$qd, quart$qh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(half$hd, half$hh, type="l", lty="dashed", lwd=0.7, col="gray66")
points(tquart$tqd, tquart$tqh, type="l", lty="dashed", lwd=0.7, col="gray66")
abline(v=(mean(cdf$pd)*max(cater$date)), lwd=0.8, lty="dashed", col="gray66")
points(days, exp(curve), type="l")

text(150, 0.0125, "55.2%", cex=0.9, col="gray40")
text(159.5, 0.0125, "44.8%", cex=0.9, col="gray40")
text(150, 0.028, "51.4%", cex=0.9, col="gray40")
text(159.5, 0.028, "48.6%", cex=0.9, col="gray40")
text(150.3, 0.044, "47.3%", cex=0.9, col="gray40")
text(159.3, 0.044, "52.7%", cex=0.9, col="gray40")

text(125, mean(quart$qh), "0.25", cex=0.9, col="gray66")
text(125.4, mean(half$hh), "0.5", cex=0.9, col="gray66")
text(125, mean(tquart$tqh), "0.75", cex=0.9, col="gray66")

arrows(x0=127.5, y0=mean(half$hh), x1=136, y1=mean(half$hh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(quart$qh), x1=130.25, y1=mean(quart$qh), length=0.1, col="gray66")
arrows(x0=127.5, y0=mean(tquart$tqh), x1=140, y1=mean(tquart$tqh), length=0.1, col="gray66")

############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(CurveShape$Sol[,1]),3)," (",
                      round(HPDinterval(CurveShape$Sol[,1])[1],3)," - ",
                      round(HPDinterval(CurveShape$Sol[,1])[2],3),")",sep=""),
                      round(effectiveSize(CurveShape$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(CurveShape$Sol[,2]),3)," (",
                      round(HPDinterval(CurveShape$Sol[,2])[1],3)," - ",
                      round(HPDinterval(CurveShape$Sol[,2])[2],3),")",sep=""),
                      round(effectiveSize(CurveShape$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(CurveShape$Sol[,3]),3)," (",
                      round(HPDinterval(CurveShape$Sol[,3])[1],3)," - ",
                      round(HPDinterval(CurveShape$Sol[,3])[2],3),")",sep=""),
                      round(effectiveSize(CurveShape$Sol[,3]))),
  c("Date³ (scaled)",paste(round(mean(CurveShape$Sol[,4]),3)," (",
                      round(HPDinterval(CurveShape$Sol[,4])[1],3)," - ",
                      round(HPDinterval(CurveShape$Sol[,4])[2],3),")",sep=""),
                      round(effectiveSize(CurveShape$Sol[,4]))))

####random 
column<-1
siteyear1<-c("SiteYear- Intercept var",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-2
siteyear2<-c("SiteYear- Intercept:Date slope covar",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-3
siteyear3<-c("SiteYear- Intercept:Date² slope covar",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-5
siteyear5<-c("SiteYear- Date slope var",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-6
siteyear6<-c("SiteYear- Date slope:Date² slope covar",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-9
siteyear9<-c("SiteYear- Date² slope var",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                              round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                              round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
                              round(effectiveSize(CurveShape$VCV[, column])))

column<-10
recorder<-c("Recorder",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape$VCV[, column])))


column<-11
siteday<-c("Site Day",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                            round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                            round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(CurveShape$VCV[, column])))


column<-12
treeID<-c("Tree ID",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                          round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                          round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(CurveShape$VCV[, column])))

column<-13
residual<-c("Residual",paste(round(posterior.mode(CurveShape$VCV[, column]),3)," (",
                             round(HPDinterval(CurveShape$VCV[, column])[1],3)," - ",
                             round(HPDinterval(CurveShape$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(CurveShape$VCV[, column])))




random<-rbind(siteyear1,siteyear2,siteyear3,siteyear5,siteyear6,siteyear9, recorder, siteday, treeID, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableCurveShape.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
