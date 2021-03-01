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

######################################################################
#### Organising data frame for caterpillar habitat paper analyses ####
######################################################################

# Categories: Alder, Ash, Birch, Beech, Elm, Hazel, Oak, Rowan, Sycamore, Willow

#### Habitat data set ####
habitat <- read.csv("~/Dropbox/2018/Habitats_2018_all.csv")
# recorded as habitat at 15m radius around each nest box (NB)
table(habitat$site) #checking all there

#### Stands, thickets and small trees combined as small trees and other deciduous combined
# X6= stand of at least 6 trunks (1=0.5 small trees)
# X21= stand of at least 21 trunks (1=1 small tree)
# z= thicket (in dataframe as no. of small trees already)
# s,m,l= small, medium, large

# combining other deciduous (anything deciduous and not one of the 10 tree taxa catergories)
habitat$X6_oth.decid <- rowSums(habitat[,c("X6_cherry", "X6_elder", "X6_holly", "X6_rose")], na.rm=TRUE)
habitat$X21_oth.decid <- rowSums(habitat[,c("X21_blackthorn", "X21_holly")], na.rm=TRUE)
habitat$s_oth.decid <- rowSums(habitat[,c("s_aspen", "s_cherry", "s_chestnut", "s_elder", "s_hawthorn", "s_holly", "s_lime", "s_whitebeam", "s_other")], na.rm=TRUE)
habitat$m_oth.decid <- rowSums(habitat[,c("m_aspen", "m_cherry", "m_chestnut", "m_holly", "m_lime", "m_whitebeam", "m_horsechestnut")], na.rm=TRUE)
habitat$l_oth.decid <- habitat$l_cherry
habitat$z_oth.decid <- rowSums(habitat[,c("z_blackthorn", "z_cherry", "z_other")], na.rm=TRUE)

# combining conifers (pine, yew, juniper and conifer)
habitat$X6_conifer_all <- habitat$X6_juniper
habitat$s_conifer_all <- rowSums(habitat[,c("s_conifer", "s_pine", "s_yew", "s_juniper")], na.rm=TRUE)
habitat$m_conifer_all <- rowSums(habitat[,c("m_conifer", "m_pine", "m_yew")], na.rm=TRUE)
habitat$l_conifer_all <- rowSums(habitat[,c("l_conifer", "l_pine")], na.rm=TRUE)

# converting X6 stands into small trees
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


#### Use the weighting function to calculate foliage as number of small trees per tree type at each NB

weighting <- function(x){return(pi*(x/(2*pi))^2)} # function for the weighting equation 
# weighted by cross sectional area of the trunk in cm^2 (min possible circumference for each category)
Sm <- weighting(40) # 40=circumference at breast height
Med <- weighting(100)
Lar <- weighting(250)
S=1 #Pointless I know but thats how I did it originally so it was consistent to look at beneath
M=Med/Sm
L=Lar/Sm

# Weighting each taxon at each size
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
Habitat_byNB[is.na(Habitat_byNB)] <- 0 # Any taxon not present as zero rather than NA

#Foliage score for each tree taxon at each NB
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

## combining NBs within each site- mean to account for different number of NBs at sites 
#so mean foliage density for each taxa within a circle of 15m radius
Habitat_Site <- Habitat_byNB[,35:46] # isolate FS columns
Habitat_Site$Site <- Habitat_byNB$Site # add site codes
Habitat_Site <- aggregate(.~Site, Habitat_Site, mean) # mean within each site

## Globally mean centre FS's = each site tree taxa FS - mean across all tree taxa at all sites
Habitat_Site$Alder_cent <- Habitat_Site$Alder_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Ash_cent <- Habitat_Site$Ash_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Beech_cent <- Habitat_Site$Beech_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Birch_cent <- Habitat_Site$Birch_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Elm_cent <- Habitat_Site$Elm_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Hazel_cent <- Habitat_Site$Hazel_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Oak_cent <- Habitat_Site$Oak_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Rowan_cent <- Habitat_Site$Rowan_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Sycamore_cent <- Habitat_Site$Sycamore_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Willow_cent <- Habitat_Site$Willow_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$Conifer_cent <- Habitat_Site$Conifer_FS-mean(as.matrix(Habitat_Site[2:13]))
Habitat_Site$OthDecid_cent <- Habitat_Site$OthDecid_FS-mean(as.matrix(Habitat_Site[2:13]))

## Globally mean centre FS's = each site tree taxa FS - mean across all tree taxa at all sites
#Habitat_Site$Alder_scal <- (Habitat_Site$Alder_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13])) #mean= 5.790518 sd=12.4424
#Habitat_Site$Ash_scal <- (Habitat_Site$Ash_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Beech_scal <- (Habitat_Site$Beech_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Birch_scal <- (Habitat_Site$Birch_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Elm_scal <- (Habitat_Site$Elm_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Hazel_scal <- (Habitat_Site$Hazel_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Oak_scal <- (Habitat_Site$Oak_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Rowan_scal <- (Habitat_Site$Rowan_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Sycamore_scal <- (Habitat_Site$Sycamore_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Willow_scal <- (Habitat_Site$Willow_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$Conifer_scal <- (Habitat_Site$Conifer_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))
#Habitat_Site$OthDecid_scal <- (Habitat_Site$OthDecid_FS-mean(as.matrix(Habitat_Site[2:13])))/sd(as.matrix(Habitat_Site[2:13]))

#TotalFS as sum of tree taxa globally mean centred FS's at each site
Habitat_Site$Total_cent <- rowSums(Habitat_Site[14:25])
#Habitat_Site$Total_scal <- rowSums(Habitat_Site[14:25])

#### TreeTaxa (beaten) categories ####
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater <- subset(cater, year!=2020) #NEW ADDITION: remove 2020 because not in analyses
# Creating other deciduous for all trees not of the focal 10 tree taxa
cater$tree.species<- revalue(cater$tree.species, c("?"="OthDecid", "Cherry"="OthDecid", "Aspen"="OthDecid", "Chestnut"="OthDecid", "Lime"="OthDecid", "Field maple"="OthDecid", "Damson"="OthDecid", "Whitebeam"="OthDecid"))


#### Organising interval censored mass ####
#anything recorded as 0.02 or less is being censored to be a min 0.001 (per caterpillar) and max 0.02 (whole sample)
cater$caterpillar.mass <- revalue(cater$caterpillar.mass, c("0"="")) #removing zeros 
cater$caterpillar.mass2 <- revalue(cater$caterpillar.mass, c("<0.01"="0.02", "0.01"="0.02")) #upper sample mass 0.02 if being interval censored
cater$caterpillar.mass2 <- as.numeric(cater$caterpillar.mass2) #needs to be numeric

#upper bound
cater$mpc2 <- cater$caterpillar.mass2/cater$caterpillars # mass per caterpillar

#lower bound
cater$mpc1 <- as.character(cater$caterpillar.mass2) #sample mass upper bound
cater$mpc1 <- revalue(cater$mpc1, c("0.02"="0.001")) #change interval censoring samples to lower bound
cater$mpc1 <- as.numeric(cater$mpc1) #make numeric again
cater$mpc1 <- cater$mpc1/cater$caterpillars #divide by abundance
cater$mpc1 <- ifelse(cater$mpc1 < 0.001,0.001,cater$mpc1) # interval censored ones back to lower bound

## Sort rows with mass.uncertain: means the mass may not be for the abundance recorded
#manually changing those cells based on the notes so data not lost
which(cater$mass.uncertain=="1") # 11694 21835 22385 25555 25743 25876 26665 27967 28137
cater[11694,] #dropped one caterpillar, uncertain whether before or after weighing so remove
cater$mpc1[11694] <- NA #removed
cater$mpc2[11694] <- NA #removed
cater[21835,] #only 1 collected so mass is for one caterpillar
cater$mpc2[21835] <- cater$caterpillar.mass2[21835] #same value
cater[22385,] #1 LOST, MASS FOR 4 so divided by 4
cater$mpc1[22385] <- cater$caterpillar.mass2[22385]/4 
cater$mpc2[22385] <- cater$caterpillar.mass2[22385]/4
cater[25555,] #MASS FOR 1 
cater$mpc2[25555] <- cater$caterpillar.mass2[25555]
cater[25743,] #MASS FOR 1 CATERPILLAR
cater$mpc2[25743] <- cater$caterpillar.mass2[25743]
cater[25876,] #MASS FOR 3 CATERPILLARS
cater$mpc1[25876] <- cater$caterpillar.mass2[25876]/3
cater$mpc2[25876] <- cater$caterpillar.mass2[25876]/3
cater[26665,] #mass for 7
cater$mpc1[26665] <- cater$caterpillar.mass2[26665]/7
cater$mpc2[26665] <- cater$caterpillar.mass2[26665]/7
cater[27967,] #MASS FOR 1 CATERPILLAR
cater$mpc1[27967] <- cater$caterpillar.mass2[27967]
cater$mpc2[27967] <- cater$caterpillar.mass2[27967]
cater[28137,] #Mass for 3
cater$mpc1[28137] <- cater$caterpillar.mass2[28137]/3
cater$mpc2[28137] <- cater$caterpillar.mass2[28137]/3

#log the final interval censored values
cater$logmpc1 <- log(cater$mpc1)
cater$logmpc2 <- log(cater$mpc2)


#### Full dataframe ####
Habitat_Site$site <- Habitat_Site$Site
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE)
cater_habitat$treeID <- paste(cater_habitat$tree, cater_habitat$site)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$date, cater_habitat$year)
cater_habitat$siteyear <- paste(cater_habitat$site, cater_habitat$year)
cater_habitat$yearday <- paste(cater_habitat$date, cater_habitat$year)
cater_habitat$weight <- as.numeric(revalue(as.character(cater_habitat$caterpillars), c("0"="1")))
cater_habitat$datescaled <- scale(cater_habitat$date) #unscale(x, center= 146.4095, scale=14.19835)   now from -2.071330 to 2.013653
mean(cater_habitat$date) # 146.4095
sd(cater_habitat$date) # 14.19835

# remove other deciduous samples, just analysing focal taxa
cater_habitat<- subset(cater_habitat, tree.species!="OthDecid")
