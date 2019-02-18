########################
#### Foliage Scores ####
########################

rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
library(tidyr)
#library(doBy)
library(lme4)
library(MCMCglmm)
library(forcats)

### Habitat data set

habitat <- read.csv("~/Dropbox/2018/Habitats_2018_all.csv")
table(habitat$site) #checking all there

######## New data frame with stands, thickets and small trees combined as small trees and other deciduous combined
#converting 6 stands into small trees and combining oth.decid and conifers (combining pine, yew, juniper and conifer)

# combining other deciduous
habitat$X6_oth.decid <- rowSums(habitat[,c("X6_cherry", "X6_elder", "X6_hazel", "X6_holly", "X6_rowan", "X6_rose")], na.rm=TRUE)
habitat$X21_oth.decid <- rowSums(habitat[,c("X21_blackthorn", "X21_hazel", "X21_holly")], na.rm=TRUE)
habitat$s_oth.decid <- rowSums(habitat[,c("s_cherry", "s_chestnut", "s_elder", "s_hawthorn", "s_hazel", "s_holly", "s_lime", "s_rowan", "s_whitebeam", "s_other")], na.rm=TRUE)
habitat$m_oth.decid <- rowSums(habitat[,c("m_cherry", "m_chestnut", "m_holly", "m_lime", "m_rowan", "m_whitebeam", "m_horsechestnut")], na.rm=TRUE)
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
habitat$syc6smt <- habitat$X6_sycamore*0.5
habitat$willow6smt <- habitat$X6_willow*0.5
habitat$oth.decid6smt <- habitat$X6_oth.decid*0.5
habitat$conifer6smt <- habitat$X6_conifer_all*0.5
# thickets already as small trees

# combining stands, thickets and small trees
habitat$Alder_S <- rowSums(habitat[,c("s_alder", "X21_alder", "alder6smt")], na.rm=TRUE)
habitat$Ash_S <- rowSums(habitat[,c("s_ash", "ash6smt")], na.rm=TRUE)
habitat$Beech_S <- rowSums(habitat[,c("s_beech", "beech6smt")], na.rm=TRUE)
habitat$Birch_S <- rowSums(habitat[,c("s_birch", "X21_birch","birch6smt")], na.rm=TRUE)
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
                           Aspen_M=(M*habitat$m_aspen),
                           Aspen_S=(S*habitat$s_aspen),
                           Beech_S=(S*habitat$Beech_S),
                           Beech_M=(M*habitat$m_beech),
                           Beech_L=(L*habitat$l_beech),
                           Birch_S=(S*habitat$Birch_S),
                           Birch_M=(M*habitat$m_birch),
                           Birch_L=(L*habitat$l_birch),
                           Elm_S=(S*habitat$s_elm),
                           Elm_M=(M*habitat$m_elm),
                           Elm_L=(L*habitat$l_elm),
                           Oak_S=(S*habitat$s_oak),
                           Oak_M=(M*habitat$m_oak),
                           Oak_L=(L*habitat$l_oak), 
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
Habitat_byNB$Aspen_FS <- rowSums(Habitat_byNB[,c("Aspen_S", "Aspen_M")], na.rm=TRUE)
Habitat_byNB$Beech_FS <- rowSums(Habitat_byNB[,c("Beech_S", "Beech_M", "Beech_L")], na.rm=TRUE)
Habitat_byNB$Birch_FS <- rowSums(Habitat_byNB[,c("Birch_S", "Birch_M", "Birch_L")], na.rm=TRUE)
Habitat_byNB$Elm_FS <- rowSums(Habitat_byNB[,c("Elm_S", "Elm_M", "Elm_L")], na.rm=TRUE)
Habitat_byNB$Oak_FS <- rowSums(Habitat_byNB[,c("Oak_S", "Oak_M", "Oak_L")], na.rm=TRUE)
Habitat_byNB$Sycamore_FS <- rowSums(Habitat_byNB[,c("Sycamore_S", "Sycamore_M", "Sycamore_L")], na.rm=TRUE)
Habitat_byNB$Willow_FS <- rowSums(Habitat_byNB[,c("Willow_S", "Willow_M", "Willow_L")], na.rm=TRUE)
Habitat_byNB$Conifer_FS <- rowSums(Habitat_byNB[,c("Conifer_S", "Conifer_M", "Conifer_L")], na.rm=TRUE)
Habitat_byNB$OthDecid_FS <- rowSums(Habitat_byNB[,c("OthDecid_S", "OthDecid_M", "OthDecid_L")], na.rm=TRUE)

## combining NB within each site- mean to account for different number of NBs at sites
site <- read.csv("Dropbox/master_data/site/site_details.csv")
Habitat_Site <- Habitat_byNB[,34:44]
Habitat_Site$Site <- Habitat_byNB$Site
Habitat_Site <- aggregate(.~Site, Habitat_Site, mean)

# getting proportions of each tree category
Habitat_Site$Total <- rowSums(Habitat_Site[2:12])
Habitat_Site$Alder_prop <- Habitat_Site$Alder_FS/Habitat_Site$Total
Habitat_Site$Ash_prop <- Habitat_Site$Ash_FS/Habitat_Site$Total
Habitat_Site$Aspen_prop <- Habitat_Site$Aspen_FS/Habitat_Site$Total
Habitat_Site$Beech_prop <- Habitat_Site$Beech_FS/Habitat_Site$Total
Habitat_Site$Birch_prop <- Habitat_Site$Birch_FS/Habitat_Site$Total
Habitat_Site$Elm_prop <- Habitat_Site$Elm_FS/Habitat_Site$Total
Habitat_Site$Oak_prop <- Habitat_Site$Oak_FS/Habitat_Site$Total
Habitat_Site$Sycamore_prop <- Habitat_Site$Sycamore_FS/Habitat_Site$Total
Habitat_Site$Willow_prop <- Habitat_Site$Willow_FS/Habitat_Site$Total
Habitat_Site$Conifer_prop <- Habitat_Site$Conifer_FS/Habitat_Site$Total
Habitat_Site$OthDecid_prop <- Habitat_Site$OthDecid_FS/Habitat_Site$Total

#checking its correct
#Habitat_Site$propadd <- rowSums(Habitat_Site_mean[,14:24]) # it is

#######################################
#### Graphs of Habitat Composition ####
#######################################

#plot of total foliage by site
# from mean
ggplot(Habitat_Site, aes(Site, Total))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))
      
#Make proportions long
Habitat_props <- Habitat_Site[,14:24]
Habitat_props$Site <- Habitat_Site$Site
Habitat_props_long <- gather(Habitat_props, key="Tree", value="Proportion", select=1:11)

#plot proportions of each tree category
ggplot(Habitat_props_long, aes(Site, Proportion))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")
  #scale_fill_brewer(palette="Set3")

#prop by latitude
pmatch(site$Site, Habitat_props_long$Site)
site$Site <- site$site
Props_long_siteinfo <- merge(Habitat_props_long, site, by="Site", duplicates.ok=TRUE)
Props_long_siteinfo <- Props_long_siteinfo[order(Props_long_siteinfo$Mean.Lat),] 
ggplot(Props_long_siteinfo, aes(fct_inorder(Site), Proportion))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")

#prop by elevation- some sites have same elevation so not quite right
Props_long_siteinfo <- Props_long_siteinfo[order(Props_long_siteinfo$Mean.Elev),] 
ggplot(Props_long_siteinfo, aes(fct_inorder(Site), Proportion))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")

#Make mean foliage scores long
Habitat_FS <- Habitat_Site[,1:12]
Habitat_FS_long <- gather(Habitat_FS, key="Tree", value="FS", select=2:12)
#plot mean foliage scores of each tree category
ggplot(Habitat_FS_long, aes(Site, FS))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")
#scale_fill_brewer(palette="Set3")

#FS graph by latitude
site$Site <- site$site
pmatch(site$Site, Habitat_FS_long$Site)
FS_long_siteinfo <- merge(Habitat_FS_long, site, by="Site", duplicates.ok=TRUE)
FS_long_siteinfo <- FS_long_siteinfo[order(FS_long_siteinfo$Mean.Lat),] 
ggplot(FS_long_siteinfo, aes(fct_inorder(Site), FS))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")

#FS graph by elevation
FS_long_siteinfo <- FS_long_siteinfo[order(FS_long_siteinfo$Mean.Elev),] 
ggplot(FS_long_siteinfo, aes(fct_inorder(Site), FS))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")

 



#########################################
#### Putting Aspen, Beech + Elm in OthDecid ####
#########################################
Habitat_Site_condensed <- Habitat_Site
Habitat_Site_condensed$OthDecid_prop <- Habitat_Site_condensed$Elm_prop+Habitat_Site_condensed$Aspen_prop+Habitat_Site_condensed$OthDecid_prop+Habitat_Site_condensed$Beech_prop
Habitat_Site_condensed$OthDecid_FS <- Habitat_Site_condensed$Elm_FS+Habitat_Site_condensed$Aspen_FS+Habitat_Site_condensed$OthDecid_FS+Habitat_Site_condensed$Beech_FS

################################################
#### Multi-membership with tree proportions ####
################################################

# Habitat_Site columns 14-24 for tree proportions
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
Habitat_Site$site <- Habitat_Site$Site
pmatch(cater$site,Habitat_Site$site,duplicates.ok=TRUE)
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE)
cater_habitat$sitetree <- paste(cater_habitat$tree, cater_habitat$site)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$date, cater_habitat$year)

Habitat_Site_condensed$site <- Habitat_Site_condensed$Site
pmatch(cater$site,Habitat_Site_condensed$site,duplicates.ok=TRUE)
cater_habitat_condensed<- merge(cater, Habitat_Site_condensed, by="site", duplicates.ok=TRUE)
cater_habitat_condensed$sitetree <- paste(cater_habitat_condensed$tree, cater_habitat_condensed$site)
cater_habitat_condensed$siteday <- paste(cater_habitat_condensed$site, cater_habitat_condensed$date, cater_habitat_condensed$year)

# Model priors
#k<-10000
#prior<-list(R=list(V=1,nu=0.002),
#            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

##### Need stronger prior --  no fixed effects
#MultiMembProp<- MCMCglmm(caterpillars~1, 
#                       random=~site+sitetree+siteday+tree.species+idv(~Alder_prop+Ash_prop+Birch_prop+Oak_prop+Sycamore_prop+Willow_prop+Conifer_prop+OthDecid_prop), 
#                       family="poisson", data=cater_habitat_condensed, prior=prior, nitt=250000, burnin=25000, pr=TRUE)
#save(MultiMembProp, file = "~/Documents/Models/MultiMembProp.RData")
load("~/Documents/Models/MultiMembProp.RData")

summary(MultiMembProp)

### finding tree species and multimemb columns
which(colnames(MultiMembProp$Sol)=="tree.species. Elm") #5040   Need to correct Elm as 2 at the moment
which(colnames(MultiMembProp$Sol)=="tree.species.Willow") #5058
which(colnames(MultiMembProp$Sol)=="Alder_prop.NA.1") #5059
which(colnames(MultiMembProp$Sol)=="OthDecid_prop.NA.1") # 5066


#dataframe for coeffs and Cis for tree species
MMtreespREcropped <- MultiMembProp$Sol[,5040:5058] # crop to just the columns wanted
MMtreesp.df <- data.frame(treesp=c(colnames(MMtreespREcropped))) #dataframe with column for yearsite 
MMtreesp.df$coeff <- apply(MMtreespREcropped,2, mean) # mean 
for(i in 1:length(MMtreesp.df$treesp)) {   # loop for CIs
  A <- HPDinterval(MMtreespREcropped[,i])
  MMtreesp.df$lowci[i] <- A["var1","lower"] 
  MMtreesp.df$upci[i] <- A["var1","upper"] 
} 
MMtreesp.df$treesp <- gsub("tree.species.","", MMtreesp.df$treesp)

#dataframe for coeffs and Cis for tree proportions
treepropREcropped <- MultiMembProp$Sol[,5059:5066] # crop to just the columns wanted
treeprop.df <- data.frame(treesp=c(colnames(treepropREcropped))) #dataframe with column for yearsite 
treeprop.df$coeff <- apply(treepropREcropped,2, mean) # mean 
for(i in 1:length(treeprop.df$treesp)) {   # loop for CIs
  A <- HPDinterval(treepropREcropped[,i])
  treeprop.df$lowci[i] <- A["var1","lower"] 
  treeprop.df$upci[i] <- A["var1","upper"] 
} 
treeprop.df$treesp <- gsub("_prop.NA.1","", treeprop.df$treesp)

#plot tree species
ggplot(MMtreesp.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))

#plot tree proportion
ggplot(treeprop.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))
