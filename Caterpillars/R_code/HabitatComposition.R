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
library(gridExtra)

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

#write.table(Habitat_Site,"Dropbox/master_data/site/Habitat_Site.csv", sep=",", col.names=T, row.names=F)
#write.table(Habitat_Site,"Dropbox/Ross Judge/data/Habitat_Site.csv", sep=",", col.names=T, row.names=F)
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
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")
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
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")

#prop by elevation- some sites have same elevation so not quite right
Props_long_siteinfo <- Props_long_siteinfo[order(Props_long_siteinfo$Mean.Elev),] 
ggplot(Props_long_siteinfo, aes(fct_inorder(Site), Proportion))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")

#Make mean foliage scores long
Habitat_FS <- Habitat_Site[,1:12]
Habitat_FS_long <- gather(Habitat_FS, key="Tree", value="FS", select=2:12)
#plot mean foliage scores of each tree category
ggplot(Habitat_FS_long, aes(Site, FS))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")
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
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")

#FS graph by elevation
FS_long_siteinfo <- FS_long_siteinfo[order(FS_long_siteinfo$Mean.Elev),] 
ggplot(FS_long_siteinfo, aes(fct_inorder(Site), FS))+
  geom_bar(aes(fill=Tree), stat="identity")+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  scale_fill_brewer(palette="Spectral")+
  xlab("Site")

 



################################################
#### Putting Aspen, Beech + Elm in OthDecid ####
################################################
Habitat_Site_condensed <- Habitat_Site
Habitat_Site_condensed$OthDecid_prop <- Habitat_Site_condensed$Elm_prop+Habitat_Site_condensed$Aspen_prop+Habitat_Site_condensed$OthDecid_prop+Habitat_Site_condensed$Beech_prop
Habitat_Site_condensed$OthDecid_FS <- Habitat_Site_condensed$Elm_FS+Habitat_Site_condensed$Aspen_FS+Habitat_Site_condensed$OthDecid_FS+Habitat_Site_condensed$Beech_FS

################################################
#### Multi-membership with tree proportions ####
################################################

# Habitat_Site columns 14-24 for tree proportions
cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
Habitat_Site$site <- Habitat_Site$Site
pmatch(cater$site,Habitat_Site$site,duplicates.ok=TRUE)
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE)
cater_habitat$sitetree <- paste(cater_habitat$tree, cater_habitat$site)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$date, cater_habitat$year)

#Habitat_Site_condensed$site <- Habitat_Site_condensed$Site
#pmatch(cater$site,Habitat_Site_condensed$site,duplicates.ok=TRUE)
#cater_habitat_condensed<- merge(cater, Habitat_Site_condensed, by="site", duplicates.ok=TRUE)
#cater_habitat_condensed$sitetree <- paste(cater_habitat_condensed$tree, cater_habitat_condensed$site)
#cater_habitat_condensed$siteday <- paste(cater_habitat_condensed$site, cater_habitat_condensed$date, cater_habitat_condensed$year)

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

##########################################################################################
#### Redo model with no fixed effects for proportion and FS not condensing otherdecid ####
##########################################################################################
# Model priors
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#MMpropNoFixed<- MCMCglmm(caterpillars~1, 
#                       random=~site+sitetree+siteday+tree.species+idv(~Alder_prop+Ash_prop+Aspen_prop+Beech_prop+Birch_prop+Elm_prop+Oak_prop+Sycamore_prop+Willow_prop+Conifer_prop+OthDecid_prop), 
#                       family="poisson", data=cater_habitat, prior=prior, nitt=200000, burnin=20000, pr=TRUE)
#save(MMpropNoFixed, file = "~/Documents/Models/MMpropNoFixed.RData")
load("~/Documents/Models/MMpropNoFixed.RData")

colnames(MMpropNoFixed$Sol[,5039:5067]) # tree sp= 5039:5056, prop= 5057:5067


#MMFSNoFixed<- MCMCglmm(caterpillars~1, 
#                         random=~site+sitetree+siteday+tree.species+idv(~Alder_FS+Ash_FS+Aspen_FS+Beech_FS+Birch_FS+Elm_FS+Oak_FS+Sycamore_FS+Willow_FS+Conifer_FS+OthDecid_FS), 
#                         family="poisson", data=cater_habitat, prior=prior, nitt=200000, burnin=20000, pr=TRUE)
#save(MMFSNoFixed, file = "~/Documents/Models/MMFSNoFixed.RData")
load("~/Documents/Models/MMFSNoFixed.RData")

#######################
#### MMpropNoFixed ####
#######################

# Dataframe for coeffs and CIs for tree species
MMproptreesp <- MMpropNoFixed$Sol[,5039:5056] # crop to just the columns wanted
MMproptreesp.df <- data.frame(treesp=c(colnames(MMproptreesp))) #dataframe with column for yearsite 
MMproptreesp.df$coeff <- apply(MMproptreesp,2, mean) # mean 
for(i in 1:length(MMproptreesp.df$treesp)) {   # loop for CIs
  A <- HPDinterval(MMproptreesp[,i])
  MMproptreesp.df$lowci[i] <- A["var1","lower"] 
  MMproptreesp.df$upci[i] <- A["var1","upper"] 
} 
MMproptreesp.df$treesp <- gsub("tree.species.","", MMproptreesp.df$treesp)

#dataframe for coeffs and CIs for tree proportions
MMproptreeprops <- MMpropNoFixed$Sol[,5057:5067] # crop to just the columns wanted
MMproptreeprops.df <- data.frame(treesp=c(colnames(MMproptreeprops))) #dataframe with column for yearsite 
MMproptreeprops.df$coeff <- apply(MMproptreeprops,2, mean) # mean 
for(i in 1:length(MMproptreeprops.df$treesp)) {   # loop for CIs
  A <- HPDinterval(MMproptreeprops[,i])
  MMproptreeprops.df$lowci[i] <- A["var1","lower"] 
  MMproptreeprops.df$upci[i] <- A["var1","upper"] 
} 
MMproptreeprops.df$treesp <- gsub("_prop.NA.1","", MMproptreeprops.df$treesp)

#plot tree species
MMpropBeaten <- ggplot(MMproptreesp.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  ggtitle("Tree Beaten, Prop")

#plot tree proportion
MMpropHabitat <- ggplot(MMproptreeprops.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  ggtitle("Habitat Composition, Prop")

MMproprow1 <- grid.arrange(MMpropHabitat)
MMproprow2 <- grid.arrange(MMpropBeaten)
MMprop.panel <- grid.arrange(MMproprow1, MMproprow2, nrow = 2, heights = c(1, 1))

#####################
#### MMFSNoFixed ####
#####################

# Dataframe for coeffs and CIs for tree species
MMFStreesp <- MMFSNoFixed$Sol[,5039:5056] # crop to just the columns wanted
MMFStreesp.df <- data.frame(treesp=c(colnames(MMFStreesp))) #dataframe with column for yearsite 
MMFStreesp.df$coeff <- apply(MMFStreesp,2, mean) # mean 
for(i in 1:length(MMFStreesp.df$treesp)) {   # loop for CIs
  A <- HPDinterval(MMFStreesp[,i])
  MMFStreesp.df$lowci[i] <- A["var1","lower"] 
  MMFStreesp.df$upci[i] <- A["var1","upper"] 
} 
MMFStreesp.df$treesp <- gsub("tree.species.","", MMFStreesp.df$treesp)

#dataframe for coeffs and CIs for tree proportions
MMFStreeprops <- MMFSNoFixed$Sol[,5057:5067] # crop to just the columns wanted
MMFStreeprops.df <- data.frame(treesp=c(colnames(MMFStreeprops))) #dataframe with column for yearsite 
MMFStreeprops.df$coeff <- apply(MMFStreeprops,2, mean) # mean 
for(i in 1:length(MMFStreeprops.df$treesp)) {   # loop for CIs
  A <- HPDinterval(MMFStreeprops[,i])
  MMFStreeprops.df$lowci[i] <- A["var1","lower"] 
  MMFStreeprops.df$upci[i] <- A["var1","upper"] 
} 
MMFStreeprops.df$treesp <- gsub("_FS.NA.1","", MMFStreeprops.df$treesp)

#plot tree species
MMFSBeaten <- ggplot(MMFStreesp.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  ggtitle("Tree Beaten, FS")

#plot tree proportion
MMFSHabitat <- ggplot(MMFStreeprops.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  ggtitle("Habitat Composition, FS")

MMFSrow1 <- grid.arrange(MMFSHabitat)
MMFSrow2 <- grid.arrange(MMFSBeaten)
MMFS.panel <- grid.arrange(MMFSrow1, MMFSrow2, nrow = 2, heights = c(1, 1))

par(mfcol=c(1,1), cex=1.5)
hist(MMFSNoFixed$VCV[,5], breaks=40)
abline(v=0, col=2,type="l", lty=2)
title(ylab="Frequency", outer=TRUE, line = 2)
title( xlab="Variance", outer=TRUE, line = 0)
legend("topright", legend="A", bty="n") 

ggplot(MMFStreeprops.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  xlab("Tree Taxa")+
  geom_hline(yintercept=0, linetype="dashed", color = "red")

## all 4 in one grid
MMrow1 <- grid.arrange(MMpropHabitat, MMFSHabitat, ncol = 2, widths = c(1, 1))
MMrow2 <- grid.arrange(MMpropBeaten, MMFSBeaten, ncol = 2, widths = c(1, 1))
MM.panel <- grid.arrange(MMrow1, MMrow2, nrow = 2, heights = c(1, 1))

# just habitat
MMrow1 <- grid.arrange(MMpropHabitat, MMFSHabitat, ncol = 1, widths = c(1, 1))
MM.panel.hab <- grid.arrange(MMrow1, nrow = 1)

#####################################
#### Tree props as fixed effects ####
#####################################
k<-10000
prior2<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#HabitatFixed<- MCMCglmm(caterpillars~Alder_prop+Ash_prop+Aspen_prop+Beech_prop+Birch_prop+Elm_prop+Oak_prop+Sycamore_prop+Willow_prop+Conifer_prop+OthDecid_prop, 
#                       random=~site+sitetree+siteday+tree.species, 
#                       family="poisson", data=cater_habitat, prior=prior2, nitt=300000, burnin=30000, pr=TRUE)
#save(HabitatFixed, file = "~/Documents/Models/HabitatFixed.RData")
load("~/Documents/Models/HabitatFixed.RData")

#Warning message:
#  In MCMCglmm(caterpillars ~ Alder_prop + Ash_prop + Aspen_prop +  :
#                some fixed effects are not estimable and have been removed. Use singular.ok=TRUE to sample these effects, but use an informative prior!

HabFixedProp <- HabitatFixed$Sol[,2:11] # crop to just the columns wanted
HabFixedProp.df <- data.frame(treesp=c(colnames(HabFixedProp))) #dataframe with column for yearsite 
HabFixedProp.df$coeff <- apply(HabFixedProp,2, mean) # mean 
for(i in 1:length(HabFixedProp.df$treesp)) {   # loop for CIs
  A <- HPDinterval(HabFixedProp[,i])
  HabFixedProp.df$lowci[i] <- A["var1","lower"] 
  HabFixedProp.df$upci[i] <- A["var1","upper"] 
} 
HabFixedProp.df$treesp <- gsub("_prop.","", HabFixedProp.df$treesp)

ggplot(HabFixedProp.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))+
  ggtitle("Habitat Composition Fixed, Prop")


######################################################
#### Random look at possible site width variation ####
######################################################

cater$sitetree <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$datecent <-cater$date-mean(cater$date)
cater$year <- as.factor(cater$year)

k<-10000
prior<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#SiteWidth<- MCMCglmm(caterpillars~datecent*year+datecent*site+I(datecent^2)*site, 
#                       random=~sitetree+siteday, 
#                       family="poisson", data=cater, prior=prior, nitt=300000, burnin=30000)
#save(SiteWidth, file = "~/Documents/Models/SiteWidth.RData")
load("~/Documents/Models/SiteWidth.RData")

summary(SiteWidth)
plot(SiteWidth$Sol)
site <- read.csv("Dropbox/master_data/site/site_details.csv")
site <- site[order(site$site),] 
SiteCurves <- data.frame(site=site$site)
siteint <- SiteWidth$Sol[,1:49]
colnames(siteint)
siteint[,2:6]= NULL

### im just wasting time...

######################################
#### Multi-memb FS with 2019 data ####
######################################

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MMFS_2019<- MCMCglmm(caterpillars~1, 
                         random=~tree.species+idv(~Alder_FS+Ash_FS+Aspen_FS+Beech_FS+Birch_FS+Elm_FS+Oak_FS+Sycamore_FS+Willow_FS+Conifer_FS+OthDecid_FS)+site+sitetree+siteday+recorder, 
                         family="poisson", data=cater_habitat, prior=prior, nitt=300000, burnin=30000, pr=TRUE)
#save(MMFS_2019, file = "~/Documents/Models/MMFS_2019.RData")
load("~/Documents/Models/MMFS_2019.RData")

MMPr_2019<- MCMCglmm(caterpillars~1, 
                     random=~tree.species+idv(~Alder_prop+Ash_prop+Aspen_prop+Beech_prop+Birch_prop+Elm_prop+Oak_prop+Sycamore_prop+Willow_prop+Conifer_prop+OthDecid_prop)+site+sitetree+siteday+recorder, 
                     family="poisson", data=cater_habitat, prior=prior, nitt=300000, burnin=30000, pr=TRUE)
#save(MMPr_2019, file = "~/Documents/Models/MMPr_2019.RData")
load("~/Documents/Models/MMPr_2019.RData")

#check if model generates sensible results
MMFS.Sim<-simulate(MMFS_2019,nsim=100)
sum(cater_habitat$caterpillars)
par(mfcol=c(1,1))
hist(apply(MMFS.Sim,2,sum), breaks=50)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(MMFS.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")

MMPr.Sim<-simulate(MMPr_2019,nsim=100)
sum(cater_habitat$caterpillars)
par(mfcol=c(1,1))
hist(apply(MMPr.Sim,2,sum))
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(MMPr.Sim,2,propzero))
abline(v=propzero(cater_habitat$caterpillars), col="red")

#dataframe for coeffs and CIs for tree proportions
MMFS19 <- MMFS_2019$Sol[,21:31] # crop to just the columns wanted
MMFS19.df <- data.frame(treetaxa=c(colnames(MMFS19))) #dataframe with column for yearsite 
MMFS19.df$coeff <- apply(MMFS19,2, mean) # mean 
for(i in 1:length(MMFS19.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(MMFS19[,i])
  MMFS19.df$lowci[i] <- A["var1","lower"] 
  MMFS19.df$upci[i] <- A["var1","upper"] 
} 
MMFS19.df$treetaxa <- gsub("_FS.NA.1","", MMFS19.df$treetaxa)


#plot tree proportion
ggplot(MMFS19.df, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  ggtitle("Habitat Composition, FS")

###############################
#### Habitat for each tree ####
###############################
treeID <- read_csv("Documents/treeID.csv") # nest boxes for each tree
treeID$site <- treeID$tree
treeID$site <- gsub("1B","", treeID$site)
treeID$site <- gsub("2B","", treeID$site)
treeID$site <- gsub("3B","", treeID$site)
treeID$site <- gsub("4B","", treeID$site)
treeID$site <- gsub("5B","", treeID$site)
treeID$site <- gsub("6B","", treeID$site)
treeID$site <- gsub("7B","", treeID$site)
treeID$site <- gsub("9B","", treeID$site)
treeID$site <- gsub("7Z","", treeID$site)
treeID$site <- gsub("9a","", treeID$site)
treeID$site <- gsub("1","", treeID$site)
treeID$site <- gsub("2","", treeID$site)
treeID$site <- gsub("3","", treeID$site)
treeID$site <- gsub("4","", treeID$site)
treeID$site <- gsub("5","", treeID$site)
treeID$site <- gsub("6","", treeID$site)
treeID$site <- gsub("7","", treeID$site)
treeID$site <- gsub("8","", treeID$site)
treeID$site <- gsub("9","", treeID$site)
treeID$site <- gsub("0","", treeID$site)
treeID$site <- gsub(" ","", treeID$site) # got site name

treeID$siteNB.1 <- paste(treeID$site, treeID$nb.1)
treeID$siteNB.2 <- paste(treeID$site, treeID$nb.2)
Habitat_byNB$siteNB <- paste(Habitat_byNB$Site, Habitat_byNB$NB)
HabitatFS_NB <- Habitat_byNB[,34:45]
treeID_long <- gather(treeID, key="sample", value="siteNB", select=5:6) #gives each tree two rows, one for each associated nest box
treeID_long<- merge(treeID_long, HabitatFS_NB, by="siteNB", duplicates.ok=TRUE) #habitat score for each box assigned

habitat_treeID <- treeID_long
habitat_treeID$siteNB <- NULL
habitat_treeID[,2:5] <- NULL

habitat_treeID <- aggregate(.~tree, habitat_treeID, mean) #mean habitat of 2 boxes (normally the same one)
habitat_treeID$sitetree <- habitat_treeID$tree
habitat_treeID$tree <- NULL

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$sitetree <- paste(cater$tree, cater$site)
cater$siteday <- paste(cater$site, cater$date, cater$year)

# match habitat for each tree to the beating data  with NAs
store <- pmatch(cater$sitetree, habitat_treeID$sitetree, duplicates.ok=TRUE)
cater_habitat_treeID <- cbind(cater,habitat_treeID[,1][store],habitat_treeID[,2][store],habitat_treeID[,3][store],habitat_treeID[,4][store],habitat_treeID[,5][store],habitat_treeID[,6][store],habitat_treeID[,7][store],habitat_treeID[,8][store],habitat_treeID[,9][store],habitat_treeID[,10][store],habitat_treeID[,11][store]) 
names(cater_habitat_treeID)[c(22:32)] <- c("Alder_FS", "Ash_FS", "Aspen_FS", "Beech_FS", "Birch_FS", "Elm_FS", "Oak_FS", "Sycamore_FS", "Willow_FS", "Conifer_FS", "OthDecid_FS")

# without NAs
caterhab_noNAs <- merge(cater, habitat_treeID, by="sitetree", duplicates.ok=TRUE )

#### FS by individual tree model ####

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


MMFS_bytree<- MCMCglmm(caterpillars~1, 
                     random=~tree.species+idv(~Alder_FS+Ash_FS+Aspen_FS+Beech_FS+Birch_FS+Elm_FS+Oak_FS+Sycamore_FS+Willow_FS+Conifer_FS+OthDecid_FS)+site+sitetree+siteday+recorder, 
                     family="poisson", data=caterhab_noNAs, prior=prior, nitt=300000, burnin=30000, pr=TRUE)
#save(MMFS_bytree, file = "~/Documents/Models/MMFS_bytree.RData")
load("~/Documents/Models/MMFS_bytree.RData")

k<-10000
prior1<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

MMFS_maintaxafixed<- MCMCglmm(caterpillars~Alder_FS+Ash_FS+Beech_FS+Birch_FS+Oak_FS+Willow_FS, 
                       random=~tree.species+site+sitetree+siteday+recorder, 
                       family="poisson", data=cater_habitat, prior=prior1, nitt=250000, burnin=25000, pr=TRUE)
save(MMFS_maintaxafixed, file = "~/Documents/Models/MMFS_maintaxafixed.RData")
load("~/Documents/Models/MMFS_maintaxafixed.RData")

MMFS_maintaxafixed.Sim<-simulate(MMFS_maintaxafixed,nsim=100)
sum(cater$caterpillars)
par(mfcol=c(1,1))
hist(apply(MMFS_maintaxafixed.Sim,2,sum), breaks=50)
abline(v=sum(cater$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(MMFS_maintaxafixed.Sim,2,propzero))
abline(v=propzero(cater$caterpillars), col="red")

### only the taxa that differ from average in tree taxa in multi memb
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

MMFS_maintaxa<- MCMCglmm(caterpillars~1, 
                     random=~tree.species+idv(~Alder_FS+Ash_FS+Beech_FS+Birch_FS+Oak_FS+Willow_FS)+site+sitetree+siteday+recorder, 
                     family="poisson", data=cater_habitat, prior=prior, nitt=250000, burnin=25000, pr=TRUE)
save(MMFS_maintaxa, file = "~/Documents/Models/MMFS_maintaxa.RData")
load("~/Documents/Models/MMFS_maintaxa.RData")

#####################################
#### tree diversity of each site ####
#####################################
habitatdiv <- read.csv("~/Dropbox/2018/Habitats_2018_all.csv")
# combining other deciduous
#habitat$X6_oth.decid <- rowSums(habitat[,c("X6_cherry", "X6_elder", "X6_hazel", "X6_holly", "X6_rowan", "X6_rose")], na.rm=TRUE)
#habitat$X21_oth.decid <- rowSums(habitat[,c("X21_blackthorn", "X21_hazel", "X21_holly")], na.rm=TRUE)
#habitat$s_oth.decid <- rowSums(habitat[,c("s_cherry", "s_chestnut", "s_elder", "s_hawthorn", "s_hazel", "s_holly", "s_lime", "s_rowan", "s_whitebeam", "s_other")], na.rm=TRUE)
#habitat$m_oth.decid <- rowSums(habitat[,c("m_cherry", "m_chestnut", "m_holly", "m_lime", "m_rowan", "m_whitebeam", "m_horsechestnut")], na.rm=TRUE)
#habitat$l_oth.decid <- habitat$l_cherry
#habitat$z_oth.decid <- rowSums(habitat[,c("z_blackthorn", "z_cherry", "z_other")], na.rm=TRUE)

# combining conifers
#habitat$X6_conifer_all <- habitat$X6_juniper
#habitat$s_conifer_all <- rowSums(habitat[,c("s_conifer", "s_pine", "s_yew", "s_juniper")], na.rm=TRUE)
#habitat$m_conifer_all <- rowSums(habitat[,c("m_conifer", "m_pine", "m_yew")], na.rm=TRUE)
#habitat$l_conifer_all <- rowSums(habitat[,c("l_conifer", "l_pine")], na.rm=TRUE)

# converting 6 stands into small trees
habitatdiv$alder6smt <- habitatdiv$X6_alder*0.5
habitatdiv$ash6smt <- habitatdiv$X6_ash*0.5
habitatdiv$birch6smt <- habitatdiv$X6_birch*0.5
habitatdiv$beech6smt <- habitatdiv$X6_beech*0.5
habitatdiv$syc6smt <- habitatdiv$X6_sycamore*0.5
habitatdiv$willow6smt <- habitatdiv$X6_willow*0.5
habitatdiv$cherry6smt <- habitatdiv$X6_cherry*0.5
habitatdiv$elder6smt <- habitatdiv$X6_elder*0.5
habitatdiv$hazel6smt <- habitatdiv$X6_hazel*0.5
habitatdiv$holly6smt <- habitatdiv$X6_holly*0.5
habitatdiv$rowan6smt <- habitatdiv$X6_rowan*0.5
habitatdiv$rose6smt <- habitatdiv$X6_rose*0.5
habitatdiv$juniper6smt <- habitatdiv$X6_juniper*0.5

# thickets already as small trees

# combining stands, thickets and small trees
habitatdiv$Alder_S <- rowSums(habitatdiv[,c("s_alder", "X21_alder", "alder6smt")], na.rm=TRUE)
habitatdiv$Ash_S <- rowSums(habitatdiv[,c("s_ash", "ash6smt")], na.rm=TRUE)
habitatdiv$Beech_S <- rowSums(habitatdiv[,c("s_beech", "beech6smt")], na.rm=TRUE)
habitatdiv$Birch_S <- rowSums(habitatdiv[,c("s_birch", "X21_birch","birch6smt")], na.rm=TRUE)
habitatdiv$Sycamore_S <- rowSums(habitatdiv[,c("s_sycamore", "syc6smt")], na.rm=TRUE)
habitatdiv$Willow_S <- rowSums(habitatdiv[,c("s_willow", "willow6smt", "X21_willow", "z_willow")], na.rm=TRUE)
habitatdiv$Cherry_S <- rowSums(habitatdiv[,c("cherry6smt", "z_cherry", "s_cherry")], na.rm=TRUE)
habitatdiv$Elder_S <- rowSums(habitatdiv[,c("elder6smt", "s_elder")], na.rm=TRUE)
habitatdiv$Hazel_S <- rowSums(habitatdiv[,c("hazel6smt", "s_hazel", "X21_hazel")], na.rm=TRUE)
habitatdiv$Holly_S <- rowSums(habitatdiv[,c("holly6smt", "s_holly", "X21_holly")], na.rm=TRUE)
habitatdiv$Rowan_S <- rowSums(habitatdiv[,c("rowan6smt", "s_rowan")], na.rm=TRUE)
habitatdiv$Rose_S <- habitatdiv$rose6smt
habitatdiv$Blackthorn_S <- rowSums(habitatdiv[,c("X21_blackthorn", "z_blackthorn")], na.rm=TRUE)
habitatdiv$Chestnut_S <- habitatdiv$s_chestnut
habitatdiv$Hawthorn_S <- habitatdiv$s_hawthorn
habitatdiv$Lime_S <- habitatdiv$s_lime
habitatdiv$Whitebeam_S <- habitatdiv$s_whitebeam
habitatdiv$Other_S <- rowSums(habitatdiv[,c("z_other", "s_other")], na.rm=TRUE)
#habitatdiv$Horsechestnut_S <- rowSums(habitatdiv[,c("")], na.rm=TRUE)
habitatdiv$Conifer_S <- habitatdiv$s_conifer
habitatdiv$Pine_S <- habitatdiv$s_pine
habitatdiv$Yew_S <- habitatdiv$s_yew
habitatdiv$Juniper_S <- rowSums(habitatdiv[,c("juniper6smt", "s_juniper")], na.rm=TRUE)

##### Use the weighting function to calculate foliage per tree type and size at each NB

## making a function for the weighting equation
weighting <- function(x){return(pi*(x/(2*pi))^2)}  ##weighted by cross secitonal area of the trunk in cm^2 (min possible size for each category)
Sm <- weighting(40)
Med <- weighting(100)
Lar <- weighting(250)
S=1
M=Med/Sm
L=Lar/Sm

habitatdiv_byNB <- data.frame(Site=habitatdiv$site, 
                           NB=habitatdiv$nestbox, 
                           Alder_S=(S*habitatdiv$Alder_S),
                           Alder_M=(M*habitatdiv$m_alder),
                           Ash_L=(L*habitatdiv$l_ash),
                           Ash_M=(M*habitatdiv$m_ash),
                           Ash_S=(S*habitatdiv$Ash_S),
                           Aspen_M=(M*habitatdiv$m_aspen),
                           Aspen_S=(S*habitatdiv$s_aspen),
                           Beech_S=(S*habitatdiv$Beech_S),
                           Beech_M=(M*habitatdiv$m_beech),
                           Beech_L=(L*habitatdiv$l_beech),
                           Birch_S=(S*habitatdiv$Birch_S),
                           Birch_M=(M*habitatdiv$m_birch),
                           Birch_L=(L*habitatdiv$l_birch),
                           Elm_S=(S*habitatdiv$s_elm),
                           Elm_M=(M*habitatdiv$m_elm),
                           Elm_L=(L*habitatdiv$l_elm),
                           Oak_S=(S*habitatdiv$s_oak),
                           Oak_M=(M*habitatdiv$m_oak),
                           Oak_L=(L*habitatdiv$l_oak), 
                           Sycamore_S=(S*habitatdiv$Sycamore_S),
                           Sycamore_M=(M*habitatdiv$m_sycamore),
                           Sycamore_L=(L*habitatdiv$l_sycamore),
                           Willow_S=(S*habitatdiv$Willow_S),
                           Willow_M=(M*habitatdiv$m_willow),
                           Willow_L=(L*habitatdiv$l_willow),
                           Conifer_S=(S*habitatdiv$Conifer_S),
                           Conifer_M=(M*habitatdiv$m_conifer),
                           Conifer_L=(L*habitatdiv$l_conifer),
                           Cherry_L=(L*habitatdiv$l_cherry),
                           Cherry_M=(M*habitatdiv$m_cherry),
                           Cherry_S=(S*habitatdiv$Cherry_S),
                           Pine_L=(L*habitatdiv$l_pine),
                           Pine_M=(M*habitatdiv$m_pine),
                           Pine_S=(S*habitatdiv$Pine_S),
                           Elder_S=(S*habitatdiv$Elder_S),
                           Hazel_S=(S*habitatdiv$Hazel_S),
                           Rose_S=(S*habitatdiv$Rose_S),
                           Blackthorn_S=(S*habitatdiv$Blackthorn_S),
                           Hawthorn_S=(S*habitatdiv$Hawthorn_S),
                           Other_S=(S*habitatdiv$Other_S),
                           Juniper_S=(S*habitatdiv$Juniper_S),
                           Holly_M=(M*habitatdiv$m_holly),
                           Holly_S=(S*habitatdiv$Holly_S),
                           Rowan_M=(M*habitatdiv$m_rowan),
                           Rowan_S=(S*habitatdiv$Rowan_S),
                           Chestnut_M=(M*habitatdiv$m_chestnut),
                           Chestnut_S=(S*habitatdiv$Chestnut_S),
                           Lime_M=(M*habitatdiv$m_lime),
                           Lime_S=(S*habitatdiv$Lime_S),
                           Whitebeam_M=(M*habitatdiv$m_whitebeam),
                           Whitebeam_S=(S*habitatdiv$Whitebeam_S),
                           Yew_M=(M*habitatdiv$m_yew),
                           Yew_S=(S*habitatdiv$Yew_S),
                           Horsechestnut_M=(M*habitatdiv$m_horsechestnut))
habitatdiv_byNB[is.na(habitatdiv_byNB)] <- 0

#Foliage score for each tree type at each NB
habitatdiv_byNB$Alder_FS <- rowSums(habitatdiv_byNB[,c("Alder_S", "Alder_M")], na.rm=TRUE)
habitatdiv_byNB$Ash_FS <- rowSums(habitatdiv_byNB[,c("Ash_S", "Ash_M", "Ash_L")], na.rm=TRUE)
habitatdiv_byNB$Aspen_FS <- rowSums(habitatdiv_byNB[,c("Aspen_S", "Aspen_M")], na.rm=TRUE)
habitatdiv_byNB$Beech_FS <- rowSums(habitatdiv_byNB[,c("Beech_S", "Beech_M", "Beech_L")], na.rm=TRUE)
habitatdiv_byNB$Birch_FS <- rowSums(habitatdiv_byNB[,c("Birch_S", "Birch_M", "Birch_L")], na.rm=TRUE)
habitatdiv_byNB$Elm_FS <- rowSums(habitatdiv_byNB[,c("Elm_S", "Elm_M", "Elm_L")], na.rm=TRUE)
habitatdiv_byNB$Oak_FS <- rowSums(habitatdiv_byNB[,c("Oak_S", "Oak_M", "Oak_L")], na.rm=TRUE)
habitatdiv_byNB$Sycamore_FS <- rowSums(habitatdiv_byNB[,c("Sycamore_S", "Sycamore_M", "Sycamore_L")], na.rm=TRUE)
habitatdiv_byNB$Willow_FS <- rowSums(habitatdiv_byNB[,c("Willow_S", "Willow_M", "Willow_L")], na.rm=TRUE)
habitatdiv_byNB$Conifer_FS <- rowSums(habitatdiv_byNB[,c("Conifer_S", "Conifer_M", "Conifer_L")], na.rm=TRUE)
habitatdiv_byNB$Cherry_FS <- rowSums(habitatdiv_byNB[,c("Cherry_S", "Cherry_M", "Cherry_L")], na.rm=TRUE)
habitatdiv_byNB$Pine_FS <- rowSums(habitatdiv_byNB[,c("Pine_S", "Pine_M", "Pine_L")], na.rm=TRUE)
habitatdiv_byNB$Elder_FS <- habitatdiv_byNB$Elder_S
habitatdiv_byNB$Hazel_FS <- habitatdiv_byNB$Hazel_S
habitatdiv_byNB$Rose_FS <- habitatdiv_byNB$Rose_S
habitatdiv_byNB$Blackthorn_FS <- habitatdiv_byNB$Blackthorn_S
habitatdiv_byNB$Hawthorn_FS <- habitatdiv_byNB$Hawthorn_S
habitatdiv_byNB$Other_FS <- habitatdiv_byNB$Other_S
habitatdiv_byNB$Juniper_FS <- habitatdiv_byNB$Juniper_S
habitatdiv_byNB$Horsechestnut_FS <- habitatdiv_byNB$Horsechestnut_S
habitatdiv_byNB$Rowan_FS <- rowSums(habitatdiv_byNB[,c("Rowan_S", "Rowan_M")], na.rm=TRUE)
habitatdiv_byNB$Chestnut_FS <- rowSums(habitatdiv_byNB[,c("Chestnut_S", "Chestnut_M")], na.rm=TRUE)
habitatdiv_byNB$Lime_FS <- rowSums(habitatdiv_byNB[,c("Lime_S", "Lime_M")], na.rm=TRUE)
habitatdiv_byNB$Whitebeam_FS <- rowSums(habitatdiv_byNB[,c("Whitebeam_S", "Whitebeam_M")], na.rm=TRUE)
habitatdiv_byNB$Yew_FS <- rowSums(habitatdiv_byNB[,c("Yew_S", "Yew_M")], na.rm=TRUE)
habitatdiv_byNB$Holly_FS <- rowSums(habitatdiv_byNB[,c("Holly_S", "Holly_M")], na.rm=TRUE)

site <- read.csv("Dropbox/master_data/site/site_details.csv")
HabitatDiv_Site <- habitatdiv_byNB[,57:81]
HabitatDiv_Site$Site <- habitatdiv_byNB$Site
HabitatDiv_Site <- aggregate(.~Site, HabitatDiv_Site, mean)
HabitatDiv_Site$Total_FS <- rowSums(HabitatDiv_Site[2:26])
HabitatDiv_Site$top <-   sum((HabitatDiv_Site[,2]*(HabitatDiv_Site[,2]-1))+
                               (HabitatDiv_Site[,3]*(HabitatDiv_Site[,3]-1))+
                               (HabitatDiv_Site[,4]*(HabitatDiv_Site[,4]-1))+
                               (HabitatDiv_Site[,5]*(HabitatDiv_Site[,5]-1))+
                               (HabitatDiv_Site[,6]*(HabitatDiv_Site[,6]-1))+
                               (HabitatDiv_Site[,7]*(HabitatDiv_Site[,7]-1))+
                               (HabitatDiv_Site[,8]*(HabitatDiv_Site[,8]-1))+
                               (HabitatDiv_Site[,9]*(HabitatDiv_Site[,9]-1))+
                               (HabitatDiv_Site[,10]*(HabitatDiv_Site[,10]-1))+
                               (HabitatDiv_Site[,11]*(HabitatDiv_Site[,11]-1))+
                               (HabitatDiv_Site[,12]*(HabitatDiv_Site[,12]-1))+
                               (HabitatDiv_Site[,13]*(HabitatDiv_Site[,13]-1))+
                               (HabitatDiv_Site[,14]*(HabitatDiv_Site[,14]-1))+
                               (HabitatDiv_Site[,15]*(HabitatDiv_Site[,15]-1))+
                               (HabitatDiv_Site[,16]*(HabitatDiv_Site[,16]-1))+
                               (HabitatDiv_Site[,17]*(HabitatDiv_Site[,17]-1))+
                               (HabitatDiv_Site[,18]*(HabitatDiv_Site[,18]-1))+
                               (HabitatDiv_Site[,19]*(HabitatDiv_Site[,19]-1))+
                               (HabitatDiv_Site[,20]*(HabitatDiv_Site[,20]-1))+
                               (HabitatDiv_Site[,21]*(HabitatDiv_Site[,21]-1))+
                               (HabitatDiv_Site[,22]*(HabitatDiv_Site[,22]-1))+
                               (HabitatDiv_Site[,23]*(HabitatDiv_Site[,23]-1))+
                               (HabitatDiv_Site[,24]*(HabitatDiv_Site[,24]-1))+
                               (HabitatDiv_Site[,25]*(HabitatDiv_Site[,25]-1))+
                               (HabitatDiv_Site[,26]*(HabitatDiv_Site[,26]-1)))
  
HabitatDiv_Site$bottom <- HabitatDiv_Site$Total_FS*(HabitatDiv_Site$Total_FS-1)
HabitatDiv_Site$SDI <- HabitatDiv_Site$top/HabitatDiv_Site$bottom # site diversity index

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)

store2 <- pmatch(cater$site, HabitatDiv_Site$Site, duplicates.ok=TRUE)
cater_habdiv <- cbind(cater,HabitatDiv_Site[,30][store2]) 
names(cater_habdiv)[20] <- c("SDI")
cater_habdiv$sitetree <- paste(cater_habdiv$site, cater_habdiv$tree)
cater_habdiv$siteday <- paste(cater_habdiv$year, cater_habdiv$site, cater_habdiv$date)
cater_habdiv$datecent <- cater_habdiv$date - mean(cater_habdiv$date)

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
treediv<- MCMCglmm(caterpillars~datecent*year+SDI*datecent+I(datecent^2)*SDI+I(SDI^2), random=~tree.species+us(1+datecent):site+sitetree+siteday+recorder, family="poisson", data=cater_habdiv, prior=prior2, nitt=200000, burnin=20000)
save(treediv, file = "~/Documents/Models/treediv.RData")
load("~/Documents/Models/treediv.RData")

# SDI from 7.3 to 158.7

# plotting curve with changing tree diversity
predday <- seq(-26.4,28.6,0.5)
div10 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*10) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*10^2) +
  mean(treediv$Sol[,"datecent:SDI"]*10)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*10)*predday^2

div30 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*30) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*30^2) +
  mean(treediv$Sol[,"datecent:SDI"]*30)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*30)*predday^2

div50 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*50) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*50^2) +
  mean(treediv$Sol[,"datecent:SDI"]*50)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*50)*predday^2

div70 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*70) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*70^2) +
  mean(treediv$Sol[,"datecent:SDI"]*70)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*70)*predday^2

div90 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*90) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*90^2) +
  mean(treediv$Sol[,"datecent:SDI"]*90)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*90)*predday^2

div110 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*110) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*110^2) +
  mean(treediv$Sol[,"datecent:SDI"]*110)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*110)*predday^2

div130 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*130) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*130^2) +
  mean(treediv$Sol[,"datecent:SDI"]*130)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*130)*predday^2

div150 <- mean(treediv$Sol[,"(Intercept)"]) + 
  mean(treediv$Sol[,"datecent"])*predday + 
  mean(treediv$Sol[,"SDI"]*150) + 
  mean(treediv$Sol[,"I(datecent^2)"])*predday^2 +
  mean(treediv$Sol[,"I(SDI^2)"]*150^2) +
  mean(treediv$Sol[,"datecent:SDI"]*150)*predday +
  mean(treediv$Sol[,"SDI:I(datecent^2)"]*150)*predday^2

predday2 <- seq(120,175,0.5)

plot(predday2, exp(div10), type="l", col="forestgreen", lwd=2)
points(predday2, exp(div30), type="l", col="springgreen4", lwd=2)
points(predday2, exp(div50), type="l", col="green4", lwd=2)
points(predday2, exp(div70), type="l", col="springgreen3", lwd=2)
points(predday2, exp(div90), type="l", col="green3", lwd=2)
points(predday2, exp(div110), type="l", col="limegreen", lwd=2)
points(predday2, exp(div130), type="l", col="green2", lwd=2)
points(predday2, exp(div150), type="l", col="green", lwd=2)

plot(predday2, div10, type="l")#, col=6, lwd=2, ylim=c(0,0.5))
points(predday2, div30, type="l")#, col=5, lwd=2)
points(predday2, div50, type="l")#, col=4, lwd=2)
points(predday2, div70, type="l")#, col=3, lwd=2)
points(predday2, div90, type="l")#, col=2, lwd=2)
points(predday2, div110, type="l")#, col=1, lwd=2)
points(predday2, div130, type="l")#, col=1, lwd=2)
points(predday2, div150, type="l")#, col=1, lwd=2)


##############################
#### MMFS using scaled FS ####
colnames(Habitat_Site)
max(Habitat_Site[,2:12]) #97.125

Habitat_Site[,25:35] <- (Habitat_Site[,2:12]/max(Habitat_Site[,2:12])) # Dividing all FS's by max FS
colnames(Habitat_Site)[25:35] <- c("Alder_Sld","Ash_Sld","Aspen_Sld","Beech_Sld","Birch_Sld","Elm_Sld","Oak_Sld","Sycamore_Sld","Willow_Sld","Conifer_Sld","OthDecid_Sld")

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
Habitat_Site$site <- Habitat_Site$Site
pmatch(cater$site,Habitat_Site$site,duplicates.ok=TRUE)
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE)
cater_habitat$sitetree <- paste(cater_habitat$tree, cater_habitat$site)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$date, cater_habitat$year)

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#MMFS_scaled<- MCMCglmm(caterpillars~1, 
#                     random=~tree.species+idv(~Alder_Sld+Ash_Sld+Aspen_Sld+Beech_Sld+Birch_Sld+Elm_Sld+Oak_Sld+Sycamore_Sld+Willow_Sld+Conifer_Sld+OthDecid_Sld)+site+sitetree+siteday+recorder, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=900000, burnin=40000, pr=TRUE)
#save(MMFS_scaled, file = "~/Documents/Models/MMFS_scaled.RData")
load("~/Documents/Models/MMFS_scaled.RData")
# removed from zero! woooo


########## Use HabitatTreeTaxaCategories code! ###########
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#TreeTaxaHabAbund<- MCMCglmm(caterpillars~1, 
#                     random=~tree.species+idv(~Alder_Scld+Ash_Scld+Beech_Scld+Birch_Scld+Elm_Scld+Hazel_Scld+Oak_Scld+Rowan_Scld+Sycamore_Scld+Willow_Scld+Conifer_Scld+OthDecid_Scld)+site+year+siteyear+sitetree+siteday+recorder, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=300000, burnin=30000, pr=TRUE, thin=50)
#save(TreeTaxaHabAbund, file = "~/Documents/Models/TreeTaxaHabAbund.RData")
load("~/Documents/Models/TreeTaxaHabAbund.RData")

TTHA.Sim<-simulate(TreeTaxaHabAbund,nsim=100)
sum(cater_habitat$caterpillars)
par(mfcol=c(1,1))
hist(apply(TTHA.Sim,2,sum), breaks=50)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(TTHA.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")


#dataframe for coeffs and CIs for beaten tree taxa
TTHA <- TreeTaxaHabAbund$Sol[,2:11] # crop to just the columns wanted
TTHA.df <- data.frame(treetaxa=c(colnames(TTHA))) #dataframe with column for beaten tree taxa
TTHA.df$coeff <- apply(TTHA,2, mean) # mean 
for(i in 1:length(TTHA.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA[,i])
  TTHA.df$lowci[i] <- A["var1","lower"] 
  TTHA.df$upci[i] <- A["var1","upper"] 
} 
TTHA.df$treetaxa <- gsub("tree.species.","", TTHA.df$treetaxa)


#plot beaten tree taxa
plot1 <- ggplot(TTHA.df, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  geom_hline(yintercept=0, linetype="dashed", colour="red")+
  ggtitle("Beaten Tree Taxa")

#dataframe for coeffs and CIs for habitat FS's
TTHA2 <- TreeTaxaHabAbund$Sol[,12:23] # crop to just the columns wanted
TTHA2.df <- data.frame(treetaxa=c(colnames(TTHA2))) #dataframe with column for beaten tree taxa
TTHA2.df$coeff <- apply(TTHA2,2, mean) # mean 
for(i in 1:length(TTHA2.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA2[,i])
  TTHA2.df$lowci[i] <- A["var1","lower"] 
  TTHA2.df$upci[i] <- A["var1","upper"] 
} 
TTHA2.df$treetaxa <- gsub("_Scld.NA.1","", TTHA2.df$treetaxa)


#plot for habitat FS's
plot2 <- ggplot(TTHA2.df, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  geom_hline(yintercept=0, linetype="dashed", colour="red")+
  ggtitle("Habitat composition FS")

row1 <- grid.arrange(plot1, ncol = 1, widths = c(1))
row2 <- grid.arrange(plot2, ncol = 1, widths = c(1))
TTHAcoeffs <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))


#### Using proportion and including total FS in fixed ####
########## Use HabitatTreeTaxaCategories code! ###########
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#TTHPA<- MCMCglmm(caterpillars~Total, 
#                     random=~tree.species+idv(~Alder_prop+Ash_prop+Beech_prop+Birch_prop+Elm_prop+Hazel_prop+Oak_prop+Rowan_prop+Sycamore_prop+Willow_prop+Conifer_prop+OthDecid_prop)+site+year+siteyear+sitetree+siteday+recorder, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=1000000, burnin=50000, pr=TRUE, thin=100)
#save(TTHPA, file = "~/Documents/Models/TTHPA.RData")
load("~/Documents/Models/TTHPA.RData")

TTHPA.Sim<-simulate(TTHPA,nsim=100)
sum(cater_habitat$caterpillars)
par(mfcol=c(1,1))
hist(apply(TTHPA.Sim,2,sum), breaks=50)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(TTHPA.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")


#dataframe for coeffs and CIs for beaten tree taxa
TTHA <- TTHPA$Sol[,3:12] # crop to just the columns wanted
TTHA.df <- data.frame(treetaxa=c(colnames(TTHA))) #dataframe with column for beaten tree taxa
TTHA.df$coeff <- apply(TTHA,2, mean) # mean 
for(i in 1:length(TTHA.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA[,i])
  TTHA.df$lowci[i] <- A["var1","lower"] 
  TTHA.df$upci[i] <- A["var1","upper"] 
} 
TTHA.df$treetaxa <- gsub("tree.species.","", TTHA.df$treetaxa)


#plot beaten tree taxa
plot1 <- ggplot(TTHA.df, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  geom_hline(yintercept=0, linetype="dashed", colour="red")+
  ggtitle("Beaten Tree Taxa")

#dataframe for coeffs and CIs for habitat FS's
TTHA2 <- TTHPA$Sol[,13:24] # crop to just the columns wanted
TTHA2.df <- data.frame(treetaxa=c(colnames(TTHA2))) #dataframe with column for beaten tree taxa
TTHA2.df$coeff <- apply(TTHA2,2, mean) # mean 
for(i in 1:length(TTHA2.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA2[,i])
  TTHA2.df$lowci[i] <- A["var1","lower"] 
  TTHA2.df$upci[i] <- A["var1","upper"] 
} 
TTHA2.df$treetaxa <- gsub("_prop.NA.1","", TTHA2.df$treetaxa)


#plot for habitat FS's
plot2 <- ggplot(TTHA2.df, aes(treetaxa, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  geom_hline(yintercept=0, linetype="dashed", colour="red")+
  ggtitle("Habitat composition prop")

row1 <- grid.arrange(plot1, ncol = 1, widths = c(1))
row2 <- grid.arrange(plot2, ncol = 1, widths = c(1))
TTHAcoeffs <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))
