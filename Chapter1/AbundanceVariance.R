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

###########################################
#### Variance in caterpillar abundance ####
###########################################

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
#cater_habitat$datescaled <- cater_habitat$date/max(cater_habitat$date) ## WRONG
cater_habitat$datescaled <- scale(cater_habitat$date) #unscale(x, center= 146.4095, scale=14.19835)   now from -2.071330 to 2.013653
mean(cater_habitat$date) # 146.4095
sd(cater_habitat$date) # 14.19835
cater_habitat$yearday <- paste(cater_habitat$date, cater_habitat$year)
cater_habitat$weight <- as.numeric(revalue(as.character(cater_habitat$caterpillars), c("0"="1")))

#!!!!!!!!!!! if removing OthDecid !!!!!!!!!!!!
cater_habitat<- subset(cater_habitat, tree.species!="OthDecid")

###############################
#### Model and simulations ####

# Model priors
#k<-10000
#prior<-list(R=list(V=1,nu=0.002),
#            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
#                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

# Model
#AbundVar4<- MCMCglmm(caterpillars~1, 
#                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=1500000, burnin=50000)
#save(AbundVar4, file = "~/Documents/Models/AbundVar4.RData")
#load("~/Documents/Models/AbundVar.RData") #recorder, year, siteyear, site, treeID, treetaxa, MMFS
#load("~/Documents/Models/AbundVar2.RData") # added siteday
#load("~/Documents/Models/AbundVar3.RData") #MM prop not FS
#load("~/Documents/Models/AbundVar4.RData") #no MM
#summary(AbundVar4)



#check if model generates sensible results
AbundVar4.Sim<-simulate(AbundVar4,nsim=500)
#sum(cater_habitat$caterpillars)
#par(mfcol=c(1,1))
hist(apply(AbundVar4.Sim,2,sum), breaks=100000, xlim=c(0,10000))
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AbundVar4.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")

##########################################################
#### Mean variance as proportion for each random term ####

mean(AbundVar$VCV[,1]/rowSums(AbundVar$VCV))

VarProp.df <- data.frame(Term=c(colnames(AbundVarFinal$VCV))) #dataframe with column for random term 

for(i in 1:length(VarProp.df$Term)) {   # loop for CIs
  VarProp.df$Proportion[i] <- mean(AbundVarFinal$VCV[,i]/rowSums(AbundVarFinal$VCV))
  } 
VarProp.df$Term <- gsub("recorder","Recorder", VarProp.df$Term)
VarProp.df$Term <- gsub("siteyear","SiteYear", VarProp.df$Term)
VarProp.df$Term <- gsub("siteday","SiteDay", VarProp.df$Term)
VarProp.df$Term <- gsub("site","Site", VarProp.df$Term)
VarProp.df$Term <- gsub("year","Year", VarProp.df$Term)
VarProp.df$Term <- gsub("treeID","TreeID", VarProp.df$Term)
VarProp.df$Term <- gsub("tree.species","TreeTaxa", VarProp.df$Term)
VarProp.df$Term <- gsub("units","Residual", VarProp.df$Term)
VarProp.df$Percentage <- VarProp.df$Proportion*100
VarProp.df$Percentage <- round(VarProp.df$Percentage, digits=2)
VarProp.df$Data <- "Variance"

#### Figure: stacked barplot ####  Need re ordering and recolouring

ggplot(VarProp.df, aes(x = Data, y = Percentage, fill = Term)) +
  geom_col() +
  geom_text(aes(label = paste0(Percentage, "%")),
            position = position_stack(vjust = 0.5)) +
  #scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 16) +
  ylab("Percentage") +
  xlab(NULL)

#### Sankey diagram ####

links <- data.frame(
  source=c("Total Variance", "Total Variance", "Total Variance", "Spatial", "Spatial", "Spatial", "Spatial", "Temporal", "Temporal", "Temporal", "Other", "Other"),
  target=c("Spatial", "Temporal", "Other", "TreeID", "TreeTaxa", "TreeTaxaFS", "Site", "SiteDay", "SiteYear", "Year", "Recorder", "Residual"),
  value=c(sum(VarProp.df[5:8,3]), sum(VarProp.df[2:4,3]), VarProp.df[1,3]+VarProp.df[9,3], VarProp.df[6,3], VarProp.df[7,3], VarProp.df[8,3], VarProp.df[5,3], VarProp.df[4,3], VarProp.df[3,3], VarProp.df[2,3], VarProp.df[1,3], VarProp.df[9,3])
)
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE)

#### Riverplot package Sankey ####
library(riverplot)
?riverplot

edges <- data.frame(
  N1=c("Total Variance ", "Total Variance", "Total Variance", "Spatial", "Spatial", "Spatial", "Temporal", "Temporal", "Temporal", "Other", "Other"),
  N2=c("Spatial", "Temporal", "Other", "Site", "Tree ID", "Tree Taxa", "Site Day", "Site Year", "Year", "Recorder", "Residual"),
  Value=c(sum(VarProp.df[5:7,3]), sum(VarProp.df[2:4,3]), VarProp.df[1,3]+VarProp.df[8,3], VarProp.df[5,3], VarProp.df[6,3], VarProp.df[7,3], VarProp.df[4,3], VarProp.df[3,3], VarProp.df[2,3], VarProp.df[1,3], VarProp.df[8,3])
)

edges$N2<-paste(edges$N2, '\n',  paste0(edges$Value, '%')) 
edges$N1<-c(rep('Total Variance', 3),
            rep(edges$N2[1], 3),
            rep(edges$N2[2], 3), 
            rep(edges$N2[3], 2)   
)


nodes <- data.frame(
  ID=c(as.character(edges$N1), 
         as.character(edges$N2)) %>% unique()
)

nodes$x=as.integer(c(1,2,2,2,3,3,3,3,3,3,3,3))
nodes$y=as.numeric(c(5.1,2.3,6.5,10,0.2,1.7,3.4,5.5,7.2,8.8,10.7,12.5))
rownames(nodes) = nodes$ID

library(RColorBrewer)
library("colorspace")
palette = paste0(rainbow_hcl(18))#, "70")
styles = lapply(nodes$y, function(n) {list(col = palette[n+1], lty = 0, textcol = "black")})
names(styles) = nodes$ID

rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
plot(rp, plot_area = 0.95, yscale=0.05, nodewidth = 4)


########################
#### Colour palette ####
########################
library("colorspace")
pal <- choose_palette()
bluegreen <- sequential_hcl(15, h=c(250,120), c.= c(65,85), l=c(20,75), power=c(1.5,1.5))

###############################################
#### Model using HabitatTreeTaxaCategories ####
###############################################

# Model priors
k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

# Model
#AbundVarFinal1<- MCMCglmm(caterpillars~1, 
#                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000)
#save(AbundVarFinal1, file = "~/Documents/Models/AbundVarFinal1.RData")
#load("~/Documents/Models/AbundVarFinal1.RData")

#### Inc day and date + date^2 ####
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

# Model
#AbundVariance<- MCMCglmm(caterpillars~datescaled+I(datescaled^2), 
#                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species+yearday, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000)
#save(AbundVariance, file = "~/Documents/Models/AbundVariance.RData")
load("~/Documents/Models/AbundVariance.RData") 

#check if model generates sensible results
AbundVariance.Sim<-simulate(AbundVariance,nsim=500)
#sum(cater_habitat$caterpillars)
#par(mfcol=c(1,1))
hist(apply(AbundVariance.Sim,2,sum), breaks=100)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AbundVariance.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")


#### Variance from fixed effects ####
# have pos as a vector indicating the position of the relevant terms (for example 2:3 if the 2nd and 3rd terms are x and x^2) and then do this:
  
V<-cov(as.matrix(AbundVariance$X[,2:3]))

R2<-apply(AbundVariance$Sol[,2:3], 1, function(x){x%*%V%*%x}) # not actual R2- variance explained by fixed effects

Variances <- data.frame(AbundVariance$VCV)
Variances$Fixed <- R2
# Variance breakdown
VarProp.df <- data.frame(Term=c(colnames(Variances))) #dataframe with column for random term 

for(i in 1:length(VarProp.df$Term)) {   # loop for CIs
  VarProp.df$Proportion[i] <- mean(Variances[,i]/rowSums(Variances))
  A <- HPDinterval(mcmc(Variances[,i]/rowSums(Variances)))
  VarProp.df$ProportionLCI[i] <- A[1]
  VarProp.df$ProportionUCI[i] <- A[2]
} 

VarProp.df$Term <- gsub("recorder","Recorder", VarProp.df$Term)
VarProp.df$Term <- gsub("yearday","Day", VarProp.df$Term)
VarProp.df$Term <- gsub("siteyear","SiteYear", VarProp.df$Term)
VarProp.df$Term <- gsub("siteday","SiteDay", VarProp.df$Term)
VarProp.df$Term <- gsub("site","Site", VarProp.df$Term)
VarProp.df$Term <- gsub("year","Year", VarProp.df$Term)
VarProp.df$Term <- gsub("treeID","TreeID", VarProp.df$Term)
VarProp.df$Term <- gsub("tree.species","TreeTaxa", VarProp.df$Term)
VarProp.df$Term <- gsub("units","Residual", VarProp.df$Term)
VarProp.df$Percentage <- VarProp.df$Proportion*100
VarProp.df$Percentage <- round(VarProp.df$Percentage, digits=2)
VarProp.df$Data <- "Variance"

#### Bar plot of prop variance with CIs ####
ggplot(VarProp.df, aes(Term, Proportion))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=ProportionUCI, ymin=ProportionLCI, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("Term")+
  ylab("Proportion of variance")
  theme_bw()+
  theme(text=element_text(size= 20))


#### Riverplot package Sankey ####
library(riverplot)
?riverplot

SpatialProp <- sum(VarProp.df[5:7,5])
TemporalProp <- sum(VarProp.df[2:4,5])+VarProp.df[8,5]+VarProp.df[10,5]
OtherProp <- VarProp.df[1,5]+VarProp.df[9,5]
RecorderProp <- VarProp.df[1,5]
YearProp <- VarProp.df[2,5]
YearSiteProp <- VarProp.df[3,5]
DaySiteYearProp <- VarProp.df[4,5]
SiteProp <- VarProp.df[5,5]
TreeProp <- VarProp.df[6,5]
TreeTaxonProp <- VarProp.df[7,5]
DayYearProp <- VarProp.df[8,5]
ResidualProp <- VarProp.df[9,5]
FixedProp <- VarProp.df[10,5]
  
edges <- data.frame( 
  N1=   c("Total Variance", "Total Variance", "Total Variance", "Spatial", "Spatial",     "Spatial", "Temporal", "Temporal",  "Temporal",   "Temporal",  "Temporal",      "Other",      "Other"),
  N2=   c("Spatial",        "Temporal",       "Other",          "\nSite",    "Tree\nTaxon",  "\nTree",    "\nYear",     "Year\nSite", "Date+\nDate²", "Day\nYear",  "Day Site\nYear", "\nRecorder",   "\nResidual"),
  Value=c(SpatialProp,      TemporalProp,     OtherProp,        SiteProp,  TreeTaxonProp, TreeProp,   YearProp,  YearSiteProp, FixedProp,   DayYearProp, DaySiteYearProp, RecorderProp, ResidualProp)
)

edges$N2<-paste(edges$N2, '\n',  paste0(edges$Value, '%')) 
edges$N1<-c(rep('Total Variance', 3),
            rep(edges$N2[1], 3),
            rep(edges$N2[2], 5), 
            rep(edges$N2[3], 2)   
)


nodes <- data.frame(
  ID=c(as.character(edges$N1), 
       as.character(edges$N2)) %>% unique()
)

nodes$x=as.integer(c(1,2,2,2,3,3,3,3,3,3,3,3,3,3))
nodes$y=as.numeric(c(7.8,2.3,9.5,15,0,2,4,6.5,8.5,10.5,12.5,14.5,17.5,19.5))
rownames(nodes) = nodes$ID

library(RColorBrewer)
library("colorspace")
palette = paste0(rainbow_hcl(n=22, c=30), "90")
styles = lapply(nodes$y, function(n) {list(col = palette[n+1], lty = 0, textcol = "black")})
names(styles) = nodes$ID

rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
par(cex=0.95)
plot(rp, plot_area = 0.95, yscale=0.08, nodewidth = 4.6) #saved as 10"x5"

#### Plot with spatiotemporal category ####

SpatialProp <- sum(VarProp.df[5:7,5])
SpatiotemporalProp <- sum(VarProp.df[3:4,5])
TemporalProp <- VarProp.df[2,5]+VarProp.df[8,5]+VarProp.df[10,5]
OtherProp <- VarProp.df[1,5]+VarProp.df[9,5]
RecorderProp <- VarProp.df[1,5]
YearProp <- VarProp.df[2,5]
YearSiteProp <- VarProp.df[3,5]
DaySiteYearProp <- VarProp.df[4,5]
SiteProp <- VarProp.df[5,5]
TreeProp <- VarProp.df[6,5]
TreeTaxonProp <- VarProp.df[7,5]
DayYearProp <- VarProp.df[8,5]
ResidualProp <- VarProp.df[9,5]
FixedProp <- VarProp.df[10,5]

edges <- data.frame( 
  N1=   c("Total Variance", "Total Variance",   "Total Variance", "Total Variance", "Spatial", "Spatial",     "Spatial", "Spatio-\ntemporal", "Spatio-\ntemporal", "Temporal", "Temporal",     "Temporal",  "Other",      "Other"),
  N2=   c("Spatial",        "Spatio-\ntemporal", "Temporal",       "Other",          "\nSite",  "Tree\nTaxon", "\nTree",  "Site\nYear",        "Day Site\nYear",    "\nYear",   "Date+\nDate²", "Day\nYear", "\nRecorder", "\nResidual"),
  Value=c(SpatialProp,     SpatiotemporalProp,   TemporalProp,     OtherProp,        SiteProp,  TreeTaxonProp, TreeProp,  YearSiteProp,        DaySiteYearProp,     YearProp,   FixedProp,      DayYearProp, RecorderProp, ResidualProp)
)

edges$N2<-paste(edges$N2, '\n',  paste0(edges$Value, '%')) 
edges$N1<-c(rep('Total Variance', 4),
            rep(edges$N2[1], 3),
            rep(edges$N2[2], 2),
            rep(edges$N2[3], 3), 
            rep(edges$N2[4], 2)   
)


nodes <- data.frame(
  ID=c(as.character(edges$N1), 
       as.character(edges$N2)) %>% unique()
)

nodes$x=as.integer(c(1,2,2,2,2,3,3,3,3,3,3,3,3,3,3))
nodes$y=as.numeric(c(9.5,2.3,7.5,13,18,0,2,4,6.5,8.5,11.5,14,16,18.5,20.5))
rownames(nodes) = nodes$ID

library(RColorBrewer)
library("colorspace")
palette = paste0(rainbow_hcl(n=22, c=30), "95")
styles = lapply(nodes$y, function(n) {list(col = palette[n+1], lty = 0, textcol = "black")})
names(styles) = nodes$ID

rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
par(cex=0.95)
plot(rp, plot_area = 0.95, yscale=0.08, nodewidth = 4.6) #saved as 10"x5"



############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(AbundVariance$Sol[,1]),3)," (",
                      round(HPDinterval(AbundVariance$Sol[,1])[1],3)," - ",
                      round(HPDinterval(AbundVariance$Sol[,1])[2],3),")",sep=""),round(effectiveSize(AbundVariance$Sol[,1]))),
  
  c("Date (scaled)",paste(round(mean(AbundVariance$Sol[,2]),3)," (",
                      round(HPDinterval(AbundVariance$Sol[,2])[1],3)," - ",
                      round(HPDinterval(AbundVariance$Sol[,2])[2],3),")",sep=""),round(effectiveSize(AbundVariance$Sol[,2]))),
  
  c("Date² (scaled)",paste(round(mean(AbundVariance$Sol[,3]),3)," (",
                      round(HPDinterval(AbundVariance$Sol[,3])[1],3)," - ",
                      round(HPDinterval(AbundVariance$Sol[,3])[2],3),")",sep=""),round(effectiveSize(AbundVariance$Sol[,3]))))

####random 
column<-1
recorder<-c("Recorder",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                                 round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                                 round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                                 round(effectiveSize(AbundVariance$VCV[, column])))

column<-2
year<-c("Year",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-3
siteyear<-c("Site Year",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-4
siteday<-c("Site Day",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-5
site<-c("Site",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-6
treeID<-c("Tree ID",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-7
treetaxa<-c("Tree Taxa",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))

column<-8
day<-c("Day",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(AbundVariance$VCV[, column])))

column<-9
residual<-c("Residual",paste(round(posterior.mode(AbundVariance$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVariance$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVariance$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVariance$VCV[, column])))




random<-rbind(site,treeID,treetaxa,siteday,day,siteyear,year,recorder,residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableAbundVariance.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
