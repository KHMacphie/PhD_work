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
#library(plotly)
library(networkD3)

############################################################
#### Variance attributed to temporal, sites and habitat ####
############################################################

# CatAbundance ~ noFixed, Random= TreeTaxa + TreeID + Site + Year + SiteYear + HabitatFSMM (+ residual)
# Staked bar chart
# Sankey plot

##############################
#### Setting up dataframe ####

#### Habitat data set- Habitat_Site ####
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

#### Habitat FS into beating data ####

cater <- read.csv("Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year)

Habitat_Site$site <- Habitat_Site$Site
pmatch(cater$site,Habitat_Site$site,duplicates.ok=TRUE)
cater_habitat<- merge(cater, Habitat_Site, by="site", duplicates.ok=TRUE) #merging habitat and beating data into one df


#### Completing data set for model ####
cater_habitat$treeID <- paste(cater_habitat$site, cater_habitat$tree)
cater_habitat$siteyear <- paste(cater_habitat$site, cater_habitat$year)
cater_habitat$siteday <- paste(cater_habitat$site, cater_habitat$year, cater_habitat$date)

#### Remove rows with treespecies= ? ####
toBeRemoved <- which(cater_habitat$tree.species=="?")
cater_habitat <- cater_habitat[-toBeRemoved,]

###############################
#### Model and simulations ####

# Model priors
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
AbundVar4<- MCMCglmm(caterpillars~1, 
                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species, 
                     family="poisson", data=cater_habitat, prior=prior, nitt=1500000, burnin=50000)
save(AbundVar4, file = "~/Documents/Models/AbundVar4.RData")
load("~/Documents/Models/AbundVar.RData") #recorder, year, siteyear, site, treeID, treetaxa, MMFS
load("~/Documents/Models/AbundVar2.RData") # added siteday
load("~/Documents/Models/AbundVar3.RData") #MM prop not FS
load("~/Documents/Models/AbundVar4.RData") #no MM
summary(AbundVar4)

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
AbundVarFinal1<- MCMCglmm(caterpillars~1, 
                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species, 
                     family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000)
save(AbundVarFinal1, file = "~/Documents/Models/AbundVarFinal1.RData")
load("~/Documents/Models/AbundVarFinal1.RData")

#check if model generates sensible results
AbundVarFinal1.Sim<-simulate(AbundVarFinal1,nsim=500)
#sum(cater_habitat$caterpillars)
#par(mfcol=c(1,1))
hist(apply(AbundVarFinal1.Sim,2,sum), breaks=100000)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(AbundVarFinal1.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")

# Variance breakdown
VarProp.df <- data.frame(Term=c(colnames(AbundVarFinal1$VCV))) #dataframe with column for random term 

for(i in 1:length(VarProp.df$Term)) {   # loop for CIs
  VarProp.df$Proportion[i] <- mean(AbundVarFinal1$VCV[,i]/rowSums(AbundVarFinal1$VCV))
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

#### Riverplot package Sankey ####
library(riverplot)
?riverplot

edges <- data.frame(
  N1=c("Total Variance ", "Total Variance", "Total Variance", "Spatial", "Spatial", "Spatial", "Temporal", "Temporal", "Temporal", "Other", "Other"),
  N2=c("Spatial", "Temporal", "Other", "Site", "Tree ID", "Tree Taxa", " Site Day", " Site Year", " Year", "Recorder", "Residual"),
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
nodes$y=as.numeric(c(5.1,2.3,6.5,10,0.2,1.7,3.4,5.5,7.4,9.0,10.7,12.5))
rownames(nodes) = nodes$ID

library(RColorBrewer)
library("colorspace")
palette = paste0(rainbow_hcl(18))#, "70")
styles = lapply(nodes$y, function(n) {list(col = palette[n+1], lty = 0, textcol = "black")})
names(styles) = nodes$ID

rp <- list(nodes = nodes, edges = edges, styles = styles)
class(rp) <- c(class(rp), "riverplot")
plot(rp, plot_area = 0.95, yscale=0.06, nodewidth = 3)


############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(AbundVarFinal1$Sol[,1]),3)," (",
                      round(HPDinterval(AbundVarFinal1$Sol[,1])[1],3)," - ",
                      round(HPDinterval(AbundVarFinal1$Sol[,1])[2],3),")",sep=""),round(effectiveSize(AbundVarFinal1$Sol[,1]))))

####random 
column<-1
recorder<-c("Recorder",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                                 round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                                 round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                                 round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-2
year<-c("Year",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-3
siteyear<-c("Site Year",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-4
siteday<-c("Site Day",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-5
site<-c("Site",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-6
treeID<-c("Tree ID",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-7
treetaxa<-c("Tree Taxa",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))

column<-8
residual<-c("Residual",paste(round(posterior.mode(AbundVarFinal1$VCV[, column]),3)," (",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(AbundVarFinal1$VCV[, column])[2],3),")",sep=""),
                             round(effectiveSize(AbundVarFinal1$VCV[, column])))




random<-rbind(site,treeID,treetaxa,siteday,siteyear,year,recorder,residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableAbundVarFinal1.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
