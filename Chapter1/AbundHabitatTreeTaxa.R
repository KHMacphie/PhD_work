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

#####################################################
#### Model: Habitat and TreeTaxa Abundance model ####
#####################################################

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


#TTHA_final1<- MCMCglmm(caterpillars~1, 
#                     random=~tree.species+idv(~Alder_Scld+Ash_Scld+Beech_Scld+Birch_Scld+Elm_Scld+Hazel_Scld+Oak_Scld+Rowan_Scld+Sycamore_Scld+Willow_Scld+Conifer_Scld+OthDecid_Scld)+site+year+siteyear+treeID+siteday+recorder, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=3000000, burnin=30000, pr=TRUE, thin=50)
#save(TTHA_final1, file = "~/Documents/Models/TTHA_final1.RData")
load("~/Documents/Models/TTHA_final1.RData")

TTHA_final1.Sim<-simulate(TTHA_final1,nsim=100)
sum(cater_habitat$caterpillars)
par(mfcol=c(1,1))
hist(apply(TTHA_final1.Sim,2,sum), breaks=50)
abline(v=sum(cater_habitat$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(TTHA_final1.Sim,2,propzero), breaks=50)
abline(v=propzero(cater_habitat$caterpillars), col="red")


#dataframe for coeffs and CIs for beaten tree taxa
TTHA <- TTHA_final1$Sol[,2:11] # crop to just the columns wanted
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
  ggtitle("a) Sampled Tree Taxon")+
  xlab("")+
  ylab("Coefficient")

#dataframe for coeffs and CIs for habitat FS's
TTHA2 <- TTHA_final1$Sol[,12:23] # crop to just the columns wanted
TTHA2.df <- data.frame(treetaxa=c(colnames(TTHA2))) #dataframe with column for beaten tree taxa
TTHA2.df$coeff <- apply(TTHA2,2, mean) # mean 
for(i in 1:length(TTHA2.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA2[,i])
  TTHA2.df$lowci[i] <- A["var1","lower"] 
  TTHA2.df$upci[i] <- A["var1","upper"] 
} 
TTHA2.df$treetaxa <- gsub("_Scld.NA.1","", TTHA2.df$treetaxa)
TTHA2.df$treetaxa <- gsub("OthDecid","Other", TTHA2.df$treetaxa)

#plot for habitat FS's
plot2 <- ggplot(TTHA2.df, aes(fct_inorder(treetaxa), coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=25),axis.text.x= element_text(angle=90))+
  geom_hline(yintercept=0, linetype="dashed", colour="red")+
  ggtitle("b) Habitat composition")+
  xlab("")+
  ylab("Coefficient")

row1 <- grid.arrange(plot1, ncol = 1, widths = c(0.5))
row2 <- grid.arrange(plot2, ncol = 1, widths = c(1))
TTHAcoeffs <- grid.arrange(row1, row2, nrow = 2, heights = c(1,1))

############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed
fixed<-rbind(
  c("Intercept",paste(round(mean(TTHA_final1$Sol[,1]),3)," (",
                      round(HPDinterval(TTHA_final1$Sol[,1])[1],3)," - ",
                      round(HPDinterval(TTHA_final1$Sol[,1])[2],3),")",sep=""),round(effectiveSize(TTHA_final1$Sol[,1]))))

####random 
column<-1
treetaxa<-c("Sampled Tree Taxa",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA_final1$VCV[, column])))

column<-2
habitat<-c("Habitat Composition",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA_final1$VCV[, column])))

column<-3
site<-c("Site",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                              round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                              round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA_final1$VCV[, column])))

column<-4
year<-c("Year",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                            round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                            round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(TTHA_final1$VCV[, column])))

column<-5
siteyear<-c("Site Year",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA_final1$VCV[, column])))

column<-6
treeID<-c("Tree ID",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                          round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                          round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(TTHA_final1$VCV[, column])))

column<-7
siteday<-c("Site Day",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                              round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                              round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA_final1$VCV[, column])))

column<-8
recorder<-c("Recorder",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA_final1$VCV[, column])))

column<-9
residual<-c("Residual",paste(round(posterior.mode(TTHA_final1$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA_final1$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA_final1$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA_final1$VCV[, column])))




random<-rbind(treetaxa, habitat, site, year, siteyear, treeID, siteday,recorder,residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableTTHA_final1.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
