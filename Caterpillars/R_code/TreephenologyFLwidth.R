rm(list=ls())
library(plyr)
library(dplyr)
library(tidyr)
library(MCMCglmm)

############################
#### Tree phen FL width ####
############################

treephen <- read.csv("~/Dropbox/master_data/trees/Tree_Phenology.csv")
treephen$year <- as.factor(treephen$year)
treephen$SY <- paste(treephen$site, treephen$year)
treephen <- data.frame(SY=treephen$SY, FL=treephen$flf)
treephen$FL <- revalue(treephen$FL, c("<123"="122", "<91"="90", "<92"="91", "XXXXX"=""))
treephen$FL <- as.numeric(as.character(treephen$FL))
treephen <- treephen[complete.cases(treephen), ]
treephenSY <- aggregate(FL ~ SY, data = treephen, max)
treephenSY <- rename(treephenSY, FLmax=FL)
treephenSY$FLmin <- tapply(treephen$FL, treephen$SY, min)
treephenSY$FLwidth <- treephenSY$FLmax-treephenSY$FLmin
SYFL <- data.frame(SY=treephenSY$SY, FLwidth=treephenSY$FLwidth)

cater <- read.csv("~/Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
cater$SY <- paste(cater$site, cater$year)
cater <- merge(cater, SYFL, by="SY", duplicates.ok=TRUE, all.x=TRUE)

cater$datecentred <- cater$date-mean(cater$date)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$sitetree <- paste(cater$site, cater$tree)
cater$year <- as.factor(cater$year)

########################################################################
#### Model with interaction between tree timing width and quadratic ####
########################################################################

k<-10000
prior<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#TreePhenQuad<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2)*FLwidth+I(FLwidth^2), random=~sitetree+site+siteday, family="poisson", data=cater, prior=prior, nitt=300000, burnin=30000)
save(TreePhenQuad, file = "~/Documents/Models/TreePhenQuad.RData")
load("~/Documents/Models/TreePhenQuad.RData")

summary(TreePhenQuad)
plot(TreePhenQuad$Sol) #fixedeffects
