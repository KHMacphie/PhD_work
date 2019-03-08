###############################
#### Tree diversity beaten ####
###############################

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

# by site year as trees beaten differs between years
# simpsons diversity index D= 1-((sum n(n-1))/N(N-1))

cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
cater$year <- as.factor(cater$year)

trees <- data.frame(year=cater$year, site=cater$site, tree=cater$tree, treesp=cater$`tree species`)
trees$SYSp <- paste(trees$site, trees$year, trees$treesp)
trees$SYT <- paste(trees$site, trees$year, trees$tree)
trees$SY <- paste(trees$site, trees$year)
# site year, species, count of each
treesunique <- data.frame(SYT=unique(trees$SYT))
match <- pmatch(treesunique$SYT, trees$SYT)
treesunique$SYSp <- trees$SYSp[match]
treesunique$SY <- trees$SY[match]
treesunique$treesp <- trees$treesp[match]

countSTSp <- data.frame(table(treesunique$SYSp)) # trees beaten per species per site per year
match2 <- pmatch(countSTSp$Var1, treesunique$SYSp)
countSTSp$treesp <- treesunique$treesp[match2]
countSTSp$SY <- treesunique$SY[match2]
countSTSp$Var1 <- NULL

#countwide <- spread(countSTSp, treesp, Freq)

SpSY <- data.frame(table(countSTSp$SY)) # species per site per year

countSTSp$npart <- countSTSp$Freq*(countSTSp$Freq-1)
divSY <- aggregate(npart~ SY, countSTSp, sum)
divSY$N <- SpSY$Freq
divSY$Npart <- divSY$N*(divSY$N-1)
divSY$D <- 1-(divSY$npart/divSY$Npart)
divSY$expD <- exp(divSY$D)
expD.SY <- data.frame(SY=divSY$SY, expD=divSY$expD)

###  using expD to remove -Inf and makes positive and smaller value variation

cater <- read.csv("~/Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
cater$SY <- paste(cater$site, cater$year)
cater <- merge(cater, expD.SY, by="SY", duplicates.ok=TRUE, all.x=TRUE)

cater$datecentred <- cater$date-mean(cater$date)
cater$siteday <- paste(cater$site, cater$date, cater$year)
cater$sitetree <- paste(cater$site, cater$tree)
cater$year <- as.factor(cater$year)

########################################################################
#### Model with interaction between expD and quadratic ####
########################################################################

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#expDivQuad<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2)*expD+I(expD^2), random=~sitetree+site+siteday, family="poisson", data=cater, prior=prior, nitt=300000, burnin=30000)
#save(expDivQuad, file = "~/Documents/Models/expDivQuad.RData")
load("~/Documents/Models/expDivQuad.RData")

summary(expDivQuad)
plot(expDivQuad$Sol) #fixedeffects
