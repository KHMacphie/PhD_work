#################################################
## Random effect output for tree species ##
#################################################


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

##### Setting up dataframe #####

site <- read.csv("Dropbox/master_data/site/site_details.csv")
colnames(site)[4:6] <- c("latitude", "longitude", "elevation")

cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating_correctingID.csv")
cater$year <- as.factor(cater$year)


pmatch(cater$site,site$site,duplicates.ok=TRUE)
all_data<- merge(cater, site, by="site", duplicates.ok=TRUE)


all_data$treespecies <- all_data$`tree species`
all_data$yearsite<- paste(all_data$site, all_data$year)
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data$siteday <- paste(all_data$site, all_data$date, all_data$year)
all_data$obs<-as.factor(seq(1,length(all_data[,1])))

########################################################################################
##### MCMCglmm for base date, year date^2 model with tree species as random effect #####

k<-1000
prior2<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))


#TreeSpeciesRE<- MCMCglmm(caterpillars~date*year+I(date^2), random=~site+sitetree+siteday+treespecies, family="poisson", data=all_data, prior=prior3, nitt=200000, burnin=20000, pr=TRUE)
#save(TreeSpeciesRE, file = "~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/TreeSpeciesRE.RData")
load("~/Documents/PhD/GitHub/R/Caterpillar analysis/Models/TreeSpeciesRE.RData")

#check if model generates sensible results
# TreeSpeciesRE
TreeSpeciesRE.Sim<-simulate(TreeSpeciesRE,nsim=100)
#sum(all_data$caterpillars)
par(mfcol=c(1,1))
hist(apply(TreeSpeciesRE.Sim,2,sum))
abline(v=sum(all_data$caterpillars),col=2)

propzero <- function(x){return(length(which(x==0))/length(x))}
hist(apply(TreeSpeciesRE.Sim,2,propzero))
abline(v=propzero(all_data$caterpillars), col="red")
## not quite as good as usual- maybe slightly over estimating zeros?

summary(TreeSpeciesRE)

###########################################################
#### Calculating the coefficient for each tree sp. with CI #### 

# TreeSpeciesRE tree species coefficients 
which(colnames(TreeSpeciesRE$Sol)=="treespecies.?") #= 5022
which(colnames(TreeSpeciesRE$Sol)=="treespecies.Willow") #= 5039

treespREcropped <- TreeSpeciesRE$Sol[,5022:5039] # crop to just the columns wanted
treesp.df <- data.frame(treesp=c(colnames(treespREcropped))) #dataframe with column for yearsite 
treesp.df$coeff <- apply(treespREcropped,2, mean) # mean 
for(i in 1:length(treesp.df$treesp)) {   # loop for CIs
  A <- HPDinterval(treespREcropped[,i])
  treesp.df$lowci[i] <- A["var1","lower"] 
  treesp.df$upci[i] <- A["var1","upper"] 
} 
treesp.df$treesp <- gsub("treespecies.","", treesp.df$treesp)

ggplot(treesp.df, aes(treesp, coeff))+
  geom_point(size=3, alpha=0.5)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(axis.text.x= element_text(angle=90))


##########################################################
#### Interaction with tree species:date and site:date ####
##########################################################


all_data$datecentred <- all_data$date-mean(all_data$date)

a<-1000
prior<-list(R=list(V=diag(1), nu=0.002), 
               G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                      G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                      G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

#TreeSpeciesDateRE<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2), random=~sitetree+us(1+datecentred):site+us(1+datecentred):treespecies, family="poisson", data=all_data, prior=prior, nitt=200000, burnin=20000, pr=TRUE)
#save(TreeSpeciesDateRE, file = "~/Documents/Models/TreeSpeciesDateRE.RData")
load("~/Documents/Models/TreeSpeciesDateRE.RData")

a<-1000
prior2<-list(R=list(V=diag(1), nu=0.002), 
            G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a),
                   G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*a)))

#TreeSpDSY<- MCMCglmm(caterpillars~datecentred*year+I(datecentred^2), random=~sitetree+us(1+datecentred):site+us(1+datecentred):treespecies+us(1+datecentred):yearsite, family="poisson", data=all_data, prior=prior2, nitt=200000, burnin=20000, pr=TRUE)
#save(TreeSpDSY, file = "~/Documents/Models/TreeSpDSY.RData")
load("~/Documents/Models/TreeSpDSY.RData")
