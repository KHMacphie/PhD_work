#load Dataframe.R script

#############################################
#### Model: Variance decomposition model ####
#############################################

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


#AbundVariance<- MCMCglmm(caterpillars~datescaled+I(datescaled^2), 
#                     random=~recorder+year+siteyear+siteday+site+treeID+tree.species+yearday, 
#                     family="poisson", data=cater_habitat, prior=prior, nitt=2000000, burnin=50000)
#save(AbundVariance, file = "~/Documents/Models/AbundVariance.RData")
load("~/Documents/Models/AbundVariance.RData") 

#### Checking model fits the data and converged ####
plot(AbundVariance) #look at fixed effect and random term trace plots 
AbundVariance.Sim<-simulate(AbundVariance,nsim=1000) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(AbundVariance.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
#normally not too many rogue values but reduce x axis a bit to see main distribution relative to observed value 
hist(apply(AbundVariance.Sim,2,sum), breaks=10000, xlim=c(0,100000))
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))}  # function for proportion of zeros
hist(apply(AbundVariance.Sim,2,propzero), breaks=50) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data

#### Variance from fixed effects ####
#from Jarrod: have pos as a vector indicating the position of the relevant terms (for example 2:3 if the 2nd and 3rd terms are x and x^2) and then do this:
V<-cov(as.matrix(AbundVariance$X[,2:3]))
R2<-apply(AbundVariance$Sol[,2:3], 1, function(x){x%*%V%*%x}) # not actual R2- variance explained by fixed effects (marginal)

Variances <- data.frame(AbundVariance$VCV)
Variances$Fixed <- R2
# Variance explained by random terms
VarProp.df <- data.frame(Term=c(colnames(Variances))) #dataframe with column for random term 

## loop for mean proportion of variance explained by each term with CIs
# mean used otherwise they don't sum to one, it is a little bit off the mode of some terms
for(i in 1:length(VarProp.df$Term)) {   
  VarProp.df$Proportion[i] <- mean(Variances[,i]/rowSums(Variances))
  A <- HPDinterval(mcmc(Variances[,i]/rowSums(Variances)))
  VarProp.df$ProportionLCI[i] <- A[1]
  VarProp.df$ProportionUCI[i] <- A[2]
} 

# correcting names
VarProp.df$Term <- gsub("recorder","Recorder", VarProp.df$Term)
VarProp.df$Term <- gsub("yearday","Day", VarProp.df$Term)
VarProp.df$Term <- gsub("siteyear","SiteYear", VarProp.df$Term)
VarProp.df$Term <- gsub("siteday","SiteDay", VarProp.df$Term)
VarProp.df$Term <- gsub("site","Site", VarProp.df$Term)
VarProp.df$Term <- gsub("year","Year", VarProp.df$Term)
VarProp.df$Term <- gsub("treeID","TreeID", VarProp.df$Term)
VarProp.df$Term <- gsub("tree.species","TreeTaxa", VarProp.df$Term)
VarProp.df$Term <- gsub("units","Residual", VarProp.df$Term)

# percentage rather than proportion
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
  ylab("Proportion of variance")+
  theme_bw()+
  theme(text=element_text(size= 20))


#### Riverplot package Sankey ####
library(riverplot)

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

####fixed effects
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

####random terms
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

#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableAbundVariance.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)


