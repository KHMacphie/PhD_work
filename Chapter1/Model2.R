#load Dataframe.R script rm(list=ls()) setwd('/Users/s1205615/')

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


#TTHA20<- MCMCglmm(caterpillars~Total_cent, 
#                  random=~tree.species+idv(~Alder_cent+Ash_cent+Beech_cent+Birch_cent+Elm_cent+Hazel_cent+Oak_cent+Rowan_cent+Sycamore_cent+Willow_cent+Conifer_cent+OthDecid_cent)+site+year+siteyear+treeID+siteday+recorder, 
#                  family="poisson", data=cater_habitat, prior=prior, nitt=4500000, burnin=500000, pr=TRUE, thin=800)
#save(TTHA20, file = "~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/TTHA20.RData") #samplesize 5000
load("~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/TTHA20.RData")
summary(TTHA20)


#### Checking model fits the data and converged ####
plot(TTHA20) #look at fixed effect and random term trace plots 
TTHA20.Sim<-simulate(TTHA20,nsim=1000) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(TTHA20.Sim,2,sum), breaks=1000) #histogram of simulation predictions for total abundance
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(TTHA20.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data


# dataframe for medians and CIs for beaten tree taxa effects
TTHA <- TTHA20$Sol[,3:12] # crop to just the columns wanted
TTHA.df <- data.frame(treetaxa=c(colnames(TTHA))) #dataframe with column for beaten tree taxa
TTHA.df$coeff <- apply(TTHA,2, median) # median 
for(i in 1:length(TTHA.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA[,i])
  TTHA.df$lowci[i] <- A["var1","lower"] 
  TTHA.df$upci[i] <- A["var1","upper"] 
} 
TTHA.df$treetaxa <- gsub("tree.species.","", TTHA.df$treetaxa) # adjust name
#write.csv(TTHA.df,'~/Documents/Models/Tables/Inc2020/TTHA.TTMetrics.csv')

#dataframe for medians and CIs for habitat FS effects
TTHA2 <- TTHA20$Sol[,13:24] # crop to just the columns wanted
TTHA2.df <- data.frame(treetaxa=c(colnames(TTHA2))) #dataframe with column for FS tree taxa
TTHA2.df$coeff <- apply(TTHA2,2, median) # median 
for(i in 1:length(TTHA2.df$treetaxa)) {   # loop for CIs
  A <- HPDinterval(TTHA2[,i])
  TTHA2.df$lowci[i] <- A["var1","lower"] 
  TTHA2.df$upci[i] <- A["var1","upper"] 
} 
TTHA2.df$treetaxa <- gsub("_cent.NA.1","", TTHA2.df$treetaxa) # adjust name
TTHA2.df$treetaxa <- gsub("OthDecid","Other", TTHA2.df$treetaxa) # adjust name
#write.csv(TTHA2.df,'~/Documents/Models/Tables/Inc2020/TTHA.FSMetrics.csv')

#plot beaten tree taxa effects
Colours <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

(plot1 <- ggplot(TTHA.df, aes(fct_rev(treetaxa), coeff, col=treetaxa))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  coord_flip()+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=Colours)+
  geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+
  xlab("Tree Taxon")+
  ylab("Coefficient (log scale)")+
  guides(linetype="none", colour="none"))

(plot1 <- ggplot(TTHA.df, aes(treetaxa, exp(coeff), col=treetaxa))+
    geom_point(size=1.5, alpha=0.9)+
    geom_errorbar(aes(ymax=exp(upci), ymin=exp(lowci), width=0.2))+
    theme_bw()+
    #coord_flip()+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    scale_colour_manual(values=Colours)+
    geom_hline(yintercept=1, linetype="dashed", colour="black", size=0.3)+
    xlab("Tree Taxon")+
    ylab("Coefficient (exponentiated)")+
    guides(linetype="none", colour="none")) #4x8

(plot1 <- ggplot(TTHA.df, aes(treetaxa, exp(coeff), col=treetaxa))+
    geom_point(size=3, alpha=0.9)+
    geom_errorbar(aes(ymax=exp(upci), ymin=exp(lowci), width=0.2), size=1)+
    theme_bw()+
    #coord_flip()+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    scale_colour_manual(values=Colours)+
    geom_hline(yintercept=1, linetype="dashed", colour="black", size=0.3)+
    xlab("Tree Taxon")+
    ylab("Coefficient (exponentiated)")+
    guides(linetype="none", colour="none")) #4x8

#plot for habitat FS's effects
Colours <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid", "gray57", "gray35")
(plot2 <- ggplot(TTHA2.df, aes(fct_rev(fct_inorder(treetaxa)), coeff, col=fct_inorder(treetaxa)))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=Colours)+
  geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+
  coord_flip()+
  xlab("Woodland composition category")+
  ylab("Coefficient (log scale)")+
  guides(linetype="none", colour="none"))

(plot2 <- ggplot(TTHA2.df, aes(fct_inorder(treetaxa), coeff, col=fct_inorder(treetaxa)))+
    geom_point(size=1.5, alpha=0.9)+
    geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.2))+
    theme_bw()+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    scale_colour_manual(values=Colours)+
    geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+
    xlab("Woodland composition category")+
    ylab("Coefficient (log scale)")+
    guides(linetype="none", colour="none"))

hist <- data.frame(var=as.numeric(TTHA20$VCV[,2]))
ggplot(hist, aes(var))+
  geom_histogram(binwidth=0.00005, col=1, fill="gray")+
  xlab("Variance")+
  ylab("Frequency")+
  theme_bw()+
  theme(text = element_text(size=15))


#Proportional difference to oak exp(each treetaxa-oak random effects [log scale])
alderOD <- (TTHA20$Sol[,3]-TTHA20$Sol[,9])
ashOD <- (TTHA20$Sol[,4]-TTHA20$Sol[,9])
beechOD <- (TTHA20$Sol[,5]-TTHA20$Sol[,9])
birchOD <- (TTHA20$Sol[,6]-TTHA20$Sol[,9])
elmOD <- (TTHA20$Sol[,7]-TTHA20$Sol[,9])
hazelOD <- (TTHA20$Sol[,8]-TTHA20$Sol[,9])
rowanOD <- (TTHA20$Sol[,10]-TTHA20$Sol[,9])
sycamoreOD <- (TTHA20$Sol[,11]-TTHA20$Sol[,9])
willowOD <- (TTHA20$Sol[,12]-TTHA20$Sol[,9])

# extract medians and CIs
TTOD <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"),
                   medianOD=c(exp(median(alderOD)), exp(median(ashOD)), exp(median(beechOD)), exp(median(birchOD)), exp(median(elmOD)), exp(median(hazelOD)), exp(median(rowanOD)), exp(median(sycamoreOD)), exp(median(willowOD))),
                   lowci=c(exp(HPDinterval(alderOD)[1]), exp(HPDinterval(ashOD)[1]),exp(HPDinterval(beechOD)[1]),exp(HPDinterval(birchOD)[1]),exp(HPDinterval(elmOD)[1]),exp(HPDinterval(hazelOD)[1]),exp(HPDinterval(rowanOD)[1]),exp(HPDinterval(sycamoreOD)[1]),exp(HPDinterval(willowOD)[1])),
                   upci=c(exp(HPDinterval(alderOD)[2]), exp(HPDinterval(ashOD)[2]),exp(HPDinterval(beechOD)[2]),exp(HPDinterval(birchOD)[2]),exp(HPDinterval(elmOD)[2]),exp(HPDinterval(hazelOD)[2]),exp(HPDinterval(rowanOD)[2]),exp(HPDinterval(sycamoreOD)[2]),exp(HPDinterval(willowOD)[2])))
#write.csv(TTOD,'~/Documents/Models/Tables/Inc2020/TTHAODMetrics.csv')


# plot of proportional differences to oak
ggplot(TTOD, aes(fct_rev(TT), medianOD))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=(upci), ymin=(lowci), width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Difference to oak (proportional)")+
  coord_flip()+ 
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1)) #saved 6"x4"

#### Plot with slope for effect on abundance for FS of each tree taxa
#Plotted for range of foliage scores of each taxon present at sites
summary(rowSums(Habitat_Site[2:13])) #real FS total scale
summary(Habitat_Site$Total_cent) #used in model
summary(as.matrix(Habitat_Site[14:25])) #taxaFS's in model

totcentfs <- seq(-44.561, 45.029, 0.001) # for slope calcs 
totfs <- seq(24.93,114.52, 0.001) # for plot in real fs scale
taxacentfs <- seq(-5.791, 91.334, 0.001)  # fortaxaFS's slope calcs

# range of FS's for each taxa
aldercent <- seq(min(Habitat_Site$Alder_cent), max(Habitat_Site$Alder_cent), 0.01)
ashcent <- seq(min(Habitat_Site$Ash_cent), max(Habitat_Site$Ash_cent), 0.01)
beechcent <- seq(min(Habitat_Site$Beech_cent), max(Habitat_Site$Beech_cent), 0.01)
birchcent <- seq(min(Habitat_Site$Birch_cent), max(Habitat_Site$Birch_cent), 0.01)
elmcent <- seq(min(Habitat_Site$Elm_cent), max(Habitat_Site$Elm_cent), 0.01)
hazelcent <- seq(min(Habitat_Site$Hazel_cent), max(Habitat_Site$Hazel_cent), 0.01)
oakcent <- seq(min(Habitat_Site$Oak_cent), max(Habitat_Site$Oak_cent), 0.01)
rowancent <- seq(min(Habitat_Site$Rowan_cent), max(Habitat_Site$Rowan_cent), 0.01)
sycamorecent <- seq(min(Habitat_Site$Sycamore_cent), max(Habitat_Site$Sycamore_cent), 0.01)
willowcent <- seq(min(Habitat_Site$Willow_cent), max(Habitat_Site$Willow_cent), 0.01)
conifercent <- seq(min(Habitat_Site$Conifer_cent), max(Habitat_Site$Conifer_cent), 0.01)
othercent <- seq(min(Habitat_Site$OthDecid_cent), max(Habitat_Site$OthDecid_cent), 0.01)

Totalslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2])*totcentfs
# each slope coefficient is the deviation from the total slope from fixed effects
Alderslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,13])*aldercent
Ashslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,14])*ashcent
Beechslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,15])*beechcent
Birchslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,16])*birchcent
Elmslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,17])*elmcent
Hazelslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,18])*hazelcent
Oakslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,19])*oakcent
Rowanslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,20])*rowancent
Sycamoreslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,21])*sycamorecent
Willowslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,22])*willowcent
Coniferslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,23])*conifercent
Otherslope <- median(TTHA20$Sol[,1])+median(TTHA20$Sol[,2]+TTHA20$Sol[,24])*othercent

## List of colours used for tree taxa (alphabetical order) throughout all figures
#AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")
#Other+Conifer <- c("gray57", "gray35")

par(mfcol=c(1,1))
plot(aldercent, exp(Alderslope), type="l", lty=6, ylim=c(0.00848,0.07467), xlim=c(-5.691,91.334), col="darkred", xlab="Deviation in FS from mean", ylab="Caterpillar Abundance")
points(ashcent, exp(Ashslope), type="l", col="firebrick3", lty=6)
points(beechcent, exp(Beechslope), type="l", col="chocolate2", lty=6)
points(birchcent, exp(Birchslope), type="l", col="goldenrod", lty=6)
points(elmcent, exp(Elmslope), type="l", col="olivedrab4", lty=6)
points(oakcent, exp(Oakslope), type="l", col="deepskyblue3")
points(sycamorecent, exp(Sycamoreslope), type="l", col="slateblue2", lty=6)
points(willowcent, exp(Willowslope), type="l", col="orchid", lty=6)
points(conifercent, exp(Coniferslope), type="l", col="gray57", lty=6)
points(othercent, exp(Otherslope), type="l", col="gray35", lty=6)
points(hazelcent, exp(Hazelslope), type="l", col="darkgreen", lty=6)
points(rowancent, exp(Rowanslope), type="l", col="royalblue4", lty=6)
legend("topleft", legend=c("Alder","Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow", "Conifer", "Other"),
       lty=c(1,1,1,1,1,1,1,1,1,1,1,1), lwd=3, 
       col=c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid","gray57","grey35"), 
       cex=0.85, seg.len=0.8, bty = "n") # saved as 7"x4.5"

## Plot in ggplot
TTFSlong <- data.frame(FS=c(aldercent,ashcent,beechcent,birchcent,elmcent,hazelcent,oakcent,rowancent,sycamorecent,willowcent,conifercent,othercent),
                       logabund=c(Alderslope,Ashslope,Beechslope,Birchslope,Elmslope,Hazelslope,Oakslope,Rowanslope,Sycamoreslope,Willowslope,Coniferslope,Otherslope),
                       TT=c(rep("Alder",length(aldercent)),rep("Ash",length(ashcent)),rep("Beech",length(beechcent)),rep("Birch",length(birchcent)),rep("Elm",length(elmcent)),
                            rep("Hazel",length(hazelcent)),rep("Oak",length(oakcent)),rep("Rowan",length(rowancent)),rep("Sycamore",length(sycamorecent)),rep("Willow",length(willowcent)),
                            rep("Conifer",length(conifercent)),rep("OtherDecid",length(othercent))))

Colours <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid", "gray57", "gray35")

(plot3 <- ggplot(TTFSlong, aes(FS, exp(logabund), col=fct_inorder(TT)))+
    geom_line(size=0.6,aes(linetype=fct_inorder(TT)))+
    theme_bw()+
    theme(text = element_text(size=15))+
    xlab("Deviation in FS from mean")+
    ylab("Caterpillar Abundance")+
    scale_colour_manual(values=Colours)+
    scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","dashed","solid","dashed","dashed","dashed","dashed","dashed"))+
    guides(color = "none", linetype="none"))
gap <- ggplot()+theme_void()
row1 <- grid.arrange(plot1,gap, plot2,gap,plot3, ncol=5, widths=c(1,0.1,1,0.1,1)) #saved as 6"x12"

(plot3 <- ggplot(TTFSlong, aes(FS, exp(logabund), col=fct_inorder(TT)))+
    geom_line(size=0.6,aes(linetype=fct_inorder(TT)))+
    theme_bw()+
    theme(text = element_text(size=15))+
    xlab("Deviation in foliage score")+
    ylab("Caterpillar Abundance")+
    scale_colour_manual(values=Colours)+
    scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","dashed","solid","dashed","dashed","dashed","dashed","dashed")))+
  guides(color = "none", linetype="none")) #5x8 with legend 6x6.5without
    

totalplot <- data.frame(fs=totfs, logabund=Totalslope)
exp(-2.464) exp(-4.677)

ggplot(totalplot, aes(fs, exp(logabund)))+
    geom_line(size=0.6)+
    theme_bw()+
    theme(text = element_text(size=15))+
    xlab("Deviation in foliage score")+
    ylab("Caterpillar Abundance")+
    ylim(exp(-4.677), exp(-2.464))+
    xlim(-5.791,91.329)
    #scale_colour_manual(values=Colours)+
    #scale_linetype_manual(values=c("dashed","dashed","dashed","dashed","dashed","dashed","solid","dashed","dashed","dashed","dashed","dashed")) #5x8


### looking at credible intervals on oak slope
oakhabitat <- data.frame(FS = seq(min(Habitat_Site$Oak_cent), max(Habitat_Site$Oak_cent), 0.1))
for(i in 1:length(oakhabitat$FS)){
  oakhabitat$median[i] <- median(TTHA20$Sol[,1]+((TTHA20$Sol[,2]+TTHA20$Sol[,19])*oakhabitat$FS[i]))
  oakhabitat$LCI[i] <- HPDinterval(TTHA20$Sol[,1]+((TTHA20$Sol[,2]+TTHA20$Sol[,19])*oakhabitat$FS[i]))[1]
  oakhabitat$UCI[i] <- HPDinterval(TTHA20$Sol[,1]+((TTHA20$Sol[,2]+TTHA20$Sol[,19])*oakhabitat$FS[i]))[2]
}

plot(oakhabitat$FS, exp(oakhabitat$median), type="l", ylim=c(0,0.3), lwd=2)
points(oakhabitat$FS, exp(oakhabitat$LCI), type="l")
points(oakhabitat$FS, exp(oakhabitat$UCI), type="l")

maxoak <- TTHA20$Sol[,1]+((TTHA20$Sol[,2]+TTHA20$Sol[,19])*max(Habitat_Site$Oak_cent))
nooak <- TTHA20$Sol[,1]+((TTHA20$Sol[,2]+TTHA20$Sol[,19])*min(Habitat_Site$Oak_cent))
exp(median(maxoak-nooak)) #median prop diff in abundance from no oak to max oak = 5.19x
exp(HPDinterval(mcmc(maxoak-nooak))) #CIs prop diff in abundance from no oak to max oak = 1.73 - 17.45x


############################
#### Model output table ####
############################

#for random terms use posterior mode and fixed terms mean
library(MCMCglmm)

####fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(TTHA20$Sol[,1]),3)," (",
                      round(HPDinterval(TTHA20$Sol[,1])[1],3)," - ",
                      round(HPDinterval(TTHA20$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(TTHA20$Sol[,1]))),
  c("Total Foliage Score",paste(round(mean(TTHA20$Sol[,2]),3)," (",
                                round(HPDinterval(TTHA20$Sol[,2])[1],3)," - ",
                                round(HPDinterval(TTHA20$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(TTHA20$Sol[,2]))))

####random terms 
column<-1
treetaxa<-c("Sampled Tree Taxa",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                                      round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA20$VCV[, column])))

column<-2
habitat<-c("Habitat Composition",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                                       round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(TTHA20$VCV[, column])))

column<-3
site<-c("Site",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA20$VCV[, column])))

column<-4
year<-c("Year",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                     round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                     round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
        round(effectiveSize(TTHA20$VCV[, column])))

column<-5
siteyear<-c("Site Year",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                              round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                              round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA20$VCV[, column])))

column<-6
treeID<-c("Tree ID",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                          round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                          round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(TTHA20$VCV[, column])))

column<-7
siteday<-c("Site Day",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                            round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                            round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(TTHA20$VCV[, column])))

column<-8
recorder<-c("Recorder",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA20$VCV[, column])))

column<-9
residual<-c("Residual",paste(round(posterior.mode(TTHA20$VCV[, column]),3)," (",
                             round(HPDinterval(TTHA20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(TTHA20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(TTHA20$VCV[, column])))




random<-rbind(treetaxa, habitat, site, year, siteyear, treeID, siteday,recorder,residual)


#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/Inc2020/TableTTHA20.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)



#### random look at site effects ####
siteeffs <- TTHA20$Sol[,25:68]
colnames(siteeffs) <- gsub("site.","", colnames(siteeffs))
site.df <- data.frame(site=c(colnames(siteeffs)))
site.df$coeff <- apply(siteeffs,2,mean)
for(i in 1:length(site.df$site)) {   # loop for CIs
  A <- HPDinterval(siteeffs[,i])
  site.df$lowci[i] <- A["var1","lower"] 
  site.df$upci[i] <- A["var1","upper"] 
}

site <- read.csv("~/Dropbox/master_data/site/site_details.csv")
colnames(site)
store <- pmatch(site.df$site, site$site)
site.df <- cbind(site.df,site$Mean.Lat[store],site$Mean.Elev[store])
site.df$elev <- site.df$'site$Mean.Elev[store]'
site.df$lat <- site.df$'site$Mean.Lat[store]'
site.df <- site.df[order(site.df$lat),] 

#plot by latitude
plotlat <- ggplot(site.df, aes(fct_inorder(site), coeff))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  coord_flip()+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+
  xlab("Site by latitude")+
  ylab("Coefficient (log scale)")

site.df <- site.df[order(site.df$elev),] 
plotelev <- ggplot(site.df, aes(fct_inorder(site), coeff))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  theme_bw()+
  coord_flip()+
  theme(text = element_text(size=15))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+
  xlab("Site by elevation")+
  ylab("Coefficient (log scale)") 
row1 <- grid.arrange(plotlat, plotelev, ncol=2, widths=c(1,1)) #saved as 8"x7"

