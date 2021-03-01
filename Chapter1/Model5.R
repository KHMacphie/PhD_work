#load Dataframe.R script rm(list=ls()) setwd('/Users/s1205615/')

# Mass only available for 2017-19
cater_habitat_1720 <- subset(cater_habitat, year!="2014")
cater_habitat_1720 <- subset(cater_habitat_1720, year!="2015")
cater_habitat_1720 <- subset(cater_habitat_1720, year!="2016")

# reordering so year 2019 is reference level
cater_habitat_1720$year <- relevel(cater_habitat_1720$year, ref="2019")

##################################################################
#### Model: Phenological distribution of mass among tree taxa ####
##################################################################

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#Mass20<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled + I(datescaled^2), 
#                  random=~us(1+datescaled):tree.species + us(1+datescaled):site + year + siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
#                  family="cengaussian", data=cater_habitat_1720, prior=prior3, pr=TRUE, nitt=5500000, burnin=500000, thin=500)
#save(Mass20, file = "~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/Mass20.RData") #sample size 10000
load("~/Dropbox/Kirsty's/Chapter1/Models/Inc2020/Mass20.RData")
summary(Mass20)

#### Checking model fits the data and converged ####
plot(Mass20) #look at fixed effect and random term trace plots 
Mass20.Sim<-simulate(Mass20,nsim=100) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(Mass20.Sim,2,mean),50) #histogram of simulation predictions for mean mass (log scale)
abline(v=mean((cater_habitat_1720$logmpc1+cater_habitat_1720$logmpc2)/2, na.rm=TRUE),col=2) # red line for mean mass (from the mean of logmpc1 and logmpc2)

# tree taxon specific coefficients for date^2 (A), date (B) and intercept (C)
AlderC <- (Mass20$Sol[,1] + Mass20$Sol[,4])
AshC <- (Mass20$Sol[,1] + Mass20$Sol[,5])
BeechC <- (Mass20$Sol[,1] + Mass20$Sol[,6])
BirchC <- (Mass20$Sol[,1] + Mass20$Sol[,7])
ElmC <- (Mass20$Sol[,1] + Mass20$Sol[,8])
HazelC <- (Mass20$Sol[,1] + Mass20$Sol[,9])
OakC <- (Mass20$Sol[,1] + Mass20$Sol[,10])
RowanC <- (Mass20$Sol[,1] + Mass20$Sol[,11])
SycamoreC <- (Mass20$Sol[,1] + Mass20$Sol[,12])
WillowC <- (Mass20$Sol[,1] + Mass20$Sol[,13])

AlderB <- (Mass20$Sol[,2] + Mass20$Sol[,14])
AshB <- (Mass20$Sol[,2] + Mass20$Sol[,15])
BeechB <- (Mass20$Sol[,2] + Mass20$Sol[,16])
BirchB <- (Mass20$Sol[,2] + Mass20$Sol[,17])
ElmB <- (Mass20$Sol[,2] + Mass20$Sol[,18])
HazelB <- (Mass20$Sol[,2] + Mass20$Sol[,19])
OakB <- (Mass20$Sol[,2] + Mass20$Sol[,20])
RowanB <- (Mass20$Sol[,2] + Mass20$Sol[,21])
SycamoreB <- (Mass20$Sol[,2] + Mass20$Sol[,22])
WillowB <- (Mass20$Sol[,2] + Mass20$Sol[,23])

AlderA <- (Mass20$Sol[,3])
AshA <- (Mass20$Sol[,3])
BeechA <- (Mass20$Sol[,3])
BirchA <- (Mass20$Sol[,3])
ElmA <- (Mass20$Sol[,3])
HazelA <- (Mass20$Sol[,3])
OakA <- (Mass20$Sol[,3])
RowanA <- (Mass20$Sol[,3])
SycamoreA <- (Mass20$Sol[,3])
WillowA <- (Mass20$Sol[,3])

# Fixed effect coefficients 
MeanA <- (Mass20$Sol[,3])
MeanB <- (Mass20$Sol[,2])
MeanC <- (Mass20$Sol[,1])

#### Plotting growth curves ####

preddayscaled <- seq(-2.1199,2.0102,0.001)
predday <- unscale(preddayscaled, center= 146.7727, scale=14.04083)$V1


meanslope <- median(MeanC)+median(MeanB)*preddayscaled+median(MeanA)*preddayscaled^2
alder <- median(AlderC)+median(AlderB)*preddayscaled+median(AlderA)*preddayscaled^2
ash <- median(AshC)+median(AshB)*preddayscaled+median(AshA)*preddayscaled^2
beech <- median(BeechC)+median(BeechB)*preddayscaled+median(BeechA)*preddayscaled^2
birch <- median(BirchC)+median(BirchB)*preddayscaled+median(BirchA)*preddayscaled^2
elm <- median(ElmC)+median(ElmB)*preddayscaled+median(ElmA)*preddayscaled^2
hazel <- median(HazelC)+median(HazelB)*preddayscaled+median(HazelA)*preddayscaled^2
oak <- median(OakC)+median(OakB)*preddayscaled+median(OakA)*preddayscaled^2
rowan <- median(RowanC)+median(RowanB)*preddayscaled+median(RowanA)*preddayscaled^2
sycamore <- median(SycamoreC)+median(SycamoreB)*preddayscaled+median(SycamoreA)*preddayscaled^2
willow <- median(WillowC)+median(WillowB)*preddayscaled+median(WillowA)*preddayscaled^2
intsamples <- subset(cater_habitat, mpc1 == 0.001, 
                     select=c(date, mpc1))
mycolblack <- rgb(100, 100, 100, max = 250, alpha = 30, names = "blacktrans")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

fixedslope <- data.frame(date=predday, mass=exp(meanslope))
fixedslope <- fixedslope[fixedslope$date < 172, ]

dataplot <- ggplot(cater_habitat_1720, aes(date, mpc2))+
  geom_point(col=mycolblack, size=0.5)+
  geom_point(data=intsamples, aes(date, mpc1), col=mycolblack, size=0.5)+
  geom_line(data=fixedslope, aes(date,mass), col=1, lty=5)+
  xlab("Ordinal Date")+
  ylab("Mass (g)")+
  coord_cartesian(xlim=c(117, 172), ylim=c(0,1))+
  scale_y_continuous(expand = c(0,0.02))+
  scale_x_continuous(breaks=seq(120,170,10))+
  theme_bw()+
  theme(text=element_text(size= 15))+
  annotate(geom="text", x=132, y=0.95, label="- - -  Fixed effect",
           color="black")

taxoncurves <- data.frame(date = predday, FixedEffect=exp(meanslope), Alder=exp(alder), Ash=exp(ash), Beech=exp(beech), Birch=exp(birch), Elm=exp(elm), Hazel=exp(hazel), Oak=exp(oak), Rowan=exp(rowan),Sycamore=exp(sycamore),Willow=exp(willow)) 
taxoncurveslong <- gather(taxoncurves, Taxon, Mass, 2:12)
taxoncurveslong$Taxon <- as.factor(taxoncurveslong$Taxon)
taxoncurveslong$Taxon <- factor(taxoncurveslong$Taxon , c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak","Rowan","Sycamore","Willow" ,"FixedEffect"))

AllCols <- c( "darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid","black")

taxonplot <- ggplot(taxoncurveslong, aes(date,Mass, col=Taxon))+
  geom_line(aes(linetype=Taxon))+
  xlab("Ordinal Date")+
  ylab("Mass (g)")+
  scale_colour_manual(values=AllCols)+
  scale_linetype_manual(values=c(1,1,1,1,1,1,1,1,1,1,5))+
  xlim(117, 170)+
  ylim(0,0.042)+
  geom_vline(xintercept=168, linetype="dotted", colour="gray45", size=0.5)+
  scale_y_continuous(expand = c(0,0.004))+
  guides(color = "none", linetype="none")+
  theme_bw()+
  theme(text=element_text(size= 15))

#### Mass on 168 ####
latedate <- scale(168,  center= 146.7727, scale=14.04083)[1,1]
Mean168 <- exp(MeanC + MeanB*latedate + MeanA*latedate^2)
Alder168 <- exp(AlderC + AlderB*latedate + AlderA*latedate^2)
Ash168 <- exp(AshC + AshB*latedate + AshA*latedate^2)
Beech168 <- exp(BeechC + BeechB*latedate + BeechA*latedate^2)
Birch168 <- exp(BirchC + BirchB*latedate + BirchA*latedate^2)
Elm168 <- exp(ElmC + ElmB*latedate + ElmA*latedate^2)
Hazel168 <- exp(HazelC + HazelB*latedate + HazelA*latedate^2)
Oak168 <- exp(OakC + OakB*latedate + OakA*latedate^2)
Rowan168 <- exp(RowanC + RowanB*latedate + RowanA*latedate^2)
Sycamore168 <- exp(SycamoreC + SycamoreB*latedate + SycamoreA*latedate^2)
Willow168 <- exp(WillowC + WillowB*latedate + WillowA*latedate^2)

#### Difference to mean on 168 (0.96) ####
Alder168dif <- Alder168-Mean168
Ash168dif <- Ash168-Mean168
Beech168dif <- Beech168-Mean168
Birch168dif <- Birch168-Mean168
Elm168dif <- Elm168-Mean168
Hazel168dif <- Hazel168-Mean168
Oak168dif <- Oak168-Mean168
Rowan168dif <- Rowan168-Mean168
Sycamore168dif <- Sycamore168-Mean168
Willow168dif <- Willow168-Mean168

#### Difference to mean on 168 (0.96) ####
Alder168prop <- Alder168/Mean168
Ash168prop <- Ash168/Mean168
Beech168prop <- Beech168/Mean168
Birch168prop <- Birch168/Mean168
Elm168prop <- Elm168/Mean168
Hazel168prop <- Hazel168/Mean168
Oak168prop <- Oak168/Mean168
Rowan168prop <- Rowan168/Mean168
Sycamore168prop <- Sycamore168/Mean168
Willow168prop <- Willow168/Mean168

# median coefficient and CIs for mass of 168 for each tree taxon
Mass168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                      median=c(median(Alder168), median(Ash168), median(Beech168), median(Birch168), median(Elm168), median(Hazel168), median(Oak168), median(Rowan168), median(Sycamore168), median(Willow168)),
                      lowci=c(HPDinterval(Alder168)[1],HPDinterval(Ash168)[1],HPDinterval(Beech168)[1],HPDinterval(Birch168)[1],HPDinterval(Elm168)[1],HPDinterval(Hazel168)[1],HPDinterval(Oak168)[1],HPDinterval(Rowan168)[1],HPDinterval(Sycamore168)[1],HPDinterval(Willow168)[1]),
                      upci=c(HPDinterval(Alder168)[2],HPDinterval(Ash168)[2],HPDinterval(Beech168)[2],HPDinterval(Birch168)[2],HPDinterval(Elm168)[2],HPDinterval(Hazel168)[2],HPDinterval(Oak168)[2],HPDinterval(Rowan168)[2],HPDinterval(Sycamore168)[2],HPDinterval(Willow168)[2]),
                      medianprop=c(median(Alder168prop), median(Ash168prop), median(Beech168prop), median(Birch168prop), median(Elm168prop), median(Hazel168prop), median(Oak168prop), median(Rowan168prop), median(Sycamore168prop), median(Willow168prop)),
                      proplowci=c(exp(HPDinterval(log(Alder168prop)))[1],exp(HPDinterval(log(Ash168prop)))[1],exp(HPDinterval(log(Beech168prop)))[1],exp(HPDinterval(log(Birch168prop)))[1],exp(HPDinterval(log(Elm168prop)))[1],exp(HPDinterval(log(Hazel168prop)))[1],exp(HPDinterval(log(Oak168prop)))[1],exp(HPDinterval(log(Rowan168prop)))[1],exp(HPDinterval(log(Sycamore168prop)))[1],exp(HPDinterval(log(Willow168prop)))[1]),
                      propupci=c(exp(HPDinterval(log(Alder168prop)))[2],exp(HPDinterval(log(Ash168prop)))[2],exp(HPDinterval(log(Beech168prop)))[2],exp(HPDinterval(log(Birch168prop)))[2],exp(HPDinterval(log(Elm168prop)))[2],exp(HPDinterval(log(Hazel168prop)))[2],exp(HPDinterval(log(Oak168prop)))[2],exp(HPDinterval(log(Rowan168prop)))[2],exp(HPDinterval(log(Sycamore168prop)))[2],exp(HPDinterval(log(Willow168prop)))[2]))
#write.csv(Mass168,'~/Documents/Models/Tables/Inc2020/Mass168Metrics.csv')

# Plot of mass difference to the fixed effect prediction
Mass168propplot <- ggplot(Mass168, aes(fct_rev(TT), medianprop, col=TT))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=propupci, ymin=proplowci, width=0.5))+
  geom_hline(yintercept=1, linetype="longdash", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Prop. difference to fixed eff.")+
  theme_bw()+
  theme(text=element_text(size= 15))+
  scale_colour_manual(values=AllTaxaCols)+
  guides(color = "none") # saved as 6"x4.5"

space <- ggplot() + theme_void()
MassPlots <- grid.arrange(dataplot,space, taxonplot,space, Mass168propplot, ncol = 5, widths = c(1,0.1,1,0.1,1.1)) #saved as 6"x12"


#### For supp mat ####
## Difference to oak on 168 (0.96)
AlderOD <- Alder168/Oak168
AshOD <- Ash168/Oak168
BeechOD <- Beech168/Oak168 
BirchOD <- Birch168/Oak168 
ElmOD <- Elm168/Oak168 
HazelOD <- Hazel168/Oak168 
RowanOD <- Rowan168/Oak168
SycamoreOD <- Sycamore168/Oak168
WillowOD <- Willow168/Oak168

# median and CIs for mass difference to oak on day 168
OD168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                    median=c(median(AlderOD), median(AshOD), median(BeechOD), median(BirchOD), median(ElmOD), median(HazelOD), median(RowanOD), median(SycamoreOD), median(WillowOD)),
                    lowci=c(exp(HPDinterval(log(AlderOD)))[1],exp(HPDinterval(log(AshOD)))[1],exp(HPDinterval(log(BeechOD)))[1],exp(HPDinterval(log(BirchOD)))[1],exp(HPDinterval(log(ElmOD)))[1],exp(HPDinterval(log(HazelOD)))[1],exp(HPDinterval(log(RowanOD)))[1],exp(HPDinterval(log(SycamoreOD)))[1],exp(HPDinterval(log(WillowOD)))[1]),
                    upci=c(exp(HPDinterval(log(AlderOD)))[2],exp(HPDinterval(log(AshOD)))[2],exp(HPDinterval(log(BeechOD)))[2],exp(HPDinterval(log(BirchOD)))[2],exp(HPDinterval(log(ElmOD)))[2],exp(HPDinterval(log(HazelOD)))[2],exp(HPDinterval(log(RowanOD)))[2],exp(HPDinterval(log(SycamoreOD)))[2],exp(HPDinterval(log(WillowOD)))[2]))
#write.csv(OD168,'~/Documents/Models/Tables/Inc2020/MassOD168Metrics.csv')

# Plot of mass difference to oak prediction
ggplot(OD168, aes(fct_rev(TT), median))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=1, linetype="dashed", colour="black", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Prop. difference to oak")+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 15),axis.text.x = element_text(angle = 45, hjust=1)) # saved as 6"x4"


meanmetrics <- data.frame("168mean"=c(exp(mean(log(Mean168)))), 
                          lci=c(exp(HPDinterval(mcmc(log(Mean168))))[1]),
                          uci=c(exp(HPDinterval(mcmc(log(Mean168))))[2]),
                          SS=c(length(which(is.na(Mean168)==FALSE))))
#write.csv(meanmetrics,'~/Documents/Models/Tables/Inc2020/MassFixedEff168Metrics.csv')

############################
#### Model output table ####   
############################

#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(Mass20$Sol[,1]),3)," (",
                                  round(HPDinterval(Mass20$Sol[,1])[1],3)," - ",
                                  round(HPDinterval(Mass20$Sol[,1])[2],3),")",sep=""),
                                  round(effectiveSize(Mass20$Sol[,1]))),
  c("Date scaled",paste(round(mean(Mass20$Sol[,2]),3)," (",
                                  round(HPDinterval(Mass20$Sol[,2])[1],3)," - ",
                                  round(HPDinterval(Mass20$Sol[,2])[2],3),")",sep=""),
                                  round(effectiveSize(Mass20$Sol[,2]))),
  c("Date² scaled",paste(round(mean(Mass20$Sol[,3]),3)," (",
                                  round(HPDinterval(Mass20$Sol[,3])[1],3)," - ",
                                  round(HPDinterval(Mass20$Sol[,3])[2],3),")",sep=""),
                                  round(effectiveSize(Mass20$Sol[,3]))))

#### random terms 
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                             round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass20$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                                          round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass20$VCV[, column])))


column<-4
treetaxa4<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                              round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass20$VCV[, column])))


column<-5
site5<-c("Site- Intercept var",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                     round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                     round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass20$VCV[, column])))

column<-6
site6<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                                  round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                                  round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass20$VCV[, column])))


column<-8
site8<-c("Site- Date slope var",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                                      round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass20$VCV[, column])))

column<-9
year<-c("Year",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                              round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                              round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass20$VCV[, column])))

column <- 10
siteyear<-c("Site-Year",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                              round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                              round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass20$VCV[, column])))


column<-13
recorder<-c("Recorder",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                             round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass20$VCV[, column])))


column<-12
siteday<-c("Site-Day",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                            round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                            round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(Mass20$VCV[, column])))


column<-11
treeID<-c("Tree ID",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                          round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                          round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(Mass20$VCV[, column])))

column<-14
weight<-c("Weighting",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                            round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                            round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(Mass20$VCV[, column])))

column<-15
residual<-c("Residual",paste(round(posterior.mode(Mass20$VCV[, column]),3)," (",
                             round(HPDinterval(Mass20$VCV[, column])[1],3)," - ",
                             round(HPDinterval(Mass20$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass20$VCV[, column])))




random<-rbind(treetaxa1,treetaxa2,treetaxa4,site5,site6,site8,year,siteyear, recorder, siteday, treeID, weight, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/Inc2020/TableMass20.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
