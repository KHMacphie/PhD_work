#load Dataframe.R script

# Mass only available for 2017-19
cater_habitat_1719 <- subset(cater_habitat, year!="2014")
cater_habitat_1719 <- subset(cater_habitat_1719, year!="2015")
cater_habitat_1719 <- subset(cater_habitat_1719, year!="2016")

# reordering so year 2019 is reference level
cater_habitat_1719$year <- relevel(cater_habitat_1719$year, ref="2019")

##################################################################
#### Model: Phenological distribution of mass among tree taxa ####
##################################################################
k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=diag(2), nu=2, alpha.mu=c(0,0), alpha.V=diag(2)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#Mass_yearint_fix<- MCMCglmm(cbind(logmpc1, logmpc2)~ datescaled*year + I(datescaled^2)*year, 
#                        random=~us(1+datescaled):tree.species + us(1+datescaled):site + siteyear + treeID + siteday + recorder + us(sqrt(1/weight)):units, 
#                        family="cengaussian", data=cater_habitat_1719, prior=prior3, pr=TRUE,nitt=5500000, burnin=100000, thin=350)
#save(Mass_yearint_fix, file = "~/Documents/Models/Mass_yearint_fix.RData")
load("~/Documents/Models/Mass_yearint_fix.RData")

#### Checking model fits the data and converged ####
plot(Mass_yearint_fix) #look at fixed effect and random term trace plots 
Mass_yearint_fix.Sim<-simulate(Mass_yearint_fix,nsim=1000) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(Mass_yearint_fix.Sim,2,mean)) #histogram of simulation predictions for mean mass (log scale)
abline(v=mean((cater_habitat_1719$logmpc1+cater_habitat_1719$logmpc2)/2, na.rm=TRUE),col=2) # red line for mean mass (from the mean of logmpc1 and logmpc2)

# tree taxon specific coefficients for date^2 (A), date (B) and intercept (C)
AlderC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,10])
AshC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,11])
BeechC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,12])
BirchC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,13])
ElmC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,14])
HazelC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,15])
OakC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,16])
RowanC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,17])
SycamoreC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,18])
WillowC <- (Mass_yearint_fix$Sol[,1] + Mass_yearint_fix$Sol[,19])

AlderB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,20])
AshB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,21])
BeechB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,22])
BirchB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,23])
ElmB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,24])
HazelB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,25])
OakB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,26])
RowanB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,27])
SycamoreB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,28])
WillowB <- (Mass_yearint_fix$Sol[,2] + Mass_yearint_fix$Sol[,29])

AlderA <- (Mass_yearint_fix$Sol[,5])
AshA <- (Mass_yearint_fix$Sol[,5])
BeechA <- (Mass_yearint_fix$Sol[,5])
BirchA <- (Mass_yearint_fix$Sol[,5])
ElmA <- (Mass_yearint_fix$Sol[,5])
HazelA <- (Mass_yearint_fix$Sol[,5])
OakA <- (Mass_yearint_fix$Sol[,5])
RowanA <- (Mass_yearint_fix$Sol[,5])
SycamoreA <- (Mass_yearint_fix$Sol[,5])
WillowA <- (Mass_yearint_fix$Sol[,5])

# Fixed effect coefficients (2019)
MeanA <- (Mass_yearint_fix$Sol[,5])
MeanB <- (Mass_yearint_fix$Sol[,2])
MeanC <- (Mass_yearint_fix$Sol[,1])

#### Plotting growth curves ####
preddayscaled <- seq(-2.071,2.014,0.001)
predday <- unscale(preddayscaled, center= 146.4095, scale=14.19835)
predday <- predday$V1

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

dataplot <- ggplot(cater_habitat, aes(date, mpc2))+
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
latedate <- scale(168, center= 146.4095, scale=14.19835)[1,1]
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

# median coefficient and CIs for mass of 168 for each tree taxon
Mass168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                      median=c(median(Alder168), median(Ash168), median(Beech168), median(Birch168), median(Elm168), median(Hazel168), median(Oak168), median(Rowan168), median(Sycamore168), median(Willow168)),
                      lowci=c(HPDinterval(Alder168)[1],HPDinterval(Ash168)[1],HPDinterval(Beech168)[1],HPDinterval(Birch168)[1],HPDinterval(Elm168)[1],HPDinterval(Hazel168)[1],HPDinterval(Oak168)[1],HPDinterval(Rowan168)[1],HPDinterval(Sycamore168)[1],HPDinterval(Willow168)[1]),
                      upci=c(HPDinterval(Alder168)[2],HPDinterval(Ash168)[2],HPDinterval(Beech168)[2],HPDinterval(Birch168)[2],HPDinterval(Elm168)[2],HPDinterval(Hazel168)[2],HPDinterval(Oak168)[2],HPDinterval(Rowan168)[2],HPDinterval(Sycamore168)[2],HPDinterval(Willow168)[2]),
                      mediandif=c(median(Alder168dif), median(Ash168dif), median(Beech168dif), median(Birch168dif), median(Elm168dif), median(Hazel168dif), median(Oak168dif), median(Rowan168dif), median(Sycamore168dif), median(Willow168dif)),
                      diflowci=c(HPDinterval(Alder168dif)[1],HPDinterval(Ash168dif)[1],HPDinterval(Beech168dif)[1],HPDinterval(Birch168dif)[1],HPDinterval(Elm168dif)[1],HPDinterval(Hazel168dif)[1],HPDinterval(Oak168dif)[1],HPDinterval(Rowan168dif)[1],HPDinterval(Sycamore168dif)[1],HPDinterval(Willow168dif)[1]),
                      difupci=c(HPDinterval(Alder168dif)[2],HPDinterval(Ash168dif)[2],HPDinterval(Beech168dif)[2],HPDinterval(Birch168dif)[2],HPDinterval(Elm168dif)[2],HPDinterval(Hazel168dif)[2],HPDinterval(Oak168dif)[2],HPDinterval(Rowan168dif)[2],HPDinterval(Sycamore168dif)[2],HPDinterval(Willow168dif)[2]))

# Plot of mass difference to the fixed effect prediction
Mass168difplot <- ggplot(Mass168, aes(fct_rev(TT), mediandif, col=TT))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=difupci, ymin=diflowci, width=0.5))+
  geom_hline(yintercept=0, linetype="longdash", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to fixed eff. (g)")+
  theme_bw()+
  theme(text=element_text(size= 15))+
  scale_colour_manual(values=AllTaxaCols)+
  guides(color = "none") # saved as 6"x4.5"

space <- ggplot() + theme_void()
MassPlots <- grid.arrange(dataplot,space, taxonplot,space, Mass168difplot, ncol = 5, widths = c(1,0.1,1,0.1,1.1)) #saved as 6"x12"


#### For supp mat ####
## Difference to oak on 168 (0.96)
AlderOD <- Alder168-Oak168
AshOD <- Ash168-Oak168
BeechOD <- Beech168-Oak168 
BirchOD <- Birch168-Oak168 
ElmOD <- Elm168-Oak168 
HazelOD <- Hazel168-Oak168 
RowanOD <- Rowan168-Oak168
SycamoreOD <- Sycamore168-Oak168
WillowOD <- Willow168-Oak168

# median and CIs for mass difference to oak on day 168
OD168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                    median=c(median(AlderOD), median(AshOD), median(BeechOD), median(BirchOD), median(ElmOD), median(HazelOD), median(RowanOD), median(SycamoreOD), median(WillowOD)),
                    lowci=c(HPDinterval(AlderOD)[1],HPDinterval(AshOD)[1],HPDinterval(BeechOD)[1],HPDinterval(BirchOD)[1],HPDinterval(ElmOD)[1],HPDinterval(HazelOD)[1],HPDinterval(RowanOD)[1],HPDinterval(SycamoreOD)[1],HPDinterval(WillowOD)[1]),
                    upci=c(HPDinterval(AlderOD)[2],HPDinterval(AshOD)[2],HPDinterval(BeechOD)[2],HPDinterval(BirchOD)[2],HPDinterval(ElmOD)[2],HPDinterval(HazelOD)[2],HPDinterval(RowanOD)[2],HPDinterval(SycamoreOD)[2],HPDinterval(WillowOD)[2]))

# Plot of mass difference to oak prediction
ggplot(OD168, aes(fct_rev(TT), median))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="black", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to oak (g)")+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 15),axis.text.x = element_text(angle = 45, hjust=1)) # saved as 6"x4"

############################
#### Model output table ####   
############################

#### fixed effects
fixed<-rbind(
  c("Intercept (Year 2019)",paste(round(mean(Mass_yearint_fix$Sol[,1]),3)," (",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,1])[1],3)," - ",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,1]))),
  c("Year 2017",paste(round(mean(Mass_yearint_fix$Sol[,3]),3)," (",
                      round(HPDinterval(Mass_yearint_fix$Sol[,3])[1],3)," - ",
                      round(HPDinterval(Mass_yearint_fix$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,3]))),
  c("Year 2018",paste(round(mean(Mass_yearint_fix$Sol[,4]),3)," (",
                      round(HPDinterval(Mass_yearint_fix$Sol[,4])[1],3)," - ",
                      round(HPDinterval(Mass_yearint_fix$Sol[,4])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,4]))),
  c("Date scaled (Year 2019)",paste(round(mean(Mass_yearint_fix$Sol[,2]),3)," (",
                                    round(HPDinterval(Mass_yearint_fix$Sol[,2])[1],3)," - ",
                                    round(HPDinterval(Mass_yearint_fix$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,2]))),
  c("Date scaled:Year 2017",paste(round(mean(Mass_yearint_fix$Sol[,6]),3)," (",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,6])[1],3)," - ",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,6])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,6]))),
  c("Date scaled:Year 2018",paste(round(mean(Mass_yearint_fix$Sol[,7]),3)," (",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,7])[1],3)," - ",
                                  round(HPDinterval(Mass_yearint_fix$Sol[,7])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,7]))),
  c("Date² scaled (Year 2019)",paste(round(mean(Mass_yearint_fix$Sol[,5]),3)," (",
                                     round(HPDinterval(Mass_yearint_fix$Sol[,5])[1],3)," - ",
                                     round(HPDinterval(Mass_yearint_fix$Sol[,5])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,5]))),
  c("Date² scaled:Year 2017",paste(round(mean(Mass_yearint_fix$Sol[,8]),3)," (",
                                   round(HPDinterval(Mass_yearint_fix$Sol[,8])[1],3)," - ",
                                   round(HPDinterval(Mass_yearint_fix$Sol[,8])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,8]))),
  c("Date² scaled:Year 2018",paste(round(mean(Mass_yearint_fix$Sol[,9]),3)," (",
                                   round(HPDinterval(Mass_yearint_fix$Sol[,9])[1],3)," - ",
                                   round(HPDinterval(Mass_yearint_fix$Sol[,9])[2],3),")",sep=""),
    round(effectiveSize(Mass_yearint_fix$Sol[,9]))))

#### random terms 
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass_yearint_fix$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                                          round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-4
treetaxa4<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                              round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-5
site5<-c("Site- Intercept var",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                     round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                     round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass_yearint_fix$VCV[, column])))

column<-6
site6<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                                  round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                                  round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-8
site8<-c("Site- Date slope var",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                                      round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
         round(effectiveSize(Mass_yearint_fix$VCV[, column])))

column<-9
siteyear<-c("Site-Year",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                              round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                              round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-12
recorder<-c("Recorder",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-11
siteday<-c("Site-Day",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                            round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                            round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(Mass_yearint_fix$VCV[, column])))


column<-10
treeID<-c("Tree ID",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                          round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                          round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(Mass_yearint_fix$VCV[, column])))

column<-13
weight<-c("Weighting",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                            round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                            round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(Mass_yearint_fix$VCV[, column])))

column<-14
residual<-c("Residual",paste(round(posterior.mode(Mass_yearint_fix$VCV[, column]),3)," (",
                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[1],3)," - ",
                             round(HPDinterval(Mass_yearint_fix$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(Mass_yearint_fix$VCV[, column])))




random<-rbind(treetaxa1,treetaxa2,treetaxa4,site5,site6,site8,siteyear, recorder, siteday, treeID, weight, residual)


write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableMass_yearint_fix.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
