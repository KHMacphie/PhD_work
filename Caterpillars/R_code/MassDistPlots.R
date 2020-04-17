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

########################################
#### Mass phenological distribution ####
########################################
load("~/Documents/Models/Mass_yearint_fix.RData")

# Also need cater_habitat data frame

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
mycolblack <- rgb(0, 0, 0, max = 250, alpha = 25, names = "blacktrans")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

par(mfrow=c(1,2), cex=1.1)
plot(intsamples$date, intsamples$mpc1, xlab="Ordinal Date", ylab="Mass (g)", pch=20, col=mycolblack, cex=0.4, xlim=c(117,172), ylim=c(0,1))
points(cater_habitat$date, exp(cater_habitat$logmpc2), col=mycolblack, pch=20, cex=0.4)
points(predday,exp(meanslope), type="l", lwd=2, col=2)
legend("topleft", legend="Fixed effect curve",
       lty=1, lwd=3, 
       col=2, 
       cex=0.8, seg.len=1.2, bty = "n")
plot(predday,exp(alder), type="l", lwd=1.5, col="darkred", xlab="Ordinal Date", ylab="Mass (g)", ylim=c(0,0.05))
points(predday,exp(ash), type="l", lwd=1.5, col="firebrick3")
points(predday,exp(beech), type="l", lwd=1.5, col="chocolate2")
points(predday,exp(birch), type="l", lwd=1.5, col="goldenrod")
points(predday,exp(elm), type="l", lwd=1.5, col="olivedrab4")
points(predday,exp(hazel), type="l", lwd=1.5, col="darkgreen")
points(predday,exp(oak), type="l", lwd=1.5, col="deepskyblue3")
points(predday,exp(rowan), type="l", lwd=1.5, col="royalblue4")
points(predday,exp(sycamore), type="l", lwd=1.5, col="slateblue2")
points(predday,exp(willow), type="l", lwd=1.5, col="orchid")
points(predday,exp(meanslope), type="l", lwd=2, col=2, lty="dotted")
points(rep(168, 101), seq(0,0.05, 0.0005), type="l", col=1, lwd=0.5, lty="dashed")
legend("topleft", legend="Fixed effect curve",
       lty=3, lwd=3, 
       col=2, 
       cex=0.8, seg.len=1.2, bty = "n") #saved as 7"x8"

legend("topleft", legend=c("Fixed effect","Alder","Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"),
       lty=c(3,1,1,1,1,1,1,1,1,1,1), lwd=2, 
       col=c(2,"darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid"), 
       cex=0.8, seg.len=0.8, bty = "n") #saved as 7"x8"

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

Mass168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Oak", "Rowan","Sycamore","Willow"),
                      median=c(median(Alder168), median(Ash168), median(Beech168), median(Birch168), median(Elm168), median(Hazel168), median(Oak168), median(Rowan168), median(Sycamore168), median(Willow168)),
                      lowci=c(HPDinterval(Alder168)[1],HPDinterval(Ash168)[1],HPDinterval(Beech168)[1],HPDinterval(Birch168)[1],HPDinterval(Elm168)[1],HPDinterval(Hazel168)[1],HPDinterval(Oak168)[1],HPDinterval(Rowan168)[1],HPDinterval(Sycamore168)[1],HPDinterval(Willow168)[1]),
                      upci=c(HPDinterval(Alder168)[2],HPDinterval(Ash168)[2],HPDinterval(Beech168)[2],HPDinterval(Birch168)[2],HPDinterval(Elm168)[2],HPDinterval(Hazel168)[2],HPDinterval(Oak168)[2],HPDinterval(Rowan168)[2],HPDinterval(Sycamore168)[2],HPDinterval(Willow168)[2]),
                      mediandif=c(median(Alder168dif), median(Ash168dif), median(Beech168dif), median(Birch168dif), median(Elm168dif), median(Hazel168dif), median(Oak168dif), median(Rowan168dif), median(Sycamore168dif), median(Willow168dif)),
                      diflowci=c(HPDinterval(Alder168dif)[1],HPDinterval(Ash168dif)[1],HPDinterval(Beech168dif)[1],HPDinterval(Birch168dif)[1],HPDinterval(Elm168dif)[1],HPDinterval(Hazel168dif)[1],HPDinterval(Oak168dif)[1],HPDinterval(Rowan168dif)[1],HPDinterval(Sycamore168dif)[1],HPDinterval(Willow168dif)[1]),
                      difupci=c(HPDinterval(Alder168dif)[2],HPDinterval(Ash168dif)[2],HPDinterval(Beech168dif)[2],HPDinterval(Birch168dif)[2],HPDinterval(Elm168dif)[2],HPDinterval(Hazel168dif)[2],HPDinterval(Oak168dif)[2],HPDinterval(Rowan168dif)[2],HPDinterval(Sycamore168dif)[2],HPDinterval(Willow168dif)[2]))

#write.csv(Mass168,'~/Documents/Models/Tables/Mass168Metrics.csv')

Mass168plot <- ggplot(Mass168, aes(fct_rev(TT), median, col=TT))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass (g)")+
  theme_bw()+
  #labs(tag = "Day168")+
  theme(text=element_text(size= 20))+
  scale_colour_manual(values=AllTaxaCols)+
  guides(color = "none")

Mass168difplot <- ggplot(Mass168, aes(fct_rev(TT), mediandif, col=TT))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=difupci, ymin=diflowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to fixed eff. (g)")+
  theme_bw()+
  #labs(tag = "Day168")+
  theme(text=element_text(size= 15))+
  scale_colour_manual(values=AllTaxaCols)+
  guides(color = "none") # saved as 6"x4.5"

#col1 <- grid.arrange(Mass168plot, nrow = 1, heights = 1)
#col2 <- grid.arrange(Mass168difplot, nrow = 1, heights = 1)
#MassDif168 <- grid.arrange(col1, col2, ncol = 2, widths = c(1,1)) #saved as 7"by 9.5"

#row1 <- grid.arrange(Mass168plot, ncol = 1)
#row2 <- grid.arrange(Mass168difplot, ncol = 1)
#MassDif168 <- grid.arrange(row1, row2, nrow = 2) #saved as 7"x5"
#### Difference to oak on 168 (0.96) #### 

AlderOD <- Alder168-Oak168
AshOD <- Ash168-Oak168
BeechOD <- Beech168-Oak168 
BirchOD <- Birch168-Oak168 
ElmOD <- Elm168-Oak168 
HazelOD <- Hazel168-Oak168 
RowanOD <- Rowan168-Oak168
SycamoreOD <- Sycamore168-Oak168
WillowOD <- Willow168-Oak168

OD168 <- data.frame(TT=c("Alder","Ash","Beech","Birch","Elm","Hazel","Rowan","Sycamore","Willow"),
                    median=c(median(AlderOD), median(AshOD), median(BeechOD), median(BirchOD), median(ElmOD), median(HazelOD), median(RowanOD), median(SycamoreOD), median(WillowOD)),
                    lowci=c(HPDinterval(AlderOD)[1],HPDinterval(AshOD)[1],HPDinterval(BeechOD)[1],HPDinterval(BirchOD)[1],HPDinterval(ElmOD)[1],HPDinterval(HazelOD)[1],HPDinterval(RowanOD)[1],HPDinterval(SycamoreOD)[1],HPDinterval(WillowOD)[1]),
                    upci=c(HPDinterval(AlderOD)[2],HPDinterval(AshOD)[2],HPDinterval(BeechOD)[2],HPDinterval(BirchOD)[2],HPDinterval(ElmOD)[2],HPDinterval(HazelOD)[2],HPDinterval(RowanOD)[2],HPDinterval(SycamoreOD)[2],HPDinterval(WillowOD)[2]))

#write.csv(OD168,'~/Documents/Models/Tables/MassOD168Metrics.csv')

MassODplot <- ggplot(OD168, aes(fct_rev(TT), median))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=upci, ymin=lowci, width=0.5))+
  geom_hline(yintercept=0, linetype="dashed", colour="red", size=0.3)+  
  coord_flip()+
  xlab("Tree Taxon")+
  ylab("Mass difference to oak (g)")+
  #annotate("text", label="*", x="Beech", y=-0.0018, size=7)+
  #annotate("text", label="*", x="Birch", y=-0.0006, size=7)+
  #annotate("text", label="*", x="Willow", y=-0.0018, size=7)+
  theme_bw()+
  labs(tag = "  ")+
  theme(text=element_text(size= 15),axis.text.x = element_text(angle = 45, hjust=1)) # saved as 6"x4"
