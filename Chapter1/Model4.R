#load Dataframe.R script

#######################################################################
#### Model: Phenological distribution of abundance among tree taxa ####
#######################################################################

k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

#ATTCyearint<- MCMCglmm(caterpillars~ datescaled + I(datescaled^2), 
#                         random=~us(1+datescaled+I(datescaled^2)):tree.species + us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + siteyear + treeID + siteday + recorder, 
#                         family="poisson", data=cater_habitat, prior=prior3, nitt=1000000, burnin=50000, pr=TRUE, thin=30)
#save(ATTCyearint, file = "~/Documents/Models/ATTCyearint.RData")
load("~/Documents/Models/ATTCyearint.RData")

#### Checking model fits the data and converged ####
plot(ATTCyearint) #look at fixed effect and random term trace plots 
ATTCyearint.Sim<-simulate(ATTCyearint,nsim=1000) #simulate 1000 times
par(mfcol=c(1,1))
hist(apply(ATTCyearint.Sim,2,sum), breaks=1000000) #histogram of simulation predictions for total abundance
#normally not too many rogue values but reduce x axis a bit to see main distribution relative to observed value 
hist(apply(ATTCyearint.Sim,2,sum), breaks=1000000, xlim=c(0,100000))
abline(v=sum(cater_habitat$caterpillars),col=2) # red line for observed value in data

propzero <- function(x){return(length(which(x==0))/length(x))} # function for proportion of zeros
hist(apply(ATTCyearint.Sim,2,propzero), breaks=100) # histogram of proportion of zeros in simulated data
abline(v=propzero(cater_habitat$caterpillars), col="red") # red line for observed proportion in data


# tree taxon specific coefficients for date^2 (A), date (B) and intercept (C)
AlderC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,4])
AshC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,5])
BeechC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,6])
BirchC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,7])
ElmC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,8])
HazelC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,9])
OakC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,10])
RowanC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,11])
SycamoreC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,12])
WillowC <- (ATTCyearint$Sol[,1] + ATTCyearint$Sol[,13])

AlderB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,14])
AshB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,15])
BeechB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,16])
BirchB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,17])
ElmB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,18])
HazelB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,19])
OakB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,20])
RowanB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,21])
SycamoreB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,22])
WillowB <- (ATTCyearint$Sol[,2] + ATTCyearint$Sol[,23])

AlderA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,24])
AshA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,25])
BeechA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,26])
BirchA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,27])
ElmA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,28])
HazelA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,29])
OakA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,30])
RowanA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,31])
SycamoreA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,32])
WillowA <- (ATTCyearint$Sol[,3] + ATTCyearint$Sol[,33])

#fixed effect coefficients
MeanA <- (ATTCyearint$Sol[,3])
MeanB <- (ATTCyearint$Sol[,2])
MeanC <- (ATTCyearint$Sol[,1])

# unscaled date
dayscal <- seq(-2.071,2.014,0.001)
days <- unscale(dayscal, center= 146.4095, scale=14.19835)
days <- days$V1

#median curve for each tree taxon and the average trend
AlderCurve <- exp(median(AlderC) +median(AlderB)*dayscal + median(AlderA)*dayscal^2)
AshCurve <- exp(median(AshC) +median(AshB)*dayscal + median(AshA)*dayscal^2)
BeechCurve <- exp(median(BeechC) +median(BeechB)*dayscal + median(BeechA)*dayscal^2)
BirchCurve <- exp(median(BirchC) +median(BirchB)*dayscal + median(BirchA)*dayscal^2)
ElmCurve <- exp(median(ElmC) +median(ElmB)*dayscal + median(ElmA)*dayscal^2)
HazelCurve <- exp(median(HazelC) +median(HazelB)*dayscal + median(HazelA)*dayscal^2)
OakCurve <- exp(median(OakC) +median(OakB)*dayscal + median(OakA)*dayscal^2)
RowanCurve <- exp(median(RowanC) +median(RowanB)*dayscal + median(RowanA)*dayscal^2)
SycamoreCurve <- exp(median(SycamoreC) +median(SycamoreB)*dayscal + median(SycamoreA)*dayscal^2)
WillowCurve <- exp(median(WillowC) +median(WillowB)*dayscal + median(WillowA)*dayscal^2)
MeanCurve <- exp(median(MeanC) +median(MeanB)*dayscal + median(MeanA)*dayscal^2)

#### Peak date per taxa ####   -b/2a
MeanPD <- unscale((-MeanB/(2*MeanA)), center= 146.4095, scale=14.19835)
MeanPD <- MeanPD$var1
AlderPD <-  unscale((-AlderB/(2*AlderA)), center= 146.4095, scale=14.19835)
AlderPD <- AlderPD$var1
AshPD <-  unscale((-AshB/(2*AshA)), center= 146.4095, scale=14.19835)
AshPD <-  AshPD$var1
BeechPD <-  unscale((-BeechB/(2*BeechA)), center= 146.4095, scale=14.19835)
BeechPD <- BeechPD$var1
BirchPD <-  unscale((-BirchB/(2*BirchA)), center= 146.4095, scale=14.19835)
BirchPD <- BirchPD$var1
ElmPD <-  unscale((-ElmB/(2*ElmA)), center= 146.4095, scale=14.19835)
ElmPD <- ElmPD$var1
HazelPD <-  unscale((-HazelB/(2*HazelA)), center= 146.4095, scale=14.19835)
HazelPD <- HazelPD$var1
OakPD <-  unscale((-OakB/(2*OakA)), center= 146.4095, scale=14.19835)
OakPD <- OakPD$var1
RowanPD <-  unscale((-RowanB/(2*RowanA)), center= 146.4095, scale=14.19835)
RowanPD <- RowanPD$var1
SycamorePD <-  unscale((-SycamoreB/(2*SycamoreA)), center= 146.4095, scale=14.19835)
SycamorePD <- SycamorePD$var1
WillowPD <-  unscale((-WillowB/(2*WillowA)), center= 146.4095, scale=14.19835)
WillowPD <- WillowPD$var1

## Peak date diff to mean
AlderPDdif <- AlderPD-MeanPD
AshPDdif <- AshPD-MeanPD
BeechPDdif <- BeechPD-MeanPD
BirchPDdif <- BirchPD-MeanPD
ElmPDdif <- ElmPD-MeanPD
HazelPDdif <- HazelPD-MeanPD
OakPDdif <- OakPD-MeanPD
RowanPDdif <- RowanPD-MeanPD
SycamorePDdif <- SycamorePD-MeanPD
WillowPDdif <- WillowPD-MeanPD

## Peak date diff to oak
AlderPDOD <- AlderPD-OakPD
AshPDOD <- AshPD-OakPD
BeechPDOD <- BeechPD-OakPD
BirchPDOD <- BirchPD-OakPD
ElmPDOD <- ElmPD-OakPD
HazelPDOD <- HazelPD-OakPD
RowanPDOD <- RowanPD-OakPD
SycamorePDOD <- SycamorePD-OakPD
WillowPDOD <- WillowPD-OakPD

#### Peak height per taxa ####
MeanPH <- exp(MeanC +MeanB*(-MeanB/(2*MeanA)) + MeanA*(-MeanB/(2*MeanA))^2)
AlderPH <- exp(AlderC +AlderB*(-AlderB/(2*AlderA)) + AlderA*(-AlderB/(2*AlderA))^2)
AshPH <- exp(AshC +AshB*(-AshB/(2*AshA)) + AshA*(-AshB/(2*AshA))^2)
BeechPH <- exp(BeechC +BeechB*(-BeechB/(2*BeechA)) + BeechA*(-BeechB/(2*BeechA))^2)
BirchPH <- exp(BirchC +BirchB*(-BirchB/(2*BirchA)) + BirchA*(-BirchB/(2*BirchA))^2)
ElmPH <- exp(ElmC +ElmB*(-ElmB/(2*ElmA)) + ElmA*(-ElmB/(2*ElmA))^2)
HazelPH <- exp(HazelC +HazelB*(-HazelB/(2*HazelA)) + HazelA*(-HazelB/(2*HazelA))^2)
OakPH <- exp(OakC +OakB*(-OakB/(2*OakA)) + OakA*(-OakB/(2*OakA))^2)
RowanPH <- exp(RowanC +RowanB*(-RowanB/(2*RowanA)) + RowanA*(-RowanB/(2*RowanA))^2)
SycamorePH <- exp(SycamoreC +SycamoreB*(-SycamoreB/(2*SycamoreA)) + SycamoreA*(-SycamoreB/(2*SycamoreA))^2)
WillowPH <- exp(WillowC +WillowB*(-WillowB/(2*WillowA)) + WillowA*(-WillowB/(2*WillowA))^2)

## Peak height diff to mean
AlderPHdif <- AlderPH-MeanPH
AshPHdif <- AshPH-MeanPH
BeechPHdif <- BeechPH-MeanPH
BirchPHdif <- BirchPH-MeanPH
ElmPHdif <- ElmPH-MeanPH
HazelPHdif <- HazelPH-MeanPH
OakPHdif <- OakPH-MeanPH
RowanPHdif <- RowanPH-MeanPH
SycamorePHdif <- SycamorePH-MeanPH
WillowPHdif <- WillowPH-MeanPH

## Peak height diff to oak
AlderPHOD <- AlderPH-OakPH
AshPHOD <- AshPH-OakPH
BeechPHOD <- BeechPH-OakPH
BirchPHOD <- BirchPH-OakPH
ElmPHOD <- ElmPH-OakPH
HazelPHOD <- HazelPH-OakPH
RowanPHOD <- RowanPH-OakPH
SycamorePHOD <- SycamorePH-OakPH
WillowPHOD <- WillowPH-OakPH

## Peak height proportion of mean
AlderPHprop <- AlderPH/MeanPH
AshPHprop <- AshPH/MeanPH
BeechPHprop <- BeechPH/MeanPH
BirchPHprop <- BirchPH/MeanPH
ElmPHprop <- ElmPH/MeanPH
HazelPHprop <- HazelPH/MeanPH
OakPHprop <- OakPH/MeanPH
RowanPHprop <- RowanPH/MeanPH
SycamorePHprop <- SycamorePH/MeanPH
WillowPHprop <- WillowPH/MeanPH

## Peak height proportion of oak
AlderPHODprop <- AlderPH/OakPH
AshPHODprop <- AshPH/OakPH
BeechPHODprop <- BeechPH/OakPH
BirchPHODprop <- BirchPH/OakPH
ElmPHODprop <- ElmPH/OakPH
HazelPHODprop <- HazelPH/OakPH 
RowanPHODprop <- RowanPH/OakPH
SycamorePHODprop <- SycamorePH/OakPH
WillowPHODprop <- WillowPH/OakPH

#### Peak width per taxa (set height) ####   x= (-b +/- sqrt(b^2 - 4ac))/2a

# set at 0.01 - roughly half the height of the lowest curve (Alder)

MeanPW1 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.01)))))/(2*MeanA)
MeanPW2 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(0.01)))))/(2*MeanA)
MeanPW <- unscale((MeanPW2), center= 146.4095, scale=14.19835)-unscale((MeanPW1), center= 146.4095, scale=14.19835)
MeanPW <- MeanPW$var1 #NaNs=206
#length(which(is.na(MeanPW)==TRUE))

AlderPW1 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.01)))))/(2*AlderA)
AlderPW2 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(0.01)))))/(2*AlderA)
AlderPW <- unscale((AlderPW2), center= 146.4095, scale=14.19835)-unscale((AlderPW1), center= 146.4095, scale=14.19835)
AlderPW <- AlderPW$var1 #NaNs=942

AshPW1 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(0.01)))))/(2*AshA)
AshPW2 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(0.01)))))/(2*AshA)
AshPW <- unscale((AshPW2), center= 146.4095, scale=14.19835)-unscale((AshPW1), center= 146.4095, scale=14.19835)
AshPW <- AshPW$var1 #NaNs=583

BeechPW1 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.01)))))/(2*BeechA)
BeechPW2 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(0.01)))))/(2*BeechA)
BeechPW <- unscale((BeechPW2), center= 146.4095, scale=14.19835)-unscale((BeechPW1), center= 146.4095, scale=14.19835)
BeechPW <- BeechPW$var1 #NaNs=273

BirchPW1 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.01)))))/(2*BirchA)
BirchPW2 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(0.01)))))/(2*BirchA)
BirchPW <- unscale((BirchPW2), center= 146.4095, scale=14.19835)-unscale((BirchPW1), center= 146.4095, scale=14.19835)
BirchPW <- BirchPW$var1 #NaNs=113

ElmPW1 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.01)))))/(2*ElmA)
ElmPW2 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(0.01)))))/(2*ElmA)
ElmPW <- unscale((ElmPW2), center= 146.4095, scale=14.19835)-unscale((ElmPW1), center= 146.4095, scale=14.19835)
ElmPW <- ElmPW$var1 #NaNs=246

HazelPW1 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.01)))))/(2*HazelA)
HazelPW2 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(0.01)))))/(2*HazelA)
HazelPW <- unscale((HazelPW2), center= 146.4095, scale=14.19835)-unscale((HazelPW1), center= 146.4095, scale=14.19835)
HazelPW <- HazelPW$var1 #NaNs=199

OakPW1 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(0.01)))))/(2*OakA)
OakPW2 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(0.01)))))/(2*OakA)
OakPW <- unscale((OakPW2), center= 146.4095, scale=14.19835)-unscale((OakPW1), center= 146.4095, scale=14.19835)
OakPW <- OakPW$var1 #NaNs=55

RowanPW1 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.01)))))/(2*RowanA)
RowanPW2 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(0.01)))))/(2*RowanA)
RowanPW <- unscale((RowanPW2), center= 146.4095, scale=14.19835)-unscale((RowanPW1), center= 146.4095, scale=14.19835)
RowanPW <- RowanPW$var1 #NaNs=335

SycamorePW1 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.01)))))/(2*SycamoreA)
SycamorePW2 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(0.01)))))/(2*SycamoreA)
SycamorePW <- unscale((SycamorePW2), center= 146.4095, scale=14.19835)-unscale((SycamorePW1), center= 146.4095, scale=14.19835)
SycamorePW <- SycamorePW$var1 #NaNs=180

WillowPW1 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.01)))))/(2*WillowA)
WillowPW2 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(0.01)))))/(2*WillowA)
WillowPW <- unscale((WillowPW2), center= 146.4095, scale=14.19835)-unscale((WillowPW1), center= 146.4095, scale=14.19835)
WillowPW <- WillowPW$var1 #NaNs=91

## Peak width diff to mean
AlderPWdif <- AlderPW-MeanPW
AshPWdif <- AshPW-MeanPW
BeechPWdif <- BeechPW-MeanPW
BirchPWdif <- BirchPW-MeanPW
ElmPWdif <- ElmPW-MeanPW
HazelPWdif <- HazelPW-MeanPW
OakPWdif <- OakPW-MeanPW
RowanPWdif <- RowanPW-MeanPW
SycamorePWdif <- SycamorePW-MeanPW
WillowPWdif <- WillowPW-MeanPW

## Peak width diff to oak
AlderPWOD <- AlderPW-OakPW
AshPWOD <- AshPW-OakPW
BeechPWOD <- BeechPW-OakPW
BirchPWOD <- BirchPW-OakPW
ElmPWOD <- ElmPW-OakPW
HazelPWOD <- HazelPW-OakPW
RowanPWOD <- RowanPW-OakPW
SycamorePWOD <- SycamorePW-OakPW
WillowPWOD <- WillowPW-OakPW

#### Peak width per taxa (half each PH) ####   x= (-b +/- sqrt(b^2 - 4ac))/2a

# at half the height of each peak

MeanPW1.5 <- (-MeanB + sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW2.5 <- (-MeanB - sqrt(MeanB^2 - (4*MeanA*(MeanC-log(MeanPH/2)))))/(2*MeanA)
MeanPW.5 <- unscale((MeanPW2.5), center= 146.4095, scale=14.19835)-unscale((MeanPW1.5), center= 146.4095, scale=14.19835)
MeanPW.5 <- MeanPW.5$var1 #NaNs=23

AlderPW1.5 <- (-AlderB + sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW2.5 <- (-AlderB - sqrt(AlderB^2 - (4*AlderA*(AlderC-log(AlderPH/2)))))/(2*AlderA)
AlderPW.5 <- unscale((AlderPW2.5), center= 146.4095, scale=14.19835)-unscale((AlderPW1.5), center= 146.4095, scale=14.19835)
AlderPW.5 <- AlderPW.5$var1 #NaNs=75

AshPW1.5 <- (-AshB + sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW2.5 <- (-AshB - sqrt(AshB^2 - (4*AshA*(AshC-log(AshPH/2)))))/(2*AshA)
AshPW.5 <- unscale((AshPW2.5), center= 146.4095, scale=14.19835)-unscale((AshPW1.5), center= 146.4095, scale=14.19835)
AshPW.5 <- AshPW.5$var1 #NaNs=25

BeechPW1.5 <- (-BeechB + sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW2.5 <- (-BeechB - sqrt(BeechB^2 - (4*BeechA*(BeechC-log(BeechPH/2)))))/(2*BeechA)
BeechPW.5 <- unscale((BeechPW2.5), center= 146.4095, scale=14.19835)-unscale((BeechPW1.5), center= 146.4095, scale=14.19835)
BeechPW.5 <- BeechPW.5$var1 #NaNs=18

BirchPW1.5 <- (-BirchB + sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW2.5 <- (-BirchB - sqrt(BirchB^2 - (4*BirchA*(BirchC-log(BirchPH/2)))))/(2*BirchA)
BirchPW.5 <- unscale((BirchPW2.5), center= 146.4095, scale=14.19835)-unscale((BirchPW1.5), center= 146.4095, scale=14.19835)
BirchPW.5 <- BirchPW.5$var1 #NaNs=31

ElmPW1.5 <- (-ElmB + sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW2.5 <- (-ElmB - sqrt(ElmB^2 - (4*ElmA*(ElmC-log(ElmPH/2)))))/(2*ElmA)
ElmPW.5 <- unscale((ElmPW2.5), center= 146.4095, scale=14.19835)-unscale((ElmPW1.5), center= 146.4095, scale=14.19835)
ElmPW.5 <- ElmPW.5$var1 #NaNs=41

HazelPW1.5 <- (-HazelB + sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW2.5 <- (-HazelB - sqrt(HazelB^2 - (4*HazelA*(HazelC-log(HazelPH/2)))))/(2*HazelA)
HazelPW.5 <- unscale((HazelPW2.5), center= 146.4095, scale=14.19835)-unscale((HazelPW1.5), center= 146.4095, scale=14.19835)
HazelPW.5 <- HazelPW.5$var1 #NaNs=197

OakPW1.5 <- (-OakB + sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW2.5 <- (-OakB - sqrt(OakB^2 - (4*OakA*(OakC-log(OakPH/2)))))/(2*OakA)
OakPW.5 <- unscale((OakPW2.5), center= 146.4095, scale=14.19835)-unscale((OakPW1.5), center= 146.4095, scale=14.19835)
OakPW.5 <- OakPW.5$var1 #NaNs=14

RowanPW1.5 <- (-RowanB + sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW2.5 <- (-RowanB - sqrt(RowanB^2 - (4*RowanA*(RowanC-log(RowanPH/2)))))/(2*RowanA)
RowanPW.5 <- unscale((RowanPW2.5), center= 146.4095, scale=14.19835)-unscale((RowanPW1.5), center= 146.4095, scale=14.19835)
RowanPW.5 <- RowanPW.5$var1 #NaNs=11

SycamorePW1.5 <- (-SycamoreB + sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW2.5 <- (-SycamoreB - sqrt(SycamoreB^2 - (4*SycamoreA*(SycamoreC-log(SycamorePH/2)))))/(2*SycamoreA)
SycamorePW.5 <- unscale((SycamorePW2.5), center= 146.4095, scale=14.19835)-unscale((SycamorePW1.5), center= 146.4095, scale=14.19835)
SycamorePW.5 <- SycamorePW.5$var1 #NaNs=16

WillowPW1.5 <- (-WillowB + sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW2.5 <- (-WillowB - sqrt(WillowB^2 - (4*WillowA*(WillowC-log(WillowPH/2)))))/(2*WillowA)
WillowPW.5 <- unscale((WillowPW2.5), center= 146.4095, scale=14.19835)-unscale((WillowPW1.5), center= 146.4095, scale=14.19835)
WillowPW.5 <- WillowPW.5$var1 #NaNs=41

## Peak width diff to mean
AlderPW.5dif <- AlderPW.5-MeanPW.5
AshPW.5dif <- AshPW.5-MeanPW.5
BeechPW.5dif <- BeechPW.5-MeanPW.5
BirchPW.5dif <- BirchPW.5-MeanPW.5
ElmPW.5dif <- ElmPW.5-MeanPW.5
HazelPW.5dif <- HazelPW.5-MeanPW.5
OakPW.5dif <- OakPW.5-MeanPW.5
RowanPW.5dif <- RowanPW.5-MeanPW.5
SycamorePW.5dif <- SycamorePW.5-MeanPW.5
WillowPW.5dif <- WillowPW.5-MeanPW.5

## Peak width diff to oak
AlderPW.5OD <- AlderPW.5-OakPW.5
AshPW.5OD <- AshPW.5-OakPW.5
BeechPW.5OD <- BeechPW.5-OakPW.5
BirchPW.5OD <- BirchPW.5-OakPW.5
ElmPW.5OD <- ElmPW.5-OakPW.5
HazelPW.5OD <- HazelPW.5-OakPW.5
RowanPW.5OD <- RowanPW.5-OakPW.5
SycamorePW.5OD <- SycamorePW.5-OakPW.5
WillowPW.5OD <- WillowPW.5-OakPW.5


#### Data frame of PDHW median and CIs ####   

AbundCurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                          PD=c(median(AlderPD),median(AshPD),median(BeechPD),median(BirchPD),median(ElmPD),median(HazelPD),median(OakPD),median(RowanPD),median(SycamorePD),median(WillowPD)),
                          PDLCI=c(HPDinterval(mcmc(AlderPD))[1], HPDinterval(mcmc(AshPD))[1], HPDinterval(mcmc(BeechPD))[1], HPDinterval(mcmc(BirchPD))[1], HPDinterval(mcmc(ElmPD))[1], HPDinterval(mcmc(HazelPD))[1], HPDinterval(mcmc(OakPD))[1], HPDinterval(mcmc(RowanPD))[1], HPDinterval(mcmc(SycamorePD))[1], HPDinterval(mcmc(WillowPD))[1]),
                          PDUCI=c(HPDinterval(mcmc(AlderPD))[2], HPDinterval(mcmc(AshPD))[2], HPDinterval(mcmc(BeechPD))[2], HPDinterval(mcmc(BirchPD))[2], HPDinterval(mcmc(ElmPD))[2], HPDinterval(mcmc(HazelPD))[2], HPDinterval(mcmc(OakPD))[2], HPDinterval(mcmc(RowanPD))[2], HPDinterval(mcmc(SycamorePD))[2], HPDinterval(mcmc(WillowPD))[2]),
                          PDdif=c(median(AlderPDdif),median(AshPDdif),median(BeechPDdif),median(BirchPDdif),median(ElmPDdif),median(HazelPDdif),median(OakPDdif),median(RowanPDdif),median(SycamorePDdif),median(WillowPDdif)),
                          PDdifLCI=c(HPDinterval(mcmc(AlderPDdif))[1], HPDinterval(mcmc(AshPDdif))[1], HPDinterval(mcmc(BeechPDdif))[1], HPDinterval(mcmc(BirchPDdif))[1], HPDinterval(mcmc(ElmPDdif))[1], HPDinterval(mcmc(HazelPDdif))[1], HPDinterval(mcmc(OakPDdif))[1], HPDinterval(mcmc(RowanPDdif))[1], HPDinterval(mcmc(SycamorePDdif))[1], HPDinterval(mcmc(WillowPDdif))[1]),
                          PDdifUCI=c(HPDinterval(mcmc(AlderPDdif))[2], HPDinterval(mcmc(AshPDdif))[2], HPDinterval(mcmc(BeechPDdif))[2], HPDinterval(mcmc(BirchPDdif))[2], HPDinterval(mcmc(ElmPDdif))[2], HPDinterval(mcmc(HazelPDdif))[2], HPDinterval(mcmc(OakPDdif))[2], HPDinterval(mcmc(RowanPDdif))[2], HPDinterval(mcmc(SycamorePDdif))[2], HPDinterval(mcmc(WillowPDdif))[2]),
                          PH=c(median(AlderPH),median(AshPH),median(BeechPH),median(BirchPH),median(ElmPH),median(HazelPH),median(OakPH),median(RowanPH),median(SycamorePH),median(WillowPH)),
                          PHLCI=c(HPDinterval(mcmc(AlderPH))[1], HPDinterval(mcmc(AshPH))[1], HPDinterval(mcmc(BeechPH))[1], HPDinterval(mcmc(BirchPH))[1], HPDinterval(mcmc(ElmPH))[1], HPDinterval(mcmc(HazelPH))[1], HPDinterval(mcmc(OakPH))[1], HPDinterval(mcmc(RowanPH))[1], HPDinterval(mcmc(SycamorePH))[1], HPDinterval(mcmc(WillowPH))[1]),
                          PHUCI=c(HPDinterval(mcmc(AlderPH))[2], HPDinterval(mcmc(AshPH))[2], HPDinterval(mcmc(BeechPH))[2], HPDinterval(mcmc(BirchPH))[2], HPDinterval(mcmc(ElmPH))[2], HPDinterval(mcmc(HazelPH))[2], HPDinterval(mcmc(OakPH))[2], HPDinterval(mcmc(RowanPH))[2], HPDinterval(mcmc(SycamorePH))[2], HPDinterval(mcmc(WillowPH))[2]),
                          PHdif=c(median(AlderPHdif),median(AshPHdif),median(BeechPHdif),median(BirchPHdif),median(ElmPHdif),median(HazelPHdif),median(OakPHdif),median(RowanPHdif),median(SycamorePHdif),median(WillowPHdif)),
                          PHdifLCI=c(HPDinterval(mcmc(AlderPHdif))[1], HPDinterval(mcmc(AshPHdif))[1], HPDinterval(mcmc(BeechPHdif))[1], HPDinterval(mcmc(BirchPHdif))[1], HPDinterval(mcmc(ElmPHdif))[1], HPDinterval(mcmc(HazelPHdif))[1], HPDinterval(mcmc(OakPHdif))[1], HPDinterval(mcmc(RowanPHdif))[1], HPDinterval(mcmc(SycamorePHdif))[1], HPDinterval(mcmc(WillowPHdif))[1]),
                          PHdifUCI=c(HPDinterval(mcmc(AlderPHdif))[2], HPDinterval(mcmc(AshPHdif))[2], HPDinterval(mcmc(BeechPHdif))[2], HPDinterval(mcmc(BirchPHdif))[2], HPDinterval(mcmc(ElmPHdif))[2], HPDinterval(mcmc(HazelPHdif))[2], HPDinterval(mcmc(OakPHdif))[2], HPDinterval(mcmc(RowanPHdif))[2], HPDinterval(mcmc(SycamorePHdif))[2], HPDinterval(mcmc(WillowPHdif))[2]),
                          PHprop=c(median(AlderPHprop),median(AshPHprop),median(BeechPHprop),median(BirchPHprop),median(ElmPHprop),median(HazelPHprop),median(OakPHprop),median(RowanPHprop),median(SycamorePHprop),median(WillowPHprop)),
                          PHpropLCI=c(HPDinterval(mcmc(AlderPHprop))[1], HPDinterval(mcmc(AshPHprop))[1], HPDinterval(mcmc(BeechPHprop))[1], HPDinterval(mcmc(BirchPHprop))[1], HPDinterval(mcmc(ElmPHprop))[1], HPDinterval(mcmc(HazelPHprop))[1], HPDinterval(mcmc(OakPHprop))[1], HPDinterval(mcmc(RowanPHprop))[1], HPDinterval(mcmc(SycamorePHprop))[1], HPDinterval(mcmc(WillowPHprop))[1]),
                          PHpropUCI=c(HPDinterval(mcmc(AlderPHprop))[2], HPDinterval(mcmc(AshPHprop))[2], HPDinterval(mcmc(BeechPHprop))[2], HPDinterval(mcmc(BirchPHprop))[2], HPDinterval(mcmc(ElmPHprop))[2], HPDinterval(mcmc(HazelPHprop))[2], HPDinterval(mcmc(OakPHprop))[2], HPDinterval(mcmc(RowanPHprop))[2], HPDinterval(mcmc(SycamorePHprop))[2], HPDinterval(mcmc(WillowPHprop))[2]),
                          PW=c(median(AlderPW, na.rm=TRUE),median(AshPW, na.rm=TRUE),median(BeechPW, na.rm=TRUE),median(BirchPW, na.rm=TRUE),median(ElmPW, na.rm=TRUE),median(HazelPW, na.rm=TRUE),median(OakPW, na.rm=TRUE),median(RowanPW, na.rm=TRUE),median(SycamorePW, na.rm=TRUE),median(WillowPW, na.rm=TRUE)),
                          PWLCI=c(HPDinterval(mcmc(AlderPW))[1], HPDinterval(mcmc(AshPW))[1], HPDinterval(mcmc(BeechPW))[1], HPDinterval(mcmc(BirchPW))[1], HPDinterval(mcmc(ElmPW))[1], HPDinterval(mcmc(HazelPW))[1], HPDinterval(mcmc(OakPW))[1], HPDinterval(mcmc(RowanPW))[1], HPDinterval(mcmc(SycamorePW))[1], HPDinterval(mcmc(WillowPW))[1]),
                          PWUCI=c(HPDinterval(mcmc(AlderPW))[2], HPDinterval(mcmc(AshPW))[2], HPDinterval(mcmc(BeechPW))[2], HPDinterval(mcmc(BirchPW))[2], HPDinterval(mcmc(ElmPW))[2], HPDinterval(mcmc(HazelPW))[2], HPDinterval(mcmc(OakPW))[2], HPDinterval(mcmc(RowanPW))[2], HPDinterval(mcmc(SycamorePW))[2], HPDinterval(mcmc(WillowPW))[2]),
                          PWdif=c(median(AlderPWdif, na.rm=TRUE),median(AshPWdif, na.rm=TRUE),median(BeechPWdif, na.rm=TRUE),median(BirchPWdif, na.rm=TRUE),median(ElmPWdif, na.rm=TRUE),median(HazelPWdif, na.rm=TRUE),median(OakPWdif, na.rm=TRUE),median(RowanPWdif, na.rm=TRUE),median(SycamorePWdif, na.rm=TRUE),median(WillowPWdif, na.rm=TRUE)),
                          PWdifLCI=c(HPDinterval(mcmc(AlderPWdif))[1], HPDinterval(mcmc(AshPWdif))[1], HPDinterval(mcmc(BeechPWdif))[1], HPDinterval(mcmc(BirchPWdif))[1], HPDinterval(mcmc(ElmPWdif))[1], HPDinterval(mcmc(HazelPWdif))[1], HPDinterval(mcmc(OakPWdif))[1], HPDinterval(mcmc(RowanPWdif))[1], HPDinterval(mcmc(SycamorePWdif))[1], HPDinterval(mcmc(WillowPWdif))[1]),
                          PWdifUCI=c(HPDinterval(mcmc(AlderPWdif))[2], HPDinterval(mcmc(AshPWdif))[2], HPDinterval(mcmc(BeechPWdif))[2], HPDinterval(mcmc(BirchPWdif))[2], HPDinterval(mcmc(ElmPWdif))[2], HPDinterval(mcmc(HazelPWdif))[2], HPDinterval(mcmc(OakPWdif))[2], HPDinterval(mcmc(RowanPWdif))[2], HPDinterval(mcmc(SycamorePWdif))[2], HPDinterval(mcmc(WillowPWdif))[2]),
                          PW.5=c(median(AlderPW.5, na.rm=TRUE),median(AshPW.5, na.rm=TRUE),median(BeechPW.5, na.rm=TRUE),median(BirchPW.5, na.rm=TRUE),median(ElmPW.5, na.rm=TRUE),median(HazelPW.5, na.rm=TRUE),median(OakPW.5, na.rm=TRUE),median(RowanPW.5, na.rm=TRUE),median(SycamorePW.5, na.rm=TRUE),median(WillowPW.5, na.rm=TRUE)),
                          PW.5LCI=c(HPDinterval(mcmc(AlderPW.5))[1], HPDinterval(mcmc(AshPW.5))[1], HPDinterval(mcmc(BeechPW.5))[1], HPDinterval(mcmc(BirchPW.5))[1], HPDinterval(mcmc(ElmPW.5))[1], HPDinterval(mcmc(HazelPW.5))[1], HPDinterval(mcmc(OakPW.5))[1], HPDinterval(mcmc(RowanPW.5))[1], HPDinterval(mcmc(SycamorePW.5))[1], HPDinterval(mcmc(WillowPW.5))[1]),
                          PW.5UCI=c(HPDinterval(mcmc(AlderPW.5))[2], HPDinterval(mcmc(AshPW.5))[2], HPDinterval(mcmc(BeechPW.5))[2], HPDinterval(mcmc(BirchPW.5))[2], HPDinterval(mcmc(ElmPW.5))[2], HPDinterval(mcmc(HazelPW.5))[2], HPDinterval(mcmc(OakPW.5))[2], HPDinterval(mcmc(RowanPW.5))[2], HPDinterval(mcmc(SycamorePW.5))[2], HPDinterval(mcmc(WillowPW.5))[2]),
                          PW.5dif=c(median(AlderPW.5dif, na.rm=TRUE),median(AshPW.5dif, na.rm=TRUE),median(BeechPW.5dif, na.rm=TRUE),median(BirchPW.5dif, na.rm=TRUE),median(ElmPW.5dif, na.rm=TRUE),median(HazelPW.5dif, na.rm=TRUE),median(OakPW.5dif, na.rm=TRUE),median(RowanPW.5dif, na.rm=TRUE),median(SycamorePW.5dif, na.rm=TRUE),median(WillowPW.5dif, na.rm=TRUE)),
                          PW.5difLCI=c(HPDinterval(mcmc(AlderPW.5dif))[1], HPDinterval(mcmc(AshPW.5dif))[1], HPDinterval(mcmc(BeechPW.5dif))[1], HPDinterval(mcmc(BirchPW.5dif))[1], HPDinterval(mcmc(ElmPW.5dif))[1], HPDinterval(mcmc(HazelPW.5dif))[1], HPDinterval(mcmc(OakPW.5dif))[1], HPDinterval(mcmc(RowanPW.5dif))[1], HPDinterval(mcmc(SycamorePW.5dif))[1], HPDinterval(mcmc(WillowPW.5dif))[1]),
                          PW.5difUCI=c(HPDinterval(mcmc(AlderPW.5dif))[2], HPDinterval(mcmc(AshPW.5dif))[2], HPDinterval(mcmc(BeechPW.5dif))[2], HPDinterval(mcmc(BirchPW.5dif))[2], HPDinterval(mcmc(ElmPW.5dif))[2], HPDinterval(mcmc(HazelPW.5dif))[2], HPDinterval(mcmc(OakPW.5dif))[2], HPDinterval(mcmc(RowanPW.5dif))[2], HPDinterval(mcmc(SycamorePW.5dif))[2], HPDinterval(mcmc(WillowPW.5dif))[2]))

#### Plotting curves in ggplot ####
Curves <- data.frame(date=days, "Fixed Effects"=MeanCurve, Alder=AlderCurve, Ash=AshCurve, Beech=BeechCurve, Birch=BirchCurve, Elm=ElmCurve, Hazel=HazelCurve, Oak=OakCurve, Rowan=RowanCurve, Sycamore=SycamoreCurve, Willow=WillowCurve)
Curveslong <- gather(Curves, key="TreeTaxa", value="Abundance", 2:12)
AllCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "black","darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")
AllTaxaCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "deepskyblue3", "royalblue4", "slateblue2", "orchid")

# Curves with legend
ggplot(Curveslong, aes(date, Abundance, col=TreeTaxa))+ #saved as 8"x9"
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))

#Curves without legend
AbundCurvesplot <- ggplot(Curveslong, aes(date, Abundance, col=TreeTaxa))+ 
  geom_line(size=0.6, aes(linetype=TreeTaxa))+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  scale_colour_manual(values=AllCols)+
  coord_cartesian(xlim=c(119,174))+
  guides(color = "none", linetype="none")+
  scale_linetype_manual(values=c("solid","solid","solid","solid","solid","dashed","solid","solid","solid","solid","solid"))

# Peak Date difference from mean
PDdifplot <- ggplot(AbundCurves, aes(fct_rev(TT), (PDdif), colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=(PDdifUCI), ymin=(PDdifLCI), width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Peak Timing (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

# Peak Height prop difference from mean
PHpropdifplot <- ggplot(AbundCurves, aes(fct_rev(TT), PHprop, colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (proportional)")+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  coord_flip()+
  scale_colour_manual(values=AllTaxaCols)

# Peak Width.5 difference from mean
PW.5difplot <- ggplot(AbundCurves, aes(fct_rev(TT), PW.5dif, colour=TT))+
  geom_point(size=2, alpha=0.5)+
  geom_errorbar(aes(ymax=PW.5difUCI, ymin=PW.5difLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Relative Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")+
  scale_colour_manual(values=AllTaxaCols)

library(gridExtra)

col5 <- grid.arrange(AbundCurvesplot, nrow = 1, widths = 1)
col7 <- grid.arrange(PDdifplot, PHpropdifplot, PW.5difplot, nrow = 3, heights = c(1,1,1))
AbundCurvesFig2 <- grid.arrange(col5, col7,  ncol = 2, widths = c(4,2)) #saved as 8"x10"

#different orientation
row5 <- grid.arrange(AbundCurvesplot, ncol = 1, heights = 1)
row7 <- grid.arrange(PDdifplot, PHpropdifplot, PW.5difplot, ncol = 3, widths = c(1,1,1))
AbundCurvesFig3 <- grid.arrange(row5, row7,  nrow = 2, heights = c(2,1)) #saved as 10"x10"


#### Supp. mat. figures ####

# Peak Width (set height) difference from mean
ggplot(AbundCurves, aes(fct_rev(TT), PWdif))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWdifUCI, ymin=PWdifLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Duration (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+#+
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
  guides(color = "none")# saved as 6"x4"



## Oak Difference PDHW 
AbundODcurves <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                            PD=c(median(AlderPDOD),median(AshPDOD),median(BeechPDOD),median(BirchPDOD),median(ElmPDOD),median(HazelPDOD),median(RowanPDOD),median(SycamorePDOD),median(WillowPDOD)),
                            PDLCI=c(HPDinterval(mcmc(AlderPDOD))[1], HPDinterval(mcmc(AshPDOD))[1], HPDinterval(mcmc(BeechPDOD))[1], HPDinterval(mcmc(BirchPDOD))[1], HPDinterval(mcmc(ElmPDOD))[1], HPDinterval(mcmc(HazelPDOD))[1], HPDinterval(mcmc(RowanPDOD))[1], HPDinterval(mcmc(SycamorePDOD))[1], HPDinterval(mcmc(WillowPDOD))[1]),
                            PDUCI=c(HPDinterval(mcmc(AlderPDOD))[2], HPDinterval(mcmc(AshPDOD))[2], HPDinterval(mcmc(BeechPDOD))[2], HPDinterval(mcmc(BirchPDOD))[2], HPDinterval(mcmc(ElmPDOD))[2], HPDinterval(mcmc(HazelPDOD))[2], HPDinterval(mcmc(RowanPDOD))[2], HPDinterval(mcmc(SycamorePDOD))[2], HPDinterval(mcmc(WillowPDOD))[2]),
                            PH=c(median(AlderPHOD),median(AshPHOD),median(BeechPHOD),median(BirchPHOD),median(ElmPHOD),median(HazelPHOD),median(RowanPHOD),median(SycamorePHOD),median(WillowPHOD)),
                            PHLCI=c(HPDinterval(mcmc(AlderPHOD))[1], HPDinterval(mcmc(AshPHOD))[1], HPDinterval(mcmc(BeechPHOD))[1], HPDinterval(mcmc(BirchPHOD))[1], HPDinterval(mcmc(ElmPHOD))[1], HPDinterval(mcmc(HazelPHOD))[1], HPDinterval(mcmc(RowanPHOD))[1], HPDinterval(mcmc(SycamorePHOD))[1], HPDinterval(mcmc(WillowPHOD))[1]),
                            PHUCI=c(HPDinterval(mcmc(AlderPHOD))[2], HPDinterval(mcmc(AshPHOD))[2], HPDinterval(mcmc(BeechPHOD))[2], HPDinterval(mcmc(BirchPHOD))[2], HPDinterval(mcmc(ElmPHOD))[2], HPDinterval(mcmc(HazelPHOD))[2], HPDinterval(mcmc(RowanPHOD))[2], HPDinterval(mcmc(SycamorePHOD))[2], HPDinterval(mcmc(WillowPHOD))[2]),
                            PHprop=c(median(AlderPHODprop),median(AshPHODprop),median(BeechPHODprop),median(BirchPHODprop),median(ElmPHODprop),median(HazelPHODprop),median(RowanPHODprop),median(SycamorePHODprop),median(WillowPHODprop)),
                            PHpropLCI=c(HPDinterval(mcmc(AlderPHODprop))[1], HPDinterval(mcmc(AshPHODprop))[1], HPDinterval(mcmc(BeechPHODprop))[1], HPDinterval(mcmc(BirchPHODprop))[1], HPDinterval(mcmc(ElmPHODprop))[1], HPDinterval(mcmc(HazelPHODprop))[1], HPDinterval(mcmc(RowanPHODprop))[1], HPDinterval(mcmc(SycamorePHODprop))[1], HPDinterval(mcmc(WillowPHODprop))[1]),
                            PHpropUCI=c(HPDinterval(mcmc(AlderPHODprop))[2], HPDinterval(mcmc(AshPHODprop))[2], HPDinterval(mcmc(BeechPHODprop))[2], HPDinterval(mcmc(BirchPHODprop))[2], HPDinterval(mcmc(ElmPHODprop))[2], HPDinterval(mcmc(HazelPHODprop))[2], HPDinterval(mcmc(RowanPHODprop))[2], HPDinterval(mcmc(SycamorePHODprop))[2], HPDinterval(mcmc(WillowPHODprop))[2]),
                            PW=c(median(AlderPWOD, na.rm=TRUE),median(AshPWOD, na.rm=TRUE),median(BeechPWOD, na.rm=TRUE),median(BirchPWOD, na.rm=TRUE),median(ElmPWOD, na.rm=TRUE),median(HazelPWOD, na.rm=TRUE),median(RowanPWOD, na.rm=TRUE),median(SycamorePWOD, na.rm=TRUE),median(WillowPWOD, na.rm=TRUE)),
                            PWLCI=c(HPDinterval(mcmc(AlderPWOD))[1], HPDinterval(mcmc(AshPWOD))[1], HPDinterval(mcmc(BeechPWOD))[1], HPDinterval(mcmc(BirchPWOD))[1], HPDinterval(mcmc(ElmPWOD))[1], HPDinterval(mcmc(HazelPWOD))[1], HPDinterval(mcmc(RowanPWOD))[1], HPDinterval(mcmc(SycamorePWOD))[1], HPDinterval(mcmc(WillowPWOD))[1]),
                            PWUCI=c(HPDinterval(mcmc(AlderPWOD))[2], HPDinterval(mcmc(AshPWOD))[2], HPDinterval(mcmc(BeechPWOD))[2], HPDinterval(mcmc(BirchPWOD))[2], HPDinterval(mcmc(ElmPWOD))[2], HPDinterval(mcmc(HazelPWOD))[2], HPDinterval(mcmc(RowanPWOD))[2], HPDinterval(mcmc(SycamorePWOD))[2], HPDinterval(mcmc(WillowPWOD))[2]),
                            PW.5=c(median(AlderPW.5OD, na.rm=TRUE),median(AshPW.5OD, na.rm=TRUE),median(BeechPW.5OD, na.rm=TRUE),median(BirchPW.5OD, na.rm=TRUE),median(ElmPW.5OD, na.rm=TRUE),median(HazelPW.5OD, na.rm=TRUE),median(RowanPW.5OD, na.rm=TRUE),median(SycamorePW.5OD, na.rm=TRUE),median(WillowPW.5OD, na.rm=TRUE)),
                            PW.5LCI=c(HPDinterval(mcmc(AlderPW.5OD))[1], HPDinterval(mcmc(AshPW.5OD))[1], HPDinterval(mcmc(BeechPW.5OD))[1], HPDinterval(mcmc(BirchPW.5OD))[1], HPDinterval(mcmc(ElmPW.5OD))[1], HPDinterval(mcmc(HazelPW.5OD))[1], HPDinterval(mcmc(RowanPW.5OD))[1], HPDinterval(mcmc(SycamorePW.5OD))[1], HPDinterval(mcmc(WillowPW.5OD))[1]),
                            PW.5UCI=c(HPDinterval(mcmc(AlderPW.5OD))[2], HPDinterval(mcmc(AshPW.5OD))[2], HPDinterval(mcmc(BeechPW.5OD))[2], HPDinterval(mcmc(BirchPW.5OD))[2], HPDinterval(mcmc(ElmPW.5OD))[2], HPDinterval(mcmc(HazelPW.5OD))[2], HPDinterval(mcmc(RowanPW.5OD))[2], HPDinterval(mcmc(SycamorePW.5OD))[2], HPDinterval(mcmc(WillowPW.5OD))[2]))

NoOakCols <- c("darkred", "firebrick3", "chocolate2", "goldenrod", "olivedrab4", "darkgreen", "royalblue4", "slateblue2", "orchid")

PDODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PD))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PDUCI, ymin=PDLCI, width=0.5))+
  theme_bw()+
  xlab("Tree Taxon")+
  ylab("Timing (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PHODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PH))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PHUCI, ymin=PHLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (abundance)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PHpropODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PHprop))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PHpropUCI, ymin=PHpropLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Height (proportional)")+
  coord_flip()+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWODplot <- ggplot(AbundODcurves, aes(fct_rev(TT), PW))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PWUCI, ymin=PWLCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Duration (days) ")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))

PWOD.5plot <- ggplot(AbundODcurves, aes(fct_rev(TT), PW.5))+
  geom_point(size=2, alpha=0.9)+
  geom_errorbar(aes(ymax=PW.5UCI, ymin=PW.5LCI, width=0.5))+
  theme_bw()+
  xlab("")+
  ylab("Width (days)")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15))


AbundODplots <- grid.arrange(PDODplot, PHpropODplot, PWOD.5plot, PWODplot, ncol = 4, widths=c(1,1,1,1)) #saved 6" by 13"

#### Posterior distribution sample size (without NaNs) ####
AbundCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Oak", "Rowan", "Sycamore", "Willow"), 
                                    PD=c(length(which(is.na(AlderPD)==FALSE)),length(which(is.na(AshPD)==FALSE)),length(which(is.na(BeechPD)==FALSE)),length(which(is.na(BirchPD)==FALSE)),length(which(is.na(ElmPD)==FALSE)),length(which(is.na(HazelPD)==FALSE)),length(which(is.na(OakPD)==FALSE)),length(which(is.na(RowanPD)==FALSE)),length(which(is.na(SycamorePD)==FALSE)),length(which(is.na(WillowPD)==FALSE))),
                                    PDdif=c(length(which(is.na(AlderPDdif)==FALSE)),length(which(is.na(AshPDdif)==FALSE)),length(which(is.na(BeechPDdif)==FALSE)),length(which(is.na(BirchPDdif)==FALSE)),length(which(is.na(ElmPDdif)==FALSE)),length(which(is.na(HazelPDdif)==FALSE)),length(which(is.na(OakPDdif)==FALSE)),length(which(is.na(RowanPDdif)==FALSE)),length(which(is.na(SycamorePDdif)==FALSE)),length(which(is.na(WillowPDdif)==FALSE))),
                                    PH=c(length(which(is.na(AlderPH)==FALSE)),length(which(is.na(AshPH)==FALSE)),length(which(is.na(BeechPH)==FALSE)),length(which(is.na(BirchPH)==FALSE)),length(which(is.na(ElmPH)==FALSE)),length(which(is.na(HazelPH)==FALSE)),length(which(is.na(OakPH)==FALSE)),length(which(is.na(RowanPH)==FALSE)),length(which(is.na(SycamorePH)==FALSE)),length(which(is.na(WillowPH)==FALSE))),
                                    PHdif=c(length(which(is.na(AlderPHdif)==FALSE)),length(which(is.na(AshPHdif)==FALSE)),length(which(is.na(BeechPHdif)==FALSE)),length(which(is.na(BirchPHdif)==FALSE)),length(which(is.na(ElmPHdif)==FALSE)),length(which(is.na(HazelPHdif)==FALSE)),length(which(is.na(OakPHdif)==FALSE)),length(which(is.na(RowanPHdif)==FALSE)),length(which(is.na(SycamorePHdif)==FALSE)),length(which(is.na(WillowPHdif)==FALSE))),
                                    PHprop=c(length(which(is.na(AlderPHprop)==FALSE)),length(which(is.na(AshPHprop)==FALSE)),length(which(is.na(BeechPHprop)==FALSE)),length(which(is.na(BirchPHprop)==FALSE)),length(which(is.na(ElmPHprop)==FALSE)),length(which(is.na(HazelPHprop)==FALSE)),length(which(is.na(OakPHprop)==FALSE)),length(which(is.na(RowanPHprop)==FALSE)),length(which(is.na(SycamorePHprop)==FALSE)),length(which(is.na(WillowPHprop)==FALSE))),
                                    PW=c(length(which(is.na(AlderPW)==FALSE)),length(which(is.na(AshPW)==FALSE)),length(which(is.na(BeechPW)==FALSE)),length(which(is.na(BirchPW)==FALSE)),length(which(is.na(ElmPW)==FALSE)),length(which(is.na(HazelPW)==FALSE)),length(which(is.na(OakPW)==FALSE)),length(which(is.na(RowanPW)==FALSE)),length(which(is.na(SycamorePW)==FALSE)),length(which(is.na(WillowPW)==FALSE))),
                                    PWdif=c(length(which(is.na(AlderPWdif)==FALSE)),length(which(is.na(AshPWdif)==FALSE)),length(which(is.na(BeechPWdif)==FALSE)),length(which(is.na(BirchPWdif)==FALSE)),length(which(is.na(ElmPWdif)==FALSE)),length(which(is.na(HazelPWdif)==FALSE)),length(which(is.na(OakPWdif)==FALSE)),length(which(is.na(RowanPWdif)==FALSE)),length(which(is.na(SycamorePWdif)==FALSE)),length(which(is.na(WillowPWdif)==FALSE))),
                                    PW.5=c(length(which(is.na(AlderPW.5)==FALSE)),length(which(is.na(AshPW.5)==FALSE)),length(which(is.na(BeechPW.5)==FALSE)),length(which(is.na(BirchPW.5)==FALSE)),length(which(is.na(ElmPW.5)==FALSE)),length(which(is.na(HazelPW.5)==FALSE)),length(which(is.na(OakPW.5)==FALSE)),length(which(is.na(RowanPW.5)==FALSE)),length(which(is.na(SycamorePW.5)==FALSE)),length(which(is.na(WillowPW.5)==FALSE))),
                                    PW.5dif=c(length(which(is.na(AlderPW.5dif)==FALSE)),length(which(is.na(AshPW.5dif)==FALSE)),length(which(is.na(BeechPW.5dif)==FALSE)),length(which(is.na(BirchPW.5dif)==FALSE)),length(which(is.na(ElmPW.5dif)==FALSE)),length(which(is.na(HazelPW.5dif)==FALSE)),length(which(is.na(OakPW.5dif)==FALSE)),length(which(is.na(RowanPW.5dif)==FALSE)),length(which(is.na(SycamorePW.5dif)==FALSE)),length(which(is.na(WillowPW.5dif)==FALSE))))
#write.csv(AbundCurvesSampleSize,'~/Documents/Models/Tables/AbundCurvesSampleSize.csv')

AbundODCurvesSampleSize <- data.frame(TT=c("Alder", "Ash", "Beech", "Birch", "Elm", "Hazel", "Rowan", "Sycamore", "Willow"), 
                                      PD=c(length(which(is.na(AlderPDOD)==FALSE)),length(which(is.na(AshPDOD)==FALSE)),length(which(is.na(BeechPDOD)==FALSE)),length(which(is.na(BirchPDOD)==FALSE)),length(which(is.na(ElmPDOD)==FALSE)),length(which(is.na(HazelPDOD)==FALSE)),length(which(is.na(RowanPDOD)==FALSE)),length(which(is.na(SycamorePDOD)==FALSE)),length(which(is.na(WillowPDOD)==FALSE))),
                                      PH=c(length(which(is.na(AlderPHOD)==FALSE)),length(which(is.na(AshPHOD)==FALSE)),length(which(is.na(BeechPHOD)==FALSE)),length(which(is.na(BirchPHOD)==FALSE)),length(which(is.na(ElmPHOD)==FALSE)),length(which(is.na(HazelPHOD)==FALSE)),length(which(is.na(RowanPHOD)==FALSE)),length(which(is.na(SycamorePHOD)==FALSE)),length(which(is.na(WillowPHOD)==FALSE))),
                                      PHprop=c(length(which(is.na(AlderPHODprop)==FALSE)),length(which(is.na(AshPHODprop)==FALSE)),length(which(is.na(BeechPHODprop)==FALSE)),length(which(is.na(BirchPHODprop)==FALSE)),length(which(is.na(ElmPHODprop)==FALSE)),length(which(is.na(HazelPHODprop)==FALSE)),length(which(is.na(RowanPHODprop)==FALSE)),length(which(is.na(SycamorePHODprop)==FALSE)),length(which(is.na(WillowPHODprop)==FALSE))),
                                      PW=c(length(which(is.na(AlderPWOD)==FALSE)),length(which(is.na(AshPWOD)==FALSE)),length(which(is.na(BeechPWOD)==FALSE)),length(which(is.na(BirchPWOD)==FALSE)),length(which(is.na(ElmPWOD)==FALSE)),length(which(is.na(HazelPWOD)==FALSE)),length(which(is.na(RowanPWOD)==FALSE)),length(which(is.na(SycamorePWOD)==FALSE)),length(which(is.na(WillowPWOD)==FALSE))),
                                      PW.5=c(length(which(is.na(AlderPW.5OD)==FALSE)),length(which(is.na(AshPW.5OD)==FALSE)),length(which(is.na(BeechPW.5OD)==FALSE)),length(which(is.na(BirchPW.5OD)==FALSE)),length(which(is.na(ElmPW.5OD)==FALSE)),length(which(is.na(HazelPW.5OD)==FALSE)),length(which(is.na(RowanPW.5OD)==FALSE)),length(which(is.na(SycamorePW.5OD)==FALSE)),length(which(is.na(WillowPW.5OD)==FALSE))))
#write.csv(AbundODCurvesSampleSize,'~/Documents/Models/Tables/AbundODCurvesSampleSize.csv')

meanmetrics <- data.frame(metric=c("PD", "PH", "PW.5", "PW.01"), 
                          median=c(median(MeanPD), median(MeanPH),median(MeanPW.5,na.rm=TRUE), median(MeanPW,na.rm=TRUE)), 
                          lci=c(HPDinterval(mcmc(MeanPD))[1], HPDinterval(mcmc(MeanPH))[1],HPDinterval(mcmc(MeanPW.5))[1], HPDinterval(mcmc(MeanPW))[1]),
                          uci=c(HPDinterval(mcmc(MeanPD))[2], HPDinterval(mcmc(MeanPH))[2],HPDinterval(mcmc(MeanPW.5))[2], HPDinterval(mcmc(MeanPW))[2]),
                          SS=c(length(which(is.na(MeanPD)==FALSE)),length(which(is.na(MeanPH)==FALSE)),length(which(is.na(MeanPW.5)==FALSE)),length(which(is.na(MeanPW)==FALSE))))
#write.csv(meanmetrics,'~/Documents/Models/Tables/AbundFixedEffMetrics.csv')


############################
#### Model output table ####  
############################


#### fixed effects
fixed<-rbind(
  c("Intercept",paste(round(mean(ATTCyearint$Sol[,1]),3)," (",
                      round(HPDinterval(ATTCyearint$Sol[,1])[1],3)," - ",
                      round(HPDinterval(ATTCyearint$Sol[,1])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,1]))),
  c("Date (scaled)",paste(round(mean(ATTCyearint$Sol[,2]),3)," (",
                          round(HPDinterval(ATTCyearint$Sol[,2])[1],3)," - ",
                          round(HPDinterval(ATTCyearint$Sol[,2])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,2]))),
  c("Date² (scaled)",paste(round(mean(ATTCyearint$Sol[,3]),3)," (",
                           round(HPDinterval(ATTCyearint$Sol[,3])[1],3)," - ",
                           round(HPDinterval(ATTCyearint$Sol[,3])[2],3),")",sep=""),
    round(effectiveSize(ATTCyearint$Sol[,3]))))

#### random terms
column<-1
treetaxa1<-c("TreeTaxa- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-2
treetaxa2<-c("TreeTaxa- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                          round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-3
treetaxa3<-c("TreeTaxa- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                           round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-5
treetaxa5<-c("TreeTaxa- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                              round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-6
treetaxa6<-c("TreeTaxa- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-9
treetaxa9<-c("TreeTaxa- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                               round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
             round(effectiveSize(ATTCyearint$VCV[, column])))

column<-10
site10<-c("Site- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                      round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-11
site11<-c("Site- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                   round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                   round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-12
site12<-c("Site- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                    round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                    round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-14
site14<-c("Site- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                       round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-15
site15<-c("Site- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                     round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                     round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-18
site18<-c("Site- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                        round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                        round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-19
year10<-c("Year- Intercept var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                      round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                      round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-20
year11<-c("Year- Intercept:Date slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                   round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                   round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-21
year12<-c("Year- Intercept:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                    round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                    round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-23
year14<-c("Year- Date slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                       round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                       round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-24
year15<-c("Year- Date slope:Date² slope covar",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                                     round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                                     round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-27
year18<-c("Year- Date² slope var",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                                        round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                                        round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-31
recorder<-c("Recorder",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))


column<-30
siteday<-c("Site-Day",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                            round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                            round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
           round(effectiveSize(ATTCyearint$VCV[, column])))

column<-28
siteyear<-c("Site-Year",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                              round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                              round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))

column<-29
treeID<-c("Tree ID",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                          round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                          round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
          round(effectiveSize(ATTCyearint$VCV[, column])))

column<-32
residual<-c("Residual",paste(round(posterior.mode(ATTCyearint$VCV[, column]),3)," (",
                             round(HPDinterval(ATTCyearint$VCV[, column])[1],3)," - ",
                             round(HPDinterval(ATTCyearint$VCV[, column])[2],3),")",sep=""),
            round(effectiveSize(ATTCyearint$VCV[, column])))

random<-rbind(treetaxa1,treetaxa2,treetaxa3,treetaxa5,treetaxa6,treetaxa9,site10,site11,site12,site14,site15,site18,year10,year11,year12,year14,year15,year18, siteyear, recorder, siteday, treeID, residual)

#write.table(rbind(c("Fixed Terms","",""),fixed,c("Random Terms","",""),random),"~/Documents/Models/Tables/TableATTCyearint.txt",sep="\t",col.names=c("","Coefficient/Variance (Mean/mode and CI)","Effective sample size"),row.names=F)
