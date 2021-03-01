library(readr)
library(lme4)
Kel <- read_csv("~/Downloads/Kel.csv")
Kelsub <- subset(Kel, species=="Montipora capitata")
Kelsub$acc <- Kelsub$`Accretion rate`
Kelsub$trmt <- Kelsub$Site
Kelsub$site <- as.character(Kelsub$Treatment)
Kelsub$site <- gsub("3","1", Kelsub$site)
Kelsub$site <- gsub("4","2", Kelsub$site)
model <- lmer(acc~site*trmt + (1|Parent_Colony), data=Kelsub)
mod <- lmer(acc~ (1|Parent_Colony), data=Kelsub)
anova <- anova(model,mod)
summary(model)

Kelsub$trmt2 <- Kelsub$trmt
Kelsub$trmt2[16:30] <- gsub("Transplant","TranTo1", Kelsub$trmt2[16:30])
Kelsub$trmt2[46:60] <- gsub("Transplant","TranTo2", Kelsub$trmt2[46:60])

model2 <- lmer(acc~site+trmt2 + (1|Parent_Colony), data=Kelsub)
Kelsub$trmt2 <- relevel(as.factor(Kelsub$trmt2), ref="TranTo1")
model3 <- lmer(acc~site+trmt2 + (1|Parent_Colony), data=Kelsub)

ggplot(Kelsub, aes(x=site, y=acc))+
  geom_boxplot(aes(colour=trmt))+
  theme_bw()

#baysmodel <- MCMCglmm(acc~site*trmt, random=~Parent_Colony, data=Kelsub)
