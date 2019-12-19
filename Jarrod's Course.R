###############################
#### Jarrod's Stats Course ####
###############################

#### Practical 1 ####

photo_long <- read.csv("~/Downloads/photo_long.csv")

m1 <- lm(y~type, data=photo_long)
summary(m1) # differece between means 1.2283
summary(m1)$coefficients[2]
grumpymean <- summary(m1)$coefficients[1]
happymean <- grumpymean+summary(m1)$coefficients[2]
sd <- summary(m1)$sigma

outputs <- data.frame(simno=seq(1,1000,1))
for(i in 1:1000){
  grumpysim <- rnorm(22, mean=grumpymean, sd=sd)
  happysim <- rnorm(22, mean=happymean, sd=sd)
  sim <- c(grumpysim, happysim)
  mood <- c(rep("grumpy", 22), rep("happy", 22))
  m2 <- lm(sim~mood)
  outputs$grumpy[i] <- summary(m2)$coefficients[1]
  outputs$diff[i] <- summary(m2)$coefficients[2]
  outputs$happy[i] <- summary(m2)$coefficients[1]+summary(m2)$coefficients[2]
  outputs$sd[i] <- summary(m2)$sigma
}

hist(outputs$grumpy, breaks=100)
abline(v=grumpymean,col=2)

hist(outputs$happy, breaks=100)
abline(v=happymean,col=2)

mean(outputs$grumpy)
grumpymean
sd(outputs$grumpy)
sd
mean(outputs$happy)
happymean
#no sd for happy
mean(outputs$diff)
sd(outputs$diff)

qqnorm(outputs$grumpy)
qqline(outputs$grumpy)
qqnorm(outputs$diff)
qqline(outputs$diff)

#grumpy expected sampling distribution from m1- mean=5.3759, sd of mean=0.2619 
#happy expected sampling distribution from m1- mean=4.1476, sd of mean cant be easily seen
#SE in table by happy is the SD of the difference in mean (like how it works for difference in mean intercept)

funcsims <- simulate(nsim=1000, m1) #simulate shortcut - categories are same structure down column as in original dataset
 

#### Practical 2 ####

collage <- read.csv("~/Downloads/collage_photo.csv")

# What is the best estimate of how many more people (as a %) do people recognise when given more time?

collage$speed <- as.factor(collage$speed)

m1 <- lm(y~speed, collage)
summary(m1)

secs3 <- summary(m1)$coefficient[1]# 0.9394 in 3 seconds
secs5 <- summary(m1)$coefficient[1]+summary(m1)$coefficient[2] #1.2727.. in 5 secs
diff <- summary(m1)$coefficient[2] # difference
percentdiff <- (diff/secs3)*100 #35.5% more people 

boxplot(y~speed,collage)
# Do people improve as they assess more arrays? no?
collage$rep <- rep(seq(1,6,1),11)
m2 <- lm(y~rep+speed,collage)
summary(m2)

ggplot(collage, aes(rep,y, col=speed))+
  geom_point()

#qqnorm(m2$residuals)
#qqline(m2$residuals) # data is integers so residulas are clumped and removed from line at ends.

plot(m2)
# Assume that the improvement is linear with respect to how many times they had tried. The expected number of IEB photos on an array was 2. 
# Does a participant on their first try and given 5 seconds recognise significantly less than this number.

collage$speed <- relevel(collage$speed, ref="5")
collage$rep <- rep(seq(0,5,1),11)

m3 <- lm(y~rep+speed,collage)
summary(m3)
# mean:1.31818    SE:0.20170
library(metRology)
pt.scaled(1.31818, mean = 2, sd = 0.20170, df = 66 - 3) # 0.0006235787   degrees of freedom is no. of datapoints - no. of location parameters (intercept, speed, rep)
pt.scaled(1.31818-2, mean = 0, sd = 0.20170, df = 66 - 3) # 0.0006235787

# Assume that the improvement is linear with respect to the time the people were given. 
# If the point estimates were the true values, how long should I have given people if I wanted them to recognise two people on their first try (on average)?

collage$speed <- as.numeric(as.character(collage$speed))
m4 <- lm(y~speed+rep, collage)
summary(m4)
# y=0.09848*speed+ intercept  speed=(2-intercept)/0.09848
recog2 <- (2-summary(m4)$coefficients[1])/summary(m4)$coefficients[2]

# Put approximate confidence intervals on this using simulation.
# y= b1*1 + b2*speed + b3*rep

m4sim <- simulate(nsim=1000,m4)
m4sim$speed <- collage$speed
m4sim$rep <- collage$rep




outputs <- data.frame(simno=seq(1,1000,1))
for(i in 1:1000){
  m5 <- lm(m4sim[,i]~m4sim$speed+m4sim$rep)
  outputs$time[i] <- (2-summary(m5)$coefficients[1])/summary(m5)$coefficients[2]
}
hist(outputs$time)   
sd(outputs$time)

#### Practical 4 ####
library(lme4)
photo_long <- read.csv("~/Downloads/photo_long.csv")
photo_wide <- read.csv("~/Downloads/photo_wide.csv")

# Using the long data estimate the effect of happy and grumpy on the scores using both lm and lmer
lm <- lm(y~type,data=photo_long) # grumpy=5.38 se=0.26, diff= -1.23 se=0.37, happy= 4.15 
lmm <- lmer(y~type + (1|person), data=photo_long) # grumpy=5.38 se=0.26, diff= -1.23 se=0.25, happy= 4.15 

# Do the standard errors change in the way you expected them to?
## Yes the intercept SE doesnt change (unless there is also a continuous variable?)

#Using the function (t.test) (and the wide data) try and get the same answer as in the linear mixed model (you will need to read the t.test help file)
t.test(photo_wide$y.grumpy, photo_wide$y.happy, paired=TRUE)
## paired t test gives the same results as the lmm



#Estimate the amount of variation due to respondent and person and try and think what this means.
full <- read.csv("~/Downloads/photo_long_full.csv")
full$respondent <- as.factor(full$respondent)
lmm2 <- lmer(y~1+(1|respondent)+(1|person), data=full)

random <- as_dataframe(ranef(lmm2))
