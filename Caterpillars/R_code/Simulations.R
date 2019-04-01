##### testing power

p <- c() # empty space for data
for(x in 1:1000){
treatment <- as.factor(c(rep(1,15), rep(2,15))) #2 treatments, 15 of each
clutches <- c(rnorm(15,8,1), rnorm(15,7,1)) #15 clutches in each treatment with means 8 and 7 and sd=1

model <- lm(clutches~treatment) # model 
p[x] <- summary(model)$coefficients[2,4] # extract p values
}
prop.sig <- length(which(p<0.05))/1000 # what proportion of the simulations are significant- chance of noticing the true difference wiht these sample sizes etc
# change to same means (no difference) to see likelihood of false positive


###### checking if a coefficient (slope here) is biased
# simulate using the slope value and then if you arent signif likely to get that value from simulation distribution its biased..?

slopes <- c()
for(x in 1:1000){
area <- runif(30,0,10) # 30 samples between 0 and 10, uniform distribution
lambda <- exp(area*0.25+2) # creating lambda for the poisson distribution of species numbers
species <- rpois(30,lambda)
### Ally did this slightly differently: before loop defined slope=0.25 and intercept=2
#Then resp<-rpois(30, exp(intercept+slope*logislandsize))
#model<-glm(resp~logislandsize,family=Poisson)

model <- glm(species~area, family = "poisson")
slopes[x] <- summary(model)$coefficients[2,1] # extract slope
}

mean.slope <- mean(slopes)
hist(slopes, breaks=30)
abline(v=mean(slopes),col=2)

length(which(slopes>0.25))
length(which(slopes<0.25))       
## to get a p value for the probability of the slop being different to the 0.25 value you put in 
# look a the proportion of samples above and beneath 0.25 and time smaller probability by 2?????

##################
### simulating fixed and random effects

#If the effect of temperature on lay date is -4, how precise an estimate would we expect with the following design?
#intercept = 150
#44 sites (variance = 20)
#5 years (variance = 40)
#8 nestboxes per site (resid var = 20).
#Temperatures for each site:year drawn at random in the range 5-10C.

slope <- -4
intercept <- 150

df <- data.frame(
  years=as.factor(c(rep(1,352), rep(2,352), rep(3,352), rep(4,352), rep(5,352))),
  sites=as.factor(c(rep(c(seq(1,44,1), seq(1,44,1), seq(1,44,1), seq(1,44,1), seq(1,44,1), seq(1,44,1), seq(1,44,1), seq(1,44,1)),5))),
  nestboxes=as.factor(c(rep(c(rep(1,44), rep(2,44), rep(3,44), rep(4,44), rep(5,44), rep(6,44), rep(7,44), rep(8,44)),5))))
df$siteyear <- as.factor(paste(df$years, df$sites))
temp <- data.frame(siteyear=unique(paste(df$years, df$sites)))
temp$temp <- runif(220,5,10)      # each siteyear random temp within 5-10 degrees
#pmatch(df$siteyear, temp$siteyear)
#matchSY <- pmatch(df$siteyear, temp$siteyear)                         
dftemp <- merge(df, temp, by="siteyear", duplicates.ok=TRUE)

store <- c()

laydate <- dftemp$temp*slope+intercept  ### need to assign the variance within each variable..?
