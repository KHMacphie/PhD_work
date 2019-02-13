#### Basic Understanding ####

rm(list=ls())
setwd('/Users/s1205615/')
library(ggplot2)
library(dplyr)
library(ggfortify)
library(readr)
#library(doBy)
library(lme4)
library(MCMCglmm)

##### Setting up dataframe #####

site <- read_csv("Dropbox/master_data/site/site_details.csv", 
                 col_types = cols(`Mean Elev` = col_double()))

#site <- read_csv("~/Dropbox/master_data/site/site_details.csv", 
#col_types = cols(`Mean Elev` = col_double()))
cater <- read_csv("Dropbox/master_data/inverts/Branch_Beating.csv", 
                  col_types = cols(year = col_factor(levels = c("2014", 
                                                                "2015", "2016", "2017", "2018"))))
#cater <- read_csv("~/Dropbox/master_data/inverts/Branch_Beating.csv", 
#col_types = cols(year = col_factor(levels = c("2014", 
#                                             "2015", "2016", "2017", "2018"))))

temp <- read_csv("Dropbox/master_data/site/temperatures.csv", 
                 col_types = cols(year = col_factor(levels = c("2014", 
                                                               "2015", "2016", "2017", "2018"))))

#cater<- mutate(cater, catbinom=caterpillars)
#cater$catbinom[which(cater$caterpillars>1)]<-1
#cater<-cater[-which(is.na(cater$catbinom)==TRUE),]
pmatch(cater$site,site$site,duplicates.ok=TRUE)
all_data<- merge(cater, site, by="site", duplicates.ok=TRUE)
all_data<- rename(all_data, latitude="Mean Lat")
all_data<- rename(all_data, longitude="Mean Long")
all_data<- rename(all_data, elevation="Mean Elev")
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data<-all_data[-which(is.na(all_data$caterpillars)==TRUE),] 

temp$yearsite<- paste(temp$site, temp$year) #nests site and year

#### Temperature data frame ####
mean_temps <- data.frame(site=temp$site, year=temp$year)
#View(mean_temps)
mean_temps <- mean_temps[!duplicated(mean_temps), ] #remove duplicated rows
mean_temps$yearsite <- paste(mean_temps$site, mean_temps$year)
pmatch(mean_temps$yearsite, temp$yearsite)
mean_temps <- mean_temps %>% arrange(site) #arrange by site to match means order

## Using mean temp through Apr

mean_temps$Apr <-   tapply(apply(temp[, 1059:1778],1, mean), temp$yearsite, mean) #calculating mean temp within time window for each logger then the mean of the two values per site

## Putting into full dataframe with all beating and site data
all_data$yearsite<- paste(all_data$site, all_data$year)
pmatch(all_data$yearsite, mean_temps$yearsite)
mean_temps <- select(mean_temps, -site, -year)
all_data<- merge(all_data, mean_temps, by="yearsite", duplicates.ok=TRUE)
all_data$sitetree <- paste(all_data$tree, all_data$site)
all_data$siteday <- paste(all_data$site, all_data$date, all_data$year)
all_data$obs<-as.factor(seq(1,length(all_data[,1])))

#### Using lme4 to quickly look at what different parts do/how the log link functions effects tihngs

basic <- glm(caterpillars~date+I(date^2)+year+Apr, family=poisson, data=all_data)
summary(basic)

pred_day <- seq(120,175,1)
pred_day2 <- seq(0,200,1)
basic14 <- basic$coef["(Intercept)"] + basic$coef["date"]*pred_day + basic$coef["I(date^2)"]*pred_day^2 + basic$coef["Apr"]*7.5
plot(pred_day, exp(basic14), type="l")
plot(pred_day, basic14, type="l")
basicnumbers <- exp(-93.51879) * exp(1.213843*pred_day) * exp(-0.00394120*pred_day^2)  * exp(-0.1308632*7.5) 
#for negative coefficients can remove - from brackets and change * to /
plot(pred_day, basicnumbers, type="l")


noquad <- glm(caterpillars~date+year+Apr, family=poisson, data=all_data)
summary(noquad)
noquad14 <- noquad$coef["(Intercept)"] + noquad$coef["date"]*pred_day +noquad$coef["Apr"]*7.5
plot(pred_day, exp(noquad14), type="l")
plot(pred_day, noquad14, type="l")

### so the classic parabola or straight line is with the logged coefficients 

# going back over how things move between log and exp() with effects on intercept and from interactions
basic15 <- basic$coef["(Intercept)"] + basic$coef["year2015"] + basic$coef["date"]*pred_day + basic$coef["I(date^2)"]*pred_day^2+basic$coef["Apr"]*7.5
plot(pred_day, exp(basic14), type="l")
points(pred_day, exp(basic15), type="l", col=2)

plot(pred_day, basic14, type="l")
points(pred_day, basic15, type="l", col=2)

noquad15 <- noquad$coef["(Intercept)"] + noquad$coef["year2015"] + noquad$coef["date"]*pred_day + noquad$coef["Apr"]*7.5
plot(pred_day, exp(noquad14), type="l", ylim=c(0,0.5))
points(pred_day, exp(noquad15), type="l", col=2)

plot(pred_day, noquad14, type="l", ylim=c(-1,-10))
points(pred_day, noquad15, type="l", col=2)

# only intercept changing as nothing interacting with day of the year so peak date is the same but height varies


interactionfactor <- glm(caterpillars~date*year+I(date^2)+Apr, family=poisson, data=all_data)
summary(interactionfactor)

interactionfactor14 <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["date"]*pred_day + interactionfactor$coef["I(date^2)"]*pred_day^2 + interactionfactor$coef["Apr"]*7.5
interactionfactor15 <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["year2015"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2015"])*pred_day + interactionfactor$coef["I(date^2)"]*pred_day^2 + interactionfactor$coef["Apr"]*7.5
interactionfactor16 <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["year2016"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2016"])*pred_day + interactionfactor$coef["I(date^2)"]*pred_day^2 + interactionfactor$coef["Apr"]*7.5
interactionfactor17 <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["year2017"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2017"])*pred_day + interactionfactor$coef["I(date^2)"]*pred_day^2 + interactionfactor$coef["Apr"]*7.5
interactionfactor18 <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["year2018"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2018"])*pred_day + interactionfactor$coef["I(date^2)"]*pred_day^2 + interactionfactor$coef["Apr"]*7.5
plot(pred_day, interactionfactor14, type="l")#, ylim=c(-1,-10))
points(pred_day, interactionfactor15, type="l", col=2)
points(pred_day, interactionfactor16, type="l", col=3)
points(pred_day, interactionfactor17, type="l", col=4)
points(pred_day, interactionfactor18, type="l", col=5)

plot(pred_day, exp(interactionfactor14), type="l") #black
points(pred_day, exp(interactionfactor15), type="l", col=2) #red
points(pred_day, exp(interactionfactor16), type="l", col=3) #green
points(pred_day, exp(interactionfactor17), type="l", col=4) #blue
points(pred_day, exp(interactionfactor18), type="l", col=5) #lightblue

#not putting in the intercept year change or quadratic to check interaction effect- different gradients as expected
interactionfactor14wrong <- interactionfactor$coef["(Intercept)"] + interactionfactor$coef["date"]*pred_day2  + interactionfactor$coef["Apr"]*7.5 + interactionfactor$coef["I(date^2)"]*pred_day2^2
interactionfactor15wrong <- interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2015"])*pred_day2 + interactionfactor$coef["Apr"]*7.5+ interactionfactor$coef["I(date^2)"]*pred_day2^2
interactionfactor16wrong <- interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2016"])*pred_day2 + interactionfactor$coef["Apr"]*7.5+ interactionfactor$coef["I(date^2)"]*pred_day2^2
interactionfactor17wrong <- interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2017"])*pred_day2 + interactionfactor$coef["Apr"]*7.5+ interactionfactor$coef["I(date^2)"]*pred_day2^2
interactionfactor18wrong <- interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2018"])*pred_day2 + interactionfactor$coef["Apr"]*7.5+ interactionfactor$coef["I(date^2)"]*pred_day2^2
plot(pred_day2, interactionfactor14wrong, type="l", ylim=c(-130,30), xlim=c(0,175)) #black
points(pred_day2, interactionfactor15wrong, type="l", col=2) #red
points(pred_day2, interactionfactor16wrong, type="l", col=3) #green
points(pred_day2, interactionfactor17wrong, type="l", col=4) #blue
points(pred_day2, interactionfactor18wrong, type="l", col=5) #lightblue

plot(pred_day2, exp(interactionfactor14wrong), type="l") #black
points(pred_day2, exp(interactionfactor15wrong), type="l", col=2) #red
points(pred_day2, exp(interactionfactor16wrong), type="l", col=3) #green
points(pred_day2, exp(interactionfactor17wrong), type="l", col=4) #blue
points(pred_day2, exp(interactionfactor18wrong), type="l", col=5) #lightblue


gradient14 <- ((interactionfactor$coef["(Intercept)"] + interactionfactor$coef["date"]*170  + interactionfactor$coef["Apr"]*7.5)-(interactionfactor$coef["(Intercept)"] + interactionfactor$coef["date"]*120  + interactionfactor$coef["Apr"]*7.5))/(170-120)
gradient15 <- ((interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2015"])*170 + interactionfactor$coef["Apr"]*7.5)-(interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2015"])*120 + interactionfactor$coef["Apr"]*7.5))/(170-120)
gradient16 <- ((interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2016"])*170 + interactionfactor$coef["Apr"]*7.5)-(interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2016"])*120 + interactionfactor$coef["Apr"]*7.5))/(170-120)
gradient17 <- ((interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2017"])*170 + interactionfactor$coef["Apr"]*7.5)-(interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2017"])*120 + interactionfactor$coef["Apr"]*7.5))/(170-120)
gradient18 <- ((interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2018"])*170 + interactionfactor$coef["Apr"]*7.5)-(interactionfactor$coef["(Intercept)"] + (interactionfactor$coef["date"]+interactionfactor$coef["date:year2018"])*120 + interactionfactor$coef["Apr"]*7.5))/(170-120)

gradient14 #1.512946 
gradient15 #1.620309
gradient16 #1.634493
gradient17 #1.507683
gradient18 #1.519909




interactioncont <- glm(caterpillars~date*Apr+I(date^2)+year, family=poisson, data=all_data)
summary(interactioncont)

interactioncont5 <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*5))*pred_day + interactioncont$coef["I(date^2)"]*pred_day^2 + interactioncont$coef["Apr"]*5
interactioncont6 <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*6))*pred_day + interactioncont$coef["I(date^2)"]*pred_day^2 + interactioncont$coef["Apr"]*6
interactioncont7 <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*7))*pred_day + interactioncont$coef["I(date^2)"]*pred_day^2 + interactioncont$coef["Apr"]*7
interactioncont8 <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*8))*pred_day + interactioncont$coef["I(date^2)"]*pred_day^2 + interactioncont$coef["Apr"]*8
interactioncont9 <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*9))*pred_day + interactioncont$coef["I(date^2)"]*pred_day^2 + interactioncont$coef["Apr"]*9

plot(pred_day, interactioncont5, type="l")#, ylim=c(-1,-10))
points(pred_day, interactioncont6, type="l", col=2)
points(pred_day, interactioncont7, type="l", col=3)
points(pred_day, interactioncont8, type="l", col=4)
points(pred_day, interactioncont9, type="l", col=5)

plot(pred_day, exp(interactioncont5), type="l") #black
points(pred_day, exp(interactioncont6), type="l", col=2) #red
points(pred_day, exp(interactioncont7), type="l", col=3) #green
points(pred_day, exp(interactioncont8), type="l", col=4) #blue
points(pred_day, exp(interactioncont9), type="l", col=5) #lightblue

#looking ta elements
pred_day2 <- seq(0,175,1)
interactioncont5wrong <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*5))*pred_day2 + interactioncont$coef["I(date^2)"]*pred_day2^2 #+ interactioncont$coef["Apr"]*5
interactioncont6wrong <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*6))*pred_day2 + interactioncont$coef["I(date^2)"]*pred_day2^2 #+ interactioncont$coef["Apr"]*6
interactioncont7wrong <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*7))*pred_day2 + interactioncont$coef["I(date^2)"]*pred_day2^2 #+ interactioncont$coef["Apr"]*7
interactioncont8wrong <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*8))*pred_day2 + interactioncont$coef["I(date^2)"]*pred_day2^2 #+ interactioncont$coef["Apr"]*8
interactioncont9wrong <- interactioncont$coef["(Intercept)"] + (interactioncont$coef["date"]+(interactioncont$coef["date:Apr"]*9))*pred_day2 + interactioncont$coef["I(date^2)"]*pred_day2^2 #+ interactioncont$coef["Apr"]*9

plot(pred_day2, interactioncont5wrong, type="l")#, ylim=c(0,-60))
points(pred_day2, interactioncont6wrong, type="l", col=2)
points(pred_day2, interactioncont7wrong, type="l", col=3)
points(pred_day2, interactioncont8wrong, type="l", col=4)
points(pred_day2, interactioncont9wrong, type="l", col=5)

plot(pred_day2, exp(interactioncont5wrong), type="l", xlim=c(120,175)) #black
points(pred_day2, exp(interactioncont6wrong), type="l", col=2) #red
points(pred_day2, exp(interactioncont7wrong), type="l", col=3) #green
points(pred_day2, exp(interactioncont8wrong), type="l", col=4) #blue
points(pred_day2, exp(interactioncont9wrong), type="l", col=5) #lightblue

##interaction affects the relationship between 2 explanatory variables and seperately influences intercept 
# (not date becasue that is x but whatever it interacts with is also included independently so affects c (y = ax^2 + bx + c))
## but its all actually y=exp(ax^2 + bx + c) so the effects look different because the basic changes are in the logarithm form 
# and once reverted what was a straight line is now curved/ normal parabola becomes the cater curve etc
## becasue natural logs use base e and e is positive e^coefficient will always be >1 so y can not go beneath 0 once reverted back from log? 

## For interpreting the coefficients:
# where negative include - in the exp() and the result will be <1 so when multiplied by the intercept it is a reduction
# e.g. when adjusting for year: times the 2014 predicted abundance by the exp(coef for other year) to see abundance given by adding the coefficient to the equation because adding in exp is multiplying if not within same brackets

## quadratic coefficient, log() +/- describes orientation of curve
# larger parameter size means narrower curve (regardless of +/- as that just determines orientation), so smaller number=broader peak

### Going over how intercept/interaction effects curve AGAIN
elevation <- glm(caterpillars~date+I(date^2)+date*year+elevation, family=poisson, data=all_data)
elevationdate <- glm(caterpillars~date+I(date^2)+date*year+elevation*date, family=poisson, data=all_data)
summary(elevation)
summary(elevationdate)

pred_day <- seq(120,175,1)
elevation14L <- elevation$coef["(Intercept)"] + elevation$coef["date"]*pred_day + elevation$coef["I(date^2)"]*pred_day^2 + elevation$coef["elevation"]*50
elevation14H <- elevation$coef["(Intercept)"] + elevation$coef["date"]*pred_day + elevation$coef["I(date^2)"]*pred_day^2 + elevation$coef["elevation"]*350

plot(pred_day, elevation14L, type="l")
points(pred_day, elevation14H, type="l", col=2)
plot(pred_day, exp(elevation14L), type="l")
points(pred_day, exp(elevation14H), type="l", col=2)   ## elevation only impacting intercept: height of peak changes, date stays the same!

elevationdate14L <- elevationdate$coef["(Intercept)"] + elevationdate$coef["date"]*pred_day + elevationdate$coef["I(date^2)"]*pred_day^2 + (elevationdate$coef["date:elevation"]*150)*pred_day + elevationdate$coef["elevation"]*150 
elevationdate14H <- elevationdate$coef["(Intercept)"] + elevationdate$coef["date"]*pred_day + elevationdate$coef["I(date^2)"]*pred_day^2 + (elevationdate$coef["date:elevation"]*350)*pred_day + elevationdate$coef["elevation"]*350 

plot(pred_day, elevationdate14H, type="l")
points(pred_day, elevationdate14L, type="l", col=2)
plot(pred_day, exp(elevationdate14H), type="l")
points(pred_day, exp(elevationdate14L), type="l", col=2)   ## elevationdate interaction affects both the height and the peak date- is this what ally was saying about with a change in gradient the intercept will also change?

peakL <- -(elevationdate$coef["date"]+elevationdate$coef["date:elevation"]*150)/(2*elevationdate$coef["I(date^2)"]) # 150.9854 
peakH <- -(elevationdate$coef["date"]+elevationdate$coef["date:elevation"]*350)/(2*elevationdate$coef["I(date^2)"]) # 155.4132 
