rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
library(lme4)
library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(forcats)
library(DHARMa)
library(metafor)
library(tidyr)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#############################################
#### Temperature sliding window analysis ####
#############################################

#Cater data
cater <- read.csv("~/Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$siteday <- paste(cater$site, cater$date ,cater$year)
cater$datescaled <- scale(cater$date)[,1] #mean 146.77, sd 14.04305
cater$resid <- 1:length(cater$site)

# Load Stan site-year peak model
post <- read.csv("covar4_post.csv")
no.chains <- 3
para <- (length(colnames(post))-1)/no.chains

# Organise + combine posterior distributions
chn1 <- data.frame(chain=rep(1, length(post$X)),iteration=post[,1]) 
chn2 <- data.frame(chain=rep(2, length(post$X)),iteration=post[,1]) 
chn3 <- data.frame(chain=rep(3, length(post$X)),iteration=post[,1]) 
#chn4 <- data.frame(chain=rep(4, length(post$X)),iteration=post[,1]) 

for(i in 1:para){
  chn1[,i+2] <- post[,i*3-1] #*4-2
  chn2[,i+2] <- post[,i*3] #*4-1
  chn3[,i+2] <- post[,i*3+1] #*4
  #chn4[,i+2] <- post[,i*4+1] 
  
  colnames(chn1)[i+2] <- colnames(post)[i*3-1]
  colnames(chn2)[i+2] <- colnames(post)[i*3]
  colnames(chn3)[i+2] <- colnames(post)[i*3+1]
  #colnames(chn4)[i+2] <- colnames(post)[i*4+1]
}

remove.start <- function(x, n){  #function to keep last n number of characters in a string
  substr(x, nchar(x)-(nchar(x)-n-1), nchar(x))
}

colnames(chn1)[c(3:ncol(chn1))] <- remove.start(colnames(chn1)[c(3:ncol(chn1))],8) #removing "chain.X."
colnames(chn2)[c(3:ncol(chn2))] <- remove.start(colnames(chn2)[c(3:ncol(chn2))],8)
colnames(chn3)[c(3:ncol(chn3))] <- remove.start(colnames(chn3)[c(3:ncol(chn3))],8)
#colnames(chn4)[c(3:ncol(chn4))] <- remove.start(colnames(chn4)[c(3:ncol(chn4))],8)

stan_post <- rbind(chn1,chn2,chn3) #,chn4

#Columns in posterior df
#stan mu: 3
#stan logsigma: 4
#stan logmax: 5
#stan site mu: 16-59
#stan site logsigma: 60-103 
#stan site logmax: 104-147 
#stan year mu: 148-154
#stan year logsigma: 155-161 
#stan year logmax: 162-168
#stan siteyear mu: 169-424 
#stan siteyear logsigma: 425-680 
#stan siteyear logmax: 681-936 

#mode and vcv matrix for each SY parameter estimate
sy.effs <- data.frame(site_id=as.numeric(as.factor(cater$site)), site=cater$site,year_id=as.numeric(as.factor(cater$year)), year=cater$year,siteyear_id=as.numeric(as.factor(cater$siteyear)), siteyear=cater$siteyear)
sy.effs <- distinct(sy.effs)
sy.effs.long <- sy.effs[rep(seq_len(nrow(sy.effs)), each = 3), ]

for (i in 1:length(sy.effs$site_id)){
  
  mu.i <- mcmc(stan_post[,3] + stan_post[,sy.effs$site_id[i]+15] + stan_post[,sy.effs$year_id[i]+147] + stan_post[,sy.effs$siteyear_id[i]+168])
  logsigma.i <- mcmc(stan_post[,4] + stan_post[,sy.effs$site_id[i]+59] + stan_post[,sy.effs$year_id[i]+154] + stan_post[,sy.effs$siteyear_id[i]+424])
  logmax.i <- mcmc(stan_post[,5] + stan_post[,sy.effs$site_id[i]+103] + stan_post[,sy.effs$year_id[i]+161] + stan_post[,sy.effs$siteyear_id[i]+680])
  
  sy.effs.long$para[i*3-2] <- "mu"
  sy.effs.long$para[i*3-1] <- "logsigma"
  sy.effs.long$para[i*3] <- "logmax"
  sy.effs.long$yi[i*3-2] <- posterior.mode(mu.i)
  sy.effs.long$yi[i*3-1] <- posterior.mode(logsigma.i)
  sy.effs.long$yi[i*3] <- posterior.mode(logmax.i)
 
  df <- data.frame(mu.i=mu.i, logsigma.i=logsigma.i, logmax.i=logmax.i)
  df.m <- as.matrix(df)
  df.cov <- cov(df.m)
  
  sy.effs.long[i*3-2,9:11] <- df.cov[1,]
  sy.effs.long[i*3-1,9:11] <- df.cov[2,]
  sy.effs.long[i*3,9:11] <- df.cov[3,]
  
  mu.i <- NULL
  logsigma.i <- NULL
  logmax.i <- NULL
  df <- NULL
  df.m <- NULL
  df.cov <- NULL
}

colnames(sy.effs.long)[9:11] <- c("vi.mu","vi.logsigma","vi.logmax")
sy.effs.long$id <- rep(seq(1,length(sy.effs$siteyear),1),1, each=3)
V <- bldiag(lapply(split(sy.effs.long[,c("vi.mu","vi.logsigma","vi.logmax")], sy.effs.long$id), as.matrix))
sy.effs.long$para.mu <- ifelse(sy.effs.long$para=="mu", 1, 0)
sy.effs.long$para.logsigma <- ifelse(sy.effs.long$para=="logsigma", 1, 0)
sy.effs.long$para.logmax <- ifelse(sy.effs.long$para=="logmax", 1, 0)

#temperature
rm(post,stan_post,chn1,chn2,chn3) #just reducing memory taken up as temp is large
#temp <- read.csv("~/Dropbox/master_data/site/temperatures.csv")
#temp$siteyear <- paste(temp$site, temp$year)
#temp.sy <- temp %>%
#            group_by(siteyear) %>%
#            summarize_all(.funs = c(mean="mean"), na.rm=T) %>% 
#            as.data.frame()
#write.csv(temp.sy, "~/Documents/Stan/tempsy.csv", row.names = F)
#temp.sy[,c("site_mean","year_mean")] <- NULL
#plot(seq(1,2808,1),temp.sy[1,2:2809],pch=20)

#NAs <- c()
#for(i in 1:2808){
#  NAs[i] <- length(which(is.na(temp.sy[,i+1])==T))
#}

#temp.sy[,which(sapply(temp.sy, function(x) sum(is.na(x)))>0)] <- NULL #removing columns with NAs

#write.csv(temp.sy, "~/Documents/Stan/tempsy.csv", row.names = F)
temp.sy <- read.csv("~/Documents/Stan/tempsy.csv")

##night-time and day-time median temp for each day: 10pm-4am and 10am-4pm
#00:00  "X58_mean"               
#01:00  "X58.0416666666667_mean" 
#02:00  "X58.0833333333333_mean" 
#03:00  "X58.125_mean"           
#04:00  "X58.1666666666667_mean" 
#05:00  "X58.2083333333333_mean"
#06:00  "X58.25_mean"            
#07:00  "X58.2916666666667_mean" 
#08:00  "X58.3333333333333_mean" 
#09:00  "X58.375_mean"           
#10:00  "X58.4166666666667_mean" 
#11:00  "X58.4583333333333_mean" 
#12:00  "X58.5_mean"            
#13:00  "X58.5416666666667_mean" 
#14:00  "X58.5833333333333_mean" 
#15:00  "X58.625_mean"           
#16:00  "X58.6666666666667_mean" 
#17:00  "X58.7083333333333_mean" 
#18:00  "X58.75_mean"            
#19:00  "X58.7916666666667_mean"
#20:00  "X58.8333333333333_mean" 
#21:00  "X58.875_mean"           
#22:00  "X58.9166666666667_mean" 
#23:00  "X58.9583333333333_mean"
temp.sy[,2:11] <- NULL #starting at day58 10am

daytemp <- data.frame(siteyear=temp.sy$siteyear)
nighttemp <- data.frame(siteyear=temp.sy$siteyear)
#daytemp$X58 <- apply(temp.sy[,2:8],1,median)
#daytemp$X59 <- apply(temp.sy[,(2+24):(8+24)],1,median)
#nighttemp$X58 <- apply(temp.sy[,14:20],1,median) #name after day night starts on
#nighttemp$X59 <- apply(temp.sy[,(14+24):(20+24)],1,median)

#105 days, 104 nights  median temp per day/night

for(i in 1:105){
  daytemp[,i+1] <- apply(temp.sy[,(2+(24*(i-1))):(8+(24*(i-1)))],1,median)
  colnames(daytemp)[i+1] <- substr(colnames(temp.sy)[(2+(24*(i-1)))], 1, 4)
} 

for(i in 1:104){
  nighttemp[,i+1] <- apply(temp.sy[,(14+(24*(i-1))):(20+(24*(i-1)))],1,median)
  colnames(nighttemp)[i+1] <- substr(colnames(temp.sy)[(14+(24*(i-1)))], 1, 4)
} 

#### Sliding window
#from 14 days-60 days 1 day incriments of increase and change in start
#105 days

### DAYTIME ###
minwin <- 15 #minimum duration
maxwin <- 60 #maximum duration
firstcol <- 2 #first column with a day's temp
dat <- daytemp #data.frame of daily temps
###############

duration <- c(seq((minwin-1),(maxwin-1),15)) #duration is -1 becasue  will be start+duration
start <- c(seq(firstcol,length(colnames(dat)),15)) #column to start in
bounds <- data.frame(start=rep(start,length(duration)))
bounds$duration <- rep(duration, 1, each=length(start))
bounds$end <- bounds$start+bounds$duration
bounds <- bounds[-c(which(bounds$end>length(colnames(dat)))),]
bounds$ID <- paste(colnames(dat)[bounds$start],"_",(bounds$duration+1))
slidwin.temp <- data.frame(siteyear=dat$siteyear)
mv.wndws <- data.frame(mu.ID=bounds$ID,logsigma.ID=bounds$ID,logmax.ID=bounds$ID)
mv.wndws <- expand.grid(mv.wndws)

i=25 #X88. _ 15 X73. _ 15 X58. _ 15
slidwin.temp$temp.mu <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$mu.ID[i])]
                                  :bounds$end[which(bounds$ID==mv.wndws$mu.ID[i])]
                                  ],1,mean)
slidwin.temp$temp.logsigma <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                  :bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                  ],1,mean)
slidwin.temp$temp.logmax <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                  :bounds$end[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                  ],1,mean)

for(i in 1:length(bounds$start)){
  #mean time window
  slidwin.temp$temp.mu <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$mu.ID[i])]
                                    :bounds$end[which(bounds$ID==mv.wndws$mu.ID[i])]
                                    ],1,mean)
  slidwin.temp$temp.logsigma <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                          :bounds$end[which(bounds$ID==mv.wndws$logsigma.ID[i])]
                                          ],1,mean)
  slidwin.temp$temp.logmax <- apply(dat[,bounds$start[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                        :bounds$end[which(bounds$ID==mv.wndws$logmax.ID[i])]
                                        ],1,mean)
  #pair with cater peak
  store <- pmatch(sy.effs.long$siteyear, slidwin.temp$siteyear, duplicates.ok = T)
  sy.effs.long$temp.mu <- slidwin.temp$temp.mu[store]
  sy.effs.long$temp.logsigma <- slidwin.temp$temp.logsigma[store]
  sy.effs.long$temp.logmax <- slidwin.temp$temp.logmax[store]
  #run model
  system.time(multivar.mod <- rma.mv(yi, V, mods=~ para.mu + para.logsigma + para.logmax + temp.mu:para.mu + temp.logsigma:para.logsigma + temp.logmax:para.logmax - 1, random=list(~para|site, ~para|year), data=sy.effs.long, struct="UN", method="ML")) #, ~para|siteyear cant have 3?!
  #ML or REML??
  #system.time1= 525.729 , 8.762mins different
  #system.time2= 542.816, 9.047mins same time
  #print(multivar.mod, digits=3)
  
  #store AIC
  mv.wndws$AIC[i] <- AIC(meta.mod)
  mv.wndws$temp.mu.beta[i] <- multivar.mod$beta[4]
  mv.wndws$temp.mu.se[i] <- multivar.mod$se[4]
  mv.wndws$temp.logsigma.beta[i] <- multivar.mod$beta[5]
  mv.wndws$temp.logsigma.se[i] <- multivar.mod$se[5]
  mv.wndws$temp.logmax.beta[i] <- multivar.mod$beta[6]
  mv.wndws$temp.logmax.se[i] <- multivar.mod$se[6]

  slidwin.temp$temp.mu <- NULL
  slidwin.temp$temp.logsigma <- NULL
  slidwin.temp$temp.logmax <- NULL
  store <- NULL
  sy.effs.long$temp.mu <- NULL
  sy.effs.long$temp.logsigma <- NULL
  sy.effs.long$temp.logmax <- NULL
  multivar.mod <- NULL
  }

### NIGHTTIME ###
minwin <- 14 #minimum duration
maxwin <- 60 #maximum duration
firstcol <- 2 #first column with a day's temp
dat <- nighttemp #data.frame of daily temps     #
#################

duration <- c(seq((minwin-1),(maxwin-1),1)) #duration is -1 becasue  will be start+duration
start <- c(seq(firstcol,length(colnames(dat)),1)) #column to start in
bounds <- data.frame(start=rep(start,length(duration)))
bounds$duration <- rep(duration, 1, each=length(start))
bounds$end <- bounds$start+bounds$duration
bounds <- bounds[-c(which(bounds$end>length(colnames(dat)))),]
bounds$ID <- paste(colnames(dat)[bounds$start],"_",(bounds$duration+1))
slidwin.temp <- data.frame(siteyear=dat$siteyear)
slidwin.metrics.night <- data.frame(ID=bounds$ID)     #

for(i in 1:length(bounds$start)){
  #mean time window
  slidwin.temp$temp <- apply(dat[,bounds$start[i]:bounds$end[i]],1,mean)
  #pair with cater peak
  store <- pmatch(sy.effs$siteyear, slidwin.temp$siteyear)
  sy.effs$temp <- slidwin.temp$temp[store]
  #run model
  mod <- lmer(mu.mode ~ temp + (1|site) + (1|year), data=sy.effs) #needs to change
  #store AIC
  #slidwin.metrics.night$warn[i] <- paste(warning())
  slidwin.metrics.night$AIC[i] <- AIC(mod)                      #
  slidwin.metrics.night$coef[i] <- summary(mod)$coeff[2,1]      #
  slidwin.metrics.night$coef.se[i] <- summary(mod)$coeff[2,2]   #
  
  slidwin.temp$temp <- NULL
  store <- NULL
  sy.effs$temp <- NULL
  mod <- NULL
}




##Meta-analysis model
mod <- lmer(mu.mode ~ temp + (1|site) + (1|year), data=sy.effs)
meta.mod <- rma.mv(yi=mu.mode, V=mu.var, mods=~temp, random=list(~1|site, ~1|year, ~1|siteyear), tdist=TRUE, data=sy.effs)
#tdist=TRUE argument specifying that test statistics and confidence intervals must be based on the t-distribution
#default settings of the rma.mv function prescribe that test statistics of individual coefficients and confidence intervals 
#are based on the normal distribution (i.e., the Z distribution). using the Z distribution in assessing the significance of 
#model coefficients and in building confidence intervals around these coefficients, may lead to an increase in the number of 
#unjustified significant results. To reduce this problem, the user can apply the Knapp and Hartung (2003) adjustment to the 
#analyses by passing the argument tdist=TRUE to the rma.mv function.

#I'm not sure why the exmaple has every effectID as a random term
#Should it be REML? it automatically is..

##store AIC, coefficient and se 
slidwin.metrics$AIC <- AIC(meta.mod)
slidwin.metrics$temp.mu.beta <- multivar.mod$beta[4]
slidwin.metrics$temp.mu.se <- multivar.mod$se[4]
slidwin.metrics$temp.logsigma.beta <- multivar.mod$beta[5]
slidwin.metrics$temp.logsigma.se <- multivar.mod$se[5]
slidwin.metrics$temp.logmax.beta <- multivar.mod$beta[6]
slidwin.metrics$temp.logmax.se <- multivar.mod$se[6]


a <- seq(-4,4,0.01)
b1 <- 1 - ((a-0.01)^2)/(2*1^2)
b2 <- 0.3 - ((a-0.01)^2)/(2*1^2)
plot(a,(b1), type="l", ylim=c(-3,1), xlab=("x"), ylab=("log y"))
points(a,(b2), type="l",col=2)
plot(a,exp(b1), type="l", ylim=c(0,2.75), xlab=("x"), ylab=("y"))
points(a,exp(b2), type="l",col=2)
