rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
library(lme4)
library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(forcats)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

######################################
#### Gaussian function peak model ####
######################################

#load data
cater <- read.csv("~/Dropbox/master_data/inverts/Branch_Beating.csv")
cater$year <- as.factor(cater$year) #year as a factor
cater$tree.species<- revalue(cater$tree.species, c("?"="OthDecid", "Cherry"="OthDecid", "Aspen"="OthDecid", "Chestnut"="OthDecid", "Lime"="OthDecid", "Field maple"="OthDecid", "Damson"="OthDecid", "Whitebeam"="OthDecid"))
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$siteday <- paste(cater$site, cater$date ,cater$year)
cater$resid <- 1:length(cater$site)


### Want model with gaussian function with random terms for site, siteday, treeID, recorder and residual
freq_mod <- glmer(caterpillars ~ datescaled + I(datescaled^2) + (1|site) + (1|siteday) + (1|treeID) + (1|recorder) + (1|resid), data = cater19, family = poisson)
summary(freq_mod)

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
mcmc_mod <- MCMCglmm(caterpillars ~ datescaled + I(datescaled^2), random=~site + siteday + treeID + recorder, data = cater19, family = "poisson", prior=prior, nitt=30000,burnin=15000, pr=TRUE)


## no QR decomposition
write("data{
      int<lower=0> N;
      vector[N] date;
      int<lower=0> N_site;
      int<lower=0> N_siteday;
      int<lower=0> N_treeID;
      int<lower=0> N_rec;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_siteday> siteday_id[N]; 
      int<lower=0,upper=N_treeID> treeID_id[N]; 
      int<lower=0,upper=N_rec> rec_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real<lower=0> sigma;
      real logmax;
      vector[N_site] site_scaled; 
      vector[N_siteday] siteday_scaled; 
      vector[N_treeID] treeID_scaled; 
      vector[N_rec] rec_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_siteday;
      real<lower=0> sigma_treeID;
      real<lower=0> sigma_rec;
      real<lower=0> sigma_resid;
      }
      
      model{
      vector[N] y_1;
      vector[N_site] site_effects; 
      vector[N_siteday] siteday_effects; 
      vector[N_treeID] treeID_effects; 
      vector[N_rec] rec_effects; 
      vector[N_resid] resid_effects;
      
      site_effects = sigma_site * site_scaled; 
      siteday_effects = sigma_siteday * siteday_scaled; 
      treeID_effects = sigma_treeID * treeID_scaled; 
      rec_effects = sigma_rec * rec_scaled; 
      resid_effects = sigma_resid * resid_scaled;
      
      y_1 = logmax - ((date-mu) .* (date-mu)) ./ (2*sigma*sigma);
      for(i in 1:N){y_1[i] += site_effects[site_id[i]] + siteday_effects[siteday_id[i]] + treeID_effects[treeID_id[i]] + rec_effects[rec_id[i]] + resid_effects[resid_id[i]];
      }
      
      y ~ poisson_log(y_1);
      
      site_scaled ~ normal(0,1); 
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      sigma ~ normal(0,10);
      logmax ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_siteday ~ cauchy(0,10);
      sigma_treeID ~ cauchy(0,10);
      sigma_rec ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);
      }",
      
      "gaussian1.stan"
)

stanc("gaussian1.stan")

stanModelgaussian <- stan_model("gaussian1.stan")

stan_data1 <-list(
  N=nrow(cater19),
  date=cater19$datescaled,
  N_site=length(unique(cater19$site)),
  N_siteday=length(unique(cater19$siteday)),
  N_treeID=length(unique(cater19$treeID)),
  N_rec=length(unique(cater19$recorder)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  site_id=as.numeric(as.factor(cater19$site)),
  siteday_id=as.numeric(as.factor(cater19$siteday)),
  treeID_id=as.numeric(as.factor(cater19$treeID)),
  rec_id=as.numeric(as.factor(cater19$recorder)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_modG1 <- sampling(object=stanModelgaussian, data=stan_data1, iter=2000, chains=4, pars=c("mu", "sigma","logmax","sigma_site","sigma_siteday","sigma_treeID","sigma_rec","sigma_resid", "site_scaled"))

summary(stan_modG1)$summary[,c(1,3,4,8,9,10)]
rstan::traceplot(stan_modG1)

plot(cater19$datescaled,cater19$caterpillars)
datameans <- data.frame(datescaled=cater19$datescaled,cater=cater19$caterpillars)
datameans <- aggregate(.~datescaled, datameans, mean)
mod_betas <- summary(stan_modG1)$summary[c(1,2,3),1]
freq_betas <- summary(freq_mod)$coefficients[c(1,2,3),1]
datameans$pred <- mod_betas[3] - ((datameans$datescaled-mod_betas[1])^2)/(2*mod_betas[2]^2)
datameans$freqpred <- freq_betas[1] + freq_betas[2]*datameans$datescaled + freq_betas[3]*datameans$datescaled^2
plot(datameans$datescaled,datameans$cater)
points(datameans$datescaled, exp(datameans$pred), type="l")
points(datameans$datescaled, exp(datameans$freqpred), type="l", col=2)

site_scaleffects <- summary(stan_modG1)$summary[9:52,1]
site_effects <- site_scaleffects*summary(stan_modG1)$summary["sigma_site",1]
cater19$sitenum <- as.numeric(as.factor(cater19$site))
sites <- data.frame(site=cater19$site,sitenum=cater19$sitenum)
sites <- aggregate(.~site, sites, mean)
sites$stanef <- site_effects
sites$freqef <- ranef(freq_mod)$site[1:44,1]
mcmc_sites <- mcmc_mod$Sol[,4:47]
for(i in 1:ncol(mcmc_sites)){sites$mcmcef[i] <- posterior.mode(mcmc_sites[,i])}
plot(sites$stanef, sites$freqef)
plot(sites$stanef, sites$mcmcef)
abline(0,1, col=2)


##random intercepts for all parameters
freq_mod2 <- glmer(caterpillars ~ datescaled + I(datescaled^2) + (1+datescaled+I(datescaled^2)|site) + (1|siteday) + (1|treeID) + (1|recorder) + (1|resid), data = cater19, family = poisson)
summary(freq_mod2)

k<-10000
prior2<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
mcmc_mod3 <- MCMCglmm(caterpillars ~ datescaled + I(datescaled^2), random=~us(1+datescaled+I(datescaled^2)):site + year + siteyear + siteday + treeID + recorder, data = cater, family = "poisson", prior=prior2, nitt=25000,burnin=15000, pr=TRUE)

write("data{
      int<lower=0> N;
      vector[N] date;
      int<lower=0> N_site;
      int<lower=0> N_siteyear;
      int<lower=0> N_year;
      int<lower=0> N_siteday;
      int<lower=0> N_treeID;
      int<lower=0> N_rec;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_siteyear> siteyear_id[N]; 
      int<lower=0,upper=N_year> year_id[N]; 
      int<lower=0,upper=N_siteday> siteday_id[N]; 
      int<lower=0,upper=N_treeID> treeID_id[N]; 
      int<lower=0,upper=N_rec> rec_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real<lower=0> logsigma;
      real logmax;
      vector[N_site] site_mu_scaled; 
      vector[N_site] site_logsigma_scaled;
      vector[N_site] site_logmax_scaled;
      vector[N_siteyear] siteyear_scaled; 
      vector[N_year] year_scaled; 
      vector[N_siteday] siteday_scaled; 
      vector[N_treeID] treeID_scaled; 
      vector[N_rec] rec_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sd_site_mu;
      real<lower=0> sd_site_logsigma;
      real<lower=0> sd_site_logmax;
      real<lower=0> sd_siteyear;
      real<lower=0> sd_year;
      real<lower=0> sd_siteday;
      real<lower=0> sd_treeID;
      real<lower=0> sd_rec;
      real<lower=0> sd_resid;
      }
      
      model{
      vector[N] y_1;
      vector[N_site] site_mu_effects = sd_site_mu * site_mu_scaled; 
      vector[N_site] site_logsigma_effects = sd_site_logsigma * site_logsigma_scaled;
      vector[N_site] site_logmax_effects = sd_site_logmax * site_logmax_scaled;
      vector[N_siteyear] siteyear_effects = sd_siteyear * siteyear_scaled; 
      vector[N_year] year_effects = sd_year * year_scaled; 
      vector[N_siteday] siteday_effects = sd_siteday * siteday_scaled; 
      vector[N_treeID] treeID_effects = sd_treeID * treeID_scaled; 
      vector[N_rec] rec_effects = sd_rec * rec_scaled; 
      vector[N_resid] resid_effects = sd_resid * resid_scaled;
      
      y_1 = logmax - ((date-(mu+site_mu_effects[site_id])) .* (date-(mu+site_mu_effects[site_id]))) ./ (2.*(exp(logsigma+site_logsigma_effects[site_id])).*(exp(logsigma+site_logsigma_effects[site_id]))) 
       + site_logmax_effects[site_id] + siteyear_effects[siteyear_id] + year_effects[year_id] + siteday_effects[siteday_id] 
       + treeID_effects[treeID_id] + rec_effects[rec_id] + resid_effects[resid_id];
      
      
      y ~ poisson_log(y_1);
      
      site_mu_scaled ~ normal(0,1); 
      site_logsigma_scaled ~ normal(0,1); 
      site_logmax_scaled ~ normal(0,1); 
      siteyear_scaled ~ normal(0,1);
      year_scaled ~ normal(0,1);
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      logsigma ~ cauchy(0,10);
      logmax ~ normal(0,10);
      sd_site_mu ~ cauchy(0,10);
      sd_site_logsigma ~ cauchy(0,10);
      sd_site_logmax ~ cauchy(0,10);
      sd_siteyear ~ cauchy(0,10);
      sd_year ~ cauchy(0,10);
      sd_siteday ~ cauchy(0,10);
      sd_treeID ~ cauchy(0,10);
      sd_rec ~ cauchy(0,10);
      sd_resid ~ cauchy(0,10);
      }",
      
      "gaussiansite1.stan"
)

stanc("gaussiansite1.stan")
stanModelGsite1 <- stan_model("gaussiansite1.stan")

stan_data3 <-list(
  N=nrow(cater),
  date=cater$datescaled,
  N_site=length(unique(cater$site)),
  N_siteyear=length(unique(cater$siteyear)),
  N_year=length(unique(cater$year)),
  N_siteday=length(unique(cater$siteday)),
  N_treeID=length(unique(cater$treeID)),
  N_rec=length(unique(cater$recorder)),
  N_resid=length(unique(cater$resid)),
  y=cater$caterpillars,
  site_id=as.numeric(as.factor(cater$site)),
  siteyear_id=as.numeric(as.factor(cater$siteyear)),
  year_id=as.numeric(as.factor(cater$year)),
  siteday_id=as.numeric(as.factor(cater$siteday)),
  treeID_id=as.numeric(as.factor(cater$treeID)),
  rec_id=as.numeric(as.factor(cater$recorder)),
  resid_id=as.numeric(as.factor(cater$resid))
)


stan_modGS1 <- sampling(object=stanModelGsite1, data=stan_data3, iter=200, chains=4, 
                       pars=c("mu", "logsigma","logmax","sd_site_mu","sd_site_logsigma","sd_site_logmax","sd_siteyear","sd_year","sd_resid","site_mu_scaled","site_logsigma_scaled",
                              "site_logmax_scaled","resid_scaled")) 

write.csv(summary(stan_modGS1)$summary, "~/Documents/Stan/GS1_sum.csv")
GS1_post <- extract(stan_modGS1, permuted=FALSE)
#could limit which columns here...
write.csv(GS1_post, "~/Documents/Stan/GS1_post.csv")

#summary(stan_modG3)$summary[,c(1,3,4,8,9,10)]
#rstan::traceplot(stan_modG3)

GS1_sum <- read.csv("GS1_sum.csv")
GS1_post <- read.csv("GS1_post.csv")

##compare mcmcglmm and stan peaks

hist(mcmc_mod3$Sol[,1])
abline(v=posterior.mode(mcmc_mod3$Sol[,1]),col=2)
abline(v=mean(mcmc_mod3$Sol[,1]),col=3)

#mcmc site int: 4-47 cols
#mcmc site date: 48-91 cols
#mcmc site date^2: 92-135 cols

cater$sitenum <- as.numeric(as.factor(cater$site))
sites <- data.frame(site=cater$site,sitenum=cater$sitenum)
sites <- aggregate(.~site, sites, mean)
for (i in 1:44){
  sites$mcmc_int[i] <- posterior.mode(mcmc_mod3$Sol[,1]+mcmc_mod3$Sol[,i+3])
  sites$mcmc_d[i] <- posterior.mode(mcmc_mod3$Sol[,2]+mcmc_mod3$Sol[,i+47])
  sites$mcmc_d2[i] <- posterior.mode(mcmc_mod3$Sol[,3]+mcmc_mod3$Sol[,i+ 91])
}

#stan site mu: 10-53 row, sd:4 
#stan site logsigma: 54-97 row, sd:5
#stan site logmax: 98-141 row, sd:6

for (i in 1:44){
  sites$stan_mu[i] <- GS1_sum[1,2]+GS1_sum[i+9,2]*GS1_sum[4,2]
  sites$stan_logsigma[i] <- GS1_sum[2,2]+GS1_sum[i+53,2]*GS1_sum[5,2]
  sites$stan_logmax[i] <- GS1_sum[3,2]+GS1_sum[i+97,2]*GS1_sum[6,2]
}

mcmcpeaks <- data.frame(dscal=seq(-2.2, 2.1, 0.01))
stanpeaks <- data.frame(dscal=seq(-2.2, 2.1, 0.01))
for(i in 1:44){
  for (j in 1:length(mcmcpeaks$dscal)){
  mcmcpeaks[j,i+1] <- sites$mcmc_int[i] + sites$mcmc_d[i]*mcmcpeaks$dscal[j] + sites$mcmc_d2[i]*mcmcpeaks$dscal[j]^2
  }
  }
colnames(mcmcpeaks)[c(2:45)] <- sites$site

for(i in 1:44){
  for (j in 1:length(stanpeaks$dscal)){
    stanpeaks[j,i+1] <- sites$stan_logmax[i] - ((stanpeaks$dscal[j]-sites$stan_mu[i])^2)/(2*exp(sites$stan_logsigma[i])^2)
  }
}
colnames(stanpeaks)[c(2:45)] <- sites$site

plot(cater$datescaled,cater$caterpillars)
datameans <- data.frame(datescaled=cater$datescaled,cater=cater$caterpillars)
datameans <- aggregate(.~datescaled, datameans, mean)
stan_betas <- GS1_sum[1:3,2]
mcmc_betas <- c(mean(mcmc_mod3$Sol[,1]),mean(mcmc_mod3$Sol[,2]),mean(mcmc_mod3$Sol[,3]))
datameans$stanpred <- stan_betas[3] - ((datameans$datescaled-stan_betas[1])^2)/(2*exp(stan_betas[2])^2)
datameans$mcmcpred <- mcmc_betas[1] + mcmc_betas[2]*datameans$datescaled + mcmc_betas[3]*datameans$datescaled^2
plot(datameans$datescaled,datameans$cater)
points(datameans$datescaled, exp(datameans$stanpred), type="l")
points(datameans$datescaled, exp(datameans$mcmcpred), type="l", col=2)

plot(mcmcpeaks$dscal, exp(mcmcpeaks$DOW), type="l", ylim=c(0,0.5))
points(stanpeaks$dscal, exp(stanpeaks$DOW), type="l", col=2)

for(i in 1:44){
  
  X <- subset(cater, site==sites$site[i])
  X <- data.frame(date=X$datescaled, cater=X$caterpillars)
  X <- aggregate(.~date, X, mean)
  
  mypath <- file.path(paste("~/Documents/Stan/ComparePeaks/",sites$site[i],".pdf"))
  pdf(file = mypath, width=8, height=6)
  
  plot(X$date, X$cater, xlim=c(-2.2, 2.1),ylim=c(0,max(X$cater)), xlab="Date scaled", ylab="Abund")
  points(mcmcpeaks$dscal, exp(mcmcpeaks[,i+1]),type="l")
  points(stanpeaks$dscal, exp(stanpeaks[,i+1]), type="l", col=2)
  text(x=-2,y=max(X$cater), "MCMCglmm", cex=0.8)
  text(x=-2,y=max(X$cater)-0.04*max(X$cater), "Stan", cex=0.8, col=2)
  title(main=sites$site[i])
  
  dev.off()
  
}


#### Site-year peaks! ####
k<-10000
prior3<-list(R=list(V=1,nu=0.002),
             G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                    G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
mcmc_mod4 <- MCMCglmm(caterpillars ~ datescaled + I(datescaled^2), random=~us(1+datescaled+I(datescaled^2)):site + us(1+datescaled+I(datescaled^2)):year + us(1+datescaled+I(datescaled^2)):siteyear + siteday + treeID + recorder, data = cater, family = "poisson", prior=prior3, nitt=50000, burnin=25000, thin=100, pr=TRUE)
save(mcmc_mod4, "~/Documents/Stan/mcmc_mod4.RData")

write("data{
      int<lower=0> N;
      vector[N] date;
      int<lower=0> N_site;
      int<lower=0> N_year;
      int<lower=0> N_siteyear;
      int<lower=0> N_siteday;
      int<lower=0> N_treeID;
      int<lower=0> N_rec;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_year> year_id[N]; 
      int<lower=0,upper=N_siteyear> siteyear_id[N]; 
      int<lower=0,upper=N_siteday> siteday_id[N]; 
      int<lower=0,upper=N_treeID> treeID_id[N]; 
      int<lower=0,upper=N_rec> rec_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real<lower=0> logsigma;
      real logmax;
      vector[N_site] site_mu_scaled; 
      vector[N_site] site_logsigma_scaled;
      vector[N_site] site_logmax_scaled;
      vector[N_year] year_mu_scaled; 
      vector[N_year] year_logsigma_scaled;
      vector[N_year] year_logmax_scaled;
      vector[N_siteyear] siteyear_mu_scaled; 
      vector[N_siteyear] siteyear_logsigma_scaled;
      vector[N_siteyear] siteyear_logmax_scaled;
      vector[N_siteday] siteday_scaled; 
      vector[N_treeID] treeID_scaled; 
      vector[N_rec] rec_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sd_site_mu;
      real<lower=0> sd_site_logsigma;
      real<lower=0> sd_site_logmax;
      real<lower=0> sd_year_mu;
      real<lower=0> sd_year_logsigma;
      real<lower=0> sd_year_logmax;
      real<lower=0> sd_siteyear_mu;
      real<lower=0> sd_siteyear_logsigma;
      real<lower=0> sd_siteyear_logmax;
      real<lower=0> sd_siteday;
      real<lower=0> sd_treeID;
      real<lower=0> sd_rec;
      real<lower=0> sd_resid;
      }
      
      model{
      vector[N] y_1;
      vector[N_site] site_mu_effects = sd_site_mu * site_mu_scaled; 
      vector[N_site] site_logsigma_effects = sd_site_logsigma * site_logsigma_scaled;
      vector[N_site] site_logmax_effects = sd_site_logmax * site_logmax_scaled;
      vector[N_year] year_mu_effects = sd_year_mu * year_mu_scaled; 
      vector[N_year] year_logsigma_effects = sd_year_logsigma * year_logsigma_scaled;
      vector[N_year] year_logmax_effects = sd_year_logmax * year_logmax_scaled;
      vector[N_siteyear] siteyear_mu_effects = sd_siteyear_mu * siteyear_mu_scaled; 
      vector[N_siteyear] siteyear_logsigma_effects = sd_siteyear_logsigma * siteyear_logsigma_scaled;
      vector[N_siteyear] siteyear_logmax_effects = sd_siteyear_logmax * siteyear_logmax_scaled;
      vector[N_siteday] siteday_effects = sd_siteday * siteday_scaled; 
      vector[N_treeID] treeID_effects = sd_treeID * treeID_scaled; 
      vector[N_rec] rec_effects = sd_rec * rec_scaled; 
      vector[N_resid] resid_effects = sd_resid * resid_scaled;
      
      y_1 = logmax - ((date-(mu+site_mu_effects[site_id]+year_mu_effects[year_id]+siteyear_mu_effects[siteyear_id])) .* 
            (date-(mu+site_mu_effects[site_id]+year_mu_effects[year_id]+siteyear_mu_effects[siteyear_id]))) ./ 
            (2.*(exp(logsigma+site_logsigma_effects[site_id]+year_logsigma_effects[year_id]+siteyear_logsigma_effects[siteyear_id])).*
            (exp(logsigma+site_logsigma_effects[site_id]+year_logsigma_effects[year_id]+siteyear_logsigma_effects[siteyear_id]))) 
            + site_logmax_effects[site_id] + year_logmax_effects[year_id] 
            + siteyear_logmax_effects[siteyear_id] + siteday_effects[siteday_id] 
            + treeID_effects[treeID_id] + rec_effects[rec_id] + resid_effects[resid_id];
      
      
      y ~ poisson_log(y_1);
      
      site_mu_scaled ~ normal(0,1); 
      site_logsigma_scaled ~ normal(0,1); 
      site_logmax_scaled ~ normal(0,1); 
      year_mu_scaled ~ normal(0,1); 
      year_logsigma_scaled ~ normal(0,1); 
      year_logmax_scaled ~ normal(0,1); 
      siteyear_mu_scaled ~ normal(0,1); 
      siteyear_logsigma_scaled ~ normal(0,1); 
      siteyear_logmax_scaled ~ normal(0,1); 
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      logsigma ~ cauchy(0,10);
      logmax ~ normal(0,10);
      sd_site_mu ~ cauchy(0,10);
      sd_site_logsigma ~ cauchy(0,10);
      sd_site_logmax ~ cauchy(0,10);
      sd_year_mu ~ cauchy(0,10);
      sd_year_logsigma ~ cauchy(0,10);
      sd_year_logmax ~ cauchy(0,10);
      sd_siteyear_mu ~ cauchy(0,10);
      sd_siteyear_logsigma ~ cauchy(0,10);
      sd_siteyear_logmax ~ cauchy(0,10);
      sd_siteday ~ cauchy(0,10);
      sd_treeID ~ cauchy(0,10);
      sd_rec ~ cauchy(0,10);
      sd_resid ~ cauchy(0,10);
      }",
      
      "gaussiansiteyear.stan"
)

stanc("gaussiansiteyear.stan")
stanModelGSY <- stan_model("gaussiansiteyear.stan")

stan_data3 <-list(
  N=nrow(cater),
  date=cater$datescaled,
  N_site=length(unique(cater$site)),
  N_siteyear=length(unique(cater$siteyear)),
  N_year=length(unique(cater$year)),
  N_siteday=length(unique(cater$siteday)),
  N_treeID=length(unique(cater$treeID)),
  N_rec=length(unique(cater$recorder)),
  N_resid=length(unique(cater$resid)),
  y=cater$caterpillars,
  site_id=as.numeric(as.factor(cater$site)),
  siteyear_id=as.numeric(as.factor(cater$siteyear)),
  year_id=as.numeric(as.factor(cater$year)),
  siteday_id=as.numeric(as.factor(cater$siteday)),
  treeID_id=as.numeric(as.factor(cater$treeID)),
  rec_id=as.numeric(as.factor(cater$recorder)),
  resid_id=as.numeric(as.factor(cater$resid))
)

#stan_modGSY was run for 2000 but at least one chain didnt do very well.. second run did better
stan_modGSY2 <- sampling(object=stanModelGSY, data=stan_data3, iter=500, chains=4, 
                        pars=c("mu", "logsigma","logmax","sd_site_mu","sd_site_logsigma","sd_site_logmax",
                               "sd_year_mu","sd_year_logsigma","sd_year_logmax",
                               "sd_siteyear_mu","sd_siteyear_logsigma","sd_siteyear_logmax",
                               "sd_resid","site_mu_scaled","site_logsigma_scaled","site_logmax_scaled",
                               "year_mu_scaled","year_logsigma_scaled","year_logmax_scaled",
                               "siteyear_mu_scaled","siteyear_logsigma_scaled","siteyear_logmax_scaled","resid_scaled")) 

write.csv(summary(stan_modGSY2)$summary, "~/Documents/Stan/GSY2_sum.csv")
GSY2_post <- extract(stan_modGSY2, permuted=FALSE)
#could limit which columns here...(each chain is different column)
write.csv(GSY2_post, "~/Documents/Stan/GSY2_post.csv")

GSY3_sum <- read.csv("GSY3_sum.csv")
GSY3_post <- read.csv("GSY3_post.csv")

##compare mcmcglmm and stan peaks ##

##Stan
#stan mu: 1
#stan logsigma: 2
#stan logmax: 3
#stan site mu: 14-57 row, sd:4 
#stan site logsigma: 58-101 row, sd:5
#stan site logmax: 102-145 row, sd:6
#stan year mu: 146-152 row, sd:7
#stan year logsigma: 153-159 row, sd:8
#stan year logmax: 160-166 row, sd:9
#stan siteyear mu: 167-422 row, sd:10
#stan siteyear logsigma: 423-678 row, sd:11
#stan siteyear logmax: 679-934 row, sd:12

sy.effs <- data.frame(site_id=as.numeric(as.factor(cater$site)), site=cater$site,year_id=as.numeric(as.factor(cater$year)), year=cater$year,siteyear_id=as.numeric(as.factor(cater$siteyear)), siteyear=cater$siteyear)
sy.effs <- distinct(sy.effs)

for (i in 1:length(sy.effs$site_id)){
  sy.effs$mu[i] <- GSY2_sum[1,2] + GSY2_sum[sy.effs$site_id[i]+13,2]*GSY2_sum[4,2] + GSY2_sum[sy.effs$year_id[i]+145,2]*GSY2_sum[7,2] + GSY2_sum[sy.effs$siteyear_id[i]+166,2]*GSY2_sum[10,2]
  sy.effs$logsigma[i] <- GSY2_sum[2,2] + GSY2_sum[sy.effs$site_id[i]+57,2]*GSY2_sum[5,2] + GSY2_sum[sy.effs$year_id[i]+152,2]*GSY2_sum[8,2] + GSY2_sum[sy.effs$siteyear_id[i]+422,2]*GSY2_sum[11,2]
  sy.effs$logmax[i] <- GSY2_sum[3,2] + GSY2_sum[sy.effs$site_id[i]+101,2]*GSY2_sum[6,2] + GSY2_sum[sy.effs$year_id[i]+159,2]*GSY2_sum[9,2] + GSY2_sum[sy.effs$siteyear_id[i]+678,2]*GSY2_sum[12,2]
}

stanpeaks <- data.frame(dscal=seq(-2.2, 2.1, 0.01))
for(i in 1:256){
  for (j in 1:length(stanpeaks$dscal)){
    stanpeaks[j,i+1] <- sy.effs$logmax[i] - ((stanpeaks$dscal[j]-sy.effs$mu[i])^2)/(2*exp(sy.effs$logsigma[i])^2)
  }
}
colnames(stanpeaks)[c(2:257)] <- sy.effs$siteyear

##MCMCglmm
mcmc_mod4 <- read.csv("~/Documents/Stan/mcmc_mod4_post.csv")
#mcmc int: 2
#mcmc date: 3
#mcmc date^2: 4
#mcmc site int: 5-48 cols
#mcmc site date: 49-92 cols
#mcmc site date^2: 93-136 cols
#mcmc year int: 137-143 cols
#mcmc year date: 144-150 cols
#mcmc year date^2: 151-157 cols
#mcmc siteyear int: 158-413 cols
#mcmc siteyear date: 414-669 cols
#mcmc siteyear date^2: 670-925 cols

for (i in 1:256){
  sy.effs$mcmc_int[i] <- posterior.mode(mcmc(mcmc_mod4[,2]+mcmc_mod4[,sy.effs$site_id[i]+4]+mcmc_mod4[,sy.effs$year_id[i]+136]+mcmc_mod4[,sy.effs$siteyear_id[i]+157]))
  sy.effs$mcmc_d[i] <- posterior.mode(mcmc(mcmc_mod4[,3]+mcmc_mod4[,sy.effs$site_id[i]+48]+mcmc_mod4[,sy.effs$year_id[i]+143]+mcmc_mod4[,sy.effs$siteyear_id[i]+413]))
  sy.effs$mcmc_d2[i] <- posterior.mode(mcmc(mcmc_mod4[,4]+mcmc_mod4[,sy.effs$site_id[i]+92]+mcmc_mod4[,sy.effs$year_id[i]+150]+mcmc_mod4[,sy.effs$siteyear_id[i]+669]))
}

mcmcpeaks <- data.frame(dscal=seq(-2.2, 2.1, 0.01))

for(i in 1:256){
  for (j in 1:length(mcmcpeaks$dscal)){
    mcmcpeaks[j,i+1] <- sy.effs$mcmc_int[i] + sy.effs$mcmc_d[i]*mcmcpeaks$dscal[j] + sy.effs$mcmc_d2[i]*mcmcpeaks$dscal[j]^2
  }
}
colnames(mcmcpeaks)[c(2:257)] <- sy.effs$siteyear

for(i in 1:256){
  
  Y <- subset(cater, siteyear==sy.effs$siteyear[i])
  X <- data.frame(date=Y$datescaled, cater=Y$caterpillars)
  X <- aggregate(.~date, X, mean)
  
  mypath <- file.path(paste("~/Documents/Stan/ComparePeaks/",sy.effs$siteyear[i],".pdf"))
  pdf(file = mypath, width=8, height=6)
  
  plot(X$date, X$cater, xlim=c(-2.2, 2.1),ylim=c(0,max(X$cater)+0.05), xlab="Date scaled", ylab="Abund")
  points(mcmcpeaks$dscal, exp(mcmcpeaks[,i+1]),type="l")
  points(stanpeaks$dscal, exp(stanpeaks[,i+1]), type="l", col=2)
  text(x=-2,y=max(X$cater)+0.05, "MCMCglmm", cex=0.8)
  text(x=-2,y=(max(X$cater)+0.05)-0.04*(max(X$cater)+0.05), "Stan", cex=0.8, col=2)
  text(x=-2,y=(max(X$cater)+0.05)-0.08*(max(X$cater)+0.05),paste("Total cater=", sum(Y$caterpillars)), cex=0.8)
  text(x=-2,y=(max(X$cater)+0.05)-0.12*(max(X$cater)+0.05),paste("Samples=", length(Y$caterpillars)), cex=0.8)
  title(main=sy.effs$siteyear[i])
  
  dev.off()
  
}

# Stan posterior distributions
chn1 <- data.frame(chain=rep(1, 250),iteration=GSY2_post[,1]) 
chn2 <- data.frame(chain=rep(2, 250),iteration=GSY2_post[,1]) 
chn3 <- data.frame(chain=rep(3, 250),iteration=GSY2_post[,1]) 
chn4 <- data.frame(chain=rep(4, 250),iteration=GSY2_post[,1]) 
for(i in 1:32150){
  chn1[,i+2] <- GSY2_post[,i*4-2]
  chn2[,i+2] <- GSY2_post[,i*4-1]
  chn3[,i+2] <- GSY2_post[,i*4]
  chn4[,i+2] <- GSY2_post[,i*4+1]
  
  colnames(chn1)[i+2] <- colnames(GSY2_post)[i*4-2]
  colnames(chn2)[i+2] <- colnames(GSY2_post)[i*4-1]
  colnames(chn3)[i+2] <- colnames(GSY2_post)[i*4]
  colnames(chn4)[i+2] <- colnames(GSY2_post)[i*4+1]
}

remove.start <- function(x, n){  #function to keep last n number of characters in a string
  substr(x, nchar(x)-(nchar(x)-n-1), nchar(x))
}

colnames(chn1)[c(3:ncol(chn1))] <- remove.start(colnames(chn1)[c(3:ncol(chn1))],8) #removing chain.X.
colnames(chn2)[c(3:ncol(chn2))] <- remove.start(colnames(chn2)[c(3:ncol(chn2))],8)
colnames(chn3)[c(3:ncol(chn3))] <- remove.start(colnames(chn3)[c(3:ncol(chn3))],8)
colnames(chn4)[c(3:ncol(chn4))] <- remove.start(colnames(chn4)[c(3:ncol(chn4))],8)

stan_post <- rbind(chn1,chn2,chn3,chn4)
#using posterior.mode from mcmcglmm

#Stan columns
#stan mu: 3
#stan logsigma: 4
#stan logmax: 5
#stan site mu: 16-59 row, sd:6 
#stan site logsigma: 60-103 row, sd:7
#stan site logmax: 104-147 row, sd:8
#stan year mu: 148-154 row, sd:9
#stan year logsigma: 155-161 row, sd:10
#stan year logmax: 162-168 row, sd:11
#stan siteyear mu: 169-424 row, sd:12
#stan siteyear logsigma: 425-680 row, sd:13
#stan siteyear logmax: 681-936 row, sd:14

for (i in 1:length(sy.effs$site_id)){
  sy.effs$mu.mode[i] <- posterior.mode(mcmc(stan_post[,3] + stan_post[,sy.effs$site_id[i]+15]*stan_post[,6] + stan_post[,sy.effs$year_id[i]+147]*stan_post[,9] + stan_post[,sy.effs$siteyear_id[i]+168]*stan_post[,12]))
  sy.effs$mu.lci[i] <- HPDinterval(mcmc(stan_post[,3] + stan_post[,sy.effs$site_id[i]+15]*stan_post[,6] + stan_post[,sy.effs$year_id[i]+147]*stan_post[,9] + stan_post[,sy.effs$siteyear_id[i]+168]*stan_post[,12]))[1]
  sy.effs$mu.uci[i] <- HPDinterval(mcmc(stan_post[,3] + stan_post[,sy.effs$site_id[i]+15]*stan_post[,6] + stan_post[,sy.effs$year_id[i]+147]*stan_post[,9] + stan_post[,sy.effs$siteyear_id[i]+168]*stan_post[,12]))[2]
  sy.effs$mu.var[i] <- var(stan_post[,3] + stan_post[,sy.effs$site_id[i]+15]*stan_post[,6] + stan_post[,sy.effs$year_id[i]+147]*stan_post[,9] + stan_post[,sy.effs$siteyear_id[i]+168]*stan_post[,12])
  
  sy.effs$logsigma.mode[i] <- posterior.mode(mcmc(stan_post[,4] + stan_post[,sy.effs$site_id[i]+59]*stan_post[,7] + stan_post[,sy.effs$year_id[i]+154]*stan_post[,10] + stan_post[,sy.effs$siteyear_id[i]+424]*stan_post[,13]))
  sy.effs$logsigma.lci[i] <- HPDinterval(mcmc(stan_post[,4] + stan_post[,sy.effs$site_id[i]+59]*stan_post[,7] + stan_post[,sy.effs$year_id[i]+154]*stan_post[,10] + stan_post[,sy.effs$siteyear_id[i]+424]*stan_post[,13]))[1]
  sy.effs$logsigma.uci[i] <- HPDinterval(mcmc(stan_post[,4] + stan_post[,sy.effs$site_id[i]+59]*stan_post[,7] + stan_post[,sy.effs$year_id[i]+154]*stan_post[,10] + stan_post[,sy.effs$siteyear_id[i]+424]*stan_post[,13]))[2]
  sy.effs$logsigma.var[i] <- var(stan_post[,4] + stan_post[,sy.effs$site_id[i]+59]*stan_post[,7] + stan_post[,sy.effs$year_id[i]+154]*stan_post[,10] + stan_post[,sy.effs$siteyear_id[i]+424]*stan_post[,13])
  
  sy.effs$logmax.mode[i] <- posterior.mode(mcmc(stan_post[,5] + stan_post[,sy.effs$site_id[i]+103]*stan_post[,8] + stan_post[,sy.effs$year_id[i]+161]*stan_post[,11] + stan_post[,sy.effs$siteyear_id[i]+680]*stan_post[,14]))
  sy.effs$logmax.lci[i] <- HPDinterval(mcmc(stan_post[,5] + stan_post[,sy.effs$site_id[i]+103]*stan_post[,8] + stan_post[,sy.effs$year_id[i]+161]*stan_post[,11] + stan_post[,sy.effs$siteyear_id[i]+680]*stan_post[,14]))[1]
  sy.effs$logmax.uci[i] <- HPDinterval(mcmc(stan_post[,5] + stan_post[,sy.effs$site_id[i]+103]*stan_post[,8] + stan_post[,sy.effs$year_id[i]+161]*stan_post[,11] + stan_post[,sy.effs$siteyear_id[i]+680]*stan_post[,14]))[2]
  sy.effs$logmax.var[i] <- var(stan_post[,5] + stan_post[,sy.effs$site_id[i]+103]*stan_post[,8] + stan_post[,sy.effs$year_id[i]+161]*stan_post[,11] + stan_post[,sy.effs$siteyear_id[i]+680]*stan_post[,14])
}

sy.effs$sigma.mode <- exp(sy.effs$logsigma.mode)
sy.effs$sigma.lci <- exp(sy.effs$logsigma.lci)
sy.effs$sigma.uci <- exp(sy.effs$logsigma.uci)
sy.effs$max.mode <- exp(sy.effs$logmax.mode)
sy.effs$max.lci <- exp(sy.effs$logmax.lci)
sy.effs$max.uci <- exp(sy.effs$logmax.uci)

#mcmc CIs  !!  NOP needs to be mu, logsigma and logmax
for (i in 1:256){
  sy.effs$mcmc_int.lci[i] <- HPDinterval(mcmc(mcmc_mod4[,2]+mcmc_mod4[,sy.effs$site_id[i]+4]+mcmc_mod4[,sy.effs$year_id[i]+136]+mcmc_mod4[,sy.effs$siteyear_id[i]+157]))[1]
  sy.effs$mcmc_int.uci[i] <- HPDinterval(mcmc(mcmc_mod4[,2]+mcmc_mod4[,sy.effs$site_id[i]+4]+mcmc_mod4[,sy.effs$year_id[i]+136]+mcmc_mod4[,sy.effs$siteyear_id[i]+157]))[2]
  sy.effs$mcmc_int.var[i] <- var(mcmc_mod4[,2]+mcmc_mod4[,sy.effs$site_id[i]+4]+mcmc_mod4[,sy.effs$year_id[i]+136]+mcmc_mod4[,sy.effs$siteyear_id[i]+157])
  
  sy.effs$mcmc_d.lci[i] <- HPDinterval(mcmc(mcmc_mod4[,3]+mcmc_mod4[,sy.effs$site_id[i]+48]+mcmc_mod4[,sy.effs$year_id[i]+143]+mcmc_mod4[,sy.effs$siteyear_id[i]+413]))[1]
  sy.effs$mcmc_d.uci[i] <- HPDinterval(mcmc(mcmc_mod4[,3]+mcmc_mod4[,sy.effs$site_id[i]+48]+mcmc_mod4[,sy.effs$year_id[i]+143]+mcmc_mod4[,sy.effs$siteyear_id[i]+413]))[2]
  sy.effs$mcmc_d.var[i] <- var(mcmc_mod4[,3]+mcmc_mod4[,sy.effs$site_id[i]+48]+mcmc_mod4[,sy.effs$year_id[i]+143]+mcmc_mod4[,sy.effs$siteyear_id[i]+413])
  
  sy.effs$mcmc_d2.lci[i] <- HPDinterval(mcmc(mcmc_mod4[,4]+mcmc_mod4[,sy.effs$site_id[i]+92]+mcmc_mod4[,sy.effs$year_id[i]+150]+mcmc_mod4[,sy.effs$siteyear_id[i]+669]))[1]
  sy.effs$mcmc_d2.uci[i] <- HPDinterval(mcmc(mcmc_mod4[,4]+mcmc_mod4[,sy.effs$site_id[i]+92]+mcmc_mod4[,sy.effs$year_id[i]+150]+mcmc_mod4[,sy.effs$siteyear_id[i]+669]))[2]
  sy.effs$mcmc_d2.var[i] <- var(mcmc_mod4[,4]+mcmc_mod4[,sy.effs$site_id[i]+92]+mcmc_mod4[,sy.effs$year_id[i]+150]+mcmc_mod4[,sy.effs$siteyear_id[i]+669])
}

forestplot <- function(d, x, y, l, u, xlab, ylab){  #function for forest plot: dataframe then column numbers
  xaxis <- d[,x]
  yaxis <- d[,y]
  lci <- d[,l]
  uci <- d[,u]
  
  ggplot(d, aes(fct_rev(fct_inorder(xaxis)), yaxis))+
    geom_point(size=2, alpha=0.5)+
    geom_errorbar(aes(ymax=uci, ymin=lci, width=0.5))+
    theme_bw()+
    xlab(xlab)+
    ylab(ylab)+
    #geom_hline(yintercept=1, linetype="dashed", color = "black")+
    theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust=1))+
    guides(color = "none")+
    coord_flip()
}
# d=dataframe
# x=xaxis (column number in d)
# y=yaxis (column number in d)
# l=lower CI (column number in d) same as y if no CIs
# u=upper CI (column number in d) same as y if no CIs
# xlab/ylab axis title

plot.mu <- forestplot(d=sy.effs, x=6, y=13, l=14, u=15, xlab="Site-year", ylab="Mu estimate")
plot.mu <- plot.mu+geom_hline(yintercept=2.01, linetype="dashed", color = "red")
plot.mu

plot.logsigma <- forestplot(d=sy.effs, x=6, y=17, l=18, u=19, xlab="Site-year", ylab="Logsigma estimate")
plot.logsigma

plot.logmax <- forestplot(d=sy.effs, x=6, y=21, l=22, u=23, xlab="Site-year", ylab="Logmax estimate")
plot.logmax

plot.sigma <- forestplot(d=sy.effs, x=6, y=25, l=26, u=27, xlab="Site-year", ylab="Sigma estimate")
plot.sigma

plot.max <- forestplot(d=sy.effs, x=6, y=28, l=29, u=30, xlab="Site-year", ylab="Max estimate")
plot.max


hist(sy.effs$mu.mode,100)
hist(sy.effs$mu.var,100)
hist(sy.effs$logsigma.mode,100)
hist(sy.effs$logsigma.var,100)
hist(sy.effs$logmax.mode,100)
hist(sy.effs$logmax.var,100)
plot(sy.effs$logsigma.mode, sy.effs$mu.mode)
cor.test(sy.effs$logsigma.mode, sy.effs$mu.mode)
plot(sy.effs$logsigma.mode, sy.effs$logmax.mode)
cor.test(sy.effs$logsigma.mode, sy.effs$logmax.mode)
plot(sy.effs$mu.mode, sy.effs$logmax.mode)
cor.test(sy.effs$mu.mode, sy.effs$logmax.mode)

plot(sy.effs$mu.mode, sy.effs$mu.var)
cor.test(sy.effs$mu.mode, sy.effs$mu.var)
plot(sy.effs$logsigma.mode, sy.effs$logsigma.var)
cor.test(sy.effs$logsigma.mode, sy.effs$logsigma.var)
plot(sy.effs$logmax.mode, sy.effs$logmax.var)
cor.test(sy.effs$logmax.mode, sy.effs$logmax.var)

plot(sy.effs$logsigma.var, sy.effs$mu.var)
cor.test(sy.effs$logsigma.var, sy.effs$mu.var)
plot(sy.effs$logsigma.var, sy.effs$logmax.var)
cor.test(sy.effs$logsigma.var, sy.effs$logmax.var)
plot(sy.effs$mu.var, sy.effs$logmax.var)
cor.test(sy.effs$mu.var, sy.effs$logmax.var)
