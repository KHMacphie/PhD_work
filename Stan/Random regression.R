rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
#library(rstanarm)
#library(brms)
library(lme4)
library(MCMCglmm)
#library(gdata)
#library(bayesplot)

## STAN IN R ##
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#data
cater <- read.csv("cater.csv")
cater$year <- as.factor(cater$year)
cater$resid <- 1:length(cater$site)
cater19 <- subset(cater, year=="2019")
cater19$datescalsqrd <- cater19$datescaled^2


## simple poisson ##

write("data{
      int<lower=0> N;
      vector[N] date;
      vector[N] datesq;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
   
      parameters{
      real a;
      real b_date;
      real b_datesq;
      vector[N_site] site_effects; 
      vector[N_resid] resid_effects;
      real<lower=0> sigma_site;
      real<lower=0> sigma_resid;
      }

      model{
      
      vector[N] mu = a + b_date*date + b_datesq*datesq + site_effects[site_id] + resid_effects[resid_id];
      
      site_effects ~ normal(0,sigma_site);
      resid_effects ~ normal(0,sigma_resid);
      a ~ normal(0,10);
      b_date ~ normal(0,10);
      b_datesq ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);
      
      y ~ poisson_log(mu);
      }",
      
      "poisson_simple.stan"
      )

stanc("poisson_simple.stan")

stanModelP <- stan_model("poisson_simple.stan")

stan_data1 <-list(
  N=nrow(cater19),
  date=cater19$datescaled,
  datesq=cater19$datescalsqrd,
  N_site=length(unique(cater19$site)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  site_id=as.numeric(as.factor(cater19$site)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_modP <- sampling(object=stanModelP, data=stan_data1, iter=4000, chains=2, pars=c("a", "b_date","b_datesq","sigma_site","sigma_resid"))

#freqentist compare
freq_mod3 <- glmer(caterpillars ~ datescaled + datescalsqrd + (1|site) + (1|resid), data = cater19, family = poisson)

summary(stan_modP)
summary(freq_mod3)


## random reg poisson ##

write("data{
      int<lower=0> N;
      vector[N] date;
      vector[N] datesq;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real a;
      real b_date;
      real b_datesq;
      vector[N_resid] resid_effects;
      real<lower=0> sigma_resid;
      vector<lower=0>[2] sigma_site; 	//subj sd (2 values for int and slope)
      cholesky_factor_corr[2] L_site; 	//Cholesky decomp of correlation matrix
      matrix[2,N_site] z_site; 			//random effects int + slope per level
      }
      
      transformed parameters{
      matrix[2,N_site] site_effects;
     
      //compute the transpose matrix product 
      site_effects <- diag_pre_multiply(sigma_site,L_site) * z_site; //random effects 
      }
      
      model{
      vector[N] mu;
            
      L_site ~ lkj_corr_cholesky(2.0);  //in brackets is eta for lkj dist
      to_vector(z_site) ~ normal(0,1); 
     

      resid_effects ~ normal(0,sigma_resid);
      a ~ normal(0,10);
      b_date ~ normal(0,10);
      b_datesq ~ normal(0,10);
      sigma_resid ~ cauchy(0,10);
      
      for(i in 1:N){mu = a + site_effects[1,site_id[i]] + resid_effects[resid_id[i]] + (site_effects[2,site_id[i]]+b_date)*date + b_datesq*datesq;
      y[i] ~ poisson_log(mu);
      }
      }",
      
      "poisson_ranreg.stan"

)

stanc("poisson_ranreg.stan")

stanModelPRR <- stan_model("poisson_ranreg.stan")

stan_data1 <-list(
  N=nrow(cater19),
  date=cater19$datescaled,
  datesq=cater19$datescalsqrd,
  N_site=length(unique(cater19$site)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  site_id=as.numeric(as.factor(cater19$site)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_modPRR <- sampling(object=stanModelPRR, data=stan_data1, iter=6000, chains=1)

#freqentist compare
freq_modPRR <- glmer(caterpillars ~ datescaled + datescalsqrd + (1+datescaled|site) + (1|resid), data = cater19, family = poisson)

summary(stan_modPRR)
summary(freq_modPRR)


