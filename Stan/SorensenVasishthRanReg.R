#### Sorensen & Vasishth: log Gaussian random intercepts and slopes ####

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

data <- read.csv("~/Documents/Stan/gibsonwu2012data.txt", sep="")

rDat <- subset(data, region=="headnoun")

write("data {
  int<lower=1> N;
  real rt[N]; 
  real<lower=-1,upper=1> so[N]; 
  int<lower=1> J;
  int<lower=1> K;
  int<lower=1, upper=J> subj[N]; 
  int<lower=1, upper=K> item[N];
}
parameters {
  vector[2] beta; 
  real<lower=0> sigma_e; 
  vector<lower=0>[2] sigma_u; 
  vector<lower=0>[2] sigma_w; 
  cholesky_factor_corr[2] L_u; 
  cholesky_factor_corr[2] L_w; 
  matrix[2,J] z_u; 
  matrix[2,K] z_w;
}
transformed parameters{
  matrix[2,J] u;
  matrix[2,K] w;
  
  u <- diag_pre_multiply(sigma_u,L_u) * z_u; //subj random effects
  w <- diag_pre_multiply(sigma_w,L_w) * z_w; //item random effects 
}
model {
  real mu;
  //priors
  L_u ~ lkj_corr_cholesky(2.0);  
  L_w ~ lkj_corr_cholesky(2.0); 
  to_vector(z_u) ~ normal(0,1); 
  to_vector(z_w) ~ normal(0,1); 

  for (i in 1:N){
    mu <- beta[1] + u[1,subj[i]] + w[1,item[i]] + (beta[2] + u[2,subj[i]] + w[2,item[i]])*so[i];
    rt[i] ~ lognormal(mu,sigma_e);
  }
}"
,

"example_ranreg.stan"
)

stanc("example_ranreg.stan")

example_ranreg <- stan_model("example_ranreg.stan")

rDat$so <- seq(-0.97,0.97,0.003546618)
stanDat<-list(subj=as.integer(factor(rDat$subj)),
              item=as.integer(factor(rDat$item)), 
              rt=rDat$rt,
              so=rDat$so,
              N=nrow(rDat), 
              J=length(unique(rDat$subj)), 
              K=length(unique(rDat$item)))


ranIntSlpFit <- sampling(object=example_ranreg, data = stanDat, iter=2000, chains = 4)
