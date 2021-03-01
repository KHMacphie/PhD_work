rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
library(rstanarm)
library(brms)
library(lme4)
library(MCMCglmm)
#library(gdata)
#library(bayesplot)

## STAN IN R ##
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#load data
MMH <- read.csv("MMHdata.csv")
MMH$year <- as.factor(MMH$year)
MMH$sy <- MMH$SY

#mcmcglmm version
freq_mod1 <- lmer(fledged ~ year + async + I(async^2) + logheight.scal + (1|site) + (1|sy), data=MMH)
summary(freq_mod1)

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))

mcmc_mod1 <- MCMCglmm(fledged ~ year + async + I(async^2) + logheight.scal, random=~site + sy, data=MMH, prior=prior)
summary(mcmc_mod1)

freq_mod2 <- lmer(fledged ~ year + async*year + I(async^2)*year + logheight.scal + (1|site) + (1|sy), data=MMH)
summary(freq_mod2)

## mod1 in RStan ##

write("data{
      int<lower=0> N;
      int<lower=0> J;
      int<lower=0> N_site;
      int<lower=0> N_sy;
      vector[N] y;
      matrix[N,J] X;
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_sy> sy_id[N];
      }
      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }
      parameters{
      vector[J] beta_tilde;
      vector[N_site] site_scaled; 
      vector[N_sy] sy_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_sy;
      real<lower=0> sigma_E;
      }
      model{
      vector[N_site] site_effects = sigma_site * site_scaled; 
      vector[N_sy] sy_effects = sigma_sy * sy_scaled;
      vector[N] mu = Q*beta_tilde + site_effects[site_id]+ sy_effects[sy_id];
      
      site_scaled ~ normal(0,1); 
      sy_scaled ~ normal(0,1);
      
      beta_tilde ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_sy ~ cauchy(0,10);
      sigma_E ~ cauchy(0,10);
      
      y ~ normal(mu, sigma_E);
      }
      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "mmh_mod1.stan")

stanc("mmh_mod1.stan")

stanModel1 <- stan_model("mmh_mod1.stan")

X<- model.matrix(fledged~year+async+I(async^2)+logheight.scal, MMH) #assigned to X as in the stan program
head(X)

stan_data1 <-list(
  N=nrow(MMH),
  J=ncol(X),
  N_site=length(unique(MMH$site)),
  N_sy=length(unique(MMH$sy)),
  y=MMH$fledged,
  X=X,
  site_id=as.numeric(as.factor(MMH$site)),
  sy_id=as.numeric(as.factor(MMH$sy))
)

stan_mod1 <- sampling(object=stanModel1, data=stan_data1, iter=5000, chains=4, pars=c("beta","sigma_site","sigma_sy","sigma_E"))

summary(stan_mod1)$summary[,c(1,3,9,10)]
summary(freq_mod1)
sumary(mcmc_mod1)

rstan::traceplot(stan_mod1, pars=c("beta","sigma_site","sigma_sy","sigma_E"))
stan_dens(stan_mod1, pars=c("beta","sigma_site","sigma_sy","sigma_E"))
stan_hist(stan_mod1, pars=c("beta","sigma_site","sigma_sy","sigma_E"))
pairs(stan_mod1, pars=c("beta","sigma_site","sigma_sy","sigma_E"))

#plot lme4 model and stan model with CIs
library(ggeffects)
freq_pred1 <- ggpredict(freq_mod1, terms=c("async[all]","year"))
ggplot(freq_pred1, aes(x,predicted, col=group))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=x, ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15)+
  theme_bw()

#extract posterior distributions
postdist1 <- extract(stan_mod1)
postdist1 <- data.frame(postdist1$beta)
colnames(postdist1) <- colnames(X) 
colnames(postdist1)[1] <- "year2014"

postdist1$year2015 <- postdist1$year2014 + postdist1$year2015
postdist1$year2016 <- postdist1$year2014 + postdist1$year2016
postdist1$year2017 <- postdist1$year2014 + postdist1$year2017
postdist1$year2018 <- postdist1$year2014 + postdist1$year2018
postdist1$year2019 <- postdist1$year2014 + postdist1$year2019

stan_pred1 <- data.frame(async=rep(seq(min(MMH$async), max(MMH$async),0.5),6), year= c(rep("year2014",93), rep("year2015",93),rep("year2016",93),rep("year2017",93),rep("year2018",93),rep("year2019",93)))

for(i in 1:length(stan_pred1$async)){
  
  A <- postdist1[,stan_pred1$year[i]]+postdist1[,7]*stan_pred1$async[i]+postdist1[,8]*stan_pred1$async[i]^2
  
  stan_pred1$mean[i] <- mean(A)
  stan_pred1$lowci[i] <- posterior_interval(as.matrix(A), prob=0.95)[1]
  stan_pred1$upci[i] <- posterior_interval(as.matrix(A), prob=0.95)[2]
}

ggplot(stan_pred1, aes(async,mean, col=year))+
  geom_line(lwd=1)+
  geom_ribbon(aes(x=async, ymin=lowci, ymax=upci, fill=year), alpha=0.15)+
  theme_bw()

## mod2- fixed interactions ##

X<- model.matrix(fledged~year+async*year+I(async^2)*year+logheight.scal, MMH) #assigned to X as in the stan program
head(X)

stan_data2 <-list(
  N=nrow(MMH),
  J=ncol(X),
  N_site=length(unique(MMH$site)),
  N_sy=length(unique(MMH$sy)),
  y=MMH$fledged,
  X=X,
  site_id=as.numeric(as.factor(MMH$site)),
  sy_id=as.numeric(as.factor(MMH$sy))
)

stan_mod2 <- sampling(object=stanModel1, data=stan_data2, iter=5000, chains=4, pars=c("beta","sigma_site","sigma_sy","sigma_E"))

summary(stan_mod2)$summary[,c(1,3,9,10)]
summary(freq_mod2)

summary(stan_mod2)$summary[1:19,1]
summary(freq_mod2)$coefficients[,1]


#### Poisson model ####
cater <- read.csv("cater.csv")
cater$year <- as.factor(cater$year)
cater$resid <- 1:length(cater$site)
cater19 <- subset(cater, year=="2019")

##Using rstanarm/brms to see code for Poisson glmms 
#in rstanarm
stan_glm1 <- stan_glmer(caterpillars ~ datescaled + I(datescaled^2) + (1|site) + (1|resid),
                      data = cater19, family = poisson,
                      chains = 4, cores = 1) 

stancode1 <- rstan::get_stancode(stan_glm1$stanfit)
cat(stancode)

#in bmrs
stan_glm2 <- brm(bf(caterpillars ~ datescaled + I(datescaled^2) + (1|site) + (1|resid),
                        family = brmsfamily('poisson')), data = cater19,
                     iter = 2000,
                     chains = 2, cores = 1)  ##bf() bypasses the tooles question that comes up somehow..

stancode2 <- stancode(stan_glm2)

#equivalent lm in bmrs
cater19$logcater <- log(cater19$caterpillars+1)
stan_lm2 <- brm(logcater ~ datescaled + I(datescaled^2) + (1|site) + (1|resid),
                    family = gaussian(), data = cater19,
                 iter = 2000,
                 chains = 2, cores = 1)

stancodelm2 <- stancode(stan_lm2)
summary(brms::standata(stan_lm2))

stan_glm2.2 <- brm(caterpillars ~ datescaled + I(datescaled^2) + (1|site) + (1|resid),
                    family = brmsfamily('poisson'), data = cater19,
                 iter = 2000,
                 chains = 2, cores = 1)

stancode2.2 <- stancode(stan_glm2.2)
summary(brms::standata(stan_glm2.2))

#### Poisson Stan model ####

freq_mod3 <- glmer(caterpillars ~ datescaled + I(datescaled^2) + (1|site) + (1|resid), data = cater19, family = poisson)
summary(freq_mod3)

## mod1 in RStan ##

write("data{
      int<lower=0> N;
      int<lower=0> J;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      matrix[N,J] X;
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }
      parameters{
      vector[J] beta_tilde;
      vector[N_site] site_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_resid;
      }
      model{
      vector[N_site] site_effects = sigma_site * site_scaled; 
      vector[N_resid] resid_effects = sigma_resid * resid_scaled;
      vector[N] mu = Q*beta_tilde + site_effects[site_id]+ resid_effects[resid_id];
      
      site_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      
      beta_tilde ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);
      
      y ~ poisson_log(mu);
      }
      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "poisson_mod1.stan"
      )

stanc("poisson_mod1.stan")

stanModel2 <- stan_model("poisson_mod1.stan")

X<- model.matrix(caterpillars~datescaled+I(datescaled^2), cater19) #assigned to X as in the stan program
head(X)

stan_data3 <-list(
  N=nrow(cater19),
  J=ncol(X),
  N_site=length(unique(cater19$site)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  X=X,
  site_id=as.numeric(as.factor(cater19$site)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_mod3 <- sampling(object=stanModel2, data=stan_data3, iter=7000, chains=3, pars=c("beta","sigma_site","sigma_resid"))

k<-10000
prior<-list(R=list(V=1,nu=0.002),
            G=list(G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k),
                   G1=list(V=1,nu=1,aplha.mu=0,alpha.V=k)))
mcmc_mod3 <- MCMCglmm(caterpillars ~ datescaled + I(datescaled^2), random=~site + resid, data = cater19, family = "poisson", prior=prior, nitt=30000,burnin=10000)

summary(stan_mod3)$summary[,c(1,3,9,10)]
summary(mcmc_mod3)
summary(freq_mod3)

rstan::traceplot(stan_mod3, pars=c("beta","sigma_site","sigma_resid"))
stan_dens(stan_mod3, pars=c("beta","sigma_site","sigma_resid"))
stan_hist(stan_mod3, pars=c("beta","sigma_site","sigma_resid"))
pairs(stan_mod3, pars=c("beta","sigma_site","sigma_resid"))



#### Gaussian function in model ####
write("data{
      int<lower=0> N;
      int<lower=0> J;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      matrix[N,J] X;
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }

      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }
      
      parameters{
      vector[J] beta_tilde;
      vector[N_site] site_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_resid;
      }
      
      model{
      vector[N_site] site_effects = sigma_site * site_scaled; 
      vector[N_resid] resid_effects = sigma_resid * resid_scaled;
      vector[N] mu = Q*beta_tilde + site_effects[site_id]+ resid_effects[resid_id];
      
      site_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      
      beta_tilde ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);
      
      y ~ poisson_log(mu);
      }
      
      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "poisson_mod1.stan"
)

## no QR decomposition
write("data{
      int<lower=0> N;
      vector[N] date;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real sigma;
      real logmax;
      vector[N_site] site_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_resid;
      }
      
      model{
      vector[N_site] site_effects = sigma_site * site_scaled; 
      vector[N_resid] resid_effects = sigma_resid * resid_scaled;
      vector[N] y_1; 
      site_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      sigma ~ normal(0,10);
      logmax ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);

      for(i in 1:N){ y_1 = (-((date[i]-mu)^2)/(2*sigma^2)) +logmax + site_effects[site_id[i]] + resid_effects[resid_id[i]];
      y[i] ~ poisson_log(y_1);
      }
      }",
      
      "gaussian_simple.stan"
)

stanc("gaussian_simple.stan")

stanModelgaussian <- stan_model("gaussian_simple.stan")

stan_data1 <-list(
  N=nrow(cater19),
  date=cater19$datescaled,
  N_site=length(unique(cater19$site)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  site_id=as.numeric(as.factor(cater19$site)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_modG <- sampling(object=stanModelgaussian, data=stan_data1, iter=4000, chains=2, pars=c("mu", "sigma","logmax","sigma_site","sigma_resid"))

## no QR decomposition
write("data{
      int<lower=0> N;
      vector[N] date;
      int<lower=0> N_site;
      int<lower=0> N_resid;
      int<lower=0> y[N];
      int<lower=0,upper=N_site> site_id[N]; 
      int<lower=0,upper=N_resid> resid_id[N];
      }
      
      parameters{
      real mu;
      real<lower=0> sigma;
      real logmax;
      vector[N_site] site_scaled; 
      vector[N_resid] resid_scaled;
      real<lower=0> sigma_site;
      real<lower=0> sigma_resid;
      }
      
      model{
      vector[N] y_1;
      vector[N_site] site_effects; 
      vector[N_resid] resid_effects;
      
      site_effects = sigma_site * site_scaled; 
      resid_effects = sigma_resid * resid_scaled;

      y_1 = logmax - ((date-mu) .* (date-mu)) ./ (2*sigma*sigma);
      for(i in 1:N){y_1[i] += site_effects[site_id[i]] + resid_effects[resid_id[i]];
      }

      y ~ poisson_log(y_1);

      site_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      sigma ~ normal(0,10);
      logmax ~ normal(0,10);
      sigma_site ~ cauchy(0,10);
      sigma_resid ~ cauchy(0,10);
      }",
      
      "gaussian_simple.stan"
)

stanc("gaussian_simple.stan")

stanModelgaussian <- stan_model("gaussian_simple.stan")

stan_data1 <-list(
  N=nrow(cater19),
  date=cater19$datescaled,
  N_site=length(unique(cater19$site)),
  N_resid=length(unique(cater19$resid)),
  y=cater19$caterpillars,
  site_id=as.numeric(as.factor(cater19$site)),
  resid_id=as.numeric(as.factor(cater19$resid))
)

stan_modG <- sampling(object=stanModelgaussian, data=stan_data1, iter=2000, chains=4)

summary(stan_modG)$summary[1:55,c(1,3,4,8,9,10)]
rstan::traceplot(stan_modG)

plot(cater19$datescaled,cater19$caterpillars)
datameans <- data.frame(datescaled=cater19$datescaled,cater=cater19$caterpillars)
datameans <- aggregate(.~datescaled, datameans, mean)
mod_betas <- summary(stan_modG)$summary[c(1,2,3),1]
datameans$pred <- mod_betas[3] - ((datameans$datescaled-mod_betas[1])^2)/(2*mod_betas[2]^2)
plot(datameans$datescaled,datameans$cater)
points(datameans$datescaled, exp(datameans$pred), type="l")
