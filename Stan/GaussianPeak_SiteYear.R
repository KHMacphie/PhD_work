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
cater$siteyear <- paste(cater$site, cater$year)
cater$treeID <- paste(cater$site, cater$tree)
cater$siteday <- paste(cater$site, cater$date ,cater$year)
cater$datescaled <- scale(cater$date)[,1] #mean 146.77, sd 14.04305
cater$resid <- 1:length(cater$site)

temp <- read.csv("~/Dropbox/master_data/site/temperatures.csv")
temp$siteyear <- paste(temp$site, temp$year)
temp.sy <- temp %>%
  group_by(siteyear) %>%
  summarize_all(.funs = c(mean="mean"), na.rm=T) %>% 
  as.data.frame()
temp.sy[,c("site_mean","year_mean")] <- NULL
aprmay.sy <- data.frame(siteyear=temp.sy$siteyear)  
aprmay.sy$temp <- apply(temp.sy[, 1081:2544],1, mean) #day 91 start to 151 end
aprmay.sy$temp.cent <- aprmay.sy$temp-mean(aprmay.sy$temp) #mean= 9.271079
store <- pmatch(cater$siteyear,aprmay.sy$siteyear, duplicates.ok = T)
cater$temp.cent <- aprmay.sy$temp.cent[store]
cater <- cater[-which(is.na(cater$temp.cent)==T),]
rm(temp, temp.sy)

#basic gaussian site-year peak model
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
      real logsigma; //removed lower bound as log scale
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


stan_modGSY3 <- sampling(object=stanModelGSY, data=stan_data3, iter=700, chains=3, 
                         pars=c("mu", "logsigma","logmax","sd_site_mu","sd_site_logsigma","sd_site_logmax",
                                "sd_year_mu","sd_year_logsigma","sd_year_logmax",
                                "sd_siteyear_mu","sd_siteyear_logsigma","sd_siteyear_logmax",
                                "sd_resid","site_mu_scaled","site_logsigma_scaled","site_logmax_scaled",
                                "year_mu_scaled","year_logsigma_scaled","year_logmax_scaled",
                                "siteyear_mu_scaled","siteyear_logsigma_scaled","siteyear_logmax_scaled")) 

write.csv(summary(stan_modGSY3)$summary, "~/Documents/Stan/GSY3_sum.csv")
GSY3_post <- extract(stan_modGSY3, permuted=FALSE)
#could limit which columns here...(each chain is different column)
write.csv(GSY3_post, "~/Documents/Stan/GSY3_post.csv")

#gaussian site-year peak model with site, year and site-year covariance from one correlation matrix
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
      real logsigma;
      real logmax;
      cholesky_factor_corr[3] L;
      matrix[3, N_site] site_scaled; 
      matrix[3, N_year] year_scaled; 
      matrix[3, N_siteyear] siteyear_scaled;
      vector[N_siteday] siteday_scaled; 
      vector[N_treeID] treeID_scaled; 
      vector[N_rec] rec_scaled; 
      vector[N_resid] resid_scaled;
      vector<lower=0>[3] sd_site;
      vector<lower=0>[3] sd_year; 
      vector<lower=0>[3] sd_siteyear;
      real<lower=0> sd_siteday;
      real<lower=0> sd_treeID;
      real<lower=0> sd_rec;
      real<lower=0> sd_resid;
      }
      
      transformed parameters{
      matrix[N_site,3] site_effs;
      matrix[N_year,3] year_effs;
      matrix[N_siteyear,3] siteyear_effs;  

      site_effs = (diag_pre_multiply(sd_site, L) * site_scaled)'; 
      year_effs = (diag_pre_multiply(sd_year, L) * year_scaled)';
      siteyear_effs = (diag_pre_multiply(sd_siteyear, L) * siteyear_scaled)';
      }
      
      model{
      vector[N] y_1;
      vector[N_siteday] siteday_effects;
      vector[N_treeID] treeID_effects;
      vector[N_rec] rec_effects;
      vector[N_resid] resid_effects;
      siteday_effects = sd_siteday * siteday_scaled; 
      treeID_effects = sd_treeID * treeID_scaled; 
      rec_effects = sd_rec * rec_scaled; 
      resid_effects = sd_resid * resid_scaled;

      y_1 = logmax - ((date-(mu+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1])) .* 
      (date-(mu+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1]))) ./ 
      (2.*(exp(logsigma+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2])).*
      (exp(logsigma+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2]))) 
      + site_effs[site_id,3] + year_effs[year_id,3] 
      + siteyear_effs[siteyear_id,3] + siteday_effects[siteday_id] 
      + treeID_effects[treeID_id] + rec_effects[rec_id] + resid_effects[resid_id];
      
      
      y ~ poisson_log(y_1);
      
      to_vector(site_scaled) ~ normal(0,1); 
      to_vector(year_scaled) ~ normal(0,1); 
      to_vector(siteyear_scaled) ~ normal(0,1); 
      siteday_scaled ~ normal(0,1); 
      treeID_scaled ~ normal(0,1); 
      rec_scaled ~ normal(0,1); 
      resid_scaled ~ normal(0,1);
      mu ~ normal(0,10);
      logsigma ~ normal(0,10);
      logmax ~ normal(0,10);
      L ~ lkj_corr_cholesky(2.0);  
      to_vector(sd_site) ~ cauchy(0,10);
      to_vector(sd_year) ~ cauchy(0,10);
      to_vector(sd_siteyear) ~ cauchy(0,10);
      sd_siteday ~ cauchy(0,10);
      sd_treeID ~ cauchy(0,10);
      sd_rec ~ cauchy(0,10);
      sd_resid ~ cauchy(0,10);
      }
      
      generated quantities {
      matrix[3, 3] omega;
      omega = L * L'; 
      }
      ",
      
      "siteyear_covar2.stan"
)
stanc("siteyear_covar2.stan")
GSY_covar2 <- stan_model("siteyear_covar2.stan")

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

GSY_covar4 <- sampling(object=GSY_covar2, data=stan_data3, iter=2500, chains=3, thin=2, 
                         pars=c("mu", "logsigma","logmax","sd_site","sd_year","sd_siteyear",
                                "sd_resid","site_effs","year_effs","siteyear_effs",
                                "site_scaled","year_scaled","siteyear_scaled",
                                "L", "omega")) 

write.csv(summary(GSY_covar4), "~/Documents/Stan/covar4_sum.csv")
covar4_post <- extract(GSY_covar4, permuted=FALSE)
#could limit which columns here...(each chain is different column)
write.csv(covar4_post, "~/Documents/Stan/covar4_post.csv")



#basic gaussian site-year peak model with temp slope from one temp variable
write("data{
      int<lower=0> N;
      vector[N] date;
      vector[N] temp;           //
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
      real logsigma;
      real logmax;
      real t_mu;            //
      real t_logsigma;      //
      real t_logmax;        //
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
      
      y_1 = logmax + t_logmax*temp - ((date-(mu+t_mu*temp+site_mu_effects[site_id]+year_mu_effects[year_id]+siteyear_mu_effects[siteyear_id])) .* 
      (date-(mu+t_mu*temp+site_mu_effects[site_id]+year_mu_effects[year_id]+siteyear_mu_effects[siteyear_id]))) ./ 
      (2.*(exp(logsigma+t_logsigma*temp+site_logsigma_effects[site_id]+year_logsigma_effects[year_id]+siteyear_logsigma_effects[siteyear_id])).*
      (exp(logsigma+t_logsigma*temp+site_logsigma_effects[site_id]+year_logsigma_effects[year_id]+siteyear_logsigma_effects[siteyear_id]))) 
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
      logsigma ~ normal(0,10);
      logmax ~ normal(0,10);
      t_mu ~ normal(0,10);           //
      t_logsigma ~ normal(0,10);     //
      t_logmax ~ normal(0,10);       //
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
      
      "gaussiantemp3.stan"
)
stanc("gaussiantemp3.stan")
stantemp3 <- stan_model("gaussiantemp3.stan")

stan_data4 <-list(
  N=nrow(cater),
  date=cater$datescaled,
  temp=cater$temp.cent,
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

stan_temp5 <- sampling(object=stantemp3, data=stan_data4, iter=2000, chains=3, thin = 2, 
                          pars=c("mu", "logsigma","logmax","t_mu", "t_logsigma","t_logmax","sd_site_mu","sd_site_logsigma","sd_site_logmax",
                                 "sd_year_mu","sd_year_logsigma","sd_year_logmax",
                                 "sd_siteyear_mu","sd_siteyear_logsigma","sd_siteyear_logmax",
                                 "sd_resid","site_mu_scaled","site_logsigma_scaled","site_logmax_scaled",
                                 "year_mu_scaled","year_logsigma_scaled","year_logmax_scaled",
                                 "siteyear_mu_scaled","siteyear_logsigma_scaled","siteyear_logmax_scaled")) 

write.csv(summary(stan_temp5), "~/Documents/Stan/stan_temp5_sum.csv")
stan_temp5_post <- extract(stan_temp5, permuted=FALSE)
#could limit which columns here...(each chain is different column)
write.csv(stan_temp5_post, "~/Documents/Stan/stan_temp5_post.csv")

date <- seq(-3,3,0.1)
mu <- 0.481424772
logsigma <- 0.070243102
logmax <- -2.518215272
t_logsigma <- -0.05
t_sigma <- exp(t_logsigma)
temp <- seq(-15,15,0.5)

test <- data.frame(date=rep(date, length(temp)), temp=rep(temp,1,each=length(date)))

test$y_1 <-  logmax  - ((test$date-(mu))^2) / (2*(exp(logsigma+t_logsigma*test$temp))^2)
test$y_2 <-  logmax  - ((test$date-(mu))^2) / (2*(exp(logsigma)+t_sigma*test$temp)^2)#+ t_logmax*temp  +t_mu*temp
y_1 <-  logmax  - ((date-(mu))^2) / (2*(exp(logsigma))^2)


ggplot(test, aes(date, exp(y_1), col=as.factor(temp)))+
  geom_line(aes(col=as.factor(temp)))+
  theme_bw()

mu <- 0.5
logsigma <- -0.2
logmax <- -2.5
t_mu <- -0.4
t_logsigma <- -0.1
t_logmax <- 0.15
date <- seq(-2.5,2.5,0.1)
temp <- 0
temp1 <- 1
tempneg1 <- -1  

y_0 <- logmax + t_logmax*temp - ((date-(mu+t_mu*temp+0+0+0)) * (date-(mu+t_mu*temp+0+0+0))) / 
        (2*(exp(logsigma+t_logsigma*temp+0+0+0))*(exp(logsigma+t_logsigma*temp+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_1 <- logmax + t_logmax*temp1 - ((date-(mu+t_mu*temp1+0+0+0)) * (date-(mu+t_mu*temp1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*temp1+0+0+0))*(exp(logsigma+t_logsigma*temp1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_neg1 <- logmax + t_logmax*tempneg1 - ((date-(mu+t_mu*tempneg1+0+0+0)) * (date-(mu+t_mu*tempneg1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempneg1+0+0+0))*(exp(logsigma+t_logsigma*tempneg1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0
