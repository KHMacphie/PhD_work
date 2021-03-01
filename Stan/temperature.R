rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
library(lme4)
library(MCMCglmm)
library(dplyr)
library(ggplot2)
library(colorspace)
library(gridExtra)
library(forcats)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

######################################################
#### Gaussian function peak model with temp slope ####
######################################################

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

####basic gaussian site-year peak model with temp slope from one temp variable ####
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
#stanc("gaussiantemp3.stan")
#stantemp3 <- stan_model("gaussiantemp3.stan")

#### gaussian site-year peak model with temp slope from one temp variable and covariance (from 1 correlation matrix) ####
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
      vector[N_siteday] siteday_effects = sd_siteday * siteday_scaled; 
      vector[N_treeID] treeID_effects = sd_treeID * treeID_scaled; 
      vector[N_rec] rec_effects = sd_rec * rec_scaled; 
      vector[N_resid] resid_effects = sd_resid * resid_scaled;
      
      y_1 = logmax + t_logmax*temp - ((date-(mu+t_mu*temp+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1])) .* 
      (date-(mu+t_mu*temp+site_effs[site_id,1]+year_effs[year_id,1]+siteyear_effs[siteyear_id,1]))) ./ 
      (2.*(exp(logsigma+t_logsigma*temp+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2])).*
      (exp(logsigma+t_logsigma*temp+site_effs[site_id,2]+year_effs[year_id,2]+siteyear_effs[siteyear_id,2]))) 
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
      t_mu ~ normal(0,10);           //
      t_logsigma ~ normal(0,10);     //
      t_logmax ~ normal(0,10);       //
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
      
      "tempcovar.stan"
)
stanc("tempcovar.stan")
tempcovar <- stan_model("tempcovar.stan")

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

tempcovar1 <- sampling(object=tempcovar, data=stan_data4, iter=1500, chains=3, 
                       pars=c("mu", "logsigma","logmax","t_mu", "t_logsigma","t_logmax",
                              "sd_site","sd_year","sd_siteyear",
                              "sd_resid","site_effs","year_effs","siteyear_effs",
                              "site_scaled","year_scaled","siteyear_scaled",
                              "L", "omega")) 

write.csv(summary(tempcovar1), "~/Documents/Stan/tempcovar1_sum.csv")
tempcovar1_post <- extract(tempcovar1, permuted=FALSE)
write.csv(tempcovar1_post, "~/Documents/Stan/tempcovar1_post.csv")

temp_sum <- read.csv("stan_temp5_sum.csv")
temp_post <- read.csv("stan_temp5_post.csv")

#### vague look ####
mu <- 0.5
logsigma <- -0.2
logmax <- -2.5
t_mu <- -0.4
t_logsigma <- -0.1
t_logmax <- 0.15
date <- seq(-2.5,2.5,0.1)
temp <- 0
temp1 <- 1
tempmax <- 1.9
tempneg1 <- -1  
tempmin <- -2.7  

y_0 <- logmax + t_logmax*temp - ((date-(mu+t_mu*temp+0+0+0)) * (date-(mu+t_mu*temp+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*temp+0+0+0))*(exp(logsigma+t_logsigma*temp+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_1 <- logmax + t_logmax*temp1 - ((date-(mu+t_mu*temp1+0+0+0)) * (date-(mu+t_mu*temp1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*temp1+0+0+0))*(exp(logsigma+t_logsigma*temp1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_max <- logmax + t_logmax*tempmax - ((date-(mu+t_mu*tempmax+0+0+0)) * (date-(mu+t_mu*tempmax+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempmax+0+0+0))*(exp(logsigma+t_logsigma*tempmax+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_neg1 <- logmax + t_logmax*tempneg1 - ((date-(mu+t_mu*tempneg1+0+0+0)) * (date-(mu+t_mu*tempneg1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempneg1+0+0+0))*(exp(logsigma+t_logsigma*tempneg1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_min <- logmax + t_logmax*tempmin - ((date-(mu+t_mu*tempmin+0+0+0)) * (date-(mu+t_mu*tempmin+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempmin+0+0+0))*(exp(logsigma+t_logsigma*tempmin+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0


plot(date, exp(y_max), type="l", col="red")
points(date, exp(y_1), type="l", col="orange")
points(date, exp(y_0), type="l")
points(date, exp(y_neg1), type="l", col="turquoise")
points(date, exp(y_min), type="l", col="darkblue")



t_mu <- 0 #too visualise impact on shape


y_0 <- logmax + t_logmax*temp - ((date-(mu+t_mu*temp+0+0+0)) * (date-(mu+t_mu*temp+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*temp+0+0+0))*(exp(logsigma+t_logsigma*temp+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_1 <- logmax + t_logmax*temp1 - ((date-(mu+t_mu*temp1+0+0+0)) * (date-(mu+t_mu*temp1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*temp1+0+0+0))*(exp(logsigma+t_logsigma*temp1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_max <- logmax + t_logmax*tempmax - ((date-(mu+t_mu*tempmax+0+0+0)) * (date-(mu+t_mu*tempmax+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempmax+0+0+0))*(exp(logsigma+t_logsigma*tempmax+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_neg1 <- logmax + t_logmax*tempneg1 - ((date-(mu+t_mu*tempneg1+0+0+0)) * (date-(mu+t_mu*tempneg1+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempneg1+0+0+0))*(exp(logsigma+t_logsigma*tempneg1+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

y_min <- logmax + t_logmax*tempmin - ((date-(mu+t_mu*tempmin+0+0+0)) * (date-(mu+t_mu*tempmin+0+0+0))) / 
  (2*(exp(logsigma+t_logsigma*tempmin+0+0+0))*(exp(logsigma+t_logsigma*tempmin+0+0+0))) + 0 + 0 + 0 + 0 + 0 + 0 + 0

plot(date, exp(y_max), type="l", col="red")
points(date, exp(y_1), type="l", col="orange")
points(date, exp(y_0), type="l")
points(date, exp(y_neg1), type="l", col="turquoise")
points(date, exp(y_min), type="l", col="darkblue")

#### Check chains and piece together ####
ggplot(temp_post, aes(chain.1.t_mu)) + 
  geom_histogram(aes(chain.1.t_mu), binwidth=0.01, fill = "red", alpha = 0.2) + 
  geom_histogram(aes(chain.2.t_mu), binwidth=0.01, fill = "blue", alpha = 0.2) +
  geom_histogram(aes(chain.3.t_mu), binwidth=0.01, fill = "yellow", alpha = 0.2) +
  theme_bw()

beta_post <- data.frame(mu=c(temp_post[,2],temp_post[,3],temp_post[,4]), logsigma=c(temp_post[,5],temp_post[,6],temp_post[,7]),logmax=c(temp_post[,8],temp_post[,9],temp_post[,10]),
                   t_mu=c(temp_post[,11],temp_post[,12],temp_post[,13]), t_logsigma=c(temp_post[,14],temp_post[,15],temp_post[,16]),t_logmax=c(temp_post[,17],temp_post[,18],temp_post[,19]))

slope_post <- data.frame(mu=c(temp_post[,11],temp_post[,12],temp_post[,13]), logsigma=c(temp_post[,14],temp_post[,15],temp_post[,16]),logmax=c(temp_post[,17],temp_post[,18],temp_post[,19]))

hist(beta_post$logmax,100) #look good
hist(slope_post$logsigma,200)
abline(v=HPDinterval(mcmc(slope_post$logsigma))[1], col=2)
abline(v=HPDinterval(mcmc(slope_post$logsigma))[2], col=2)

hist(temp_post[,14],200)
abline(v=HPDinterval(mcmc(temp_post[,14]))[1], col=2)
abline(v=HPDinterval(mcmc(temp_post[,14]))[2], col=2)
hist(temp_post[,15],200)
abline(v=HPDinterval(mcmc(temp_post[,15]))[1], col=2)
abline(v=HPDinterval(mcmc(temp_post[,15]))[2], col=2)
hist(temp_post[,16],200)
abline(v=HPDinterval(mcmc(temp_post[,16]))[1], col=2)
abline(v=HPDinterval(mcmc(temp_post[,16]))[2], col=2)

beta_sum <- data.frame(beta=colnames(beta_post))
for(i in 1:length(beta_sum$beta)){
  beta_sum$mean[i] <- mean(beta_post[,i])
  beta_sum$lowci[i] <- HPDinterval(mcmc(beta_post[,i]))[1]
  beta_sum$upci[i] <- HPDinterval(mcmc(beta_post[,i]))[2]
}

slope_sum <- data.frame(beta=colnames(slope_post))
for(i in 1:length(slope_sum$beta)){
  slope_sum$mean[i] <- mean(slope_post[,i])
  slope_sum$lowci[i] <- HPDinterval(mcmc(slope_post[,i]))[1]
  slope_sum$upci[i] <- HPDinterval(mcmc(slope_post[,i]))[2]
}

ggplot(beta_sum, aes(beta, mean))+
    geom_point(size=1)+
    geom_errorbar(aes(ymax=(lowci), ymin=(upci), width=0.1))+
    theme_bw()+
    xlab("Parameter")+
    ylab("Coefficient (link scale)")+
    coord_flip()+
    geom_hline(yintercept=0, linetype="dashed", color = "black")+
    theme(text = element_text(size=15), axis.text.x = element_text(angle = 45, hjust=1))+
    guides(color = "none")

(coeffs <- ggplot(slope_sum, aes(beta, mean))+
  geom_point(size=1.5)+
  geom_errorbar(aes(ymax=(lowci), ymin=(upci), width=0.1))+
  theme_bw()+
  xlab("")+
  ylab("Slope coefficient")+
  coord_flip()+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(text = element_text(size=15), axis.title.y=element_blank())+
  guides(color = "none")) #4x3.5"

#temp cent: -2.77642 - 1.86131  #mean= 9.271079    6.494659 - 11.13239
#date scaled: -2.1199128 - 2.0102438, 117-175    #mean 146.77, sd 14.04305
peaks <- data.frame(date=rep(seq(117,175,0.5),length(seq(6.5,11,0.5))), temp=rep(seq(6.5,11,0.5),1,each=length(seq(117,175,0.5))))
peaks$datescal <- (peaks$date-146.77)/14.04305
peaks$tempcent <- peaks$temp-9.271079 

for(i in 1:length(peaks$date)){
  peaks$mean[i] <- mean(beta_post$logmax + beta_post$t_logmax*peaks$tempcent[i] - ((peaks$datescal[i]-(beta_post$mu+beta_post$t_mu*peaks$tempcent[i]))^2) / 
                                    (2*(exp(beta_post$logsigma+beta_post$t_logsigma*peaks$tempcent[i]))^2))
  peaks$lowci[i] <- HPDinterval(mcmc(beta_post$logmax + beta_post$t_logmax*peaks$tempcent[i] - ((peaks$datescal[i]-(beta_post$mu+beta_post$t_mu*peaks$tempcent[i]))^2) / 
                                    (2*(exp(beta_post$logsigma+beta_post$t_logsigma*peaks$tempcent[i]))^2)))[1]
  peaks$upci[i] <- HPDinterval(mcmc(beta_post$logmax + beta_post$t_logmax*peaks$tempcent[i] - ((peaks$datescal[i]-(beta_post$mu+beta_post$t_mu*peaks$tempcent[i]))^2) / 
                                    (2*(exp(beta_post$logsigma+beta_post$t_logsigma*peaks$tempcent[i]))^2)))[2]
  peaks$mean2[i] <- mean(beta_post$logmax + beta_post$t_logmax*peaks$tempcent[i] - ((peaks$datescal[i]-beta_post$mu)^2) / 
                          (2*(exp(beta_post$logsigma+beta_post$t_logsigma*peaks$tempcent[i]))^2))
}

peaks$Temperature <- as.factor(peaks$temp)
PinkBlue <- c("#27408B", "#3A5FCD", "#436EEE", "#4876FF", "#9F79EE", "#AB82FF", "#E066FF", "#EE7AE9", "#FF83FA", "#F29EEE")
(peakplot <- ggplot(peaks, aes(date, exp(mean), group=temp, col=Temperature))+
  geom_line()+
  theme_bw()+
  xlab("Ordinal Date")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  #coord_cartesian(ylim = c(0,0.1))+
  guides(color = "none")+
  scale_colour_manual(values=PinkBlue)) #7x7"
peaks$date2 <- peaks$date-(mean(beta_post$mu)*14.04305 + 146.77)
(centpeakplot <- ggplot(peaks, aes(date2, exp(mean2), group=temp, col=Temperature))+
  geom_line()+
  theme_bw()+
  xlab("Days from mu")+
  ylab("Abundance")+
  theme(text = element_text(size=15))+
  #coord_cartesian(ylim = c(0,0.1))+
  guides(color = "none")+
  scale_colour_manual(values=PinkBlue))
peaks.hilow <- subset(peaks, temp==11 | temp==6.5)


#Plot timing, height and sd by temp

slopes <- data.frame(temp=seq(6.5,11,0.5))
slopes$tempcent <- slopes$temp-9.271079 
for(i in 1:length(slopes$temp)){
  slopes$mean.mu[i] <- mean(beta_post$mu+beta_post$t_mu*slopes$tempcent[i])*14.04305 + 146.77
  slopes$lowci.mu[i] <- HPDinterval(mcmc(beta_post$mu+beta_post$t_mu*slopes$tempcent[i]))[1]*14.04305 + 146.77
  slopes$upci.mu[i] <- HPDinterval(mcmc(beta_post$mu+beta_post$t_mu*slopes$tempcent[i]))[2]*14.04305 + 146.77
  
  slopes$mean.logsigma[i] <- mean(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i])
  slopes$lowci.logsigma[i] <- HPDinterval(mcmc(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i]))[1]
  slopes$upci.logsigma[i] <- HPDinterval(mcmc(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i]))[2]
  
  slopes$mean.sigma[i] <- exp(mean(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i]))*14.04305
  slopes$lowci.sigma[i] <- exp(HPDinterval(mcmc(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i]))[1])*14.04305
  slopes$upci.sigma[i] <- exp(HPDinterval(mcmc(beta_post$logsigma+beta_post$t_logsigma*slopes$tempcent[i]))[2])*14.04305
  
  slopes$mean.logmax[i] <- mean(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i])
  slopes$lowci.logmax[i] <- HPDinterval(mcmc(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i]))[1]
  slopes$upci.logmax[i] <- HPDinterval(mcmc(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i]))[2]
  
  slopes$mean.max[i] <- exp(mean(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i]))
  slopes$lowci.max[i] <- exp(HPDinterval(mcmc(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i]))[1])
  slopes$upci.max[i] <- exp(HPDinterval(mcmc(beta_post$logmax+beta_post$t_logmax*slopes$tempcent[i]))[2])
}

(muslope <- ggplot(data=slopes, aes(temp, mean.mu))+
  geom_line(lwd=1, col=1)+
  geom_ribbon(aes(x=temp, ymin=lowci.mu, ymax=upci.mu), alpha=0.25, fill="gray")+
  xlab("")+
  ylab("Mu (date)")+
  theme_bw()+
  theme(text=element_text(size= 15),axis.text.x=element_blank()))

ggplot(data=slopes, aes(temp, mean.logsigma))+
  geom_line(lwd=1, col=1)+
  geom_ribbon(aes(x=temp, ymin=lowci.logsigma, ymax=upci.logsigma), alpha=0.25, fill="gray")+
  xlab("Temperature (°C)")+
  ylab("log σ")+
  theme_bw()+
  theme(text=element_text(size= 15))

(sigmaslope <- ggplot(data=slopes, aes(temp, mean.sigma))+
  geom_line(lwd=1, col=1)+
  geom_ribbon(aes(x=temp, ymin=lowci.sigma, ymax=upci.sigma), alpha=0.25, fill="gray")+
  xlab("")+
  ylab("Sigma (days)")+
  theme_bw()+
  theme(text=element_text(size= 15),axis.text.x=element_blank()))

#ggplot(data=slopes, aes(temp, mean.logmax))+
#  geom_line(lwd=1, col=1)+
#  geom_ribbon(aes(x=temp, ymin=lowci.logmax, ymax=upci.logmax), alpha=0.25, fill="gray")+
#  xlab("Temperature (°C)")+
#  ylab("log max")+
#  theme_bw()+
#  theme(text=element_text(size= 15))

(maxslope <- ggplot(data=slopes, aes(temp, mean.max))+
  geom_line(lwd=1, col=1)+
  geom_ribbon(aes(x=temp, ymin=lowci.max, ymax=upci.max), alpha=0.25, fill="gray")+
  xlab("Temperature (°C)")+
  ylab("Max. abund. ")+
  theme_bw()+
  theme(text=element_text(size= 15)))

col1 <- grid.arrange(muslope,sigmaslope,maxslope, nrow = 3, heights = c(1,1,1))
col2 <- grid.arrange(peakplot, nrow = 1)
col3 <- grid.arrange(centpeakplot, nrow = 1)
space <- ggplot() + theme_void()
SlopePeak <- grid.arrange(col1, space, col2, ncol = 3, widths = c(1,0.1,2))#7x9
SlopePeakCent <- grid.arrange(col1, space, col3, ncol = 3, widths = c(1,0.1,2))#7x9


d <- peaks$date
dscal <- peaks$datescal
scal <- exp(-2.560951 - ((dscal - 0.8278605)^2)/(2*0.828601^2))
notscal <- exp(-2.560951 - ((d - 158.3957)^2)/(2*11.63609^2))

par(mfrow=c(1,2))
plot(dscal, scal)
plot(d,notscal)

#1 sd is 34.1% of the area under the curve
area <- data.frame(temp=slopes$temp, tempcent=slopes$tempcent)
for(i in 1:length(area$tempcent)){
area$logarea[i] <- mean(beta_post$logsigma+beta_post$logmax+(beta_post$t_logsigma+beta_post$t_logmax)*area$tempcent[i]+log(sqrt(2*pi))) 
area$logarea.lci[i] <- HPDinterval(mcmc(beta_post$logsigma+beta_post$logmax+(beta_post$t_logsigma+beta_post$t_logmax)*area$tempcent[i]+log(sqrt(2*pi))))[1] 
area$logarea.uci[i] <- HPDinterval(mcmc(beta_post$logsigma+beta_post$logmax+(beta_post$t_logsigma+beta_post$t_logmax)*area$tempcent[i]+log(sqrt(2*pi))))[2] 
}

ggplot(data=area, aes(temp, exp(logarea)))+
    geom_line(lwd=1, col=1)+
    geom_ribbon(aes(x=temp, ymin=exp(logarea.lci), ymax=exp(logarea.uci)), alpha=0.25, fill="gray")+
    xlab("Temperature (°C)")+
    ylab("Area under peak")+
    theme_bw()+
    theme(text=element_text(size= 15))
