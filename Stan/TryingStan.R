###############
#### RStan ####
###############

rm(list=ls())
setwd('/Users/s1205615/Documents/Stan/')
library(rstan)
library(gdata)
library(bayesplot)

## STAN IN R ##
#library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#The first one means that rstan detects the number of cores on your computer, 
#and so if you run multiple chains, it runs them in parallel.
#The second means that when we compile a model (see below) it saves that model 
#so we do not have to re-compile it.

#Dimensions are declared in square brackets, e.g. vector[N].

#All lines of code in stan have to end with ;

#Comments can be made with \\

#A stan program has three required “blocks”:
  #data
  #parameters
  #model
#There are also four optional blocks:
  #functions
  #transformed data
  #transformed parameters
  #generated quantities

#think of a stan model as writing a model out in maths. 
#For example, a linear regression can be written as:
  #y=α+βx+ϵ
  #ϵ∼N(0,σϵ)

## DATA BLOCK ##
#declare all the data that you are inputting 
  #the data type, their dimensions, 
  #any restrictions (i.e. upper= or lower= , which allows Stan to sanity check your inputs),
  #and their names

data{
  int<lower=0> N; 
  vector[N] y;
  vector[N] x;
}

## PARAMETERS BLOCK ##
#In this block you indicate the parameter’s dimensions, restrictions and names. 
#For the simple linear regression above, we will want to model the intercept (a), 
#slopes (b) and the standard deviation of the residuals (sigma).

parameters{
  real a;
  real b;
  real<lower=0> sigma;
}

## MODEL BLOCK ##
#This is where you include any sampling statements, including the “likelihood” (model) 
#you are using and any prior distributions you want to include for your parameters. 
#If no prior is defined, Stan uses default priors of uniform(-infinity, +infinity).
#every parameter will have a prior whther you define it or not. 
#You can restrict priors using upper or lower when declaring the parameters 
#(e.g. sigma above has <lower=0> to make sure it is positive). 
#Sampling is indicated by the ∼ symbol, and Stan already includes many common distributions 
#as vectorized functions. The manual contains a comprehensive list.
#link for prioi recommendations: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

#The final sampling statement has to be in terms of data. Our linear regression can also be written as
  #y∼N(α+βx+ϵ,σϵ)

model{
  vector[N] mu = a + b*x;
  
  a ~ normal(0,10);
  b ~ normal(0,10);
  sigma ~ cauchy(0,10); #because I specified <lower=0> for sigma in the parameters block, this is a half-cauchy prior
  
  y ~ normal(mu, sigma)
}

#I tend to write the stan program in a separate .stan file, then call it from inside R. 
#The text editor sublime (and probably others) provides syntax highlighting for stan

write("data{
      int<lower=0> N; 
      vector[N] y;
      vector[N] x;
      }
      parameters{
      real a;
      real b;
      real<lower=0> sigma;
      }
      model{
      vector[N] mu = a + b*x;
  
      a ~ normal(0,10);
      b ~ normal(0,10);
      sigma ~ cauchy(0,10);
  
      y ~ normal(mu, sigma);
      }",

      "simpleregression.stan")

stanc("simpleregression.stan")
#stan_model2 <- "simpleregression.stan" #this just stores the file path

# simulate some data 
a <- 5
b <- sqrt(0.5)
sigma <- sqrt(0.5)
x <- rnorm(100)
y <- rnorm(100,a + b*x, sigma)

# input data as a list (with the names of the items matching those in the stan program)
stan_data <- list(
  N=length(y),
  y = y,
  x = x)

#compile our model
stanModel <- stan_model("simpleregression.stan") 

#run the model
mod_stan <- sampling(object=stanModel, data=stan_data, chains=1)
#Note that in the same folder as the .stan file a .rds file by the same name has appeared. This is the compiled model.
#The defaults to 2000 iterations and a warmup period of 1000.  
#The warmup period is analogous to (but not he same as) the burnin period of Gibbs Sampling (mcmcglmm).
#if warmupnot defined it is half the iterations
# can manipulate: e.g. add warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1

#look at the output
summary(mod_stan)$summary

#n_eff is a measure of effective sample size, which looks good. 
#Rhat is a measure of convergence and is best interpreted for multiple chains 
#a Rhat of 1 is perfect, >1.1 is bad, and signifies that the chains have converged to different places (or haven’t converged). 
#lp__ is the model likelihood.

# Look at the traceplots
rstan::traceplot(mod_stan)

#look at pairs plot to look for posterior correlations between variables
pairs(mod_stan)

#Look at parameter estimates and compare to simulated data and a freq lm
summary(mod_stan)$summary[,c(1,3)]
c(a,b,sigma)
lm_mod <- lm(y~x)
rbind(summary(lm_mod)$coef[,c(1,2)],sigma=c(summary(lm_mod)$sigma,NA))

#look at full posterior
#extract() puts the posterior estimates for each parameter into a list.
posterior <- extract(mod_stan)
str(posterior)
hist(posterior$a,100)
plot(seq(1,1000,1), posterior$a)
plot(density(posterior$a), main = "Intercept")
abline(v = summary(mod_stan)$summary[1,1], col = 4, lty = 2)
#or
stan_dens(mod_stan)
stan_hist(mod_stan)

#plot parameter mean and CIs
plot(mod_stan, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")

#### Blue tit example ####
library(MCMCglmm)
data(BTdata)

#Make and run a model to test whether chick tarsus length is affected by sex and hatch date.
head(BTdata)
levels(BTdata$sex) # male female and unknown

#make the stan program
write("data{
      int<lower=0> N; // data type and sample size
      vector[N] tarsus;
      vector[N] hatch;
      vector[N] sexM;
      vector[N] sexU;
      }
      parameters{
      real a;  // intercept female
      real aM; // intercept dif male
      real aU; // intercept dif unknown
      real bH; // slope for hatch date
      real<lower=0> sigma; // standard deviation of residuals, has to be positive
      }
      model{
      vector[N] mu =  a + aM*sexM + aU*sexU + bH*hatch; //model structions for sex factor and hatch dtae slope
      
      a ~ normal(0,10); //priors
      aM ~ normal(0,10);
      aU ~ normal(0,10);
      bH ~ normal(0,10);
      sigma ~ cauchy(0,10); //half cauchy because parameters defines as >0
      
      tarsus ~ normal(mu, sigma); //normally distributed residuals
      }",

      "simple_regression.stan")

stanc("simple_regression.stan")  #Translates Stan model specification to C++ code
#stan_model2 <- "simple_regression.stan" #this just stores the file path


#Make data
stan_data2 <-list(
  N=nrow(BTdata), 
  tarsus=BTdata$tarsus,
  sexM = as.numeric(BTdata$sex=="Male"),
  sexU = as.numeric(BTdata$sex=="UNK"),
  hatch = BTdata$hatchdate
)

#compile the model
stanModel2 <- stan_model("simple_regression.stan") 

#run the model
mod_stan2 <- sampling(object=stanModel2, data=stan_data2, chains=1)

#look at the output
summary(mod_stan2)
summary(mod_stan2)$summary[,c(1,3,9)]

#compare with freq lm
lm_mod_BT <- lm(tarsus~sex+hatchdate, BTdata)
rbind(summary(lm_mod_BT)$coef[,c(1,2)],sigma=c(summary(lm_mod_BT)$sigma,NA))


#### Simple optimising for linear models ####
#linear model can also be written as:
  #y=Xβ+ϵ  X is the design matrix, and β is a vector or fixed effects
  #ϵ∼N(0,σϵ)

write("data{
      int<lower=0> N;
      int<lower=0> J; // all fixed parameters
      vector[N] y;
      matrix[N,J] X; // matrix of length N and J wide
      }
      parameters{
      vector[J] beta;
      real<lower=0> sigma;
      }
      model{
      vector[N] mu = X*beta; // model as a matrix
  
      beta ~ normal(0,10); // priors same for all fixed terms
      sigma ~ cauchy(0,10);
  
      y ~ normal(mu, sigma);
      }",

      "linear_regression.stan")

stanc("linear_regression.stan")

#Get the design matrix in r using the model.matrix() function
X<- model.matrix(tarsus~sex+hatchdate, BTdata) #assigned to X as in the stan program
head(X)

stan_data3 <-list(
  N=nrow(BTdata),
  J=ncol(X),
  y=BTdata$tarsus,
  X=X
)

stanModel3 <- stan_model("linear_regression.stan")

mod_stan3 <- sampling(object=stanModel3, data=stan_data3, chains=1)

summary(mod_stan3)$summary[,c(1,3,9)] #same as previous


#### QR decomposition ####
#You can do some fancy matrix algebra to uncorrelate your variables, so your beta estimates (beta_tilde), 
#associated with the transformed design matrix (Q), are more easily estimated. 
#Then you can back transform the estimated betas, and recover the actual estimates. 

#The transformation of the design matrix can be done in the transformed data block. 
#The calculations in this block are only done once, before the iterations are done, and not every iteration.
transformed data{
  // Compute, thin, and then scale QR decomposition
  matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
  matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
  matrix[J, J] R_inv = inverse(R);
}

#We then estimate beta_tilde instead of beta in the parameters block
parameters{
  vector[J] beta_tilde;
  real<lower=0> sigma;
}

#then Q and beta_tilde are multiplied instead of X and beta in the model block
model{
  vector[N] mu = Q*beta_tilde;
  
  beta_tilde ~ normal(0,10);
  sigma ~ cauchy(0,10);
  
  y ~ normal(mu, sigma);
}

#finally betas are back calculated
generated quantities {
  vector[J] beta = R_inv * beta_tilde; // coefficients on x
}

#Stan program:
write("data{
      int<lower=0> N;
      int<lower=0> J; 
      vector[N] y;
      matrix[N,J] X;
      }

      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }

      parameters{
      vector[J] beta_tilde;
      real<lower=0> sigma;
      }

      model{
      vector[N] mu = Q*beta_tilde;
  
      beta_tilde ~ normal(0,10);
      sigma ~ cauchy(0,10);
  
      y ~ normal(mu, sigma);
      }

      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "QRdecomp_lm.stan")

stanc("QRdecomp_lm.stan")

stanModel4 <- stan_model("QRdecomp_lm.stan")

mod_stan4 <- sampling(object=stanModel4, data=stan_data3, chains=1) #same data as before

#look at output and compare to previous version (not QR decomposed)
summary(mod_stan3)$summary[,c(1,3,9)] #previous
summary(mod_stan4)$summary[,c(1,3,9)] #new model 


#### Linear mixed model ####

#Mixed model written as:
  #yi=Xβi+zjϵi
  #z∼N(0,σz)
  #ϵ∼N(0,σϵ)

#Indexing is a key part of making a LMM in stan. Lets think of this in terms of chick in nests. 
#We can read this equation as chick i with weight yi and predictors Xi in nest j. 
#z the random nest effects are also parameters that we are estimating (i.e. they are a a latent, unobserved variable).

#so needs more info
data{
  int<lower=0> N;
  int<lower=0> J;
  int<lower=0> N_nest;
  vector[N] y;
  matrix[N,J] X;
  int<lower=0> nest_id[N];
}
#Here N_nest is the number of nests (i.e. number of levels in this factor). 
#nest_id is an index, that should go from 1 to N_nest. 
#Within the model we are going to create a vector of nest effects, 1 effect for each nest, 
#and this index refers to the position in that vector. You can’t give stan a r factor or something like this, 
#you have to create this indexing yourself.

#We need to estimate more parameters, one for the random effect SD, and a vector of random effects.
parameters{
  vector[J] beta_tilde;
  vector[N_nest] nest_effects;
  real<lower=0> sigma_nest;
  real<lower=0> sigma_E;
}

#Then in our model we can state that our random nest effects come from a normal distribution with SD sigma_nest, 
#and we can add these nest effects to our fixed effects, using the indexes we created.

model{
  vector[N] mu = Q*beta_tilde + nest_effects[nest_id];
  
  nest_effects ~ normal(0,sigma_nest);
  
  beta_tilde ~ normal(0,10);
  sigma_nest ~ cauchy(0,10);
  sigma_E ~ cauchy(0,10);
  
  y ~ normal(mu, sigma_E);
}


write("data{
      int<lower=0> N;
      int<lower=0> J;
      int<lower=0> N_nest;
      vector[N] y;
      matrix[N,J] X;
      int<lower=0> nest_id[N];
      }

      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }

      parameters{
      vector[J] beta_tilde;
      vector[N_nest] nest_effects;
      real<lower=0> sigma_nest;
      real<lower=0> sigma_E;
      }

      model{
      vector[N] mu = Q*beta_tilde + nest_effects[nest_id];
  
      nest_effects ~ normal(0,sigma_nest);
  
      beta_tilde ~ normal(0,10);
      sigma_nest ~ cauchy(0,10);
      sigma_E ~ cauchy(0,10);
  
      y ~ normal(mu, sigma_E);
      }

      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "mixedmodel.stan")

stanc("mixedmodel.stan")

#Run this on simulated data

N <- 1000
N_nest <-100
J <- 5 #number of fixed parameters
beta <- rnorm(J)
X <- matrix(rnorm(N*J),N,J) #when rnorm has no mean and sd specified it defaults to 0 and 1, the matrix is of N*J values with N rows and J columns

sigma_nest <- sqrt(0.7)
sigma_E <- sqrt(0.3)

z <- rnorm(N_nest,0,sigma_nest) #random effects
nest_id <- rep(1:N_nest,N/N_nest) #nest id's
mu <- X %*% beta + z[nest_id]
y <- rnorm(N, mu, sigma_E)

#as data for stan
stan_data5 <- list(
  N=N,
  J=J,
  N_nest=N_nest,
  y=y,
  X=X,
  nest_id=nest_id
)

stanModel5 <- stan_model("mixedmodel.stan")

mod_stan5 <- sampling(object=stanModel5, data=stan_data5, chains=1, pars=c("beta","sigma_nest","sigma_E"))
#The pars argument specifies what we want stan to output - otherwise we get all the nest random effects in the output, 
#which we most likely don’t want (we can also use this to get rid of the beta_tilde estimates, which also are not interesting).

summary(mod_stan5)$summary[,c(1,3,9)]

#compare with lme4
library(lme4)

lmm1 <- lmer(y ~ X + (1|nest_id))
summary(lmm1)


#### Run a LMM with the blue tit data, fitting random effects of foster nest and dam ####

head(BTdata)
summary(BTdata)

#Stan program:
write("data{
      int<lower=0> N;
      int<lower=0> J; 
      int<lower=0> N_fnest;
      int<lower=0> N_dam;
      vector[N] y;
      matrix[N,J] X;
      int<lower=0> fnest_id[N];
      int<lower=0> dam_id[N];
      }

      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }
      
  parameters{
      vector[J] beta_tilde;
      vector[N_fnest] fnest_effects;
      vector[N_dam] dam_effects;
      real<lower=0> sigma_fnest;
      real<lower=0> sigma_dam;
      real<lower=0> sigma_E;
      }

      model{
      vector[N] mu = Q*beta_tilde + fnest_effects[fnest_id] + dam_effects[dam_id];
      
      fnest_effects ~ normal(0,sigma_fnest);
      dam_effects ~ normal(0,sigma_dam);

      beta_tilde ~ normal(0,10);
      sigma_fnest ~ cauchy(0,10);
      sigma_dam ~ cauchy(0,10);
      sigma_E ~ cauchy(0,10);
      
      y ~ normal(mu, sigma_E);
      }
      
      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "bluetit_lmm.stan")

stanc("bluetit_lmm.stan")

stanModel6 <- stan_model("bluetit_lmm.stan")

X<- model.matrix(tarsus~sex+hatchdate, BTdata) #assigned to X as in the stan program
head(X)

stan_data6 <-list(
  N=nrow(BTdata),
  J=ncol(X),
  N_fnest=length(unique(BTdata$fosternest)),
  N_dam=length(unique(BTdata$dam)),
  y=BTdata$tarsus,
  X=X,
  fnest_id=as.integer(BTdata$fosternest),
  dam_id=as.integer(BTdata$dam)
)


mod_stan6 <- sampling(object=stanModel6, data=stan_data6, iter=5000, chains=4, pars=c("beta","sigma_fnest","sigma_dam","sigma_E"))


summary(mod_stan6)$summary[,c(1,3,9,10)]
rstan::traceplot(mod_stan6)

#compare with lme4
library(lme4)

lmm2 <- lmer(tarsus ~ sex+hatchdate + (1|fosternest) + (1|dam), data=BTdata)
summary(lmm2)


#### Non-centred parameterisation ####
#can also write a mixed model as:
  #yi=Xβi+σzzjϵi (alternate characters are subscript)
  #z∼N(0,1)
  #ϵ∼N(0,σϵ)
#Stan likes this more, I think essentially because it is easier for it to estimate the random effects as unit normal,
#and then scale them, than to jointly estimate the effects and their variance.

write("data{
      int<lower=0> N;
      int<lower=0> J;
      int<lower=0> N_dam;
      int<lower=0> N_fnest;
      vector[N] y;
      matrix[N,J] X;
      int<lower=0,upper=N_dam> dam_id[N]; //upper bound is new
      int<lower=0,upper=N_fnest> fnest_id[N];
      }
      transformed data{
      // Compute, thin, and then scale QR decomposition
      matrix[N, J] Q = qr_Q(X)[, 1:J] * sqrt(N - 1.0);
      matrix[J, J] R = qr_R(X)[1:J, ] / sqrt(N - 1.0);
      matrix[J, J] R_inv = inverse(R);
      }
      parameters{
      vector[J] beta_tilde;
      vector[N_dam] dam_scaled; //random effects are now scaled as sd=1 lower down so needs conversion
      real<lower=0> sigma_dam;
      vector[N_fnest] fnest_scaled;
      real<lower=0> sigma_fnest;
      real<lower=0> sigma_E;
      }
      model{
      vector[N_dam] dam_effects = sigma_dam * dam_scaled; //unscales the random effects
      vector[N_fnest] fnest_effects = sigma_fnest * fnest_scaled;
      vector[N] mu = Q*beta_tilde + dam_effects[dam_id]+ fnest_effects[fnest_id];
      
      dam_scaled ~ normal(0,1); //random effects have sd of 1 rather than the sigma
      fnest_scaled ~ normal(0,1);
      
      beta_tilde ~ normal(0,10);
      sigma_dam ~ cauchy(0,10);
      sigma_fnest ~ cauchy(0,10);
      sigma_E ~ cauchy(0,10);
      
      y ~ normal(mu, sigma_E);
      }
      generated quantities {
      // recalculate coefficients on x
      vector[J] beta = R_inv * beta_tilde; 
      }",
      
      "ncparameterisation.stan")

stanc("ncparameterisation.stan")

stanModel7 <- stan_model("ncparameterisation.stan")

mod_stan7 <- sampling(object=stanModel7, data=stan_data6, iter=5000, chains=4, pars=c("beta","sigma_dam","sigma_fnest","sigma_E"))

summary(mod_stan7)$summary[,c(1,3,9)]
summary(mod_stan6)$summary[,c(1,3,9)] #same, with non-centred paramaterisation has higher ess
