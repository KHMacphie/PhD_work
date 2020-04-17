###############
#### RStan ####
###############
## https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started

#If you ever need to change anything with your C++ toolchain configuration, you can execute
#M <- file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
#file.edit(M)

seaice <- read.csv("Downloads/CC-Stan-intro-master/seaice.csv")
plot(extent_north ~ year, pch = 20, data = seaice)

lm1 <- lm(extent_north ~ year, data = seaice)
summary(lm1)
abline(lm1, col = 2, lty = 2, lw = 3)

# equation for linear model: y = α + β∗x + error

#renaming variables and index years from 1 to 39
x <- I(seaice$year - 1978)
y <- seaice$extent_north
N <- length(seaice$year)

#rerun lm
lm1 <- lm(y ~ x)
summary(lm1)

# extract summary stats
lm_alpha <- summary(lm1)$coeff[1]  # the intercept
lm_beta <- summary(lm1)$coeff[2]  # the slope
lm_sigma <- sigma(lm1)  # the residual error

# Data passed to Stan needs to be a list of named objects. 
# The names given here need to match the variable names used in the models
stan_data <- list(N = N, x = x, y = y)

library(rstan)
library(gdata)
library(bayesplot)

#Comments are indicated by // in Stan
#The write("model code", "file_name") bit allows us to write the Stan model 
#in our R script and output the file to the working directory

write("// Stan model for simple linear regression

      data {
      int < lower = 1 > N; // Sample size
      vector[N] x; // Predictor
      vector[N] y; // Outcome
      }
      
      parameters {
      real alpha; // Intercept
      real beta; // Slope (regression coefficients)
      real < lower = 0 > sigma; // Error SD
      }
      
      model {
      y ~ normal(alpha + x * beta , sigma);
      }
      
      generated quantities {
      } // The posterior predictive distribution",

      "stan_model1.stan")

stanc("stan_model1.stan")
stan_model1 <- "stan_model1.stan"

fit <- stan(file = stan_model1, data = stan_data, warmup = 500, iter = 1000, chains = 4, cores = 2, thin = 1)
fit
#assess model convergence by looking at the Rhat values 
#at or near 1, the chains have converged.

#extract() puts the posterior estimates for each parameter into a list.
posterior <- extract(fit)
str(posterior)

# compare to previous lm
plot(y ~ x, pch = 20)

abline(lm1, col = 2, lty = 2, lw = 3)
abline( mean(posterior$alpha), mean(posterior$beta), col = 6, lw = 2)

# One way to visualize the variability in our estimation of the 
# regression line is to plot multiple estimates from the posterior.

plot(y ~ x, pch = 20)

for (i in 1:500) {
  abline(posterior$alpha[i], posterior$beta[i], col = "gray", lty = 1)
}

abline(mean(posterior$alpha), mean(posterior$beta), col = 6, lw = 2)

# convergence diagnostics
#n_eff is a crude measure of the effective sample size. 
#You usually only need to worry is this number is less than 
#1/100th or 1/1000th of your number of iterations.
#'Anything over an `n_eff` of 100 is usually "fine"' - Bob Carpenter

plot(posterior$alpha, type = "l")
plot(posterior$beta, type = "l")
plot(posterior$sigma, type = "l")

#Parameter summaries
par(mfrow = c(1,3))

plot(density(posterior$alpha), main = "Alpha")
abline(v = lm_alpha, col = 4, lty = 2)

plot(density(posterior$beta), main = "Beta")
abline(v = lm_beta, col = 4, lty = 2)

plot(density(posterior$sigma), main = "Sigma")
abline(v = lm_sigma, col = 4, lty = 2)

# directly calculate the probability of any parameter being over or under
#a certain value of interest.

#Probablility that beta is >0:
  sum(posterior$beta>0)/length(posterior$beta)

#diagnostic plots
  traceplot(fit)
  
#We can also look at the posterior densities & histograms.
  stan_dens(fit)
  stan_hist(fit)

#generate plots which indicate the mean parameter estimates and any credible intervals 
#we may be interested in. Note that the 95% credible intervals for the beta and 
#sigma parameters are very small, thus you only see the dots  
  plot(fit, show_density = FALSE, ci_level = 0.5, outer_level = 0.95, fill_color = "salmon")

#Generated quantities block  
#Stan can use random number generators to generate predicted values for each data point, at each iteration.
  write("// Stan model for simple linear regression

        data {
        int < lower = 1 > N; // Sample size
        vector[N] x; // Predictor
        vector[N] y; // Outcome
        }
        
        parameters {
        real alpha; // Intercept
        real beta; // Slope (regression coefficients)
        real < lower = 0 > sigma; // Error SD
        }
        
        model {
        y ~ normal(x * beta + alpha, sigma);
        }
        
        generated quantities {
        real y_rep[N];
        
        for (n in 1:N) {
        y_rep[n] = normal_rng(x[n] * beta + alpha, sigma);
        }
        
        }",

        "stan_model2_GQ.stan")
  
stan_model2_GQ <- "stan_model2_GQ.stan"

fit3 <- stan(stan_model2_GQ, data = stan_data, iter = 1000, chains = 4, cores = 2, thin = 1)  

#Extracting the y_rep values from posterior.
y_rep <- as.matrix(fit3, pars = "y_rep")
dim(y_rep)
#Each row is an iteration (single posterior estimate) from the model.

##bayesplot package to make some prettier looking plots. 
#This package is a wrapper for many common ggplot2 plots, and has a lot of built-in 
#functions to work with posterior predictions.

#Comparing density of y with densities of y over 200 posterior draws.
ppc_dens_overlay(y, y_rep[1:200, ])
#We can also use this to compare estimates of summary statistics.
ppc_stat(y = y, yrep = y_rep, stat = "mean")


###############################
#### Joel's RStan tutorial ####
###############################

library(MCMCglmm)
data("BTdata")
library(rstan)

# Make and run a model to test whether chick tarsus length is affected by sex and hatch date

write("data{
  int<lower=0> N; 
  vector[N] tarsus;
  vector[N] hatch;
  vector[N] sexM;
  vector[N] sexU;
}
parameters{
  real a;
  real aM;
  real aU;
  real bH;
  real<lower=0> sigma;
}
model{
  vector[N] mu =  a + aM*sexM + aU*sexU + bH*hatch;
  
  a ~ normal(0,10);
  aM ~ normal(0,10);
  aU ~ normal(0,10);
  bH ~ normal(0,10);
  sigma ~ cauchy(0,10);
  
  tarsus ~ normal(mu, sigma);
}",

"simple_regression.stan")

stanc("simple_regression.stan")
stan_model2 <- "simple_regression.stan"


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_data2 <-list(
  N=nrow(BTdata),
  tarsus=BTdata$tarsus,
  sexM = as.numeric(BTdata$sex=="Male"),
  sexU = as.numeric(BTdata$sex=="UNK"),
  hatch = BTdata$hatchdate
)

stanModel2 <- stan_model("simple_regression.stan") #compiling the model

mod_stan2 <- sampling(object=stanModel2, data=stan_data2, chains=1)

summary(mod_stan2)
