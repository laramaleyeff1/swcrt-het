#-----------------------------------------------------#
#              Model params                           #
#-----------------------------------------------------#

#-----------------------------------------------------
# Loading/defining all important functions
#-----------------------------------------------------
start_time = Sys.time()

library(foreach)
library(doParallel)
library(lme4)
library(dplyr)
library(MASS)
library(e1071)

source("../helpers/helper.R")

params <- commandArgs(trailingOnly=TRUE)

# Number of bootstrap resamples, B=0 implies no bootstrap
B <- as.numeric(params[1])

# Number of individuals per cluster per time
n_per <- as.numeric(params[2])

# Scenario that we are looking at - remains fixed throughout study
scenario_fixed <- as.numeric(params[3])

# Parameters fixed for all simulations
t_max = as.numeric(params[4])

# Number of clusters that cross over at each time period
each <- as.numeric(params[5])

# Index for bash jobs
index <- as.numeric(params[6])

# In all simulations, we assume the first crossover occurs at the 2nd time
# period and an equal number of clusters crossover at each subsequent time 
# period. For more customizable code, see ../../code_for_trial_planning
k_max <- (t_max-1)*each
expt_max <- t_max - 1
overall_eff = 0.173314
sd_expt = 0.324037
sigma_alpha_sq = 0.01716

logit <- function(x) {
  return(log(x/(1-x)))
}

rFleish <- function(n,b,c,d) {
  X <- rnorm(n)
  return(b*X + c*X^2 + d*X^3)
}

scale2 <- function(x) {
  return(sd_expt*(x-mean(x))/sd(x) + overall_eff)
}

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

#########################################
#   Set seed for one-time generation    #
#   of random effects (normal and non). #
#   Then, set seed to unique simulation #
#   iteration (to ensure variability)   #
#########################################
set.seed(123, kind = "L'Ecuyer-CMRG")
seeds <- runif(125, 1,10000)
set.seed(123, kind = "L'Ecuyer-CMRG")
expt_ranef_normal <- rnorm(t_max-1)
expt_ranef_fleish <-  rFleish(t_max-1,0.22,0.25,0.1)
Sigma1 = ar1_cor(t_max-1,0.2)
expt_ranef_ar1 <- mvrnorm(1,rep(0,t_max-1),
                                 Sigma1)

# Reset seed (different seed for each job array) so that
# simulations are reproducible
set.seed(seeds[index],kind = "L'Ecuyer-CMRG")

one_scen <- function(scenario) {
  alphas <- rnorm(k_max, 0, sqrt(sigma_alpha_sq))
  data = expand.grid(k = 1:k_max, t = 1:t_max, id = 1:n_per) 
  
  crossover_df = data.frame(k = 1:k_max)
  crossover_df$crossover_time = rep(2:t_max, each = each)
  data = merge(data, crossover_df, by = c("k"), all.x = TRUE) %>%
    arrange(k,t)
  data$expt <- data$t - data$crossover_time + 1
  data[data$expt<0,]$expt = 0
  data$trt = ifelse(data$expt==0,0,1)
  
  data$expt_factor_1 = as.numeric(data$expt == 1)
  data$expt_factor_2 = as.numeric(data$expt > 1)
  
  data$expt_adj = data$expt - 1
  data[data$expt_adj<0,]$expt_adj = 0
  
  t.data <- c(0, 
              -0.013356, 
              0.166609, 
              0.036514,
              -0.066010,
              -0.191512,
              -0.188063,
              0.005889,
              -0.013356, 
              0.166609, 
              0.036514,
              -0.066010,
              -0.191512,
              -0.188063,
              0.005889)
  ###########################################################
  #   Intercept, cluster intercept, and background time     #
  #           (Scenarios 1-24)                              #
  ###########################################################
  
  data$g <- 0.773631 +
    apply(data,
          1,
          function(x) {
            # Random effect for k and fixed effect
            # for time (oscillating)
            alphas[as.numeric(x["k"])]+
            t.data[as.numeric(x["t"])]
          })
  

  
  if (scenario == 1) {
    true <- rep(overall_eff,t_max-1)
  }
  
  # Increasing treatment effect
  if (scenario == 2) {
    true <- 1:(t_max-1)
  }
  
  # Decreasing treatment effect
  if (scenario == 3) {
    true <- -1*1:(t_max-1) 
  }
  
  # Delayed treatment effect
  if (scenario == 4) {
    true <- c(log(1),
              rep( ( (t_max-1)/(t_max-2) )*log(1.2), t_max-2))
  }
  
  # Initial treatment effect
  if (scenario == 5) {
    coef_initial = (t_max-1)*log(1.2)-(t_max-2)*log(1.1)
    true <- c(coef_initial,
              rep( log(1.1), t_max-2))
  }
  
  
  e = 1:(t_max-1)
  # Increasing, then slight decrease
  if (scenario == 6) {
    true <- sin(0.8*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2)))
  }
  
  # Decreasing, then slight increase
  if (scenario == 7) {
    true <- -1*sin(0.8*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2)))+2*log(1.2)
  }
  
  
  # Oscillating - increasing first
  if (scenario == 8) {
    true <- 0.5*sin(2*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2))) + log(1.2)
    
  }
  
  # Oscillating - decreasing first
  if (scenario == 9) {
    true <- -0.5*sin(2*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2))) + log(1.2)
    
  }
  
  # Normal random effects
  if (scenario == 10) {
    true <- expt_ranef_normal
  }
  
  # Non-normal random effects
  if (scenario == 11) {
    true <- expt_ranef_fleish
  }
  
  
  # Ar1 correlated random effects
  if (scenario == 12) {
    true <- expt_ranef_ar1
  }
  
  if (scenario > 1) {
    true = scale2(true)
  }
  
  expt_eff = c(0,true)
  data$g <- data$g + 
    apply(data,
          1,
          function(x) {
            expt_eff[as.numeric(x["expt"])+1]
          })
  
  data$prob <- exp(data$g)/(1+exp(data$g))
  data$Y <- rbinom(nrow(data),1,data$prob)
  
  if (B > 0) {
    registerDoParallel(cores=10)
    boots <- foreach(iter = 1:B, .combine=rbind) %dopar% one_boot(data, bin = TRUE)
    boot_ses = apply(boots, 2, sd)
  }
  
  if (B == 0) {
    boot_ses = rep(-1, 6)
  }
  

  ##########################
  #     Fit models 1-5     #
  ##########################
  

  true_all = c(mean(true), true)
  returned = cbind( 
              data.frame(
                 t_max = t_max,
                 n_per = n_per,
                 k_max = k_max, 
                 expt_max = expt_max,
                 sigma_alpha_sq = sigma_alpha_sq,
                 sd_expt = sd_expt,
                 B = B,
                 scenario = scenario,
                 param = true_all
                 ),
                 rbind(
                   fit_model_1(data, t_max, boot_ses[1], boot_ses[2], TRUE),
                   fit_model_4(data, t_max, boot_ses[3], boot_ses[4], TRUE),
                   fit_model_5(data, t_max, boot_ses[5], boot_ses[6], TRUE)
                   )
              )
               
  return(returned)
}


registerDoParallel(cores=10)
out <- foreach(iter = 1:10, .combine=rbind) %dopar% one_scen(scenario_fixed)


end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")

out.name <- paste0(paste("bin_param", 
                         t_max,
                         expt_max,
                         k_max,
                         n_per,
                         B,
                         scenario_fixed,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

