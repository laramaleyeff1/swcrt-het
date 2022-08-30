#-----------------------------------------------------#
#              Estimation of overall treatment        #
#   effect for continuous outcome simulated SW-CRTs   #
#-----------------------------------------------------#

#-----------------------------------------------------
# Loading/defining all important functions
#-----------------------------------------------------
start_time = Sys.time()

library(foreach)
library(doParallel)
library(lme4)
library(dplyr)

source("../helpers/helper.R")

params <- commandArgs(trailingOnly=TRUE)

# Number of bootstrap resamples, B=0 implies no bootstrap
B <- as.numeric(params[1])

# Number of individuals per cluster per time
n_per <- as.numeric(params[2])

# Parameters fixed for all simulations
t_max = as.numeric(params[3])

# Number of clusters that cross over at each 
# time period
each <- as.numeric(params[4])

# standard deviation of exposure time het
# i.e. sigma_delta - this is varied in Table S2
sd_expt <- as.numeric(params[5])


# Array index (for internal use)
index <- as.numeric(params[6])

# Number of clusters
k_max <- (t_max-1)*each
# Maximum number of exposure times observed
expt_max <- t_max - 1
# sd of individual level heterogeneity 
sd_epsilon = 1
# overall treatment effect, or average
# treatment effect
overall_eff = 2
# cluster-level heterogeneity
sigma_alpha_sq = 0.02

# takes x and returns logit(x)
logit <- function(x) {
  return(log(x/(1-x)))
}

# scales vector to have mean "overall_eff"
# and sd "sd_expt"
scale2 <- function(x) {
  return(sd_expt*(x-mean(x))/sd(x) + overall_eff)
}

# Generate n draws from a Fleishman distribution
rFleish <- function(n,b,c,d) {
  X <- rnorm(n)
  return(b*X + c*X^2 + d*X^3)
}

# Generate an n x n ar-1 correlation matrix
# with correlation rho
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
  # cluster-specific intercepts
  alphas <- rnorm(k_max, 0, sqrt(sigma_alpha_sq))
  
  # Create simulated SW-CRT with k_max clusters, t_max time periods
  # and n_per individuals per cluster per time period
  data = expand.grid(k = 1:k_max, t = 1:t_max, id = 1:n_per) 
  
  # Create dataset containing the crossover time of 
  # each cluster, and merge with"data"
  crossover_df = data.frame(k = 1:k_max)
  crossover_df$crossover_time = rep(2:t_max, each = each)
  data = merge(data, crossover_df, by = c("k"), all.x = TRUE) %>%
    arrange(k,t)
  
  # Exposure time for treated cluster-time periods
  # equals time period minus crossover time 
  data$expt <- data$t - data$crossover_time 
  # Exposure time for untreated cluster-time periods
  # is always 0
  data[data$expt<0,]$expt = 0
  # Treatment indicator
  data$trt = ifelse(data$expt==0,0,1)
  
  # Create helper variables for fitting Models 2 and 3
  data$expt_factor_1 = as.numeric(data$expt == 1)
  data$expt_factor_2 = as.numeric(data$expt > 1)
  data$expt_adj = data$expt - 1
  data[data$expt_adj<0,]$expt_adj = 0
  
  # Background calendar time trend
  t_fixef_osc <- 0.5*sin(pi*2*(1:t_max-1)/(t_max-1))

  ###########################################################
  #   Intercept, cluster intercept, and background time     #
  #           (Scenarios 1-12)                              #
  ###########################################################
  
  data$Y <- 14 +
    apply(data,
          1,
          function(x) {
            # Random effect for k and fixed effect
            # for time (oscillating)
            alphas[as.numeric(x["k"])]+
            t_fixef_osc[as.numeric(x["t"])]
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
  
  # If there exists heterogeneity, scale it to be
  # sd_expt
  if (scenario > 1) {
    true = scale2(true)
  }
  
  expt_eff = c(0,true)
  
  # Generate continuous outcomes
  data$Y <- data$Y + 
    apply(data,
          1,
          function(x) {
            # Find corresponding entry in expt_eff
            # since expt is lower bounded by 0, add 1
            expt_eff[as.numeric(x["expt"])+1]
          }) + 
    rnorm(nrow(data),0,sd_epsilon)
  
  
  if (B > 0) {
    boots <- foreach(iter = 1:B, .combine=rbind) %do% one_boot(data, bin = FALSE)
    boot_ses = apply(boots, 2, sd)
  }
  
  true_all = c(mean(true), true)
  returned = cbind( 
              data.frame(
                t_max = t_max,
                n_per = n_per,
                k_max = k_max, 
                expt_max = expt_max,
                sigma_alpha_sq = sigma_alpha_sq,
                sd_expt = sd_expt,
                sd_epsilon = sd_epsilon,
                B = B,
                scenario = scenario,
                param = true_all
                ),
                 rbind(
                   fit_model_1(data, t_max, boot_ses[1], boot_ses[2]),
                   fit_model_2(data, t_max, boot_ses[3], boot_ses[4]),
                   fit_model_3(data, t_max, boot_ses[5], boot_ses[6]),
                   fit_model_4(data, t_max, boot_ses[7], boot_ses[8]),
                   fit_model_5(data, t_max, boot_ses[9], boot_ses[10])
                   )
              )
  
               
  return(returned)
}



one_run <- function() {
  ret <- foreach(scenario = c(1,2,4,10), .combine=rbind) %do% one_scen(scenario)
  # for varying sd_expt, uncomment line below
  # ret <- foreach(scenario = c(10), .combine=rbind) %do% one_scen(scenario)
  return(ret)
}

registerDoParallel(cores=10)
out <- foreach(iter = 1:10, .combine=rbind) %dopar% one_run()

end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")
out.name <- paste0(paste("cont_param", 
                         t_max,
                         expt_max,
                         k_max,
                         n_per,
                         B,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

