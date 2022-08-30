#-----------------------------------------------------#
#             Code to simulate                        #
#    SW-CRTs with a binary outcome, and assess        #
# estimation of exposure-time-specific effects        #
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

B <- as.numeric(params[1])

# Number of time points observed
t_max <- as.numeric(params[2])

# Number of individuals per cluster per time
n_per <- as.numeric(params[3])

# Index for bash jobs
index <- as.numeric(params[4])

# Number of clusters that cross over at each time period
each <- 1
# cluster-level heterogeneity
sigma_alpha_sq = 0.05
# Number of clusters
k_max <- (t_max-1)*each
expt_max <- t_max - 1

# Helper functions
logit <- function(x) {
  return(log(x/(1-x)))
}

rFleish <- function(n,b,c,d) {
  X <- rnorm(n)
  return(b*X + c*X^2 + d*X^3)
}

scale2 <- function(x) {
  return(0.1*(x-mean(x))/sd(x))
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
expt.ranef.normal <- c(0,
               rnorm(t_max-1))
expt.ranef.fleish <-  c(0,
                        rFleish(t_max-1,0.22,0.25,0.1))
Sigma1 = ar1_cor(t_max-1,0.2)
expt.ranef.ar1.rho1 <- c(0,mvrnorm(1,
                                   rep(0,t_max-1),
                                   Sigma1))

# Here we scale the non-zero exposure times to have mean 0 and
# variance 0.1
expt.ranef.normal[-1] <- scale2(expt.ranef.normal[-1])
expt.ranef.fleish[-1] <- scale2(expt.ranef.fleish[-1])
expt.ranef.ar1.rho1[-1] <- scale2(expt.ranef.ar1.rho1[-1])
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
  
  # oscillating
  t.fixef.osc <- 0.5*sin(pi*2*(1:t_max-1)/(t_max-1))
  # increasing
  t.fixef.inc <- 1.5*(sigmoid(1:t_max-1)-0.75)
  
  ###########################################################
  #   Intercept, cluster intercept, and background time     #
  #           (Scenarios 1-24)                              #
  ###########################################################
  
  data$g <- logit(0.7) +
    apply(data,
          1,
          function(x) {
            # Random effect for group and fixed effect
            # for time (oscillating)
            alphas[as.numeric(x["k"])]+
            (scenario <= 12)*t.fixef.osc[as.numeric(x["t"])]+
            (scenario > 12)*t.fixef.inc[as.numeric(x["t"])]
          })
  

  
  # Constant treatment effect of log(1.2)
  if (scenario == 1 | scenario == 13) {
    data$g <- data$g + log(1.2)*data$trt
    true <- rep(log(1.2),t_max-1)
  }
  
  # Increasing treatment effect
  if (scenario == 2 | scenario == 14) {
    slope = log(1.2)/mean(1:(t_max-1))
    data$g <- data$g + slope*data$expt
    true <- slope*1:(t_max-1)
  }
  
  # Decreasing treatment effect
  if (scenario == 3 | scenario == 15) {
    slope = -log(1.2)/mean(1:(t_max-1))
    data$g <- data$g + slope*data$expt 
    true <- slope*1:(t_max-1) 
  }
  
  # Delayed treatment effect
  if (scenario == 4 | scenario == 16) {
    data$g <- data$g + 
      log(1)*data$expt.factor.1 +
      ((t_max-1)/(t_max-2))*log(1.2)*data$expt.factor.2
    true <- c(log(1),
              rep( ( (t_max-1)/(t_max-2) )*log(1.2), t_max-2))
  }
  
  # Initial treatment effect
  if (scenario == 5 | scenario == 17) {
    coef_initial = (t_max-1)*log(1.2)-(t_max-2)*log(1.1)
    data$g <- data$g + 
      coef_initial*data$expt.factor.1 +
      log(1.1)*data$expt.factor.2
    true <- c(coef_initial,
              rep( log(1.1), t_max-2))
  }
  
  
  e = 1:(t_max-1)
  # Increasing, then slight decrease
  if (scenario == 6 | scenario == 18) {
    expt.fixef <- c(0, 
                    sin(0.8*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2))))
    data$g <- data$g + 
      apply(data,
            1,
            function(x) {
              # Find corresponding entry in expt.fixef
              # since expt is lower bounded by 0, add 1
              expt.fixef[as.numeric(x["expt"])+1]
            })
    true <- expt.fixef[-1]
    
  }
  
  # Decreasing, then slight increase
  if (scenario == 7 | scenario == 19) {
    expt.fixef <- c(0, 
                    -1*sin(0.8*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2)))+2*log(1.2))
    data$g <- data$g + 
      apply(data,
            1,
            function(x) {
              # Find corresponding entry in expt.fixef
              # since expt is lower bounded by 0, add 1
              expt.fixef[as.numeric(x["expt"])+1]
            })
    true <- expt.fixef[-1]
    
  }
  
  
  # Oscillating - increasing first
  if (scenario == 8 | scenario == 20) {
    expt.fixef <- c(0,
                    0.5*sin(2*pi*(e-1)/(t_max-2))*log(1.2)/
                      mean(sin(0.8*pi*(e-1)/(t_max-2))) + log(1.2)
    )
    data$g <- data$g + 
      apply(data,
            1,
            function(x) {
              # Find corresponding entry in expt.fixef
              # since expt is lower bounded by 0, add 1
              expt.fixef[as.numeric(x["expt"])+1]
            })
    true <- expt.fixef[-1]
  }
  
  # Oscillating - decreasing first
  if (scenario == 9 | scenario == 21) {
    expt.fixef <- c(0,
                  -0.5*sin(2*pi*(e-1)/(t_max-2))*log(1.2)/
                    mean(sin(0.8*pi*(e-1)/(t_max-2))) + log(1.2)
                  )
    data$g <- data$g + 
      apply(data,
            1,
            function(x) {
              # Find corresponding entry in expt.fixef
              # since expt is lower bounded by 0, add 1
              expt.fixef[as.numeric(x["expt"])+1]
            })
    true <- expt.fixef[-1]
  }
  
  # Normal random effects
  if (scenario == 10 | scenario == 22) {
    data$g <- data$g + 
      log(1.2)*data$trt +
      apply(data,
            1,
            function(x) {
              expt.ranef.normal[as.numeric(x["expt"])+1]
            })
    true <- log(1.2) + expt.ranef.normal[-1]
  }
  # Non-normal random effects
  if (scenario == 11 | scenario == 23) {
    data$g <- data$g + 
      log(1.2)*data$trt +
      apply(data,
            1,
            function(x) {
              expt.ranef.fleish[as.numeric(x["expt"])+1]
            })
    true <- log(1.2) + expt.ranef.fleish[-1]
  }
  
  
  # Ar1 correlated random effects
  if (scenario == 12 | scenario == 24) {
    data$g <- data$g + 
      log(1.2)*data$trt +
      apply(data,
            1,
            function(x) {
              expt.ranef.ar1.rho1[as.numeric(x["expt"])+1]
            })
    true <- log(1.2) + expt.ranef.ar1.rho1[-1]
  }
  
  
  data$prob <- exp(data$g)/(1+exp(data$g))
  data$Y <- rbinom(nrow(data),1,data$prob)
  
  
  if (B > 0) {
    boots <- foreach(iter = 1:B, .combine=rbind) %do% one_boot_expt(data)
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
                 scenario = scenario,
                 sigma_alpha_sq = sqrt(sigma_alpha_sq),
                 sd_expt = 0.1,
                 true = true_all
                 ),
                 rbind(
                   fit_model_1(data, t_max, boot_ses[1:t_max], 0, TRUE),
                   fit_model_4(data, t_max, boot_ses[(t_max + 1):(2*(t_max))], 0, TRUE),
                   fit_model_5(data, t_max, boot_ses[(2*t_max+1):(3*t_max)], 0, TRUE)
                   )
              )
               
  return(returned)
}


one_run <- function() {
  ret <- foreach(scenario = c(1:9,
                              11:21,
                              23,24), .combine=rbind) %do% one_scen(scenario)
  return(ret)
}

one_run_boot <- function() {
  ret <- foreach(scenario = c(10,22), .combine=rbind) %do% one_scen(scenario)
  return(ret)
}



if (B == 0) {
  registerDoParallel(cores=10)
  out <- foreach(iter = 1:10, .combine=rbind) %dopar% one_run()
} else {
  registerDoParallel(cores=10)
  out <- foreach(iter = 1:10, .combine=rbind) %dopar% one_run_boot()
}

end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")
out.name <- paste0(paste("param", 
                         t_max,
                         expt_max,
                         k_max,
                         n_per,
                         B,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)

