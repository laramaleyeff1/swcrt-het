######################################
#    Code to simulate a SW-CRT with  #
#     binary outcome, as             #
# described in Section 5.1 "Trial    #
#           Planning"                #
######################################
library(doParallel)
library(lme4)
library(dplyr)
library(tidyr)
library(emdbook)

## Params used in Web Appendix F.2.2
## Results are saved as:
nsims = 1000
intercept = 0.773639
average_eff = log(1.19)
sigma_cluster_sq = 0.01716
sd_expt = 0.324037
n_per = c(59,61,51,17,38,14,36,43,29,32,32,
          10,6,39,54,66,41,22,35,20,43,42,
          40,34,27,27,11,65,51,47,40,21,41,
          11,42,33,34,27,27,12,15,63,69,62,
          40,22,31,13,37,35,27,30,22,21,14,
          95,47,48,42,22,42,14,37,50,44,25,
          19,25,15,96,60,58,34,25,31,12,35,
          29,56,32,29,17,10,77,52,47,49,19,
          35,12,37,38,57,25,29,9,9,53,46,51,
          28,18,37,16,43,24,45,29,22,6,12,71)
t_max = 8
k_max = 14
crossover_t = rep(2:t_max, each = 2)
expt_max = 7
betas = c(-0.0134, 0.1667, 0.0365,
          -0.0660,-0.1915,-0.1881, 0.0059)
expt_eff <- c(-0.375,
              -0.252,  0.412,
              -0.141, -0.119,
              0.471,  0.004)

#
# Simulates a SW-CRT based on Model 5, and fits Models 4 and 5
# @param {intercept} Value of mu (overall intercept)
# @param {betas} Vector of time-specific effect estimates
# @param {average_eff} Value of phi
# @param {sigma_cluster_sq} Value of the random cluster intercept
# (sigma_alpha^2)
# @param {sigma_expt_sq} Value of the random exposure time slope
# (sigma_delta^2)
# @param {n_per} A vector of length k_max*t_max of cluster-period sizes,
# or a constant cluster-period size
# @param {t_max} Number of time periods
# @param {k_max} [default] Number of clusters 
# @param {crossover_t} Vector with time of cross-over
# for each cluster 
# 
# @returns {ret} A data.frame with the average treatment 
# effect and standard error estimate from Model 4 and the
# average treatment effect estimate from Model 5
sim.binary <- function(ate, 
                       t_max, 
                       k_max, 
                       n_per, 
                       crossover_t, 
                       intercept,  
                       betas, 
                       expt_eff, 
                       sigma_alpha_sq
) {
  
  betas <- c(0,betas)
  expt_eff <- c(0,expt_eff + average_eff)
  if (length(crossover_t) != k_max) {
    print("Error: crossover_t must be of length k_max")
    return()
  }
  
  if (length(n_per) == 1) {
    count = rep(n_per, t_max*k_max)
  } else {
    count = n_per
  }
  
  freq_table = expand.grid(k = rep(1:k_max),
                           t = rep(1:t_max)
  )
  freq_table$count = count
  
  data = freq_table[rep(1:nrow(freq_table), freq_table[["count"]]), ]
  N = nrow(data)
  
  data = data %>%
    mutate(id = 1:N,
           cross_t = k, 
           cross_t = as.numeric(as.character(factor(k,
                      labels=crossover_t))))

  data$expt <- data$t - data$cross_t + 1
  
  data[data$expt<0,]$expt = 0
  data$trt = ifelse(data$expt==0,0,1)
  
  cluster_ranef <- rnorm(k_max, 0, sqrt(sigma_cluster_sq))
  
  data$g <- intercept +
    apply(data,
          1,
          function(x) {
            # cluster intercept 
            # Background calendar time
            # exposure time random effect
            cluster_ranef[as.numeric(x["k"])]+
              betas[as.numeric(x["t"])] +
              expt_eff[as.numeric(x["expt"]+1)]
          })
  
  # Prob(Y=1) = expit(g)
  data$prob <- exp(data$g)/(1+exp(data$g))
  # Generate 0/1 outcome Y
  data$Y <- rbinom(nrow(data),1,data$prob)
  
  model_4 <- glmer(Y ~ as.factor(expt) +
                              as.factor(t) + 
                              (1 | k),
                            data, family=binomial)
  model_5 <- glmer(Y ~ trt + as.factor(t) + 
                              (1 | k)+
                              (0 + trt | expt), 
                            data, family=binomial)
  vcov_params_4 = vcov(model_4)[2:(expt_max + 1),2:(expt_max + 1)]
  C = rep(1,(expt_max))/(expt_max)
  if(is.null(model_4@optinfo$conv$lme4$code)) {
    se_overall_model_4 = as.numeric(sqrt(C %*% vcov_params_4 %*% C))
  } else {
    se_overall_model_4 = NA
  }
  
  ret = data.frame(
    mod_4_param = mean(fixef(model_4)[2:(expt_max + 1)]),
    mod_4_se = se_overall_model_4,
    mod_5_param = fixef(model_5)[2]
  )
  return(ret)
}


ATE.sim.power.binary <- function(ate, 
                                 t_max, 
                                 k_max, 
                                 n_per, 
                                 crossover_t, 
                                 intercept,  
                                 betas, 
                                 expt_eff, 
                                 sigma_alpha_sq, 
                                 nsims) {
  
  sims.binary <- foreach(iter = 1:nsims, .combine=rbind) %do% sim.binary(ate, 
                                                                 t_max, 
                                                                 k_max, 
                                                                 n_per, 
                                                                 crossover_t, 
                                                                 intercept,  
                                                                 betas, 
                                                                 expt_eff, 
                                                                 sigma_alpha_sq)
  sigma_sq = c(sd(sims.binary$mod_4_param),
               sd(sims.binary$mod_5_param))
  return((pnorm(sqrt(1/sigma_sq) * abs(ate) - qnorm(0.975)) + 
            pnorm(-1 * sqrt(1/sigma_sq) * abs(ate) - qnorm(0.975)))*100)
}

set.seed(123)
ATE.sim.power.binary(ate = log(3), t_max = 8, k_max = 14, n_per = 34, 
                     crossover_t = rep(2:8,each=2), intercept = 0,  betas = rnorm(7, 0, 1), 
                     expt_eff = c(-0.38, -0.25,  0.41, -0.14, -0.12, 0.47, 0.004),
                     sigma_alpha_sq = 0.34, nsims = 1000)
