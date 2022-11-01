######################################
#    Code to simulate a SW-CRT       #
# with continuous outcomes, as       #
# described in Section 5.1 "Trial    #
#           Planning"                #
######################################
library(doParallel)
library(lme4)
library(dplyr)
library(tidyr)
library(emdbook)


#
# Simulates a SW-CRT based on Model 5, and fits Model 5
# @param {ate} Value of average treatment effect
# @param {t_max} Number of time periods
# @param {k_max} [default] Number of clusters (defaults to t_max-1)
# @param {n_per} A vector of length k_max*t_max of cluster-period sizes,
# or a constant cluster-period size
# @param {crossover_t} Vector with time of cross-over
# for each cluster 
# @param {intercept} Value of mu (overall intercept)
# @param {betas} Vector of time-specific effect estimates
# @param {expt_eff} Exposure-time specific 
# effects centered at 0 (\theta_1, \dots, \theta_E)
# @param {sigma_alpha_sq} Value of the random cluster intercept
# (sigma_alpha^2)
# @param {sigma_epsilon_sq} Value of the individual level heterogeneity
# 
# @returns {ret} A data.frame with the estimated ATE from Model 4, 
# estimated ATE standard error from Model 4, and the estimated 
# ATE from Model 5
sim.cont <- function(ate, 
                           t_max, 
                           k_max, 
                           n_per, 
                           crossover_t, 
                           intercept,  
                           betas, 
                           expt_eff, 
                           sigma_alpha_sq,
                           sigma_epsilon_sq
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
           cross_t = as.numeric(
                      as.character(
                        factor(cross_t,
                            labels=crossover_t))))
  
  data$expt <- data$t - data$cross_t + 1
  data[data$expt<0,]$expt = 0
  data$trt = ifelse(data$expt==0,0,1)
  
  cluster_ranef <- rnorm(k_max, 0, sqrt(sigma_cluster_sq))
 
  data$Y <- intercept +
    ate*data$trt + 
    apply(data,
          1,
          function(x) {
            # cluster intercept 
            # Background calendar time
            # exposure time random effect
            cluster_ranef[as.numeric(x["k"])]+
            betas[as.numeric(x["t"])] +
            expt_eff[as.numeric(x["expt"]+1)]
          }) +
    rnorm(nrow(data),0,sigma_epsilon_sq)
  

  model_4 <- lmer(Y ~ as.factor(expt) +
                              as.factor(t) + 
                              (1 | k),
                            data)
  model_5 <- lmer(Y ~ trt +
                              as.factor(t) + 
                              (1 | k)+
                              (0 + trt | expt), 
                            data)
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


ATE.sim.power.cont <- function(ate, 
                                 t_max, 
                                 k_max, 
                                 n_per, 
                                 crossover_t, 
                                 intercept,  
                                 betas, 
                                 expt_eff, 
                                 sigma_alpha_sq, 
                                 sigma_epsilon_sq,
                                 nsims) {
  
  sims.cont <- foreach(iter = 1:nsims, .combine=rbind) %do% sim.cont(ate, 
                                                                         t_max, 
                                                                         k_max, 
                                                                         n_per, 
                                                                         crossover_t, 
                                                                         intercept,  
                                                                         betas, 
                                                                         expt_eff, 
                                                                         sigma_alpha_sq,
                                                                         sigma_epsilon_sq)
  sigma_sq = c(sd(sims.cont$mod_4_param),
               sd(sims.cont$mod_5_param))
  return((pnorm(sqrt(1/sigma_sq) * abs(ate) - qnorm(0.975)) + 
            pnorm(-1 * sqrt(1/sigma_sq) * abs(ate) - qnorm(0.975)))*100)
}


set.seed(123)
ATE.sim.power.cont(ate = 1, t_max = 8, k_max = 7, n_per = 30, 
                     crossover_t = 2:8, intercept = 0,  betas = rnorm(7, 0, 1), 
                     expt_eff = c(-0.38, -0.25,  0.41, -0.14, -0.12, 0.47, 0.004),
                     sigma_alpha_sq = 0.34, sigma_epsilon_sq = 1, nsims = 1000)


data_het = data_het %>%
  select(group, tau, Y)
