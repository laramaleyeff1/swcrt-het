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

## Sample Parameter values:

# nsims is number of simulations to be run
nsims = 200
# R is the number of permutations in permuation test
R = 2000
intercept = 0.3
average_eff = 1
sigma_cluster_sq = 0.1
sigma_expt_sq = 0.2
n_per = 30
t_max = 8
k_max = 7
# crossover_t must have length k_max
crossover_t = 2:t_max
expt_max = t_max - min(crossover_t) + 1
betas = c(0,rnorm(t_max-1))
expt_eff <- c(0,
              rnorm(expt_max, 0, sqrt(sigma_expt_sq)))



#
# Simulates a SW-CRT based on Model 5, and fits Model 5
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
# @returns {ret} A data.frame with the p-values from the Model 4
# LR test, permutation test, and Model 5 LR test, respectively.
sim <- function(intercept, 
                betas, 
                average_eff,
                sigma_cluster_sq,
                sigma_expt_sq,
                n_per,
                t_max,
                k_max,
                crossover_t
) {

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
    average_eff*data$trt + 
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
  model_1 <- glmer(Y ~ trt +
                           as.factor(t) + 
                           (1 | k),
                         data, family=binomial)
  q_obs = anova(model_5,model_1)$Chisq[2]
  
  onePermute <- function() {
    data = data %>%
      filter(expt > 0) %>%
      group_by(k) %>%
      mutate(expt_perm = sample(expt)) %>%
      select(id, expt_perm) %>%
      right_join(data, by=c("id","k")) %>%
      mutate(expt_perm = replace_na(expt_perm, 0))
    
    model_5_perm = glmer(Y ~ trt + 
                              as.factor(t) +
                              (0 + trt | expt_perm)+
                              (1 | k),
                            family="binomial",
                            data=data)
    
    return(c(anova(model_5_perm,model_1)$Chisq[2]))
  }
  
  perm_dist = unlist(lapply(1:R, function(i) {
    onePermute()
  }))
  
  ret = data.frame(
    mod_4_lr_pval = anova(model_4,model_1)[2,8],
    permutation_pval = length(perm_dist[perm_dist>=q_obs])/R,
    mod_5_lr_pval = 1-pchibarsq(q_obs,2)
  )
  return(ret)
}

out <- foreach(iter = 1:nsims, .combine=rbind) %do% sim(intercept, 
                                                    betas, 
                                                    average_eff,
                                                    sigma_cluster_sq,
                                                    sigma_expt_sq,
                                                    n_per,
                                                    t_max,
                                                    k_max,
                                                    crossover_t)
# To compute power
colMeans(out < 0.05)
