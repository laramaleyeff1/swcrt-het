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

## Sample Parameter values:

# nsims is number of simulations to be run
nsims = 200
# R is the number of permutations in permuation test
R = 2000
intercept = 0.3
overall_eff = 1
sigma_cluster_sq = 0.1
sigma_expt_sq = 0.2
sigma_epsilon_sq = 1
n_per = 30
sd_n_per = 0
t_max = 8
k_max = 7
# crossover_t must have length k_max
crossover_t = 2:t_max
expt_max = t_max - min(crossover_t) + 1
betas = rnorm(t_max)
expt_eff <- c(0,
              rnorm(expt_max, 0, sqrt(sigma_expt_sq)))


#
# Simulates a SW-CRT based on Model 5, and fits Model 5
# @param {intercept} Value of mu (overall intercept)
# @param {betas} Vector of time-specific effect estimates
# @param {overall_eff} Value of phi
# @param {sigma_cluster_sq} Value of the random cluster intercept
# (sigma_alpha^2)
# @param {sigma_expt_sq} Value of the random exposure time slope
# (sigma_delta^2)
# @param {sigma_epsilon_sq} Value of the individual level heterogeneity
# @param {n_per} The average number of individuals per cluster-time 
# period
# @param {sd_n_per} Standard deviation of number of individuals per cluster-time 
# period (0=all equal)
# @param {t_max} Number of time periods
# @param {k_max} [default] Number of clusters (defaults to t_max-1)
# @param {crossover_t} [default] Vector with time of cross-over
# for each cluster (defaults to c(2,3,4,...,t_max))
# 
# @returns {ret} A data.frame with the p-values from the Model 4
# LR test, permutation test, and Model 5 LR test, respectively.
sim <- function(intercept, 
                betas, 
                overall_eff,
                sigma_cluster_sq,
                sigma_expt_sq,
                sigma_epsilon_sq,
                n_per,
                sd_n_per,
                t_max,
                k_max,
                crossover_t
) {

  if (length(crossover_t) != k_max) {
    print("Error: crossover_t must be of length k_max")
    return()
  }
  
  count = round(rnorm(t_max*k_max, n_per, sd_n_per))
  count[count < 2] = 2
  freq_table = expand.grid(t = rep(1:t_max),
                           k = rep(1:k_max))
  freq_table$count = count
  
  data = freq_table[rep(1:nrow(freq_table), freq_table[["count"]]), ]
  N = nrow(data)
  rownames(data) = 1:N
  data$id = 1:N
  
  data = data %>%
    mutate(id = 1:N,
           cross_t = k, 
           cross_t = as.numeric(
                      as.character(
                        factor(cross_t,
                            labels=crossover_t))))
  
  data$expt <- data$t - data$cross_t
  data[data$expt<0,]$expt = 0
  data$trt = ifelse(data$expt==0,0,1)
  
  cluster_ranef <- rnorm(k_max, 0, sqrt(sigma_cluster_sq))
 
  data$Y <- intercept +
    overall_eff*data$trt + 
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
  

  model_full_fixef <- lmer(Y ~ as.factor(expt) +
                              as.factor(t) + 
                              (1 | k),
                            data)
  model_full_ranef <- lmer(Y ~ trt +
                              as.factor(t) + 
                              (1 | k)+
                              (0 + trt | expt), 
                            data)
  model_reduced <- lmer(Y ~ trt +
                           as.factor(t) + 
                           (1 | k),
                         data)
  q_obs = anova(model_full_ranef,model_reduced)$Chisq[2]
  
  onePermute <- function() {
    data = data %>%
      filter(expt > 0) %>%
      group_by(k) %>%
      mutate(expt.perm = sample(expt)) %>%
      select(id, expt.perm) %>%
      right_join(data, by=c("id","k")) %>%
      mutate(expt.perm = replace_na(expt.perm, 0))
    
    model_full_perm = lmer(Y ~ trt + 
                              as.factor(t) +
                              (0 + trt | expt.perm)+
                              (1 | k),
                            data=data)
    
    return(c(anova(model_full_perm,model_reduced)$Chisq[2]))
  }
  
  perm_dist = unlist(lapply(1:R, function(i) {
    onePermute()
  }))
  
  ret = data.frame(
    mod_4_lr_pval = anova(model_full_fixef,model_reduced)[2,8],
    permutation_pval = length(perm_dist[perm_dist>=q_obs])/R,
    mod_5_lr_pval = 1-pchibarsq(q_obs,2)
  )
  return(ret)
}

out <- foreach(iter = 1:nsims, .combine=rbind) %do% sim(intercept, 
                                                    betas, 
                                                    overall_eff,
                                                    sigma_cluster_sq,
                                                    sigma_expt_sq,
                                                    sigma_epsilon_sq,
                                                    n_per,
                                                    sd_n_per,
                                                    t_max,
                                                    k_max,
                                                    crossover_t)
# To compute power
colMeans(out < 0.05)
