##############################################################
#    Code to compare variances of treatment effect parameter #
# based on Model 1 (no treatment effect het.) vs. Model 4    #
#         (accounts for treatment effect het.)               #
#         Note: for continuous outcomes only                 #
#         Used to create data for Figure S7                  #
##############################################################
library(ggplot2)
library(Matrix)
library(dplyr)
# Compute variance of treatment effect parameter from 
# Model 1
compute_model1_var <- function(t_max,
                               k_max,
                               expt_max,
                               n_per,
                               crossover_t,
                               icc,
                               sigma_epsilon_sq
                               ) {
  sigma_alpha_sq = sigma_epsilon_sq*(icc/(1-icc))
  
  data = expand.grid(t = 1:t_max,
                     k = 1:k_max,
                     i = 1:n_per)
  data = data %>%
    mutate(cross_t = k, 
           cross_t = as.numeric(
             as.character(
               factor(cross_t,
                      labels=crossover_t))))
  
  data$expt <- data$t - data$cross_t 
  if (t_max > 2) {
    data[data$expt<0,]$expt = 0
  }
  data$trt = ifelse(data$expt==0,0,1)
  
  data = data[order(data$k),]
  XWX_k = lapply(1:k_max, function(x) {
    X <- model.matrix(~as.factor(t) + trt,data=data[data$k == x,])
    V <- diag(rep(sigma_epsilon_sq,nrow(data[data$k == x,]))) + sigma_alpha_sq 
    return(t(X) %*% solve(V) %*% X)
  })
  
  return(solve(Reduce('+', XWX_k))[t_max+1,t_max+1])
}


# 
# Compute variance for average treatment effect estimator
# based on Model 4, derivation in Web Appendix
#
compute_model4_var <- function(t_max,
                               k_max,
                               expt_max,
                               n_per,
                               crossover_t,
                               icc,
                               sigma_epsilon_sq) {
  sigma_alpha_sq = sigma_epsilon_sq*(icc/(1-icc))
  data = expand.grid(t = 1:t_max,
                     k = 1:k_max,
                     i = 1:n_per)
  data = data %>%
    mutate(cross_t = k, 
           cross_t = as.numeric(
             as.character(
               factor(cross_t,
                      labels=crossover_t))))
  
  data$expt <- data$t - data$cross_t 
  if (t_max > 2) {
    data[data$expt<0,]$expt = 0
  }
  data$trt = ifelse(data$expt==0,0,1)
  
  
  V_star <-  diag(rep(sigma_epsilon_sq,nrow(data[data$k == 1,]))) + sigma_alpha_sq
  firstterm = 0
  V_star_inv = solve(V_star)
  
  for (i in 1:k_max) {
    Z_k <- rbind(matrix(0,ncol=expt_max,nrow=i),
                 diag(1,nrow=expt_max)[1:(t_max-i),])
    D_k <- kronecker(Z_k, rep(1,n_per))
    firstterm = firstterm + t(D_k) %*% V_star_inv %*% D_k
  }
  
  B = 0 
  for (i in 1:k_max) {
    Z_k <- rbind(matrix(0,ncol=expt_max,nrow=i),
                 diag(1,nrow=expt_max)[1:(t_max-i),])
    B = B + kronecker(Z_k, rep(1,n_per))
  }
  
  G = kronecker(diag(1,nrow=t_max), rep(1,n_per))
  second_first = t(B) %*% V_star_inv %*% G
  second_middle = solve(t(G) %*% V_star_inv %*% G)
  second_last = t(G) %*% V_star_inv %*% B
  
  secondterm = (1/k_max) *(second_first %*% second_middle %*% second_last)
  
  
  C = matrix(rep(1,expt_max)/expt_max,ncol=1)
  return( t(C) %*% solve(firstterm - secondterm) %*% C )
}

#
# Example of how to use code for a SW-CRT with 8 time 
# periods, 7 clusters, first crossover in the second
# time period, 20 individuals per cluster-step an 
# icc of 0.01 and individual level heterogeneity of 1
#

compute_model1_var(
                t_max = 8,
                k_max = 7,
                expt_max = 7,
                n_per = 20,
                crossover_t = 2:8,
                icc = 0.01,
                sigma_epsilon_sq = 1)

compute_model4_var(
  t_max = 8,
  k_max = 7,
  expt_max = 7,
  n_per = 20,
  crossover_t = 2:8,
  icc = 0.01,
  sigma_epsilon_sq = 1)

