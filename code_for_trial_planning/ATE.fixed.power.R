##############################################################
#    Code to compare variances of treatment effect parameter #
# based on Model 4 (accounts for treatment effect het.)      #
#         
#         Note: for continuous outcomes only                 #
#         Used in Web Appendix F.2                           #
##############################################################
#
# Computes variance for average treatment effect based 
# on Model 4
# @param {t_max} Number of time periods
# @param {k_max} Number of clusters 
# @param {expt_max} Number of exposure time periods 
# @param {n_per} Constant cluster-period size
# @param {crossover_t} Vector with time of cross-over
# for each cluster 
# @param {rho} ICC value
# @param {sigma_epsilon_sq} Individual level heterogeneity
# 
# @returns {sigma_sq} Value of variance for average 
# treatment effect
compute_model4_var <- function(t_max,
                               k_max,
                               n_per,
                               crossover_t,
                               rho,
                               sigma_epsilon_sq) {
  expt_max = t_max - min(crossover_t) + 1
  sigma_alpha_sq = sigma_epsilon_sq*(rho/(1-rho))
  lambda_1 = 1 - rho
  lambda_2  = 1 + (t_max*n_per - 1)*rho
  
  Zk_list <- lapply(1:k_max, function(k) {
    rbind(matrix(0,ncol=expt_max,nrow=(crossover_t[k]-1)),
          diag(1,nrow=expt_max)[1:(t_max-crossover_t[k]+1),])})
  ones_t <- matrix(rep(1,t_max),ncol=1)
  
  U_1 = Reduce("+",lapply(Zk_list, function(x) {t(x) %*% ones_t}))
  U_2 = Reduce("+",lapply(Zk_list, function(x) {t(x) %*% x}))
  W_1 = t(Reduce("+",Zk_list)) %*% Reduce("+",Zk_list)
  W_2 = Reduce("+",lapply(Zk_list, function(x) {t(x) %*% (ones_t %*% t(ones_t)) %*% x}))
  
  ones_e = matrix(rep(1,expt_max),nrow = 1)
  sigma_sq =  ((sigma_alpha_sq + sigma_epsilon_sq) * 
                 (k_max*t_max * lambda_1 * lambda_2) / 
                 (n_per * expt_max^2)) * 
    ones_e %*% 
    (solve(lambda_2 * (U_1 %*% t(U_1) + 
                         t_max * k_max * U_2 - 
                         t_max * W_1 - k_max * W_2) + 
             lambda_1 * (k_max * W_2 - 
                           U_1 %*% t(U_1)))) %*% t(ones_e) 
  return(sigma_sq)
}

compute_power_linear <- function(average_eff,
                                 t_max,
                                 k_max,
                                 expt_max,
                                 n_per,
                                 crossover_t,
                                 rho,
                                 sigma_epsilon_sq) {
  sigma_sq = compute_model4_var(
    t_max,
    k_max,
    expt_max,
    n_per,
    crossover_t,
    rho,
    sigma_epsilon_sq)
  return((pnorm(sqrt(1/sigma_sq) * abs(average_eff) - qnorm(0.975)) + 
            pnorm(-1 * sqrt(1/sigma_sq) * abs(average_eff) - qnorm(0.975)))*100 )
}

compute_detectable_effect <- function(t_max,
                                 k_max,
                                 expt_max,
                                 n_per,
                                 crossover_t,
                                 rho,
                                 sigma_epsilon_sq) {
  sigma_sq = compute_model4_var(
    t_max,
    k_max,
    expt_max,
    n_per,
    crossover_t,
    rho,
    sigma_epsilon_sq)
  return((qnorm(0.975)-qnorm(0.2))*sqrt(sigma_sq))
}

sizes = data.frame(rho=c(0,0.01,0.05,0.1,0.2))
sizes$size = sapply(sizes$rho, function(rho) {
  return(compute_detectable_effect(t_max = 8,
                                   k_max = 14,
                                   expt_max = 7,
                                   n_per = 34,
                                   crossover_t = rep(2:8,each=2),
                                   rho,
                                   sigma_epsilon_sq = 1))
})

sizes



