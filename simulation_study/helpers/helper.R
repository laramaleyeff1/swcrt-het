# 
# Runs one bootstrap resample and returns estimated parameter
# values for the overall effect and exposure-time-specific
# effects from Models 1, 4, and 5
# Bootstrap is within cluster
# @param {data} A simulated SW-CRT data.frame
# @returns A vector with the estimated average treatment effect, 
# followed by the estimated exposure-time-specific effects for
# each model
# - {model_1_boots_k} vector containing model 1 parameters 
# - {model_4_boots_k} vector containing model 4 parameters
# - {model_5_boots_k} vector containing model 5 parameters
one_boot_expt <- function(data) {
  # Take a sample with replacement within cluster
  data_boot_k = data %>%
    dplyr::group_by(k) %>%
    sample_frac(replace=TRUE)
  t_max = max(data$t)
  # Fit three models 
  model_1_k <- glmer(Y ~ trt +
                       as.factor(t) +
                       (1 | k),
                     data_boot_k,
                     family=binomial)
  
  
  model_4_k <- glmer(Y ~ as.factor(expt) +
                       as.factor(t) +
                       (1 | k),
                     data_boot_k,
                     family=binomial)
  
  model_5_k <- glmer(Y ~  trt +
                       as.factor(t) +
                       (0 + trt | expt) +
                       (1 |k),
                     data_boot_k,
                     family=binomial)
    
  return(c(
    model_1_boots_k = rep(fixef(model_1_k)[2],t_max),
    model_4_boots_k = c(mean(fixef(model_4_k)[2:t_max]),
                        fixef(model_4_k)[2:t_max]),
    model_5_boots_k = c(fixef(model_5_k)[2],
                        fixef(model_5_k)[2] + ranef(model_5_k)$expt$trt[-1])
  ))
}

# 
# Runs one bootstrap resample for within cluster boostrap
# and returns estimated parameter values for the average treatment effect 
# from models 1, 4, and 5
# 
# @param {data} A simulated SW-CRT data.frame
# @param {bin} Boolean indicating if outcome is binary or not
# @returns A vector with the estimated average treatment effect, 
# for each model
# - {model_1_boots_k} estimated average treatment effect from Model 1 
# in within cluster boostrap
# - {model_1_boots_k_t} hardcoded to 0
# - {model_4_boots_k} estimated average treatment effect from Model 4 
# in within cluster boostrap
# - {model_4_boots_k_t} hardcoded to 0
# - {model_5_boots_k} estimated average treatment effect from Model 5 
# in within cluster boostrap
# - {model_5_boots_k_t} hardcoded to 0
one_boot_cluster_only <- function(data, bin = FALSE) {
  # Take a sample with replacement within cluster
  data_boot_k = data %>%
    dplyr::group_by(k) %>%
    sample_frac(replace=TRUE)

  if (bin) {
    # For binary outcomes use glmer with 
    # family=binomial
    model_1_k <- glmer(Y ~ trt +
                         as.factor(t) +
                         (1 | k),
                       data_boot_k,
                       family=binomial)
    
    
    model_4_k <- glmer(Y ~ as.factor(expt) +
                         as.factor(t) +
                         (1 | k),
                       data_boot_k,
                       family=binomial)
   
    model_5_k <- glmer(Y ~  trt +
                         as.factor(t) +
                         (0 + trt | expt) +
                         (1 |k),
                       data_boot_k,
                       family=binomial)
    
 
  } else {
    # For continuous outcomes use lmer
    model_1_k <- lmer(Y ~ trt +
                        as.factor(t) +
                        (1 |k),
                      data_boot_k)
    
  
    
    model_4_k <- lmer(Y ~ as.factor(expt) +
                        as.factor(t) +
                        (1 | k),
                      data_boot_k)
    
    
  
    model_5_k <- lmer(Y ~  trt +
                        as.factor(t) +
                        (0 + trt | expt) +
                        (1 |k),
                      data_boot_k)
    
  }
  
  return(c(
    model_1_boots_k = fixef(model_1_k)[2],
    model_1_boots_k_t = 0,
    model_4_boots_k = mean(fixef(model_4_k)[2:t_max]),
    model_4_boots_k_t = 0,
    model_5_boots_k = fixef(model_5_k)[2],
    model_5_boots_k_t = 0
  ))
}

# 
# Runs one bootstrap resample for each type (within cluster
# and within cluster-step) and returns estimated parameter
# values for the average treatment effect from model 1-5 for
# continuous outcomes and models 1, 4, and 5 for binary outcomes
# @param {data} A simulated SW-CRT data.frame
# @param {bin} Boolean indicating if outcome is binary or not
# @returns A vector with the estimated average treatment effect, 
# for each model and each boostrap type
# - {model_1_boots_k} estimated average treatment effect from Model 1 
# in within cluster boostrap
# - {model_1_boots_k_t} estimated average treatment effect from Model 1 
# in within-step cluster boostrap
# - if (!bin) {model_2_boots_k} estimated average treatment effect from Model 2
# in within cluster boostrap
# - if (!bin) {model_2_boots_k_t} estimated average treatment effect from Model 2
# in within cluster-step boostrap
# - if (!bin) {model_3_boots_k} estimated average treatment effect from Model 3 
# in within cluster boostrap
# - if (!bin) {model_3_boots_k_t} estimated average treatment effect from Model 3 
# in within cluster-step boostrap
# - {model_4_boots_k} estimated average treatment effect from Model 4 
# in within cluster boostrap
# - {model_4_boots_k_t} estimated average treatment effect from Model 4 
# in within cluster-step boostrap
# - {model_5_boots_k} estimated average treatment effect from Model 5 
# in within cluster boostrap
# - {model_5_boots_k_t} estimated average treatment effect from Model 5 
# in within cluster-step boostrap
one_boot <- function(data, bin = FALSE) {
  data_boot_k = data %>%
    dplyr::group_by(k) %>%
    sample_frac(replace=TRUE) %>%
    arrange(k,t)
  data_boot_k_t = data %>%
    dplyr::group_by(k,t) %>%
    sample_frac(replace=TRUE) %>%
    arrange(k,t)
  
  if (bin) {
    model_1_k <- glmer(Y ~ trt +
                         as.factor(t) +
                         (1 |k),
                       data_boot_k,
                       family=binomial)
    
    model_1_k_t <- glmer(Y ~ trt +
                           as.factor(t) +
                           (1 |k),
                         data_boot_k_t,
                         family=binomial)
    
    
    model_4_k <- glmer(Y ~ as.factor(expt) +
                         as.factor(t) +
                         (1 | k),
                       data_boot_k,
                       family=binomial)
    
    
    model_4_k_t <- glmer(Y ~ as.factor(expt) +
                           as.factor(t) +
                           (1 | k),
                         data_boot_k_t,
                         family=binomial)
    
    
    model_5_k <- glmer(Y ~  trt +
                         as.factor(t) +
                         (0 + trt | expt) +
                         (1 |k),
                       data_boot_k,
                       family=binomial)
    
    model_5_k_t <- glmer(Y ~  trt +
                           as.factor(t) +
                           (0 + trt | expt) +
                           (1 | k),
                         data_boot_k_t,
                         family=binomial)
    return(c(
      model_1_boots_k = fixef(model_1_k)[2],
      model_1_boots_k_t = fixef(model_1_k_t)[2],
      model_4_boots_k = mean(fixef(model_4_k)[2:t_max]),
      model_4_boots_k_t = mean(fixef(model_4_k_t)[2:t_max]),
      model_5_boots_k = fixef(model_5_k)[2],
      model_5_boots_k_t = fixef(model_5_k_t)[2]
    ))
  } else {
    model_1_k <- lmer(Y ~ trt +
                         as.factor(t) +
                         (1 |k),
                       data_boot_k)
    
    model_1_k_t <- lmer(Y ~ trt +
                           as.factor(t) +
                           (1 |k),
                         data_boot_k_t)
    
    model_2_k <- lmer(Y ~ trt + 
                        expt_adj +
                        as.factor(t) +
                        (1 | k),
                      data_boot_k)
    
    model_2_k_t <- lmer(Y ~ trt + 
                          expt_adj +
                          as.factor(t) +
                          (1 | k),
                        data_boot_k_t)
    
    model_3_k <-lmer(Y ~ expt_factor_1 +
                                 expt_factor_2 +
                                 as.factor(t) +
                                 (1 | k),
                               data_boot_k)
    
    model_3_k_t <- lmer(Y ~ expt_factor_1 +
                          expt_factor_2 +
                          as.factor(t) +
                          (1 | k),
                        data_boot_k_t)
    
    
    model_4_k <- lmer(Y ~ as.factor(expt) +
                         as.factor(t) +
                         (1 | k),
                       data_boot_k)
    
    
    model_4_k_t <- lmer(Y ~ as.factor(expt) +
                           as.factor(t) +
                           (1 | k),
                         data_boot_k_t)
    
    
    model_5_k <- lmer(Y ~  trt +
                         as.factor(t) +
                         (0 + trt | expt) +
                         (1 |k),
                       data_boot_k)
    
    model_5_k_t <- lmer(Y ~  trt +
                           as.factor(t) +
                           (0 + trt | expt) +
                           (1 | k),
                         data_boot_k_t)

    return(c(
      model_1_boots_k = fixef(model_1_k)[2],
      model_1_boots_k_t = fixef(model_1_k_t)[2],
      model_2_boots_k = mean(fixef(model_2_k)[2] + 0:(t_max-2)*fixef(model_2_k)[3]),
      model_2_boots_k_t = mean(fixef(model_2_k_t)[2] + 0:(t_max-2)*fixef(model_2_k_t)[3]),
      model_3_boots_k = mean(c(fixef(model_3_k)[2], rep(fixef(model_3_k)[3], t_max - 2))),
      model_3_boots_k_t = mean(c(fixef(model_3_k_t)[2], rep(fixef(model_3_k_t)[3], t_max - 2))),
      model_4_boots_k = mean(fixef(model_4_k)[2:t_max]),
      model_4_boots_k_t = mean(fixef(model_4_k_t)[2:t_max]),
      model_5_boots_k = fixef(model_5_k)[2],
      model_5_boots_k_t = fixef(model_5_k_t)[2]
    ))
  }
  
  
}
# 
# The following five functions each create a data.frame 
# with the results of fitting Models 1-5 for a SW-CRT
# within either a continuous or binary outcome
# @param {data} A simulated SW-CRT data.frame
# @param {t_max} The number of time periods
# @param {se_boots_k} Within-cluster boostrap standard error, -1 if boostrap
# not run
# @param {se_boots_k_t} Within-cluster-step boostrap standard error
# @param {bin} TRUE if outcome is binary, FALSE otherwise
# @returns {ret} A data.frame with 
# - {model} Numeric value of model
# - {expt} Indexes the exposure time effect (0=overall)
# - {est_param} Value of the exposure time effect
# - {est_se} Standard deviation of {est_param}
# - {lower} Lower bound of 95% CI for {est_param}
# - {upper} Upper bound of 95% CI for {est_param}
# - {est_sigma_alpha_sq} Estimated value of sigma_alpha_sq
# - {est_sigma_delta_sq} Estimated value of sigma_delta_sq
# (only for Model 5)
# if se_boots_k != -1 (i.e. if we ran boostraps) return {
# - {se_boots_k} Within-cluster bootstrap standard error
# - {lower_boots_normal_k} Lower bound for within-cluster
# boostrap 95% CI for {est_param}
# - {upper_boots_normal_k} Upper bound for within-cluster
# boostrap 95% CI for {est_param}
# - {se_boots_k_t} Within-cluster-step bootstrap standard error
# - {lower_boots_normal_k_t} Lower bound for within-cluster-step
# boostrap 95% CI for {est_param}
# - {upper_boots_normal_k_t} Upper bound for within-cluster-step
# boostrap 95% CI for {est_param}
# }
# - {conv} Indicator for whether or not the model
#   converged (1=yes, 0=no)
fit_model_1 <- function(data, t_max, se_boots_k, se_boots_k_t, bin = FALSE) {
  if (bin) {
    model <- glmer(Y ~ trt +
                    as.factor(t) + 
                    (1 | k), 
                  data,
                  family=binomial)
  } else {
    model <- lmer(Y ~ trt +
                    as.factor(t) + 
                    (1 | k), 
                  data)
  }
  
  param = fixef(model)[2]
  se_param = summary(model)$coefficients[2,2]
  
  if (se_boots_k == -1) {
    ret = data.frame(
      model = 1,
      expt = 0:(t_max-1),
      est_param = param,
      est_se  = se_param,
      lower = param -  
        qnorm(0.975)*se_param,
      upper = param +  
        qnorm(0.975)*se_param,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0)
    )
  } else {
    ret = data.frame(
      model = 1,
      expt = 0:(t_max-1),
      est_param = param,
      est_se  = se_param,
      lower = param -  
        qnorm(0.975)*se_param,
      upper = param +  
        qnorm(0.975)*se_param,
      se_boots_k = se_boots_k,
      lower_boots_normal_k = param -  
        qnorm(0.975)*se_boots_k,
      upper_boots_normal_k = param +  
        qnorm(0.975)*se_boots_k,
      se_boots_k_t = se_boots_k_t,
      lower_boots_normal_k_t = param -  
        qnorm(0.975)*se_boots_k_t,
      upper_boots_normal_k_t = param +  
        qnorm(0.975)*se_boots_k_t,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0)
    )
  }
  
  return(ret)
}

fit_model_2 <- function(data, t_max, se_boots_k, se_boots_k_t, bin = FALSE) {
  if (bin) {
    model <- glmer(Y ~ trt + 
                    expt_adj +
                    as.factor(t) +
                    (1 | k),
                  data,
                  family=binomial)
  } else {
    model <- lmer(Y ~ trt + 
                    expt_adj +
                    as.factor(t) +
                    (1 | k),
                  data)
  }
  
  params = fixef(model)[2:3]
  vcov_params = vcov(model)[2:3,2:3]
  
  expt_effs = params[1] + 0:(t_max-2)*params[2]
  overall_eff = mean(expt_effs)
  
  se_expt_effs = sqrt(
    vcov_params[1,1] + 
      (0:(t_max-2))^2*vcov_params[2,2] + 
      2*0:(t_max-2)*vcov_params[1,2]
  )
  
  const = mean(0:(t_max-2))
  se_overall_eff = sqrt(
    vcov_params[1,1] + const^2*vcov_params[2,2] + 
      2*const*vcov_params[1,2]
  )
  effs = c(overall_eff, expt_effs)
  se_effs = c(se_overall_eff, se_expt_effs)
  
  if (se_boots_k == -1) {
    ret = data.frame(
      model = 2,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
    
  } else {
    ret = data.frame(
      model = 2,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      se_boots_k = se_boots_k,
      lower_boots_normal_k = effs -  
        qnorm(0.975)*se_boots_k,
      upper_boots_normal_k = effs +  
        qnorm(0.975)*se_boots_k,
      se_boots_k_t = se_boots_k_t,
      lower_boots_normal_k_t = effs -  
        qnorm(0.975)*se_boots_k_t,
      upper_boots_normal_k_t = effs +  
        qnorm(0.975)*se_boots_k_t,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  }
  
  return(ret)
}

fit_model_3 <- function(data, t_max, se_boots_k, se_boots_k_t, bin = FALSE) {
  if (bin) {
    model <- glmer(Y ~ expt_factor_1 +
                    expt_factor_2 +
                    as.factor(t) +
                    (1 | k),
                  data,
                  family=binomial)
  } else {
    model <- lmer(Y ~ expt_factor_1 +
                    expt_factor_2 +
                    as.factor(t) +
                    (1 | k),
                  data)
  }
  
  params = fixef(model)[2:3]
  vcov_params = vcov(model)[2:3,2:3]
  
  expt_effs = c(params[1], rep(params[2], t_max - 2))
  overall_eff = mean(expt_effs)
  
  se_expt_effs = sqrt(c(vcov_params[1,1], 
                   rep(vcov_params[2,2], t_max - 2)))
  C = c(1, (t_max-2))/(t_max-1)
  se_overall_eff = as.numeric(sqrt(C %*% vcov_params %*% C))
  
  effs = c(overall_eff, expt_effs)
  se_effs = c(se_overall_eff, se_expt_effs)
  
  if (se_boots_k == -1) {
    ret = data.frame(
      model = 3,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  } else {
    ret = data.frame(
      model = 3,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      se_boots_k = se_boots_k,
      lower_boots_normal_k = effs -  
        qnorm(0.975)*se_boots_k,
      upper_boots_normal_k = effs +  
        qnorm(0.975)*se_boots_k,
      se_boots_k_t = se_boots_k_t,
      lower_boots_normal_k_t = effs -  
        qnorm(0.975)*se_boots_k_t,
      upper_boots_normal_k_t = effs +  
        qnorm(0.975)*se_boots_k_t,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  }
  
  return(ret)
}

fit_model_4 <- function(data, t_max, se_boots_k, se_boots_k_t, bin = FALSE) {
  if (bin) {
    model <- glmer(Y ~ as.factor(expt) +
                    as.factor(t) +
                    (1 | k), 
                  data,
                  family=binomial)
  } else {
    model <- lmer(Y ~ as.factor(expt) +
                    as.factor(t) +
                    (1 | k), 
                  data)
  }
  
  expt_effs = fixef(model)[2:t_max]
  vcov_params = vcov(model)[2:t_max,2:t_max]
  
  overall_eff = mean(expt_effs)

  se_expt_effs = sqrt(diag(vcov_params))
  C = rep(1,(t_max-1))/(t_max-1)
  se_overall_eff = as.numeric(sqrt(C %*% vcov_params %*% C))
  
  effs = c(overall_eff, expt_effs)
  se_effs = c(se_overall_eff, se_expt_effs)
  
  if (se_boots_k == -1) {
    ret = data.frame(
      model = 4,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  } else {
    ret = data.frame(
      model = 4,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      se_boots_k = se_boots_k,
      lower_boots_normal_k = effs -  
        qnorm(0.975)*se_boots_k,
      upper_boots_normal_k = effs +  
        qnorm(0.975)*se_boots_k,
      se_boots_k_t = se_boots_k_t,
      lower_boots_normal_k_t = effs -  
        qnorm(0.975)*se_boots_k_t,
      upper_boots_normal_k_t = effs +  
        qnorm(0.975)*se_boots_k_t,
      est_sigma_alpha_sq = as.data.frame(VarCorr(model))[1,4],
      est_sigma_delta_sq = NA,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  }
  
  return(ret)
}

fit_model_5 <- function(data, t_max, se_boots_k, se_boots_k_t, bin = FALSE) {
  if (bin) {
    model <- glmer(Y ~ trt + 
                    as.factor(t) + 
                    (0 + trt | expt) +
                    (1 | k), 
                  data,
                  family = binomial)
  } else {
    model <- lmer(Y ~ trt + 
                    as.factor(t) + 
                    (0 + trt | expt) +
                    (1 | k), 
                  data)
  }
  
  overall_eff = fixef(model)[2]
  expt_effs = fixef(model)[2] + ranef(model)$expt$trt[-1]
  se_overall_eff = summary(model)$coefficients[2,2]
  
  coef = ranef(model, condVar=TRUE)
  se_expt_effs = unlist(c(sqrt(se_overall_eff^2+
                      attr (coef$expt, "postVar"))))[-1]
  effs = c(overall_eff, expt_effs)
  se_effs = c(se_overall_eff, se_expt_effs)

  df_varcorr = as.data.frame(VarCorr(model))
  if (se_boots_k == -1) {
    ret = data.frame(
      model = 5,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      est_sigma_alpha_sq = df_varcorr[df_varcorr$grp == "k",]$vcov,
      est_sigma_delta_sq = df_varcorr[df_varcorr$grp == "expt",]$vcov,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  } else {
    ret = data.frame(
      model = 5,
      expt = 0:(t_max-1),
      est_param = effs,
      est_se = se_effs,
      lower = effs - qnorm(0.975)*se_effs,
      upper = effs + qnorm(0.975)*se_effs,
      se_boots_k = se_boots_k,
      lower_boots_normal_k = effs -  
        qnorm(0.975)*se_boots_k,
      upper_boots_normal_k = effs +  
        qnorm(0.975)*se_boots_k,
      se_boots_k_t = se_boots_k_t,
      lower_boots_normal_k_t = effs -  
        qnorm(0.975)*se_boots_k_t,
      upper_boots_normal_k_t = effs +  
        qnorm(0.975)*se_boots_k_t,
      est_sigma_alpha_sq = df_varcorr[df_varcorr$grp == "k",]$vcov,
      est_sigma_delta_sq = df_varcorr[df_varcorr$grp == "expt",]$vcov,
      conv = ifelse(is.null(
        model@optinfo$conv$lme4$code), 1, 0))
  }
 
  return(ret)
}

