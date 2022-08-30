library(lme4)
library(emdbook)
library(dplyr)
library(tidyr)

# Example data.frame (simulated)
load("./data.Rda")

# B = number of boostrap resamples
B = 100
# R = number of permutations 
R = 2000
t_max = max(data$t)
expt_max = max(data$expt)

model_5 <- glmer(Y ~ trt+
                   as.factor(t)+
                   (0+trt|expt)+
                   (1|k),
                 data=data, 
                 family=binomial)

model_5_params <- data.frame(
  model = 5,
  number = 0:(t_max-1),
  log_or = c(fixef(model_5)[2],
            fixef(model_5)[2] + ranef(model_5)$expt$trt[-1]))

model_5_boots = matrix(nrow=t_max,ncol=B)

for (b in 1:B) {
  # Cluster-level bootstrap, for cluster-time boostrap use group_by(k,t)
  data_boot = data %>%
    group_by(k) %>%
    sample_frac(replace=TRUE)

  model_5_boot <- glmer(Y ~  trt +
                     as.factor(t) + 
                     (0 + trt | expt) +
                     (1 |k), 
                   data_boot,
                   family=binomial)
  model_5_boots[,b] = c(fixef(model_5_boot)[2], 
                        fixef(model_5_boot)[2] + ranef(model_5_boot)$expt$trt[-1])
}

model_5_params$se_boot_or = apply(model_5_boots,1,function(x) {return(sd(exp(x)))})
model_5_params$se_boot = apply(model_5_boots,1,sd)

model_5_params %>%
  mutate(format_or_se = paste(format(round(exp(log_or),2),nsmall=2),
                              " (",
                              format(round(se_boot_or,2),nsmall=2)
                              ,")", sep="")
  )

model_4 <- glmer(Y ~ as.factor(expt) +
                   as.factor(t) + 
                   (1 | k),
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

tests = data.frame(
  mod_4_lr_pval = anova(model_4,model_1)[2,8],
  permutation_pval = length(perm_dist[perm_dist>=q_obs])/R,
  mod_5_lr_pval = 1-pchibarsq(q_obs,2)
)

tests


