#-----------------------------------------------------#
#           Power script: Permutation test,           #
#           asymptotic chi-squared                    #
#-----------------------------------------------------#

#-----------------------------------------------------
# Loading/defining all important functions
#-----------------------------------------------------
start_time = Sys.time()

library(foreach)
library(doParallel)
library(lme4)
library(emdbook)
library(dplyr)
library(tidyr)

# Reading in parameter vector from the command line
params <- commandArgs(trailingOnly=TRUE)

# Number of time points observed
n.t <- as.numeric(params[1])

# SD of exposure time random effect
sd.expt.ranef <- as.numeric(params[2])

# SD of cluster heterogeneity
sd.n.cl <- as.numeric(params[3])

# Index
index <- as.numeric(params[4])

# Number of resamples
R = 1000
# Number of individuals per cluster per time point
n.per.bar <- 30
# SD of group random effect
sd.group.ranef <- 0.1

n.cl <- (n.t-1)


###########################
#    Helper functions     #
###########################

scale2 <- function(x, sd.expt.ranef) {
  return(sd.expt.ranef*(x-mean(x))/sd(x))
}

#########################################
#   Set seed for one-time generation    #
#   of random effects (normal and non). #
#   Then, set seed to unique simulation #
#   iteration (to ensure variability)   #
#########################################
set.seed(123)
seeds <- runif(125, 1,10000)

expt.ranef <- c(0,
                rnorm(n.t-1))

expt.ranef[-1] <- scale2(expt.ranef[-1],sd.expt.ranef)
count = round(rnorm(n.t*n.cl, n.per.bar, sd.n.cl))
count[count < 2] = 2

# Reset seed (different seed for each job array) so that
# simulations are reproducible
set.seed(seeds[index],kind = "L'Ecuyer-CMRG")

one_run = function(n.t, 
                   n.per.bar,
                   expt.ranef,
                   sd.group.ranef,
                   R) {
  # Create stepped wedge structure, where there is one 
  # cluster per time trajectory
  # As of now, clusters all of an equal number of individuals
  freq_table = expand.grid(t = rep(1:n.t),
                           group= rep(1:n.cl))
  freq_table$count = count
  
  data = freq_table[rep(1:nrow(freq_table), freq_table[["count"]]), ]
  N = nrow(data)
  rownames(data) = 1:N
  data$id = 1:N

  data$tau <- data$t-data$group
  data[data$tau<0,]$tau = 0
  data$trt = ifelse(data$tau==0,0,1)
  
  # Fixed effect for background calendar time
  t.fixef <- sin(pi*(1:n.t-1)/(2*(n.t-1)))
  
  # Random effect for cluster/group
  group.ranef <- rnorm(n.cl, 0, sd.group.ranef)
  
  # Model on link-scale
  data$g <- log(0.7/(1-0.7)) +
    log(1.2)*data$trt + 
    apply(data,
          1,
          function(x) {
              group.ranef[as.numeric(x["group"])]+
              t.fixef[as.numeric(x["t"])]+
              expt.ranef[as.numeric(x["tau"])+1]
          })
  
  # Convert to binary outcome
  data$Y <- rbinom(N, 1, exp(data$g)/(1+exp(data$g)))
  
  model.full.fixef <- glmer(Y ~ as.factor(tau) +
                              as.factor(t) + 
                              (1 | group),
                            data, family=binomial)
  model.full.ranef <- glmer(Y ~ trt +
                        as.factor(t) + 
                        (1 | group)+
                        (0 + trt | tau), 
                      data, family=binomial)
  model.reduced <- glmer(Y ~ trt +
                           as.factor(t) + 
                           (1 | group),
                         data, family=binomial)
  q_obs = anova(model.full.ranef,model.reduced)$Chisq[2]
  
  onePermute <- function() {
    data = data %>%
      filter(tau > 0) %>%
      group_by(group) %>%
      mutate(tau.perm = sample(tau)) %>%
      select(id, tau.perm) %>%
      right_join(data, by=c("id","group")) %>%
      mutate(tau.perm = replace_na(tau.perm, 0))

    model.full.perm = glmer(Y ~ trt + 
                              as.factor(t) +
                              (0 + trt | tau.perm)+
                              (1 | group),
                            family="binomial",
                            data=data)
    
    return(c(anova(model.full.perm,model.reduced)$Chisq[2]))
  }
  
  perm_dist = unlist(lapply(1:R, function(i) {
    onePermute()
  }))
  
  ret = data.frame(
    n.t = n.t,
    n.per = n.per.bar,
    sd.expt.ranef = sd.expt.ranef,
    sd.group.ranef = sd.group.ranef,
    sd.n.cl = sd.n.cl,
    fixef = anova(model.full.fixef,model.reduced)[2,8],
    perm = length(perm_dist[perm_dist>=q_obs])/R,
    chisq = 1-pchibarsq(q_obs,2)
  )

  return(ret)
}


# Running the simulation 10 times in parallel
registerDoParallel(cores=10)
out <- foreach(n = 1:10, .combine=rbind) %dopar% one_run( n.t, 
                                                          n.per.bar,
                                                          expt.ranef,
                                                          sd.group.ranef,
                                                          R
                                                          )
end_time = Sys.time()
out$time = difftime(end_time, start_time, units = "mins")

# Creating and writing output dataframe
setwd("results")
out.name <- paste0(paste("power", 
                         n.t, 
                         sd.expt.ranef,
                         sd.n.cl,
                         R,
                         index, sep="_"), ".csv")
write.csv(out, out.name, row.names=FALSE)
