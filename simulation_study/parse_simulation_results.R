library(dplyr)
library(plyr)
library(readr)
library(xtable)
library(ggplot2)
library(e1071)
library(cowplot)
library(tidyr)

# set working directory to source file location

################################################
# Power simulations, presented in Section 4.2  #
################################################

# Uncomment to parse raw simulation results
# power <- list.files(path = "./tests/results",     
#                      pattern = "power_", 
#                      full.names = TRUE) %>% 
#   lapply(read_csv) %>%                                          
#   bind_rows  
# 
# save(power,
#      file="./tests/results.Rda")

# load simulation results that have been pre-parsed with 
# the above code
load("./tests/power.Rda")

power %>% 
  group_by(n.t, sd.expt.ranef, sd.n.cl) %>%
  dplyr::summarize( perm=mean(perm < 0.05), 
                    chisq = mean(chisq < 0.05),
                    fixef = mean(fixef < 0.05),
                    n=n(),
                    time=mean(time)) %>%
  View()

# 
# Create Figure 1
#
power %>% 
  group_by(n.t, sd.expt.ranef) %>%
  filter(sd.n.cl == 0) %>%
  dplyr::summarize( perm=mean(perm < 0.05), 
                    chisq = mean(chisq < 0.05),
                    fixef = mean(fixef < 0.05)) %>%
  gather(test, power, perm, chisq, fixef) %>%
  mutate(test = factor(test, levels = c("fixef", "chisq", "perm"),
                       labels = c("LR Model 4", "LR Model 5", "Permutation")),
         n_expt = factor(n.t,
                         labels = c("E=3","E=5","E=7","E=9")),
         model = factor(test != "LR Model 4", labels=c("Model 4", 
                                                       "Model 5")),
         type_test = factor(test != "Permutation", labels = c("Permutation",
                                                              "LR"
         )),
         sd.power = sqrt(power*(1-power))/sqrt(500),
         width = case_when(as.numeric(n_expt) == 1 ~ 0.8/15,
                           as.numeric(n_expt) == 2 ~ 0.3/15,
                           as.numeric(n_expt) == 3 ~ 0.15/15,
                           TRUE ~ 0.1/15)) %>%
  ggplot(aes(x=sd.expt.ranef^2,y=power*100,col=test,linetype=test,shape=test)) +
  geom_hline(yintercept = 80, col="black")+
  geom_line(alpha=0.8) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 5, col="red")+
  geom_errorbar(aes(ymin=(power - 1.96*sd.power)*100, ymax=(power + 1.96*sd.power)*100,
                    width = width)) + 
  scale_y_continuous(breaks=c(5,80)) +
  scale_colour_manual(name = "Test",
                      labels = c("LR Model 4", "LR Model 5", "Permutation"),
                      values = c("#0072B2","#E69F00", "#E69F00")) +
  scale_linetype_manual(name = "Test",
                        labels = c("LR Model 4", "LR Model 5", "Permutation"),
                        values = c("twodash", "dotted", "solid")) +
  scale_shape_manual(name = "Test",
                     labels = c("LR Model 4", "LR Model 5", "Permutation"),
                     values = c(0, 17, 16)) +
  labs(x=expression(paste("Heterogeneity of exposure time effect (",
                          sigma[delta]^2,")",sep="")), y="Empirical power") +
  facet_wrap(~n_expt, scales = "free",nrow=2
  ) 

ggsave(file = "./tests/power.eps",
       width=6,
       height=5,
       device=cairo_ps)


####################################
# Continuous sims presented in     #
# Section 4.2 and Web Appendix     #
####################################
# Uncomment to parse raw simulation results
# cont_all <- list.files(path = "./cont_overall_estimation/results",
#                              pattern = "cont_param", full.names = TRUE) %>%
#   lapply(read_csv) %>%
#   bind_rows %>%
#   mutate(bias = (est_param-param),
#          width = upper-lower,
#          coverage_model = ifelse(param > lower & param < upper,1,0),
#          coverage_boots_normal_k = ifelse(param > lower_boots_normal_k & param < upper_boots_normal_k,1,0),
#          coverage_boots_normal_k_t = ifelse(param > lower_boots_normal_k_t & param < upper_boots_normal_k_t,1,0)
#   ) 
# 
# cont_all[cont_all$scenario == 1,]$sd_expt = 0
# save(cont_all, file="./cont_overall_estimation/cont_all.Rda")

load("./cont_overall_estimation/cont_all.Rda")

cont_all %>%
  filter(expt == 0) %>%
  group_by(expt, model, scenario, n_per, sd_expt, t_max) %>%
  dplyr::summarize(n=n(), conv=mean(conv),
                   param=mean(param)) %>%
  View()

#
# Table 3
#
cont_8_7_100 = cont_all %>%
  filter(conv == 1, B > 0, expt == 0, t_max == 8, n_per == 100) %>%
  mutate(scenario = factor(scenario, levels = c(1,10,2,4))) %>%
  group_by(scenario, model) %>%
  dplyr::summarize( 
    mean_param = mean(est_param),
    mean_sigma_alpha = mean(sqrt(est_sigma_alpha_sq)),
    mean_sigma_delta = mean(sqrt(est_sigma_delta_sq),na.rm=TRUE),
    empirical_se = sd(est_param),
    model_se = mean(est_se),
    boot_se_k = mean(se_boots_k),
    coverage_boots_normal_k = mean(coverage_boots_normal_k)*100,
    boot_se_k_t = mean(se_boots_k_t),
    coverage_boots_normal_k_t = mean(coverage_boots_normal_k_t)*100
  ) %>%
  dplyr::mutate(mean_sigma_delta = ifelse(is.na(mean_sigma_delta), '-',
                                          format(round(mean_sigma_delta,3),
                                                 nsmall = 3))) 
print(xtable(cont_8_7_100[,c(-1)],digits = c(0,0,3,3,3,3,3,3,1,3,1)
), include.rownames = FALSE)

#
# Table S1
#
cont_30_29_100 = cont_all %>%
  filter(conv == 1, B > 0, expt == 0, 
         t_max == 30, n_per == 100,
         (sd_expt == 2 | sd_expt == 0)) %>%
  mutate(scenario = factor(scenario, levels = c(1,10,2,4))) %>%
  group_by(scenario, model) %>%
  dplyr::summarize( 
    mean_param = mean(est_param),
    mean_sigma_alpha = mean(sqrt(est_sigma_alpha_sq)),
    mean_sigma_delta = mean(sqrt(est_sigma_delta_sq),na.rm=TRUE),
    empirical_se = sd(est_param),
    model_se = mean(est_se),
    # coverage_model=mean(coverage_model)*100,
    boot_se_k = mean(se_boots_k),
    coverage_boots_normal_k = mean(coverage_boots_normal_k)*100,
    boot_se_k_t = mean(se_boots_k_t),
    coverage_boots_normal_k_t = mean(coverage_boots_normal_k_t)*100
  ) %>%
  dplyr::mutate(mean_sigma_delta = ifelse(is.na(mean_sigma_delta), '-',
                                          format(round(mean_sigma_delta,3),
                                                 nsmall = 3)))
print(xtable(cont_30_29_100[,c(-1)],digits = c(0,0,3,3,3,3,3,3,1,3,1)
), include.rownames = FALSE)

#
# Table S2
#
cont_vary_sd_expt = cont_all %>%
  filter(conv == 1, 
         B > 0, 
         expt == 0, 
         t_max == 30, 
         (scenario %in% c(1,10))) %>%
  group_by(n_per, sd_expt, model) %>%
  dplyr::summarize( 
    mean_param = mean(est_param),
    mean_sigma_alpha = mean(sqrt(est_sigma_alpha_sq)),
    mean_sigma_delta = mean(sqrt(est_sigma_delta_sq),na.rm=TRUE),
    empirical_se = sd(est_param),
    model_se = mean(est_se),
    # coverage_model=mean(coverage_model)*100,
    boot_se_k = mean(se_boots_k),
    coverage_boots_normal_k = mean(coverage_boots_normal_k)*100,
    boot_se_k_t = mean(se_boots_k_t),
    coverage_boots_normal_k_t = mean(coverage_boots_normal_k_t)*100
  ) %>%
  dplyr::mutate(mean_sigma_delta = ifelse(is.na(mean_sigma_delta), '-',
                                          format(round(mean_sigma_delta,3),
                                                 nsmall = 3)))

print(xtable(cont_vary_sd_expt,digits = c(0,0,1,0,3,3,3,3,3,3,1,3,1)
), include.rownames = FALSE)

################################################
# Binary outcome simulations:                  #
# Average treatment effect, then               #
# exposure-time-specific treatment effects     #
################################################

# bin_vary <- list.files(path = "./binary_overall_estimation/results_temp",
#                        pattern = "from_data", full.names = TRUE) %>%
#   lapply(read_csv) %>%
#   bind_rows %>%
#   mutate(bias = (est_param-param),
#          width = upper-lower,
#          coverage_model = ifelse(param > lower & param < upper,1,0),
#          coverage_boots_normal_k = ifelse(param > lower_boots_normal_k & param < upper_boots_normal_k,1,0),
#          coverage_boots_normal_k_t = ifelse(param > lower_boots_normal_k_t & param < upper_boots_normal_k_t,1,0)
#   ) 
# 
# bin_42 <- list.files(path = "./binary_overall_estimation/results",
#                      pattern = "bin_param_8_7_42", full.names = TRUE) %>%
#   lapply(read_csv) %>%
#   bind_rows %>%
#   mutate(bias = (est_param-param),
#          width = upper-lower,
#          coverage_model = ifelse(param > lower & param < upper,1,0),
#          coverage_boots_normal_k = ifelse(param > lower_boots_normal_k & param < upper_boots_normal_k,1,0),
#          coverage_boots_normal_k_t = ifelse(param > lower_boots_normal_k_t & param < upper_boots_normal_k_t,1,0)
#   ) 
# 
# bin_all <- rbind(bin_vary,bin_42)
# bin_all %>%
#   filter(expt == 0) %>%
#   group_by(model, scenario, n_per) %>%
#   dplyr::summarize(n=n(), conv=mean(conv),
#                    param=mean(param)) %>%
#   View()
# 
# set.seed(123)
# bin_all_less = bin_all %>%
#   filter(expt == 0) %>%
#   group_by(model, scenario, n_per) %>%
#   sample_n(500)
# 
# save(bin_all_less, file="./binary_overall_estimation/bin_all_less.Rda")

load("./binary_overall_estimation/bin_all_less.Rda")
bin_all_less %>%
  filter(expt == 0) %>%
  group_by(model, scenario, n_per) %>%
  dplyr::summarize(n=n(), conv=mean(conv),
                   param=mean(param)) %>%
  View()

#
# Create Table 4: assessing estimation
# of overall treatment effect
bin_a_boot_table = bin_all_less %>%
  filter(conv == 1,B > 0, expt == 0) %>%
  mutate(scenario = factor(scenario, levels = c(1,10,2,4)),
         n_per = factor(n_per,
                        levels=c("from_data", 100, 500),
                        labels=c("34 [6-96]", 100, 500))) %>%
  group_by(scenario, k_max, n_per, model) %>%
  dplyr::summarize( 
    mean_param = mean(est_param),
    mean_sigma_alpha = mean(sqrt(est_sigma_alpha_sq)),
    mean_sigma_delta = mean(sqrt(est_sigma_delta_sq),na.rm=TRUE),
    empirical_se = sd(est_param),
    model_se = mean(est_se),
    boot_se_k = mean(se_boots_k),
    coverage_boots_normal_k = mean(coverage_boots_normal_k)*100,
    boot_se_k_t = mean(se_boots_k_t),
    coverage_boots_normal_k_t = mean(coverage_boots_normal_k_t)*100
  ) %>%
  dplyr::mutate(mean_sigma_delta = ifelse(is.na(mean_sigma_delta), '-',
                                          round(mean_sigma_delta,3))) 

print(xtable(bin_a_boot_table[c(-1)],digits = c(0,0,0,0,3,3,3,3,3,3,1,3,1)
), include.rownames = FALSE)

#
# Create MSE plots Figure 2 and Figures
#
# param_no_boots <- list.files(path = "./binary_expt_estimation/results_no_boots",
#                         pattern = "param_8_30", full.names = TRUE) %>%
#    lapply(read_csv) %>%
#    bind_rows %>%
#    mutate(bias = param-true,
#          percent_bias = (param-true)*100/abs(true),
#          width = exp(upper)-exp(lower),
#          coverage = ifelse(true > lower & true < upper,1,0))
# 
# set.seed(123)
# param_no_boots = param_no_boots %>%
#      group_by(scenario,model,number) %>%
#      sample_n(2000)
# save(param_no_boots, file="./binary_expt_estimation/param_no_boots.Rda")

load("./binary_expt_estimation/param_no_boots.Rda")

#
# Create Table S3: Background calendar time trends
#
n.t <- 8
time <- 1:n.t
# oscillating
t.fixef.osc <- 0.5*sin(pi*2*(time-1)/(n.t-1))
# increasing
t.fixef.inc <- 1.5*(sigmoid(time-1)-0.75)

logit <- function(x) { 
  return(log(x/(1-x)))
}

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}
p1 <- data.frame(time = time,
                 osc = t.fixef.osc,
                 inc = t.fixef.inc
) %>%
  gather(type, value, -time) %>%
  mutate(type = factor(type, 
                       levels = c("osc", "inc"),
                       labels = c("Pattern A", "Pattern B"))) %>%
  mutate(value = expit(logit(0.7) + value)) %>%
  ggplot(aes(x=time, y=value)) +
  geom_line() + 
  facet_wrap(~type)+
  ylab("Baseline prevalence")+
  xlab("t")

p2 <- param_no_boots %>%
  filter(number != 0, scenario %in% 1:12) %>%
  group_by(scenario, number) %>%
  dplyr::summarize(true = unique(true)) %>%
  mutate(scenario = factor(scenario, 
                           levels = 1:12,
                           labels = paste("Scenario", 1:12))) %>%
  ggplot(aes(x=number, y=true)) +
  geom_line() + 
  geom_point() +
  facet_wrap(~scenario)+
  labs(x=expression(e), y=expression(g(e))) 

plot_grid(p1, p2, labels=c("A", "B"),nrow=2)
# ggsave(file = "./binary_expt_estimation/etime_pattern.eps",
#        width=8,
#        height=12)

param_no_boots_summ = param_no_boots %>%
  filter(conv==1) %>%
  group_by(scenario,model,number) %>%
  dplyr::summarize(empirical.sd = sd(param),
                   mean_param= mean(param),
                   coverage=mean(coverage),
                   bias=mean(bias), 
                   width=mean(width),
                   est.sd = mean(sd),
                   true = mean(true),
                   lower= mean(lower),
                   upper=mean(upper)
  ) 

#
# Create Figure 2
#
param_no_boots_summ %>%
  filter(number != 0,
         model %in% 4:5) %>% 
  mutate(mse = bias^2 + empirical.sd^2,
         scenario = factor(scenario, levels=1:24, 
                           labels=c(paste("Scenario A",1:12,sep=""),
                                    paste("Scenario B",1:12,sep="")))) %>%
  ggplot(aes(x=number, y=mse, col=factor(model))) +
  geom_point(alpha=.6)+
  geom_line(aes(group=factor(model)),alpha=0.6)+
  labs(x="Exposure time", y="Mean squared error")+
  scale_colour_manual(name = "",
                      labels = c(
                        "Model 4", "Model 5"),
                      values = c(
                        "#0072B2", "#E69F00"))+
  facet_wrap(~scenario, scale="free_y", ncol=4)

# ggsave(file = "./binary_expt_estimation/expt_mse_45.eps",
#        width=8,
#        height=6,
#        device=cairo_ps)

#
# Create Supp Figures
# 
param_no_boots_summ %>%
  filter(number != 0, scenario <= 12
  ) %>% 
  mutate(mse = bias^2 + empirical.sd^2,
         scenario = factor(scenario, levels=1:24, 
                           labels=c(paste("Scenario A",1:12,sep=""),
                                    paste("Scenario B",1:12,sep="")))) %>%
  ggplot(aes(x=number, y=mse, col=factor(model))) +
  geom_point(alpha=.6)+
  geom_line(aes(group=factor(model)),alpha=.6)+
  coord_trans(y="log")+
  scale_y_continuous(breaks=c(0.1,0.2,0.5))+
  labs(x="Exposure time", y="Mean squared error")+
  scale_colour_manual(name = "",
                      labels = c("Model 1",
                                 "Model 2", "Model 3",
                                 "Model 4", "Model 5"),
                      values = c("#CC79A7",
                                 "#56B4E9", "#009E73",
                                 "#0072B2", "#E69F00"))+
  facet_wrap(~scenario, scale="free_y", ncol=4)

# ggsave(file = "./binary_expt_estimation/expt_mse_all_A.eps",
#        width=8,
#        height=6,
#        device=cairo_ps)

param_no_boots_summ %>%
  filter(number != 0, scenario > 12
  ) %>% 
  mutate(mse = bias^2 + empirical.sd^2,
         scenario = factor(scenario, levels=1:24, 
                           labels=c(paste("Scenario A",1:12,sep=""),
                                    paste("Scenario B",1:12,sep="")))) %>%
  ggplot(aes(x=number, y=mse, col=factor(model))) +
  geom_point(alpha=.6)+
  geom_line(aes(group=factor(model)),alpha=.6)+
  coord_trans(y="log")+
  scale_y_continuous(breaks=c(0.1,0.2,0.5))+
  labs(x="Exposure time", y="Mean squared error")+
  scale_colour_manual(name = "",
                      labels = c("Model 1",
                                 "Model 2", "Model 3",
                                 "Model 4", "Model 5"),
                      values = c("#CC79A7",
                                 "#56B4E9", "#009E73",
                                 "#0072B2", "#E69F00"))+
  facet_wrap(~scenario, scale="free_y", ncol=4)

# ggsave(file = "./binary_expt_estimation/expt_mse_all_B.eps",
#        width=8,
#        height=6,
#        device=cairo_ps)

#
# Create scenario 10 plot
#


# bin_boots_10 <- list.files(path = "./binary_expt_estimation/results",
#                            pattern = "param_8_7_7_30_100", full.names = TRUE) %>%
#   lapply(read_csv) %>%
#   bind_rows %>%
#   mutate(bias = (est_param-true),
#          width = upper-lower,
#          coverage_model = ifelse(true > lower & true < upper,1,0),
#          coverage_boots_normal_k = ifelse(true > lower_boots_normal_k & true < upper_boots_normal_k,1,0)
#   )
# 
# save(bin_boots_10, file="./binary_expt_estimation/bin_boots_10.Rda")
load("./binary_expt_estimation/bin_boots_10.Rda")

bin_boots_10 %>%
  group_by(model,expt,scenario) %>%
  dplyr::summarize(n=n(), conv=mean(conv)) 

bin_boots_summ = bin_boots_10 %>%
  filter(conv==1) %>%
  group_by(scenario,model,expt) %>%
  dplyr::summarize(est_param= mean(est_param),
                   coverage_boots_normal_k=mean(coverage_boots_normal_k),
                   bias=mean(bias), 
                   width=mean(width)
  ) 


bin_boots_summ %>% 
  filter(expt != 0, scenario %in% c(10,22)) %>% 
  group_by(model) %>%
  dplyr::summarize(mean(coverage_boots_normal_k))

bin_boots_summ %>% 
  filter(expt != 0, model %in% c(1:5),scenario %in% c(10,22)) %>% 
  mutate(p_coverage=coverage_boots_normal_k,
          coverage=100*coverage_boots_normal_k,
         scenario = factor(scenario, levels=1:24,
                           labels=c(paste("Scenario A",1:12,sep=""),
                                    paste("Scenario B",1:12,sep="")))) %>%
  ggplot(aes(x=expt, y=coverage, col=factor(model))) +
  geom_errorbar(aes(ymin = 100*(p_coverage - 1.96*sqrt(p_coverage*(1-p_coverage)/1000)), 
                    ymax = 100*(p_coverage + 1.96*sqrt(p_coverage*(1-p_coverage)/1000)))) +
  geom_hline(yintercept = 95)+
  geom_hline(yintercept = 96.4, linetype="dotted")+
  geom_hline(yintercept = 93.6, linetype="dotted")+
  geom_point(aes(col=factor(model),size=width), alpha=0.6)+
  scale_colour_manual(name = "",
                      labels = c("Model 1",
                                 # "Model 2", "Model 3",
                                 "Model 4", "Model 5"),
                      values = c("#CC79A7",
                                 # "#56B4E9", "#009E73",
                                 "#0072B2", "#E69F00"))+
  labs(x="Exposure time", y="Coverage", size="CI width")+
  facet_wrap(~scenario)

# ggsave(file = "./binary_expt_estimation/scen_10.eps",
#        width=8,
#        height=6,
#        device=cairo_ps)


