load("./binary_sim_results_or1.19.Rda")
load("./binary_sim_results_or3.Rda")

get_power <- function(sigma_sq, average_eff) {
  return((pnorm(sqrt(1/sigma_sq) * abs(average_eff) - qnorm(0.975)) + 
     pnorm(-1 * sqrt(1/sigma_sq) * abs(average_eff) - qnorm(0.975)))*100)
}

get_power(sd(binary_sim_results_or3$mod_5_param), log(3))
get_power(sd(binary_sim_results_or3$mod_4_param), log(3))
get_power(mean(binary_sim_results_or3$mod_4_se),log(3))

get_power(sd(binary_sim_results_or1.19$mod_5_param), log(1.19))
get_power(sd(binary_sim_results_or1.19$mod_4_param), log(1.19))
get_power(mean(binary_sim_results_or1.19$mod_4_se),log(1.19))
