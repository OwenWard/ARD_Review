library(MASS)
par(mfrow = c(2,3))
ego_names_verbose <- c("Male 18-24","Male 25-64","Male 65+","Female 18-24","Female 25-64","Female 65+")

# Simulate and Fit Data from Names Model

names_data_sim <- simulate_mixing(degree$NameSexAge, omega_names_sex_age, M_names_sex_age, beta_names[,-rm_names], ego_sex_age)
plot_comparison(names_data[,-rm_names], names_data_sim, ego_sex_age, ego_names_verbose)

mcmc_data_names_sim <- list(E = length(unique(ego_sex_age)), A = nrow(beta_names), K = ncol(beta_names[,-rm_names]), 
                        N = nrow(names_data_sim), y = names_data_sim, ego = ego_sex_age, 
                        Beta = beta_names[,-rm_names], theta_d = c(6.2,.5), theta_o = c(3,2), 
                        alpha = rep(1,nrow(beta_names))/nrow(beta_names))
fit_names_sim <- sampling(degree_mixing_fit, data = mcmc_data_names_sim, iter = 1000, chains = 4)

# Simulate and Fit Data from Occupation Model

occ_data_sim <- simulate_mixing(degree$OccSexAge, omega_occ_sex_age, M_occ_sex_age, beta_occ_sex_age[,-rm_occ], ego_sex_age)
plot_comparison(occ_data[,-rm_occ], occ_data_sim, ego_sex_age, ego_names_verbose)

mcmc_data_occ_sim <- list(E = length(unique(ego_sex_age)), A = nrow(beta_names), K = ncol(beta_occ_sex_age[,-rm_occ]), 
                            N = nrow(occ_data_sim), y = occ_data_sim, ego = ego_sex_age, 
                            Beta = beta_occ_sex_age[,-rm_occ], theta_d = c(6.2,.5), theta_o = c(3,2), 
                            alpha = rep(1,nrow(beta_occ_sex_age))/nrow(beta_occ_sex_age), p = 1)
fit_occ_sim <- sampling(degree_mixing_fit, data = mcmc_data_occ_sim, iter = 1000, chains = 4)

# Simulate and Fit Data from Combined Model

names_data_sim_comb <- simulate_mixing(degree$CombSexAge, omega_comb_names_sex_age, M_comb_sex_age, beta_names[,-rm_names], ego_sex_age)
occ_data_sim_comb <- simulate_mixing(degree$CombSexAge, omega_comb_occ_sex_age, M_comb_sex_age, beta_occ_sex_age[,-rm_occ], ego_sex_age)
plot_comparison(names_data[,-rm_names], names_data_sim_comb, ego_sex_age, ego_names_verbose)
plot_comparison(occ_data[,-rm_occ], occ_data_sim_comb, ego_sex_age, ego_names_verbose)

mcmc_data_comb_sim <- list(E = length(unique(ego_sex_age)), A = nrow(beta_names), J = ncol(beta_names[,-rm_names]), K = ncol(beta_occ_sex_age[,-rm_occ]), 
                       N = nrow(names_data), y = names_data_sim, z = occ_data_sim, ego = ego_sex_age, BetaY = beta_names[,-rm_names], 
                       BetaZ = beta_occ_sex_age[,-rm_occ], theta_d = c(6.2,.5), theta_o = c(3,2), alpha = rep(1,nrow(beta_names))/nrow(beta_names), p = 1)
fit_comb_sim <- sampling(degree_combined_mixing_fit, data = mcmc_data_comb_sim, iter = 1000, chains = 4)