# Simulate and then fit with fewer and fewer names
library(gtools);
library(rstan);
load("data/mix_sim/mix_mat.Rdata");
load(here("Sahai_et_al_2018/", "Swupnil_Code", "data",
          "mix_sim", "mix_mat.Rdata"))

# Simulate fake data
n_names = 14;

omega_sim <- runif(14, 0.85, 1);
degree_sim <- degree$CombSexAge;
## this parameter is missing it seems


beta_sim <- matrix(rexp(14*8, 1/mean(beta_names[beta_names>0])), nrow = 8, ncol = 14);
colnames(beta_sim) <- paste("Name", 1:ncol(beta_sim));
rownames(beta_sim) <- paste("Alter", 1:nrow(beta_sim));

ego_sim <- ego_sex_age;

# 1. Test across a different number of names
M_sim <- list();
for(i in 1) {
  # simulate the data
  M_sim[[i]] <- matrix(c(rdirichlet(1, mixing.sim[1,]), rdirichlet(1, mixing.sim[2,]), rdirichlet(1, mixing.sim[3,]), rdirichlet(1, mixing.sim[4,]), rdirichlet(1, mixing.sim[5,]), rdirichlet(1, mixing.sim[6,])),
                  nrow = 6, ncol = 8, byrow = TRUE);
  rownames(M_sim[[i]]) <- paste("Ego", 1:nrow(M_sim[[i]]));
  colnames(M_sim[[i]]) <- paste("Alter", 1:ncol(M_sim[[i]]));
  
  names_data_sim <- simulate_mixing(degree_sim, omega_sim, M_sim[[i]], beta_sim, ego_sim);
  
  # Plot and fit matrix results
  par(mfrow=c(2,3));
  for(k in c(4,6,8,10,12,14)) {
    cat("\n Currently on Loop: ",i," !!\n");
    mcmc_data_sim_k <- list(E = 6, A = 8, K = k, N = nrow(names_data_sim), y = names_data_sim[,c(1:k)], ego = ego_sim, Beta = beta_sim[,c(1:k)], theta_d = c(6.2,.5), theta_o = c(3,2), alpha = rep(1,8)/8);
    fit_comb_k <- sampling(degree_mixing_fit, data = mcmc_data_sim_k, iter = 1000, chains = 4);
    plot_mix_comp(fit_comb_k, M_sim[[i]], paste(k,"Names"));
  }
}

# 2. Test Across 10 names but mixing around the names
# simulate the data
names_data_sim <- simulate_mixing(degree_sim, omega_sim, M_sim, beta_sim, ego_sim);

# Plot and fit matrix results
n_names <- 10;
index <- fit_comb_10_test <- list();
par(mfrow=c(2,5));
for(i in 1:9) {
  cat("\n Currently on Loop: ",i," !!\n");
  index[[i]] <- sample(1:14)[1:n_names];
  mcmc_data_sim_k <- list(E = 6, A = 8, K = n_names, N = nrow(names_data_sim), y = names_data_sim[,index[[i]]], ego = ego_sim, Beta = beta_sim[,index[[i]]], theta_d = c(6.2,.5), theta_o = c(3,2), alpha = rep(1,8)/8);
  fit_comb_10_test[[i]] <- sampling(degree_mixing_fit, data = mcmc_data_sim_k, iter = 1000, chains = 4);
  plot_mix_comp(fit_comb_10_test[[i]], M_sim, paste(index[[i]], collapse = ' '));
}

# Compare with the averaged value
plot_mix_comp(fit_comb_10_test, M_sim, "Averaged");
