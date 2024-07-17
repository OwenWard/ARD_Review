# ---------------------------------------------------------------------
# ---------- Estimating Degree from Names Using Mixing Model ----------
# ---------------------------------------------------------------------

require("rstan");
options(mc.cores = parallel::detectCores());

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# ... Run the Model in Stan (10 Times with Random Names) ....
fit_degree_names <- rand_names <- list();
for(i in 1:10) {
  cat("\n Currently on Loop: ",i," !!\n");
  curr_names <- sample(1:12)[1:10];
  mcmc_data_names <- list(E = length(unique(ego_sex_age)), 
                          A = nrow(beta_names), 
                          K = ncol(beta_names[,curr_names]), 
                          N = nrow(names_data), 
                          y = names_data[,curr_names], 
                          ego = ego_sex_age, 
                          Beta = beta_names[,curr_names], 
                          theta_d = c(6.2,0.5), 
                          theta_o = c(10,2), 
                          alpha = rep(1,nrow(beta_names))/nrow(beta_names));
  rand_names[[i]] <- curr_names;
  fit_degree_names[[i]] <- sampling(degree_mixing_fit, 
                                    data = mcmc_data_names, 
                                    iter = 1000, 
                                    chains = 4);
}

# --------------------------------
# ---------- Assessment ----------
# --------------------------------

degree$NameSexAge <- round(colMeans(extract(fit_degree_names[[1]])$d));
omega_names_sex_age <- colMeans(extract(fit_degree_names[[1]])$inv_omega);

# ... Plot Degree Distributions ....

par(mfrow = c(2,3));
plot_degrees(degree$NameSexAge, ego_sex_age, ego_names_verbose);
plot_degrees(degree$NameSexAge, ego_sex_age_race, ego_names_verbose_race);

#
# ... Plot All Mixing Matrices ...
#

par(mfrow=c(3,3));
M_names_sex_age <- apply(extract(fit_degree_names[[1]])$M, c(2,3), mean);
plot_mixing(M_names_sex_age, 32, ego_names, paste(rand_names[[1]], collapse = ' '));
for(i in 2:8) {
  M_names_sex_age_i <- apply(extract(fit_degree_names[[i]])$M, c(2,3), mean);
  plot_mixing(M_names_sex_age_i, 32, ego_names, paste(rand_names[[i]], collapse = ' '));
  M_names_sex_age <- M_names_sex_age + M_names_sex_age_i;
  rm(M_names_sex_age_i);
}
M_names_sex_age <- M_names_sex_age/8;
plot_mixing(M_names_sex_age, 32, ego_names, "Averaged");

#
# ... Plot Mixing Matrices ...
#

par(mfrow=c(1,1));
M_names_sex_age <- apply(extract(fit_degree_names[[1]])$M, c(2,3), mean);
M_names_weighted <- M_names_sex_age;
for(i in 1:6) { M_names_weighted[i,] <- M_names_weighted[i,] * pop_age_sex[i]; }

# ... 1. Mixing Matrix for Sex and Age ... 

plot_mixing(M_names_sex_age, 32, ego_names);

# ... 2. Mixing Matrix for Sex ... 

M_names_sex <- matrix(0, nrow = 2, ncol = 8);
M_names_sex[1,] <- colSums(M_names_weighted[1:3,])/sum(pop_age_sex[1:3]);
M_names_sex[2,] <- colSums(M_names_weighted[4:6,])/sum(pop_age_sex[4:6]);

plot_mixing(M_names_sex, 11, c("Male", "Female"));

# ... 3. Mixing Matrix for Age ... 

M_names_age <- matrix(0, nrow = 3, ncol = 8);
M_names_age[1,] <- colSums(M_names_weighted[c(1,4),])/sum(pop_age_sex[c(1,4)]);
M_names_age[2,] <- colSums(M_names_weighted[c(2,5),])/sum(pop_age_sex[c(2,5)]);
M_names_age[3,] <- colSums(M_names_weighted[c(3,6),])/sum(pop_age_sex[c(3,6)]);

plot_mixing(M_names_age, 16, c("18-24", "25-64", "65+"));

