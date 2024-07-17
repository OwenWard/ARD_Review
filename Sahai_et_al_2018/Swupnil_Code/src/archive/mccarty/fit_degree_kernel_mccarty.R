require("rstan");
options(mc.cores = parallel::detectCores());

# ---------------------------------
# ---------- Stan Models ----------
# ---------------------------------

kernel_continuous_fit <- stan_model(file = "Code_Stan/Degree_Kernel_Continuous.stan", 
                                    model_name = "Degree Kernel Continuous NonConstant");

# ----------------------------------
# ---------- Derived Data ----------
# ----------------------------------

age_max_2001 <- 84;
age_mean_2001 <- sum(c(0:age_max_2001) * (pop_raw_mccarty$Total)/(sum(pop_raw_mccarty$Total)));

g_k_mccarty <- c(rep(2, 6), rep(1, 6));
x_axis <- 2001-c(1917:2001);
mu_k_mccarty <- sigma_k_mccarty <- sum_prob_k_mccarty <- c();

for(k in 1:nrow(beta_names_alive_mccarty)) { 
  if(g_k[k] == 1) {
    # Male Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_mccarty[k] <- sum(beta_names_alive_mccarty[k,]/rev(pop_raw_mccarty$Female));
    y_axis <- (beta_names_alive_mccarty[k,]/rev(pop_raw_mccarty$Female))/sum_prob_k[k];
  }
  else {
    # Female Alters: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k_mccarty[k] <- sum(beta_names_alive_mccarty[k,]/rev(pop_raw_mccarty$Male));
    y_axis <- (beta_names_alive_mccarty[k,]/rev(pop_raw_mccarty$Male))/sum_prob_k[k];
  }
  # Ego Perspective: mu and sigma for p(j in G_k | a_j, g_j)/[sum_(a_j) p(j in G_k | a_j, g_j)]
  mu_k_mccarty[k] <- x_axis%*%y_axis/sum(y_axis);
  sigma_k_mccarty[k] <- sqrt((x_axis-mu_k_ego[k])^2 %*% y_axis/sum(y_axis));
}

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# Initial values
init_data_mccarty <- list();
for(i in 1:4) {
  init_data_mccarty[[i]] <- list(log_d = rep(6, nrow(names_data_mccarty)),
                                 inv_omega = rep(5/6, ncol(names_data_mccarty)),
                                 log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                                 rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                                 beta = rep(0, 3),
                                 log_eta_d = 0);
}

# Run the Model in Stan
mcmc_data_kernel_mccarty <- list(K = ncol(names_data_mccarty), 
                                 N = nrow(names_data_mccarty),  
                                 age_mean = age_mean_2001,
                                 degree_mean = degree_mean_sim, 
                                 y = names_data_mccarty, 
                                 age = mccarty$Age, 
                                 g_n = mccarty$Sex, 
                                 g_k = g_k_mccarty, 
                                 mu_k = mu_k_mccarty, 
                                 sigma_k = sigma_k_mccarty, 
                                 prop_k = sum_prob_k_mccarty, 
                                 mu_d = 6, 
                                 sigma_d = 0.6, 
                                 alpha_omega = 4.5, 
                                 beta_omega = 0.5, 
                                 mu_lambda = log(100), 
                                 sigma_lambda = 0.5, 
                                 alpha_rho = c(5,5), 
                                 recallPower = 0, 
                                 isDegreeRegression = 0);
fit_names_kernel_mccarty <- sampling(kernel_continuous_fit, 
                                     data = mcmc_data_kernel_mccarty, 
                                     iter = 1000, 
                                     chains = 4,
                                     init = init_data_mccarty);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

source("Code_Kernel/Kernel_Helper.R");

# ... Plot Degree Distributions ....

degree_mccarty <- round(colMeans(extract(fit_names_kernel_mccarty)$d));
omega_names_mccarty <- colMeans(extract(fit_names_kernel_mccarty)$inv_omega);

par(mfrow=c(2,3));
plot_degrees(degree_mccarty, ego_sex_age_mccarty, ego_names_verbose);

lambda_names_mccarty <- apply(extract(fit_names_kernel_mccarty)$lambda, c(2,3), mean);
rho_names_mccarty <- apply(extract(fit_names_kernel_mccarty)$rho, c(2,3), mean);
omega_names_mccarty <- apply(extract(fit_names_kernel_mccarty)$omega, 2, median);
beta_names_mccarty <- colMeans(extract(fit_names_kernel_mccarty)$beta);
eta_names_mccarty <- mean(extract(fit_names_kernel_mccarty)$eta);

# ... 2. Plot Kernel ...

par(mfrow=c(1,2))
plot_kernel(rev(lambda_names_mccarty[2,]), 
            rev(rho_names_mccarty[2,]), 
            y_cent = age_mean_2001,
            y_min = 0,
            y_max = 85,
            main = "Female Ego");
plot_kernel(rev(lambda_names_mccarty[1,]), 
            rev(rho_names_mccarty[1,]), 
            y_cent = age_mean_2001,
            y_min = 0,
            y_max = 85,
            main = "Male Ego");

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# ... Run the Model in Stan (Without Michael Robert David) ...
rm_names_mccarty <- c(1,5,6);

# Initial values
init_data_mccarty <- list();
for(i in 1:4) {
  init_data_mccarty[[i]] <- list(log_d = rep(6, nrow(names_data_mccarty[,-rm_names_mccarty])),
                                 inv_omega = rep(5/6, ncol(names_data_mccarty[,-rm_names_mccarty])),
                                 log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                                 rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                                 beta = rep(0, 3),
                                 log_eta_d = 0);
}

# ... Run the Model in Stan (Without Michael Robert David) ...
mcmc_data_kernel_mccarty_recalled <- list(K = ncol(names_data_mccarty[,-rm_names_mccarty]), 
                                          N = nrow(names_data_mccarty[,-rm_names_mccarty]),  
                                          age_mean = age_mean_2001, 
                                          degree_mean = degree_mean_sim, 
                                          y = names_data_mccarty[-rm_names_mccarty], 
                                          age = mccarty$Age, 
                                          g_n = mccarty$Sex, 
                                          g_k = g_k_mccarty[-rm_names_mccarty], 
                                          mu_k = mu_k_mccarty[-rm_names_mccarty], 
                                          sigma_k = sigma_k_mccarty[-rm_names_mccarty],
                                          prop_k = sum_prob_k_mccarty[-rm_names_mccarty],
                                          mu_d = 6, sigma_d = 0.6, 
                                          alpha_omega = 4.5, beta_omega = 0.5, 
                                          mu_lambda = log(100), sigma_lambda = 0.5, 
                                          alpha_rho = c(5,5), 
                                          recallPower = 0, 
                                          isDegreeRegression = 1);
fit_names_kernel_mccarty_recalled <- sampling(kernel_continuous_fit, 
                                              data = mcmc_data_kernel_mccarty_recalled, 
                                              iter = 1000, 
                                              chains = 4,
                                              init = init_data_mccarty);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

source("Code_Kernel/Kernel_Helper.R");

# ... Plot Degree Distributions ....

degree_mccarty_recalled <- round(colMeans(extract(fit_names_kernel_mccarty_recalled)$d));
omega_names_mccarty_recalled <- colMeans(extract(fit_names_kernel_mccarty_recalled)$inv_omega);

par(mfrow=c(2,3))
plot_degrees(degree_mccarty_recalled, ego_sex_age_mccarty, ego_names_verbose);

lambda_names_mccarty_recalled <- apply(extract(fit_names_kernel_mccarty_recalled)$lambda, c(2,3), median);
rho_names_mccarty_recalled <- apply(extract(fit_names_kernel_mccarty_recalled)$rho, c(2,3), median);
omega_names_mccarty_recalled <- apply(extract(fit_names_kernel_mccarty_recalled)$omega, 2, median);
