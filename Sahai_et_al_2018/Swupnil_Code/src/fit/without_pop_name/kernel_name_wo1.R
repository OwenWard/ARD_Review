# ---------------------------------------------------------------------
# ---------- Estimating Degree from Names Using Kernel Model ----------
# ---------------------------------------------------------------------

require("rstan")
options(mc.cores = parallel::detectCores())

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# Initial Values
init_data <- list();
for(i in 1:4) {
  init_data[[i]] <- list(log_d = rep(6, nrow(names_data[,-7])),
                         inv_omega = rep(5/6, ncol(names_data[,-7])),
                         log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                         rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                         beta = c(6, 0.1, -3),
                         log_eta = 0);
}

# Run the Model in Stan
mcmc_data <- list(K_1 = ncol(names_data[,-7]),
                  K_2 = 0,
                  N = nrow(names_data[,-7]),  
                  age_mean = age_mean_2014, 
                  Y = as.matrix(names_data[,-7]),  
                  w = data$wgt,
                  age = data$age, 
                  g_n = as.numeric(data$Gender), 
                  g_k = c(rep(1, 6), rep(2, 6))[-7], 
                  mu_k = mu_k_name[-7], 
                  sigma_k = sigma_k_name[-7], 
                  sum_prob_k = sum_prob_k_name[-7], 
                  mu_d = 6, sigma_d = 0.6, 
                  alpha_omega = 4.5, beta_omega = 0.5, 
                  mu_lambda = log(100), sigma_lambda = 0.5, 
                  alpha_rho = c(5,5), 
                  mu_beta = c(6, 0.1, -3), sigma_beta = rep(1, 3),
                  recall_power = 0, 
                  degree_regression = 1);
fit_names_kernel_wo1 <- sampling(model_kernel, 
                             data = mcmc_data, 
                             iter = 2000,
                             init = init_data);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

# ... 1. Plot Degree Distributions ....

degree$NameSexAgeKernelWo1 <- round(colMeans(extract(fit_names_kernel_wo1)$d));

par(mfrow=c(2,3));
plot_degrees(degree$NameSexAgeKernelWo1, ego_sex_age, ego_names_verbose, data$wgt);
plot_degrees(degree$NameSexAgeKernelWo1, ego_sex_age_race, ego_names_verbose_race, data$wgt);

lambda_names_kernel_wo1 <- apply(extract(fit_names_kernel_wo1)$lambda, c(2,3), mean);
rho_names_kernel_wo1 <- apply(extract(fit_names_kernel_wo1)$rho, c(2,3), mean);
beta_names_kernel_wo1 <- colMeans(extract(fit_names_kernel_wo1)$beta);

# ... 2. Plot Kernel ...

par(mfrow=c(1,2))
plot_kernel(lambda_names_kernel_wo1[1,], 
            rho_names_kernel_wo1[1,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Female Ego");
plot_kernel(lambda_names_kernel_wo1[2,], 
            rho_names_kernel_wo1[2,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Male Ego");

# ... 3. Plot Degree Difference ....

par(mfrow=c(2,3));
plot_degrees(degree$NameSexAgeKernelWo1 - degree$NameSexAgeKernel, ego_sex_age, ego_names_verbose, data$wgt, c(-100,200), xlab = "Degree Difference");
