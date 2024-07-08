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
  init_data[[i]] <- list(log_d = rep(6, nrow(names_data)),
                         inv_omega = rep(5/6, ncol(names_data)),
                         log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                         rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                         beta = c(6, 0.1, -3),
                         log_eta = 0);
}

# Run the Model in Stan
mcmc_data <- list(K_1 = ncol(names_data),
                  K_2 = 0,
                  N = nrow(names_data),  
                  age_mean = age_mean_2014, 
                  Y = as.matrix(names_data),  
                  w = data$wgt,
                  age = data$age, 
                  g_n = as.numeric(data$Gender), ## these should be 1/2?
                  g_k = c(rep(1, 6), rep(2, 6)), 
                  mu_k = mu_k_name, sigma_k = sigma_k_name, 
                  sum_prob_k = sum_prob_k_name, 
                  mu_d = 6, sigma_d = 0.6, 
                  alpha_omega = 4.5, beta_omega = 0.5, 
                  mu_lambda = log(100), sigma_lambda = 0.5, 
                  alpha_rho = c(5,5), 
                  mu_beta = c(6, 0.1, -3), sigma_beta = rep(1, 3),
                  recall_power = 0, 
                  degree_regression = 1);
fit_names_kernel <- sampling(model_kernel, 
                             data = mcmc_data, 
                             iter = 2000,
                             init = init_data);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

# ... 1. Plot Degree Distributions ....

degree$NameSexAgeKernel <- round(colMeans(extract(fit_names_kernel)$d));
omega_names_kernel <- colMeans(extract(fit_names_kernel)$inv_omega);

par(mfrow=c(2,3));
plot_degrees(degree$NameSexAgeKernel, ego_sex_age, ego_names_verbose, data$wgt);
plot_degrees(degree$NameSexAgeKernel, ego_sex_age_race, ego_names_verbose_race, data$wgt);

lambda_names_kernel <- apply(extract(fit_names_kernel)$lambda, c(2,3), mean);
rho_names_kernel <- apply(extract(fit_names_kernel)$rho, c(2,3), mean);
beta_names_kernel <- colMeans(extract(fit_names_kernel)$beta);
eta_names_kernel <- mean(extract(fit_names_kernel)$eta);

# ... 2. Plot Kernel ...

par(mfrow=c(1,2))
plot_kernel(lambda_names_kernel[1,], 
            rho_names_kernel[1,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Female Ego");
plot_kernel(lambda_names_kernel[2,], 
            rho_names_kernel[2,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Male Ego");

# ... 3. Plot Residuals ....

names_data_pred <- simulate_kernel_continuous(data$age, 
                                              as.numeric(data$Gender),
                                              ## these turn to NA, causes the 
                                              ## issues below
                                              g_k, 
                                              rho_names_kernel, 
                                              lambda_names_kernel, 
                                              omega_names_kernel,
                                              degree$NameSexAgeKernel,
                                              mu_k_name, 
                                              sigma_k_name,
                                              sum_prob_k_name);
colnames(names_data_pred) <- colnames(names_data);
resid_names <- (names_data - names_data_pred) / ifelse(names_data_pred == 0, mean(names_data_pred), sqrt(names_data_pred));

par(mfrow = c(1,2))
for(j in 1:ncol(cancer_data)) {
    for(g in 1:2) {
      respondents <- which((cancer_data[,j] != -1) & (data$Gender==ifelse(g==1,"Female","Male")));
      plot_title <- paste(ifelse(g==1, "Female", "Male"), "Ego");
      plot(as.factor(cancer_data[respondents, j]), 
           rowMeans(resid_names[respondents, 1:6]) - rowMeans(resid_names[respondents, 7:12]), 
           ylab = "Residual Differential (Female-Male)",
           xlab = paste(colnames(cancer_data)[j], " Known"), 
           main = plot_title, bty = 'l', col = 'gray',
           ylim = c(-3,4));
      rm(respondents, plot_title);
    }
}
