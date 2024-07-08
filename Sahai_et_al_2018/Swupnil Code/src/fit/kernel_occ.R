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
  init_data[[i]] <- list(log_d = rep(6, nrow(occ_data)),
                         inv_omega = rep(5/6, ncol(occ_data)),
                         log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                         rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                         beta = c(6, 0.1, -3),
                         log_eta = 0);
}

# Run the Model in Stan
mcmc_data <- list(K_1 = 0,
                  K_2 = ncol(occ_data),
                  N = nrow(occ_data),   
                  age_mean = age_mean_2014, 
                  Y = as.matrix(occ_data),  
                  w = data$wgt,
                  age = data$age, 
                  g_n = as.numeric(data$Gender), ## does this need to be 1/2?
                  g_k = g_k_occ, 
                  mu_k = mu_k_occ, 
                  sigma_k = sigma_k_occ, 
                  sum_prob_k = sum_prob_k_occ, 
                  mu_d = 6, 
                  sigma_d = 0.6, 
                  alpha_omega = 4.5, 
                  beta_omega = 0.5, 
                  mu_lambda = log(100), 
                  sigma_lambda = 0.5, 
                  alpha_rho = c(5,5), 
                  mu_beta = c(6, 0.1, -3), 
                  sigma_beta = rep(1, 3),
                  recall_power = 0, 
                  degree_regression = 1);
fit_occs_kernel <- sampling(model_kernel, 
                            data = mcmc_data, 
                            iter = 2000,
                            init = init_data);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

# ... 1. Plot Degree Distributions ....

degree$OccSexAgeKernel <- round(colMeans(extract(fit_occs_kernel)$d));
omega_occs_kernel <- colMeans(extract(fit_occs_kernel)$inv_omega);

par(mfrow=c(2,3));
plot_degrees(degree$OccSexAgeKernel, ego_sex_age, ego_names_verbose, data$wgt);
plot_degrees(degree$OccSexAgeKernel, ego_sex_age_race, ego_names_verbose_race, data$wgt);

lambda_occs_kernel <- apply(extract(fit_occs_kernel)$lambda, c(2,3), mean);
rho_occs_kernel <- apply(extract(fit_occs_kernel)$rho, c(2,3), mean);
beta_occs_kernel <- colMeans(extract(fit_occs_kernel)$beta);
eta_occs_kernel <- mean(extract(fit_occs_kernel)$eta);

# ... 2. Plot Kernel ...

par(mfrow=c(1,2))
plot_kernel(lambda_occs_kernel[1,], 
            rho_occs_kernel[1,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Female Ego");
plot_kernel(lambda_occs_kernel[2,], 
            rho_occs_kernel[2,], 
            y_cent = age_mean_2014,
            y_min = 0,
            y_max = 100,
            main = "Male Ego");

# ... 3. Plot Residuals ....

occ_data_pred <- simulate_kernel_continuous(data$age, 
                                              as.numeric(data$Gender), 
                                              g_k, 
                                              rho_occs_kernel, 
                                              lambda_occs_kernel, 
                                              omega_occs_kernel,
                                              degree$OccSexAgeKernel,
                                              mu_k_name, 
                                              sigma_k_name,
                                              sum_prob_k_name);
colnames(occ_data_pred) <- colnames(occ_data);
resid_occs <- (occ_data - occ_data_pred) / ifelse(occ_data_pred == 0, mean(occ_data_pred), sqrt(occ_data_pred));

par(mfrow = c(1,2))
for(j in 1:ncol(cancer_data)) {
    for(g in 1:2) {
      respondents <- which((cancer_data[,j] != -1) & (data$Gender==ifelse(g==1,"Female","Male")));
      plot_title <- paste(ifelse(g==1, "Female", "Male"), "Ego");
      plot(as.factor(cancer_data[respondents, j]), 
           rowMeans(resid_occs[respondents, 1:6]) - rowMeans(resid_occs[respondents, 7:8]), 
           ## changed last part from 7:12 to 7:8, need to check this
           ylab = "Residual Differential (Female-Male)",
           xlab = paste(colnames(cancer_data)[j], " Known"), 
           main = plot_title, bty = 'l', col = 'gray',
           ylim = c(-3,4));
      rm(respondents, plot_title);
    }
}
