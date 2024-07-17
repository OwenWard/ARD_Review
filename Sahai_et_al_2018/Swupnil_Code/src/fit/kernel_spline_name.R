# ----------------------------------------------------------------------------
# ---------- Estimating Degree from Names Using Kernel Spline Model ----------
# ----------------------------------------------------------------------------

require("rstan")
require("splines")
options(mc.cores = parallel::detectCores())

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

age_grid <- seq(min(data$age), max(data$age), 1);
knots <- quantile(age_grid, probs = seq(0, 1, .1));
N_K <- length(knots);
D <- 3;
N_B <- N_K + D - 1;

# Initial Values
init_data <- list();
for(i in 1:4) {
  init_data[[i]] <- list(log_d = rep(6, nrow(names_data)),
                         inv_omega = rep(5/6, ncol(names_data)),
                         a_raw = array(0, dim = c(2,2,N_B)),
                         a0 = matrix(0, nrow = 2, ncol = 2),
                         tau = matrix(0.5, nrow = 2, ncol = 2),
                         beta = c(6, 0.1, -3),
                         log_eta = 0);
}

# Run the Model in Stan

mcmc_data <- list(K_1 = ncol(names_data),
                  K_2 = 0,
                  N = nrow(names_data),  
                  KG = ncol(names_data),
                  N_K = N_K,
                  N_S = length(age_grid),
                  D = D,
                  age_mean = age_mean_2014, 
                  Y = as.matrix(names_data),  
                  w = data$wgt,
                  age = data$age, 
                  g_n = as.numeric(data$Gender), 
                  g_k = c(rep(1, 6), rep(2, 6)), 
                  mu_k = mu_k_name, 
                  sigma_k = sigma_k_name, 
                  sum_prob_k = sum_prob_k_name, 
                  knots = knots,
                  X = age_grid,
                  mu_d = 6, sigma_d = 0.6, 
                  alpha_omega = 4.5, beta_omega = 0.5, 
                  alpha_rho = c(5,5), 
                  mu_beta = c(6, 0.1, -3), sigma_beta = rep(1, 3),
                  recall_power = 0, 
                  degree_regression = 1);
fit_names_kernel_spline <- sampling(model_kernel_spline, 
                                    data = mcmc_data, 
                                    iter = 2000,
                                    init = init_data);

# -------------------------------
# ---------- Assesment ----------
# -------------------------------

# ... 1. Plot Degree Distributions ....

degree$NameSexAgeKernelSpline <- apply(extract(fit_names_kernel_spline)$d,2,median);

par(mfrow=c(2,3));
plot_degrees(degree$NameSexAgeKernelSpline, ego_sex_age, ego_names_verbose, data$wgt, rm = rm_na_names);
plot_degrees(degree$NameSexAgeKernelSpline, ego_sex_age_race, ego_names_verbose_race, data$wgt, rm = rm_na_names);

rho_names_kernel_spline <- apply(extract(fit_names_kernel_spline)$rho, c(2,3), mean);
beta_names_kernel_spline <- colMeans(extract(fit_names_kernel_spline)$beta);
eta_names_kernel_spline <- mean(extract(fit_names_kernel_spline)$eta);

# ... 2. Plot Kernel Bandwidth ...

plot_kernel_spline(age_grid, fit_names_kernel_spline, overlay = T, ylim = c(18,55));
plot_kernel_spline(age_grid, fit_names_kernel_spline, ylim = c(18,55));

# ... 3. Plot Kernel ...

par(mfrow=c(3,2));
plot_kernels(fit_ = fit_names_kernel_spline, 
             rho_ = rho_names_kernel_spline,
             y_cent = 70,
             y_index = floor(70)-17,
             y_min = 0,
             y_max = 100);
plot_kernels(fit_ = fit_names_kernel_spline, 
             rho_ = rho_names_kernel_spline,
             y_cent = age_mean_2014,
             y_index = floor(age_mean_2014)-17,
             y_min = 0,
             y_max = 100);
plot_kernels(fit_ = fit_names_kernel_spline, 
             rho_ = rho_names_kernel_spline,
             y_cent = 21,
             y_index = floor(21)-17,
             y_min = 0,
             y_max = 100);

# ... 4. Plot Residuals ....

names_data_pred <- simulate_kernel_continuous(data$age, 
                                              as.numeric(data$Gender), 
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
