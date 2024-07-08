load("Omni.RData");
require("rstan");
options(mc.cores = parallel::detectCores());

# ---------------------------------
# ---------- Stan Models ----------
# ---------------------------------

model_kernel <- stan_model(file = "stan/degree_kernel.stan");

# ----------------------------------
# ---------- Derived Data ----------
# ----------------------------------

age_mean_2014 <- sum(c(1:nrow(pop_raw)) * (pop_raw$Male+pop_raw$Female)/(sum(pop_raw$Male)+sum(pop_raw$Female)));

g_k <- c(rep(1, 6), rep(2, 6));
x_axis <- 2016-c(1916:2014);
mu_k <- sigma_k <- mu_k_ego <- sigma_k_ego <- sum_prob_k <- c();

# Alter Perspective: p(j in G_k, g_j)
prob_k <- rowSums(beta_names_alive)/(sum(pop_raw$Female) + sum(pop_raw$Male));

for(k in 1:nrow(beta_names_alive)) {
  # Alter Perspective: mu and sigma for p(a_j | g_j, j in G_k)
  y_axis <- beta_names_alive[k,];
  mu_k[k] <- x_axis%*%y_axis/sum(y_axis);
  sigma_k[k] <- sqrt((x_axis-mu_k[k])^2 %*% y_axis/sum(y_axis));
  
  # Ego Perspective: mu and sigma for p(j in G_k | a_j, g_j)/[sum_(a_j) p(j in G_k | a_j, g_j)]
  if(g_k[k] == 1) {
    # Ego Perspective: sum_(a_j) p(j in G_k | a_j, g_j)
    sum_prob_k[k] <- sum(beta_names_alive[k,]/rev(pop_raw$Female));
    y_axis <- (beta_names_alive[k,]/rev(pop_raw$Female))/sum_prob_k[k];
  }
  else {
    sum_prob_k[k] <- sum(beta_names_alive[k,]/rev(pop_raw$Male));
    y_axis <- (beta_names_alive[k,]/rev(pop_raw$Male))/sum_prob_k[k];
  }
  mu_k_ego[k] <- x_axis%*%y_axis/sum(y_axis);
  sigma_k_ego[k] <- sqrt((x_axis-mu_k_ego[k])^2 %*% y_axis/sum(y_axis));
}

# --------------------------------
# ---------- Simulation ----------
# --------------------------------

source("Code_Kernel/Kernel_Helper.R");

# Simulate Kernel Bandwidth, Gender Mixing, Overdispersion, and Random Degrees
rho_kernel_sim <- matrix(c(0.6, 0.4, 0.45, 0.55), nrow = 2, ncol = 2, byrow = T);
lambda_kernel_sim <- matrix(c(225, 144, 100, 256), nrow = 2, ncol = 2, byrow = T);
omega_kernel_sim <- 1/(1/omega_names_sex_age-1);
# degree_kernel_sim <- degree$NameSexAge;
degree_kernel_sim <- exp(rnorm(length(degree$NameSexAge), 6.17, 0.54));
hist(log(degree_kernel_sim));

# Simulate Regression-Based Degrees
beta_sim <- c(6.25, 0.12, -3);
eta_sim <- 0.7;
degree_kernel_reg_sim <- c();
for(n in 1:length(data$age)) { 
  mu_degree <- ex_log_degree_reg(beta_sim, data$age[n], as.numeric(data$Gender[n]), age_mean_sim);
  log_degree_kernel_reg_sim_n <- rnorm(1, mu_degree, eta_sim);
  degree_kernel_reg_sim[n] <- exp(log_degree_kernel_reg_sim_n);
  rm(mu_degree, log_degree_kernel_reg_sim_n);
}
degree_mean_sim <- mean(degree_kernel_reg_sim);
hist(log(degree_kernel_reg_sim));

# Without Degree Regression (Doesn't Work)
names_data_kernel <- simulate_kernel_continuous(data$age, 
                                                as.numeric(data$Gender), 
                                                g_k, 
                                                rho_kernel_sim, 
                                                lambda_kernel_sim, 
                                                omega_kernel_sim, 
                                                degree_kernel_reg_sim, 
                                                mu_k_ego, 
                                                sigma_k_ego, 
                                                sum_prob_k);
# With Degree Regression (Works)
names_data_kernel_reg <- simulate_kernel_continuous(data$age, 
                                                    as.numeric(data$Gender), 
                                                    g_k, 
                                                    rho_kernel_sim, 
                                                    lambda_kernel_sim, 
                                                    omega_kernel_sim, 
                                                    degree_kernel_reg_sim, 
                                                    mu_k_ego, 
                                                    sigma_k_ego, 
                                                    sum_prob_k, 
                                                    age_mean_2014, 
                                                    degree_mean_sim, 
                                                    beta_sim, 
                                                    eta_sim, 
                                                    T);
colnames(names_data_kernel) <- colnames(names_data_kernel_reg) <- colnames(names_data);

# -----------------------------
# ---------- Fitting ----------
# -----------------------------

# Initial values
init_data <- list();
for(i in 1:4) {
  init_data[[i]] <- list(log_d = rep(6, nrow(names_data)),
                        inv_omega = rep(5/6, ncol(names_data)),
                        log_lambda = matrix(rep(log(100), 4), nrow = 2, ncol = 2),
                        rho = matrix(rep(0.5, 4), nrow = 2, ncol = 2),
                        beta = rep(6, 0.1, -3),
                        log_eta_d = 0);
}

# Prepare Stan Data
mcmc_data_kernel_sim <- list(K = ncol(names_data_kernel), 
                             N = nrow(names_data_kernel), 
                             age_mean = age_mean_2014, 
                             degree_mean = degree_mean_sim, 
                             y = names_data_kernel, 
                             w = rep(1, nrow(names_data_kernel)),
                             age = data$age, 
                             g_n = as.numeric(data$Gender), 
                             g_k = g_k, 
                             mu_k = mu_k_ego, 
                             sigma_k = sigma_k_ego, 
                             sum_prob_k = sum_prob_k, 
                             mu_d = 6, sigma_d = 0.6, 
                             alpha_omega = 4.5, beta_omega = 0.5, 
                             mu_lambda = log(100), sigma_lambda = 0.5, 
                             alpha_rho = c(5,5), 
                             mu_beta = c(6, 0.1, -3), sigma_beta = rep(1, 3),
                             recall_power = 0, 
                             degree_regression = 1);
mcmc_data_kernel_sim_reg <- mcmc_data_kernel_sim;
mcmc_data_kernel_sim_reg$y = names_data_kernel_reg;

# Fit kernel model
fit_kernel_sim <- sampling(model_kernel, 
                           data = mcmc_data_kernel_sim, 
                           iter = 1000,
                           init = init_data);
# Fit kernel regression model
fit_kernel_sim_reg <- sampling(model_kernel, 
                               data = mcmc_data_kernel_sim_reg, 
                               iter = 1000,
                               init = init_data);
rm(mcmc_data_kernel_sim, mcmc_data_kernel_sim_reg, init_data);

# ------------------------------------
# ---------- Model Checking ----------
# ------------------------------------

source("Code_Kernel/Kernel_Helper.R");

# Check the Fits
# Degree...
plot(degree_kernel_reg_sim, 
     colMeans(extract(fit_kernel_sim)$d), 
     main = "Kernel Degree Fit", 
     ylim = c(0,4000), 
     xlim = c(0,4000), 
     xlab = "Simulated Degree with Regression", 
     ylab = "Fitted Degree", 
     col = "gray", 
     bty = 'l');
abline(a = 0, b= 1, lty = 2, col = 1);

plot(degree_kernel_reg_sim, 
     colMeans(extract(fit_kernel_sim_reg)$d), 
     main = "Kernel Degree Fit (Regression)", 
     ylim = c(0,4000),  
     xlim = c(0,4000), 
     xlab = "Simulated Degree with Regression", 
     ylab = "Fitted Degree with Regression", 
     col = "gray", 
     bty = 'l');
abline(a = 0, b= 1, lty = 2, col = 1);

# Overdispersion...
plot_kernel_comp(x_fit = extract(fit_kernel_sim)$omega, 
                 x_act = omega_kernel_sim, 
                 main = "Kernel Overdispersion Fit", 
                 xlab = "Simulated Overdispersion", 
                 ylab = "Fitted Overdispersion", 
                 func = median, 
                 xoff = 2, 
                 yoff = -70,
                 labs = colnames(names_data_kernel_continuous),
                 y_max = 100);
plot_kernel_comp(x_fit = extract(fit_kernel_sim_reg)$omega, 
                 x_act = omega_kernel_sim, 
                 main = "Kernel Overdispersion Fit (Degree Regression)", 
                 xlab = "Simulated Overdispersion", 
                 ylab = "Fitted Overdispersion", 
                 func = median, 
                 xoff = 2, 
                 yoff = -70,
                 labs = colnames(names_data_kernel_continuous),
                 y_max = 100);

#Beta
plot_kernel_comp(x_fit = extract(fit_kernel_sim_reg)$beta, 
                 x_act = beta_sim, 
                 main = "Kernel Beta Fit (Degree Regression)", 
                 xlab = "Simulated Beta", 
                 ylab = "Fitted Beta", 
                 func = median, 
                 xoff = 1, 
                 yoff = -1,
                 y_min = -7, 
                 y_max = 7,
                 labs = c("1","2","3"),
                 type = "beta");

# Lambda...
plot_kernel_comp_mat(x_fit = extract(fit_kernel_sim)$lambda, 
                     x_act = lambda_kernel_sim, 
                     main = "Kernel Lambda Fit", 
                     xlab = "Simulated Lambda", 
                     ylab = "Fitted Lambda", 
                     xoff = 5,
                     type = "lambda");
plot_kernel_comp_mat(x_fit = extract(fit_kernel_sim_reg)$lambda, 
                     x_act = lambda_kernel_sim, 
                     main = "Kernel Lambda Fit (Degree Regression)", 
                     xlab = "Simulated Lambda", 
                     ylab = "Fitted Lambda", 
                     xoff = 5,
                     type = "lambda");

# Gender Mixing...
plot_kernel_comp_mat(x_fit = extract(fit_kernel_sim)$rho, 
                     x_act = rho_kernel_sim, 
                     main = "Kernel Mixing Fit", 
                     xlab = "Simulated Mixing", 
                     ylab = "Fitted Mixing", 
                     laboff = 0, 
                     xoff = 0.02,
                     type = "rho");
plot_kernel_comp_mat(x_fit = extract(fit_kernel_sim_reg)$rho, 
                     x_act = rho_kernel_sim, 
                     main = "Kernel Mixing Fit (Degree Regression)", 
                     xlab = "Simulated Mixing", 
                     ylab = "Fitted Mixing", 
                     laboff = 0, 
                     xoff = 0.02,
                     type = "rho");
