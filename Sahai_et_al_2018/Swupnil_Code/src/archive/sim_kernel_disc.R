require("rstan")
options(mc.cores = parallel::detectCores())

# Initialize Parameter Values That We'll Try to Recover
p_gg <- matrix(c(0.6, 0.4, 0.45, 0.55), nrow = 2, ncol = 2, byrow = T);
lambda_kernel <- matrix(c(225, 144, 100, 256), nrow = 2, ncol = 2, byrow = T);
omega_kernel <- 1/(1/rbeta(12, 10, 2)-1);
degree_kernel <- exp(rnorm(nrow(data),6.2,0.5));

# Simulate Responses from kernel model
names_data_kernel <- simulate_kernel(data, beta_names_raw, p_gg, lambda_kernel, omega_kernel, degree_kernel);
colnames(names_data_kernel) <- colnames(names_data);

mcmc_data_names_kernel <- list(K = ncol(names_data_kernel), N = nrow(names_data_kernel), y = names_data_kernel, J = 99, 
                            age = data$age, sex = as.numeric(data$Gender), Beta = beta_names_raw, s = c(100:2),
                            theta_d = c(6.2,0.5), theta_o = c(10,2), theta_l = c(log(100), 0.5), theta_m = c(5,5));
init_data <- list();
for(i in 1:4) { 
  init_data[[i]] <- list(log_d = rep(6.2, nrow(data)), inv_omega = rep(10/12, ncol(names_data)), 
                       log_lambda = matrix(log(100), nrow=2, ncol=2), M = matrix(0.5, nrow=2, ncol=2)); 
}
fit_names_kernel <- sampling(degree_kernel_fit, data = mcmc_data_names_kernel, iter = 1000, chains = 4, init = init_data);

# Check the Fits
# Degree...
plot(degree_kernel, apply(extract(fit_names_kernel)$d, 2, mean), main = "Kernel Degree Fit", 
     ylim = c(0,4000), xlim = c(0,4000), xlab = "Simulated Degree", ylab="Fitted Degree", 
     col = "gray", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

# Overdispersion...
omega_kernel_fit <- apply(extract(fit_names_kernel)$omega, 2, median);
omega_25 <- apply(extract(fit_names_kernel)$omega, 2, function(x) quantile(x,0.25));
omega_975 <- apply(extract(fit_names_kernel)$omega, 2, function(x) quantile(x,0.975));
plot(omega_kernel, omega_kernel_fit, main = "Kernel Overdispersion Fit", 
     xlim=c(1, 20), ylim = c(0, 20), xlab = "Simulated Overdispersion", ylab="Fitted Overdispersion", 
     col = "red", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

for (i in 1:length(omega_kernel)){
    arrows(omega_kernel[i], omega_25[i], omega_kernel[i],
           omega_975[i], length = 0, col = 2);
    points(omega_kernel[i], omega_kernel_fit[i], pch = 16, col = 1);
    text(omega_kernel[i], omega_25[i] - 0.25, pos = 4,
         bquote(omega[.(colnames(names_data_kernel)[i], sep = "")]));
}

# Lambda...
lambda_kernel_fit <- apply(extract(fit_names_kernel)$lambda, c(2,3), mean);
lambda_25 <- apply(extract(fit_names_kernel)$lambda, c(2,3), function(x) quantile(x,0.025));
lambda_975 <- apply(extract(fit_names_kernel)$lambda, c(2,3), function(x) quantile(x,0.975));
plot(lambda_kernel, lambda_kernel_fit, main = "Kernel Lambda Fit", 
     xlim=c(90, 270), ylim = c(100, 17000), xlab = "Simulated Lambda", ylab="Fitted Lambda", 
     col = "red", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

for (i in 1:2){
  for(j in 1:2){
    arrows(lambda_kernel[i,j], lambda_25[i,j], lambda_kernel[i,j],
           lambda_975[i,j], length = 0, col = 2);
    points(lambda_kernel[i,j], lambda_kernel_fit[i,j], pch = 16, col = 1);
    text(lambda_kernel[i,j], lambda_25[i,j], pos = 4,
         bquote(lambda[
               .(matrix(c("FF", "MF", "FM", "MM"), nrow = 2, ncol = 2)[i,j], sep = "")]));
  }
}
