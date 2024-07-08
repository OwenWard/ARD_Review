# Initialize Parameter Values That We'll Try to Recover
x_axis <- 2001-c(1900:2001);
mu_k_mccarty <- c();
sigma_k_mccarty <- c();
for(k in 1:nrow(beta_names_raw_mccarty)) {
  y_axis <- beta_names_raw_mccarty[k,];
  mu_k_mccarty <- c(mu_k_mccarty, x_axis%*%y_axis/sum(y_axis));
  sigma_k_mccarty <- c(sigma_k_mccarty, sqrt((x_axis-mu_k_mccarty[k])^2 %*% y_axis/sum(y_axis)));
}

# Simulate Responses from kernel model
names_data_kernel_mccarty <- simulate_kernel_continuous(mccarty, mu_k_mccarty, sigma_k_mccarty, rho_names_mccarty, 
                                     lambda_names_mccarty, omega_names_mccarty, degree_mccarty);
colnames(names_data_kernel_mccarty) <- colnames(names_data_mccarty);

g_k <- c(rep(1, 6), rep(2, 6));
mcmc_data_names_kernel_mccarty <- list(K = ncol(names_data_kernel_mccarty), N = nrow(names_data_kernel_mccarty), y = names_data_kernel_mccarty,
                               age = mccarty$Age, g_n = mccarty$Sex, g_k = g_k, mu_k = mu_k_mccarty, sigma_k = sigma_k_mccarty,
                               mu_d = 6.2, sigma_d = 0.5, alpha_omega = 10, beta_omega = 2, mu_lambda = log(100), 
                               sigma_lambda = 0.5, alpha_rho = c(5,5), recall = 0);
fit_names_kernel_mccarty <- sampling(kernel_continuous_fit, data = mcmc_data_names_kernel_mccarty, iter = 1000, chains = 4);

# Check the Fits
# Degree...
plot(degree_mccarty, apply(extract(fit_names_kernel_mccarty)$d, 2, mean), main = "Kernel Degree Fit", 
     ylim = c(0,4000), xlim = c(0,4000), xlab = "Simulated Degree", ylab="Fitted Degree", 
     col = "gray", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

# Overdispersion...
omega_kernel_fit_mccarty <- apply(extract(fit_names_kernel_mccarty)$omega, 2, median);
omega_25_mccarty <- apply(extract(fit_names_kernel_mccarty)$omega, 2, function(x) quantile(x,0.025));
omega_975_mccarty <- apply(extract(fit_names_kernel_mccarty)$omega, 2, function(x) quantile(x,0.975));
plot(omega_names_mccarty, omega_kernel_fit_mccarty, main = "Kernel Overdispersion Fit", 
     xlim=c(0, 5), ylim = c(0, 5), xlab = "Simulated Overdispersion", ylab="Fitted Overdispersion", 
     col = "red", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

for (i in 1:length(omega_names_mccarty)){
  arrows(omega_names_mccarty[i], omega_25_mccarty[i], omega_names_mccarty[i],
         omega_975_mccarty[i], length = 0, col = 2);
  points(omega_names_mccarty[i], omega_kernel_fit_mccarty[i], pch = 16, col = 1);
  text(omega_names_mccarty[i], omega_25_mccarty[i] - 0.25, pos = 4,
       bquote(omega[.(colnames(names_data_mccarty)[i], sep = "")]));
}

# Lambda...
lambda_kernel_fit_mccarty <- apply(extract(fit_names_kernel_mccarty)$lambda, c(2,3), mean);
lambda_25_mccarty <- apply(extract(fit_names_kernel_mccarty)$lambda, c(2,3), function(x) quantile(x,0.025));
lambda_975_mccarty <- apply(extract(fit_names_kernel_mccarty)$lambda, c(2,3), function(x) quantile(x,0.975));
plot(lambda_names_mccarty, lambda_kernel_fit_mccarty, main = "Kernel Lambda Fit", 
     xlim=c(500, 2000), ylim = c(500, 2000), xlab = "Simulated Lambda", ylab="Fitted Lambda", 
     col = "red", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

for (i in 1:2){
  for(j in 1:2){
    arrows(lambda_names_mccarty[i,j], lambda_25_mccarty[i,j], lambda_names_mccarty[i,j],
           lambda_975_mccarty[i,j], length = 0, col = 2);
    points(lambda_names_mccarty[i,j], lambda_kernel_fit_mccarty[i,j], pch = 16, col = 1);
    text(lambda_names_mccarty[i,j], lambda_25_mccarty[i,j], pos = 4,
         bquote(lambda[
           .(matrix(c("FF", "MF", "FM", "MM"), nrow = 2, ncol = 2)[i,j], sep = "")]));
  }
}