# Plot Histograms to Assess Normality
path <- file.path("Option1_Histograms.png");
png(filename = path, width = 1200, height = 1300);
par(mfrow=c(4,3));
x_axis <- 2016-c(1916:2014);
for(k in 1:nrow(beta_names_raw)) {
  y_axis <- beta_names_raw[k,];
  plot(x_axis, y_axis, type = 'l', xlab = "Age", ylab = "Density", main = rownames(beta_names_raw)[k]);
}
dev.off()

path <- file.path("Option2_Histograms.png");
png(filename = path, width = 1200, height = 1600);
par(mfrow=c(5,3));
x_axis <- 2016-c(1916:2014);
for(k in 1:nrow(beta_names_raw)) {
  y_axis <- beta_names_alive[k,]/sum(beta_names_alive[k,]);
  plot(x_axis, y_axis, type = 'l', xlab = "Age", ylab = "Density", main = rownames(beta_names_raw)[k]);
}
plot(x_axis, pop_raw$Male/sum(pop_raw$Male), type = 'l', xlab = "Age", ylab = "Density", main = "Male");
plot(x_axis, pop_raw$Female/sum(pop_raw$Female), type = 'l', xlab = "Age", ylab = "Density", main = "Female");
dev.off()

# Initialize Parameter Values That We'll Try to Recover
x_axis <- 2016-c(1916:2014);
mu_k <- c();
sigma_k <- c();
for(k in 1:nrow(beta_names_raw)) {
  y_axis <- beta_names_raw[k,];
  mu_k <- c(mu_k, x_axis%*%y_axis/sum(y_axis));
  sigma_k <- c(sigma_k, sqrt((x_axis-mu_k[k])^2 %*% y_axis/sum(y_axis)));
}

rho_kernel <- matrix(c(0.6, 0.4, 0.45, 0.55), nrow = 2, ncol = 2);
lambda_kernel <- matrix(c(225, 144, 100, 256), nrow = 2, ncol = 2);
omega_kernel <- 1/(1/rbeta(12, 10, 2)-1);
degree_kernel <- exp(rnorm(nrow(data),6.2,0.5));

# Simulate Responses from kernel model
names_data_kernel <- simulate_kernel_continuous(data, mu_k, sigma_k, rho_kernel, lambda_kernel, omega_kernel, degree_kernel);
colnames(names_data_kernel) <- colnames(names_data);

mcmc_data_names_kernel <- list(K = ncol(names_data_kernel), N = nrow(names_data_kernel), y = names_data_kernel,
                            age = data$age, sex = as.numeric(data$Gender), mu = mu_k, sigma = sigma_k,
                            theta_d = c(6.2,0.5), theta_o = c(10,2), theta_l = c(log(100), 0.5));
fit_names_kernel <- sampling(kernel_continuous_fit, data = mcmc_data_names_kernel, iter = 1000, chains = 4);

# Check the Fits
# Degree...
plot(degree_kernel, apply(extract(fit_names_kernel)$d, 2, mean), main = "Kernel Degree Fit", 
     ylim = c(0,4000), xlim = c(0,4000), xlab = "Simulated Degree", ylab="Fitted Degree", 
     col = "gray", bty = 'l')
abline(a = 0, b= 1, lty = 2, col = 1)

# Overdispersion...
omega_kernel_fit <- apply(extract(fit_names_kernel)$omega, 2, median);
omega_25 <- apply(extract(fit_names_kernel)$omega, 2, function(x) quantile(x,0.025));
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
     xlim=c(90, 270), ylim = c(90, 270), xlab = "Simulated Lambda", ylab="Fitted Lambda", 
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
