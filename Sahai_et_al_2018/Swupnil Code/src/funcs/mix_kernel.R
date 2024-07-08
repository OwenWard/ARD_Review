library(MASS);
source(paste0(code_path, "src/funcs/mix_matrix.R"))

# --------------------------------------
# ---------- Helper Functions ----------
# --------------------------------------

# Returns the value of the discrete kernel evaluated at a1 and a2
kernel <- function(a1, a2, lambda) {
  return(exp(-0.5*(a1 - a2)^2/lambda))
}

# Returns the integration of the product of two normals
norm_prod_int <- function(mu1, mu2, var1, var2) {
  y <- 1/sqrt(2 * pi * (var1 + var2)) * exp(-(mu1 - mu2)^2/(2 * (var1 + var2)));
  return(y);
}

# Returns the expectation of the product of two normals
norm_prod_ex <- function(mu1, mu2, var1, var2) {
  y <- (mu1 * var2 + mu2 * var1)/(var1 + var2);
  return(y);
}

# Returns the variance of the product of two normals
norm_prod_var <- function(mu1, mu2, var1, var2) {
  y <- (var1 * var2)/(var1 + var2);
  return(y);
}

# Expectation of alter degree kernel over all alter ages.
kernel_mean <- function(a, mu, sigma, lambda) {
  ex_aj <- norm_prod_int(a, mu, lambda, sigma^2) * norm_prod_ex(a, mu, lambda, sigma^2);
  return(ex_aj);
}

# Expectation of alter degree squared kernel over all alter ages.
kernel_squared_mean <- function(a, mu, sigma, lambda) {
  ex_aj_2 <- norm_prod_var(a, mu, lambda, sigma^2) + norm_prod_ex(a, mu, lambda, sigma^2)^2;
  ex_aj_2 <- ex_aj_2 * norm_prod_int(a, mu, lambda, sigma^2);
  return(ex_aj_2);
}

# Expectation of the log degree under the regression assumption
ex_log_degree_reg <- function(beta, age, gender, age_mean) {
  ex_log_degree <- beta[1] + beta[2] * (gender - 1) - exp(beta[3]) * ((age - age_mean)/age_mean)^2;
  return(ex_log_degree);
}

# ------------------------------------------
# ---------- Simulation Functions ----------
# ------------------------------------------

# Simulates responses from the kernel model with discrete age.
simulate_kernel_discrete <- function(age, gender, gender_k, rho, lambda, omega, degree, beta_k) {
  y <- matrix(nrow = length(d), ncol = length(omega));
  
  for(n in 1:nrow(y)) {
    a_n <- age[n];
    for(k in 1:ncol(y)) {
      gn <- gender[n];
      gk <- gender_k[k];
      
      kernel_nk <- kernel(a_n, c(100:2), lambda[gn,gk]);
      C_nk <- sum(kernel_nk);
      
      mu_nk <- rho[gn,gk] * degree[n] / C_nk * (beta_k[k,] %*% kernel_nk);
      y[n,k] <- rnegbin(1, mu_nk, omega[k] * mu_nk);
    }
  }
  
  return(y);
}

# Simulates responses from the kernel model with continuous age, nonconstant degree.
simulate_kernel_continuous <- function(age, gender, gender_k, rho, lambda, omega, degree, mu_k, sigma_k, prop_k, 
                                       age_mean = NA, degree_mean = NA, beta = NA, eta = NA, is_degree_varying = FALSE) {
  y <- matrix(nrow = length(degree), ncol = length(omega));
  
  for(n in 1:nrow(y)) {
    a_n <- age[n];
    gn <- gender[n];
    for(k in 1:ncol(y)) {
      gk <- gender_k[k];
      mu_nk <- prop_k[k];
      mu_nk <- mu_nk * rho[gn,gk];
      mu_nk <- mu_nk * norm_prod_int(a_n, mu_k[k], lambda[gn,gk], sigma_k[k]^2);
      mu_nk <- mu_nk * degree[n];
      y[n,k] <- rnegbin(1, mu_nk, omega[k] * mu_nk);
    }
  }
  
  return(y);
}

# ----------------------------------------
# ---------- Plotting Functions ----------
# ----------------------------------------

# Plots the kernel model results against the true data
plot_kernel_comp <- function(x_fit, 
                             x_act, 
                             main, 
                             xlab, 
                             ylab, 
                             func = mean,
                             xoff = 0, 
                             yoff = 0, 
                             labs = c(""),
                             laboff = -0.25, 
                             type = "omega", 
                             y_min = 0, 
                             y_max = 1) {
  # Generate quantiles
  x_mid <- apply(x_fit, 2, func);
  x_25 <- apply(x_fit, 2, function(x) quantile(x,0.25));
  x_975 <- apply(x_fit, 2, function(x) quantile(x,0.975));
  
  # Plot the data
  plot(x = x_act, 
       y = x_mid, 
       main = main, 
       xlim = c(min(x_act), max(x_act) + xoff), 
       ylim = c(y_min, y_max), 
       xlab = xlab, 
       ylab = ylab, 
       col = "red", 
       bty = 'l');
  abline(a = 0, 
         b = 1, 
         lty = 2, 
         col = 1);
  
  # Add the arrows
  for (i in 1:length(x_act)){
    if(type == "beta") { 
      labels <- bquote(beta[.(labs[i], sep = "")]); 
    }
    else { 
      labels <- bquote(omega[.(labs[i], sep = "")]); 
    }
    
    # Plot the arrows, points, and annotations
    arrows(x_act[i], 
           x_25[i], 
           x_act[i], 
           x_975[i], 
           length = 0, 
           col = 2);
    points(x_act[i], 
           x_mid[i], 
           pch = 16, 
           col = 1);
    text(x_act[i], 
         x_25[i] + laboff, 
         pos = 4, 
         labels);
  }
}

# Compares the fitted kernel lengthscales / gender mixing matrix values against their true values
plot_kernel_comp_mat <- function(x_fit, 
                                 x_act, 
                                 main, 
                                 xlab,
                                 ylab, 
                                 xoff = 0, 
                                 yoff = 0, 
                                 laboff = -0.25, 
                                 type = "lambda") {
  # Generate quantiles
  x_mid <- colMeans(x_fit);
  x_25 <- apply(x_fit, c(2,3), function(x) quantile(x,0.25));
  x_975 <- apply(x_fit, c(2,3), function(x) quantile(x,0.975));
  
  # Plot the data
  plot(x = x_act, 
       y = x_mid, 
       main = main, 
       xlim = c(min(x_act), max(x_act) + xoff), 
       ylim = c(min(x_25), max(x_975) + yoff), 
       xlab = xlab, 
       ylab = ylab, 
       col = "red", 
       bty = 'l');
  abline(a = 0, 
         b = 1, 
         lty = 2, 
         col = 1);
  
  # Add the arrows
  for (i in 1:nrow(x_act)){
    for(j in 1:ncol(x_act)){
      if(type == "lambda") {
        labels <- bquote(lambda[.(matrix(c("FF", "MF", "FM", "MM"), nrow = 2, ncol = 2)[i,j], sep = "")]);
      }
      else if(type == "rho") {
        labels <- bquote(rho[.(matrix(c("FF", "MF", "FM", "MM"), nrow = 2, ncol = 2)[i,j], sep = "")]);
      }
      
      # Plot the arrows, points, and annotations
      arrows(x_act[i,j], 
             x_25[i,j], 
             x_act[i,j], 
             x_975[i,j], 
             length = 0, 
             col = 2);
      points(x = x_act[i,j], 
             y = x_mid[i,j], 
             pch = 16, 
             col = 1);
      text(x = x_act[i,j], 
           y = x_25[i,j] + laboff, 
           pos = 4, 
           labels);
    }
  }
}

# Plot the fitted kernels adjusted by mixing rate.
plot_kernel <- function(lambda, rho = c(.5, .5),  y_cent = 50, y_min = 0, y_max = 100,
                        main = "Male Ego Network", ylab = "Age",
                        labels = c("Female Alter", "Male Alter"),
                        top = T) {
  grid <- c(y_min:y_max);
  sd <- sqrt(lambda);
  left <- dnorm(grid, y_cent, sd[1]);
  right <- dnorm(grid, y_cent, sd[2]);
  
  left <- left/sum(left) * rho[1];
  right <- right/sum(right) * rho[2];
  
  # Initialize the plot
  x_lim <- c(-max(left), max(right));
  plot(c(), 
       c(), 
       ylim = c(y_min, y_max), 
       xlim = c(-0.012,0.012),
       xlab = "Density",
       ylab = ylab,
       main = paste(round(y_cent), "Year Old", main),
       bty = 'l');
  
  # Fill in the kernels
  poly_x_left <- c(0, -left, 0);
  poly_x_right <- c(0, right, 0);
  poly_y <- c(y_min, grid, y_max);
  polygon(poly_x_left, poly_y, col = adjustcolor(2, alpha.f = 0.1));
  polygon(poly_x_right, poly_y, col = adjustcolor(4, alpha.f = 0.1));
  
  # Fill in 1 SD
  left_1sd <- c(max(floor(y_cent - sd[1]), y_min):min(ceiling(y_cent + sd[1]), y_max)) - y_min + 1;
  right_1sd <- c(max(floor(y_cent - sd[2]), y_min):min(ceiling(y_cent + sd[2]), y_max)) - y_min + 1;
  poly_x_left_1sd <- c(0, -left[left_1sd], 0);
  poly_x_right_1sd <- c(0, right[right_1sd], 0);
  poly_y_left_1sd <- grid[c(left_1sd[1], left_1sd, tail(left_1sd, 1))];
  poly_y_right_1sd <- grid[c(right_1sd[1], right_1sd, tail(right_1sd, 1))];
  polygon(poly_x_left_1sd, poly_y_left_1sd, col = adjustcolor(2, alpha.f = 0.4));
  polygon(poly_x_right_1sd, poly_y_right_1sd, col = adjustcolor(4, alpha.f = 0.4));
  
  # Add lines
  lines(-left, grid);
  lines(right, grid);
  lines(x_lim, rep(y_cent, 2), lty = 2);
  
  # Plot the text
  labels[1] <- paste(labels[1], "\n", round(rho[1] * 100, 1), "%", sep = "");
  labels[2] <- paste(labels[2], "\n", round(rho[2] * 100, 1), "%", sep = "");
  text(-0.012, y_max * ifelse(top,0.95,0.05), labels[1], pos = 4);
  text(0.012, y_max * ifelse(top,0.95,0.05), labels[2], pos = 2);
}

# Plots female and male kernels.
plot_kernels <- function(fit_,  rho_,  y_cent = 50, y_index = 50,
                         y_min = 0, y_max = 100, top = T) {
  lambda_mean_ <- apply(extract(fit_)$lambda, c(2,3,4), mean);
  plot_kernel(c(lambda_mean_[1,1,y_index],
                lambda_mean_[1,2,y_index]), 
              rho = rho_[1,], 
              y_cent = y_cent,
              y_min = y_min,
              y_max = y_max,
              main = "Female Ego",
              top = top);
  plot_kernel(c(lambda_mean_[2,1,y_index],
                lambda_mean_[2,2,y_index]),
              rho = rho_[2,], 
              y_cent = y_cent,
              y_min = y_min,
              y_max = y_max,
              main = "Male Ego",
              top = top);
}

# Plots the kernel spline.
plot_kernel_spline <- function(age_grid_, fit_, lambda_const_ = NA, overlay = F,
                               ci_samps = 400, ci = T, ci_alpha = 0.01, ylim = c(20,50)) {
  if (overlay) {
    labels <- matrix(c("Female Ego", "Female Ego", "Male Ego", "Male Ego"), 
                     nrow = 2, ncol = 2, byrow = TRUE);
    par(mfrow=c(2,1));
  }
  else {
    labels <- matrix(c("Female to Female", "Female to Male", "Male to Female", "Male to Male"), 
                     nrow = 2, ncol = 2, byrow = TRUE);
    par(mfrow=c(2,2));
  }
  lambda_ <- extract(fit_)$lambda;
  N <- dim(lambda_)[1];
  K <- dim(lambda_)[4];
  idx <- sample(1:N, ci_samps);
  for(i in 1:2) {
    for(j in 1:2) {
      lambda_med_ij <- c();
      for(k in 1:K) { lambda_med_ij[k] <- median(lambda_[,i,j,k]); }
      if ((overlay) & (j == 2)) {
        lines(age_grid, 
             sqrt(lambda_med_ij),
             lwd = 2,
             col = 2 * j);
      }
      else {
        plot(age_grid, 
             sqrt(lambda_med_ij), 
             ylab = "Mixing Kernel SD", 
             xlab = "Ego Age", 
             xlim = c(18,70),
             ylim = ylim, 
             main = labels[i,j], 
             type = 'l', 
             bty = 'l',
             lwd = 2,
             col = 2 * j);
      }
      if (ci) {
        for(n in idx) {
          lines(age_grid, 
                sqrt(lambda_[n,i,j,]), 
                col = adjustcolor(2*j, 0.1));
        }
      }
      if(!is.na(lambda_const_)) {
        abline(h = sqrt(lambda_const_[i,j]), col = 2);
      }
    }
  }
}
