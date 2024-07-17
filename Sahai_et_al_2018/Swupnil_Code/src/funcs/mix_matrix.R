# Plots the degree distributions
plot_degrees <- function(d, ego_, ego_labels, w = c(), 
                         xlim = c(0, 3000), xlab = "Degree",
                         dens = F, rm_rows = c()) {
  if (length(w) == 0) {
    w <- rep(1, length(d))
  }
  if (length(rm_rows) > 0) {
    d <- d[-rm_rows]
    ego_ <- ego_[-rm_rows]
    w <- w[-rm_rows]
  }
  d_cut <- d[(d < xlim[2]) & (d > xlim[1])]
  for (i in 1:length(unique(ego_))) {
    subset_ego_i <- which((ego_ == i) & (d < xlim[2]) & (d > xlim[1]))
    d_i <- d[subset_ego_i]
    w_i <- w[subset_ego_i]
    mean_d_i <- weighted.mean(d_i, w_i)
    median_d_i <- median(d_i)
    if (dens) {
      d_dens <- density(d_i)
      plot(d_dens$x,
        d_dens$y * length(subset_ego_i),
        xlab = xlab,
        xlim = xlim,
        main = ego_labels[i],
        type = "l",
        bty = "l"
      )
      text(median_d_i, max(d_dens$y) * 0.7, paste("Median = ", 
                                                  round(median_d_i), sep = ""),
        col = "blue", pos = 4
      )
    } else {
      h_i <- hist(d_i,
        xlab = xlab,
        xlim = xlim,
        main = ego_labels[i],
        breaks = seq(min(d_cut), max(d_cut), length.out = 20)
      )
      text(median_d_i, max(h_i$counts) * 0.7, paste("Median = ", 
                                                    round(median_d_i), sep = ""),
        col = "blue", pos = 4
      )
    }
    abline(v = median_d_i, lty = 2, col = "blue")
  }
}

# Plots a mixing matrix
plot_mixing <- function(M, y_max, ego_labels, title = "") {
  barplot(t(M[, 1:4]), horiz = T, beside = T, xlim = c(-0.6, 0.6), ylim = c(0, y_max), main = title)
  barplot(t(-M[, 4 + 1:4]), horiz = T, beside = T, add = T)
  for (i in 1:nrow(M)) {
    text(-0.5, 3 + 5 * (i - 1), ego_labels[i])
    text(0.5, 2.25 + 5 * (i - 1) - 1, "1-17")
    text(0.5, 2.25 + 5 * (i - 1) + .125, "18-24", col = "darkgray")
    text(0.5, 2.25 + 5 * (i - 1) + 1.125, "25-64", col = "gray")
    text(0.5, 2.25 + 5 * (i - 1) + 2.25, "65+", col = "lightgray")
    lines(c(-.42, -.42), 3 + 5 * (i - 1) + c(-1.5, 1.5))
  }
  text(-.5, y_max, "Ego Groups")
  text(0.2, y_max, "Male Alters")
  text(-0.2, y_max, "Female Alters")
  text(.5, y_max, "Alter Ages")
}

# Simulates responses from parameters of a mixing stan model
simulate_mixing <- function(d, omega, M, beta, ego) {
  library(MASS)
  y <- matrix(nrow = length(d), ncol = ncol(beta))

  for (i in 1:nrow(y)) {
    for (j in 1:ncol(y)) {
      mu_ij <- d[i] * M[ego[i], ] %*% beta[, j]
      y[i, j] <- rnegbin(1, mu_ij, omega[j] * mu_ij)
    }
  }

  return(y)
}

# Plots a comparison between actual and simulated responses
plot_comparison <- function(data1, data2, ego, ego_names) {
  for (e in 1:length(unique(ego))) {
    mu1 <- apply(data1[ego == e, ], 2, function(x) mean(x[which(x >= 0)]))
    mu2 <- apply(data2[ego == e, ], 2, function(x) mean(x[which(x >= 0)]))
    plot(mu1, mu2, pch = 21, bg = "grey", main = ego_names[e], xlim = c(min(mu1) - 0.5, max(mu1) + 0.5), ylim = c(min(mu2) - 0.2, max(mu2) + 0.2), xlab = "Mean Responses (Actual)", ylab = "Mean Responses (Simulated)")
    abline(a = 0, b = 1, lty = 2)
    text(mu1, mu2, pos = ifelse(mu1 > mu2, 1, 3), colnames(data1))
  }
}

# Plots bias and variance of actual and fitted matrix
plot_mix_comp <- function(fit, M_true, title = "Hello") {
  if (class(fit) == "list") {
    n <- length(fit)
    M_mu <- apply(extract(fit[[1]])$M, c(2, 3), mean) / n
    M_sd <- apply(extract(fit[[1]])$M, c(2, 3), sd) / n^2
    for (i in 2:n) {
      M_mu <- M_mu + apply(extract(fit[[i]])$M, c(2, 3), mean) / n
      M_sd <- M_sd + apply(extract(fit[[i]])$M, c(2, 3), sd) / n^2
    }
  } else {
    M_mu <- apply(extract(fit)$M, c(2, 3), mean)
    M_sd <- apply(extract(fit)$M, c(2, 3), sd)
  }

  plot(x = c(), y = c(), xlim = c(0, ncol(M_mu) + 1), ylim = c(0, nrow(M_mu) + 1), xlab = "Alter", ylab = "Ego", main = title)
  for (i in 1:nrow(M_mu)) {
    for (j in 1:ncol(M_mu)) {
      points(j, i, pch = 20, cex = abs(M_mu[i, j] - M_true[i, j]) / M_true[i, j], col = adjustcolor("black", alpha.f = .5))
      points(j, i, pch = 20, cex = M_sd[i, j] / M_true[i, j], col = adjustcolor("red", alpha.f = .5))
    }
  }
}
